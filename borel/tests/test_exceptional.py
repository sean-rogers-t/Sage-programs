"""Exact verification of the exceptional generator catalogue.

Run from ``borel/`` with ``PYTHONPATH=src sage -python -m unittest
tests.test_exceptional -v``.  Every row is recomputed, including all of psi12.
``BOREL_EXCEPTIONAL_REPORT`` optionally names a JSON verification report.

The frozen fixture comes from the earlier, separately implemented BGG table
audit.  It is not a second mathematical algorithm.  A further comparison uses
the root/coroot Chevalley formula below, without calling the library's divisor
operations or divided differences.  Runtime package data contains inputs only.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import json
import hashlib
import os
from pathlib import Path
import platform
import time
import unittest

from sage.all import PolynomialRing, QQ, WeylGroup, prod
from sage.env import SAGE_VERSION

from sage_borel import (
    exceptional_borel, exceptional_generators, exceptional_image,
    exceptional_polynomial,
)


FIXTURE_PATH = Path(__file__).parent / "data" / "exceptional_expected.json"
CORRECTED_IDS = {
    "F4/P2/row5", "E6/P3/row6", "E6/P6/row2", "E6/P6/row3",
}


def independent_polynomial(row, ring):
    """Reconstruct the audited input, without runtime generator builders.

    Elementary functions are literal subset sums, unlike the runtime dynamic
    programme.  Simple roots are obtained from ambient dot products.  The psi
    formula is copied from the mathematical audit specification.
    """
    x = tuple(ring.gens())
    spec = row["input"]
    ys = {int(i): sum(QQ(c) * a for c, a in zip(vector, x))
          for i, vector in row.get("y_definitions", {}).items()}
    operation = spec["operation"]
    if operation == "linear":
        return ring(ys[spec["variable"]])
    if operation == "omega":
        return x[spec["index"] - 1]
    if operation == "elementary":
        values = [ys[i] ** (2 if spec.get("squared") else 1)
                  for i in spec["variables"]]
        return ring(sum(prod(subset) for subset in
                        combinations(values, spec["degree"])))
    if operation == "g2quadratic":
        return (x[1] - QQ(3) / 2 * x[0]) ** 2
    if operation == "g2product":
        return x[0] * x[1] - x[0] ** 2
    if operation == "psi":
        W = WeylGroup(row["kind"], prefix="input_s")
        roots = W.domain().simple_roots()
        a = {i: sum(2 * roots[i].dot_product(roots[j]) /
                    roots[j].dot_product(roots[j]) * x[j - 1]
                    for j in range(1, len(x) + 1))
             for i in range(1, len(x) + 1)}
        xx = [5*a[1]+4*a[3]+3*a[4]+2*a[5]+a[6],
              -a[1]+4*a[3]+3*a[4]+2*a[5]+a[6],
              -a[1]-2*a[3]+3*a[4]+2*a[5]+a[6],
              -a[1]-2*a[3]-3*a[4]+2*a[5]+a[6],
              -a[1]-2*a[3]-3*a[4]-4*a[5]+a[6],
              -a[1]-2*a[3]-3*a[4]-4*a[5]-5*a[6]]
        z = -3*(a[1]+2*a[2]+2*a[3]+3*a[4]+2*a[5]+a[6])
        alphabet = [v+z for v in xx] + [v-z for v in xx]
        alphabet += [-u-v for u, v in combinations(xx, 2)]
        return ring(sum(v ** spec["degree"] for v in alphabet))
    raise ValueError("Unrecognized independently audited input: %s" % operation)


def independent_reflection(kind, polynomial, node):
    """Use ambient simple-root dot products, outside BorelMap.reflect."""
    W = WeylGroup(kind, prefix="input_s")
    roots = W.domain().simple_roots()
    x = tuple(polynomial.parent().gens())
    alpha = sum(2 * roots[node].dot_product(roots[j]) /
                roots[j].dot_product(roots[j]) * x[j - 1]
                for j in range(1, len(x) + 1))
    images = list(x)
    images[node - 1] -= alpha
    return polynomial(*images)


def term_dictionary(W, terms):
    result = defaultdict(QQ)
    for term in terms:
        result[W.from_reduced_word(term["word"])] += QQ(term["coefficient"])
    return {w: c for w, c in result.items() if c}


class RootChevalleyOracle:
    """A genuinely different algorithm for selected complete expansions.

    No BorelMap/OperationsMixin methods are used here.  Positive-root
    reflections, 2*alpha/(alpha,alpha), and Bruhat length give divisor products.
    Monomials are then accumulated by exact linearity in the full flag ring.
    """

    def __init__(self, kind):
        self.W = WeylGroup(kind, prefix="chevalley_s")
        self.fundamental = self.W.domain().fundamental_weights()
        reflections = self.W.reflections()
        self.reflections = tuple(
            (reflections[root], 2 * root / root.dot_product(root))
            for root in reflections.keys())
        self.rank = int(self.W.cartan_type().rank())

    @lru_cache(None)
    def covers(self, element):
        length = element.length() + 1
        answer = []
        for reflection, coroot in self.reflections:
            target = element * reflection
            if target.length() == length:
                answer.append((target, tuple(QQ(self.fundamental[i].dot_product(coroot))
                                             for i in range(1, self.rank + 1))))
        return tuple(answer)

    @lru_cache(None)
    def monomial(self, powers):
        if not any(powers):
            return ((self.W.one(), QQ.one()),)
        position = next(i for i, exponent in enumerate(powers) if exponent)
        lower = list(powers)
        lower[position] -= 1
        answer = defaultdict(QQ)
        for w, coefficient in self.monomial(tuple(lower)):
            for target, scalars in self.covers(w):
                if scalars[position]:
                    answer[target] += coefficient * scalars[position]
        return tuple((w, c) for w, c in answer.items() if c)

    def expand(self, polynomial):
        answer = defaultdict(QQ)
        for exponent, coefficient in polynomial.dict().items():
            for w, scalar in self.monomial(tuple(exponent)):
                answer[w] += coefficient * scalar
        return {w: c for w, c in answer.items() if c}


class ExceptionalCatalogueTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = time.monotonic()
        cls.fixture = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
        cls.rows = cls.fixture["rows"]
        cls.maps = {}
        cls.results = {}
        cls.row_reports = []
        cls.chevalley_reports = []
        cls.input_polynomials_verified = 0
        cls.levi_reflections_verified = 0

    @classmethod
    def tearDownClass(cls):
        report = {
            "schema_version": 1,
            "sage_version": SAGE_VERSION,
            "python_version": platform.python_version(),
            "seconds": round(time.monotonic() - cls.started, 3),
            "catalogue_rows": len(cls.rows),
            "frozen_fixture_sha256_normalized_lf": hashlib.sha256(
                FIXTURE_PATH.read_text(encoding="utf-8").encode("utf-8")).hexdigest(),
            "freshly_recomputed_rows": len(cls.row_reports),
            "input_polynomials_verified": cls.input_polynomials_verified,
            "levi_reflections_verified": cls.levi_reflections_verified,
            "coefficient_positions_checked": sum(r["coefficient_positions_checked"]
                                                 for r in cls.row_reports),
            "all_frozen_coefficient_comparisons_pass": all(
                r["coefficients_match"] for r in cls.row_reports) and
                len(cls.row_reports) == len(cls.rows),
            "all_independent_chevalley_comparisons_pass": len(cls.chevalley_reports) == 8 and
                all(r["coefficients_match"] for r in cls.chevalley_reports),
            "evidence": {
                "frozen_comparison": "Earlier separately implemented BGG audit; same mathematical algorithm, independently written implementation.",
                "input_check": "Literal elementary subset sums; ambient dot-product simple roots; independent psi power-sum construction.",
                "second_algorithm": "Root/coroot Chevalley multiplication of polynomial monomials; no library divisor or divided-difference methods.",
                "e8_scope": "Nine displayed dissertation rows; partial catalogue, not a complete E8 invariant-generator list.",
                "timing_scope": "Per-row times measure fresh generator_image expansion; input polynomial construction and independent validation precede that timer.",
            },
            "rows": cls.row_reports,
            "independent_chevalley": cls.chevalley_reports,
        }
        print("\nExceptional verification: %d fresh rows; %d second-algorithm comparisons; %.3fs" %
              (len(cls.row_reports), len(cls.chevalley_reports), report["seconds"]), flush=True)
        if os.environ.get("BOREL_EXCEPTIONAL_REPORT"):
            Path(os.environ["BOREL_EXCEPTIONAL_REPORT"]).write_text(
                json.dumps(report, indent=2) + "\n", encoding="utf-8")

    def borel(self, row):
        key = (row["kind"], row["parabolic"])
        if key not in self.maps:
            self.maps[key] = exceptional_borel(*key)
        return self.maps[key]

    def result(self, row):
        source_id = row["source_id"]
        if source_id not in self.results:
            b = self.borel(row)
            start = time.monotonic()
            answer = b.generator_image(row["name"])
            expected = term_dictionary(b.W, row["terms"])
            self.results[source_id] = answer
            self.row_reports.append({
                "source_id": source_id, "kind": row["kind"],
                "parabolic": row["parabolic"], "name": row["name"],
                "coefficients_match": dict(answer) == expected,
                "output_terms": len(answer),
                "coefficient_positions_checked": len(set(answer) | set(expected)),
                "seconds": round(time.monotonic() - start, 3),
                "algorithm_stats": dict(b.last_stats),
            })
            print("  %s (%s): %d terms, %.3fs" %
                  (source_id, row["name"], len(answer), self.row_reports[-1]["seconds"]), flush=True)
        return self.results[source_id]

    def test_10_complete_catalogue_recomputation(self):
        self.assertEqual(len(self.rows), 128)
        self.assertEqual(len({r["source_id"] for r in self.rows}), 128)
        self.assertEqual(sum(r["kind"] == "E8" for r in self.rows), 9)
        for row in self.rows:
            with self.subTest(source_id=row["source_id"]):
                b = self.borel(row)
                info = b.generator_info(row["name"])
                polynomial = b.generator_polynomial(row["name"])
                expected_polynomial = independent_polynomial(row, b.ring)
                self.assertEqual(polynomial, expected_polynomial)
                type(self).input_polynomials_verified += 1
                self.assertEqual(polynomial.total_degree(), info.polynomial_degree)
                self.assertEqual(info.cohomological_degree, 2 * info.polynomial_degree)
                for node in range(1, b.rank + 1):
                    if node != row["parabolic"]:
                        self.assertEqual(independent_reflection(row["kind"], polynomial, node), polynomial)
                        type(self).levi_reflections_verified += 1
                answer = self.result(row)
                self.assertEqual(dict(answer), term_dictionary(b.W, row["terms"]))
                for w, coefficient in answer.items():
                    self.assertEqual(w.length(), info.polynomial_degree)
                    self.assertIn(coefficient, QQ)
                    self.assertTrue(all(not w.has_descent(node, side="right")
                                        for node in b.retained_nodes))

    def test_20_named_corrections_and_complete_psi12(self):
        rows = {row["source_id"]: row for row in self.rows}
        self.assertTrue(CORRECTED_IDS <= set(rows))
        for source_id in sorted(CORRECTED_IDS):
            row = rows[source_id]
            with self.subTest(source_id=source_id):
                self.assertTrue(self.borel(row).generator_info(row["name"]).corrections)
                self.assertEqual(dict(self.result(row)), term_dictionary(self.borel(row).W, row["terms"]))
        psi12 = [row for row in self.rows if row["kind"] == "E7" and
                 row["input"]["operation"] == "psi" and row["input"]["degree"] == 12]
        self.assertEqual(len(psi12), 1)
        row, = psi12
        self.assertEqual(len(self.result(row)), len(row["terms"]))
        self.assertEqual(len(row["terms"]), 3)
        euler = [row for row in self.rows if row["source_id"] == "E6/P6/omitted_euler_product"]
        self.assertEqual(len(euler), 1)
        self.assertEqual(euler[0]["name"], "euler")

    def test_30_independent_root_coroot_chevalley(self):
        selected = [row for row in self.rows if row["source_id"] in CORRECTED_IDS or
                    row["source_id"] == "E6/P6/omitted_euler_product"]
        for kind in ("G2", "E7", "E8"):
            selected.append(next(row for row in self.rows if row["kind"] == kind and
                                 self.borel(row).generator_info(row["name"]).polynomial_degree == 2))
        self.assertEqual(len(selected), 8)
        oracles = {}
        for row in selected:
            with self.subTest(source_id=row["source_id"]):
                start = time.monotonic()
                oracle = oracles.setdefault(row["kind"], None)
                if oracle is None:
                    oracle = oracles[row["kind"]] = RootChevalleyOracle(row["kind"])
                b = self.borel(row)
                polynomial = independent_polynomial(row, b.ring)
                expected = oracle.expand(polynomial)
                actual = term_dictionary(oracle.W, self.result(row).to_terms())
                matches = expected == actual
                self.chevalley_reports.append({
                    "source_id": row["source_id"], "coefficients_match": matches,
                    "output_terms": len(expected),
                    "seconds": round(time.monotonic() - start, 3),
                })
                print("  Chevalley %s: %d terms, %.3fs" %
                      (row["source_id"], len(expected), self.chevalley_reports[-1]["seconds"]), flush=True)
                self.assertEqual(actual, expected)

    def test_40_aliases_wrappers_and_immutable_metadata(self):
        b = exceptional_borel("F4", 1)
        descriptors = b.generators()
        self.assertIsInstance(descriptors, tuple)
        self.assertEqual(tuple(info.name for info in exceptional_generators("F4", 1)),
                         tuple(info.name for info in descriptors))
        info = b.generator_info("e2")
        expected = b.generator_polynomial("e2")
        self.assertEqual(exceptional_polynomial("F4", 1, "e2"), expected)
        self.assertEqual(exceptional_image("F4", parabolic=1, generator="e2"), b.generator_image("e2"))
        for alias in info.aliases:
            self.assertEqual(b.generator_info(alias).name, "e2")
            self.assertEqual(b.generator_polynomial(alias), expected)
        with self.assertRaises((AttributeError, TypeError)):
            info.name = "changed"
        with self.assertRaises(TypeError):
            info.spec["degree"] = 0
        with self.assertRaises(TypeError):
            info.source["unexpected"] = True
        for row in self.rows:
            info = self.borel(row).generator_info(row["name"])
            self.assertTrue(info.source)
            self.assertTrue(info.normalization)

    def test_50_validation_and_partial_e8_coverage(self):
        for kind in ("A3", "B3", "C3", "D4"):
            with self.subTest(kind=kind):
                with self.assertRaises(ValueError):
                    exceptional_borel(kind, 1)
        for node in (0, 5, True, QQ(3) / 2):
            with self.subTest(node=node):
                with self.assertRaises((ValueError, TypeError)):
                    exceptional_borel("F4", node)
        b = exceptional_borel("F4", 2)
        with self.assertRaises((KeyError, ValueError)):
            b.generator_info("not-a-generator")
        with self.assertRaises((KeyError, ValueError)):
            b.generator_info("e2")  # two distinct elementary alphabets
        for p in (1, 2):
            context = exceptional_borel("E8", p)
            self.assertTrue(context.generators())
            self.assertFalse(context.catalogue_coverage["complete_generator_lists"])
        # Unsupported named catalogues must never silently imply completeness.
        for p in (3, 4, 5, 6, 7, 8):
            with self.subTest(parabolic=p):
                with self.assertRaises((KeyError, ValueError, NotImplementedError)):
                    exceptional_borel("E8", p)


if __name__ == "__main__":
    unittest.main()
