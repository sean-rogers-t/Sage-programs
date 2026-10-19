"""Independent partition, recorded lattice, and fresh BGG checks of the API."""

from collections import Counter
from dataclasses import FrozenInstanceError
from hashlib import sha256
import json
import os
from pathlib import Path
from time import perf_counter
import unittest
from unittest.mock import patch

from sage.all import Partitions, QQ, ZZ

from sage_borel import GeneratorAlgebra, classical_borel, orthogonal_smith


COUNTS = Counter()
ROOT = Path(__file__).resolve().parents[1]


class OrthogonalFormulaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.started = perf_counter()
        cls.methods_run = 0

    def setUp(self):
        type(self).methods_run += 1
        type(self).test_result = self._outcome.result

    @classmethod
    def tearDownClass(cls):
        failed = any(isinstance(test, cls)
                     for test, _ in cls.test_result.failures + cls.test_result.errors)
        files = ("src/sage_borel/orthogonal.py", "src/sage_borel/operations.py",
                 "tests/test_orthogonal.py",
                 "experiments/orthogonal_smith_search_results.json",
                 "experiments/type_b_multiplicity_search_results.json",
                 "experiments/type_d_even_gluing_verification.json")
        report = {"suite_result": "failed" if failed else "passed",
                  "test_methods_run": cls.methods_run,
                  "command": "PYTHONPATH=src sage -python -m unittest tests.test_orthogonal -v",
                  "scope": "finite implementation checks, not an arbitrary-rank proof",
                  "source_sha256": {name: sha256((ROOT / name).read_bytes()).hexdigest()
                                    for name in files},
                  "independent_checks": dict(COUNTS),
                  "seconds": round(perf_counter() - cls.started, 6)}
        print("\nOrthogonal formula checks: " + json.dumps(report, sort_keys=True), flush=True)
        if os.environ.get("BOREL_ORTHOGONAL_REPORT"):
            Path(os.environ["BOREL_ORTHOGONAL_REPORT"]).write_text(
                json.dumps(report, indent=2) + "\n")

    def test_independent_partition_statistic_every_degree(self):
        # The implementation uses Gaussian products.  This oracle uses
        # the different beta-number rule on individual rectangle partitions.
        for a in range(1, 7):
            for k in range(1, 7):
                n, length, width = a + k, min(a, k), max(a, k)
                rank = 0
                for m in range(a * k + 1):
                    expected = Counter()
                    for part in Partitions(m, max_part=width, max_length=length):
                        padded = list(part) + [0] * (length - len(part))
                        beta = [int(value) + length - i - 1
                                for i, value in enumerate(padded)]
                        low = sum(value < n // 2 for value in beta)
                        exponent = sum(value >= (n + 1) // 2 + low for value in beta)
                        expected[exponent] += 1
                        COUNTS["independent_partitions"] += 1
                    actual = orthogonal_smith("B", n, k, 2 * m)
                    self.assertEqual(dict(actual.exponent_multiplicities), expected,
                                     (n, k, m))
                    rank += actual.rank
                    COUNTS["partition_degree_profiles"] += 1
                self.assertEqual(rank, ZZ(n).binomial(k))

    def test_recorded_type_b_integer_matrices(self):
        for filename in ("orthogonal_smith_search_results.json",
                         "type_b_multiplicity_search_results.json"):
            data = json.loads((ROOT / "experiments" / filename).read_text())
            for case in data["cases"]:
                if case.get("family", "B") != "B":
                    continue
                for row in case["degrees"]:
                    if row["status"] != "computed":
                        continue
                    actual = orthogonal_smith("B", case["n"], case["k"],
                                              row["polynomial_degree"])
                    self.assertEqual(actual.smith_factors(),
                                     tuple(map(int, row["smith_factors"])))
                    COUNTS["recorded_B_smith_lists"] += 1

    def test_recorded_type_d_full_and_sector_lattices(self):
        data = json.loads((ROOT / "experiments" /
                           "type_d_even_gluing_verification.json").read_text())
        rows = data["three_head_components"] + data["higher_Cartesian_components"]
        for row in rows:
            for sector, field in (("full", "full_exponents"),
                                  ("square", "square_exponents"),
                                  ("euler", "Euler_exponents")):
                result = orthogonal_smith("D", row["n"], row["k"],
                                          row["polynomial_degree"], sector=sector)
                self.assertEqual(dict(result.exponent_multiplicities),
                                 Counter(row[field]), (row["n"], row["k"], sector))
                if sector == "full":
                    self.assertEqual(result.unit_gluing_dimension,
                                     row["FULL_saturated_gluing_dimension"])
                else:
                    self.assertEqual(result.unit_gluing_dimension, 0)
                COUNTS["recorded_D_sector_or_full_lists"] += 1

    def test_fresh_bgg_matrices_and_character_normalization(self):
        for family, n, k, degrees in (
                ("B", 4, 2, (0, 2, 4, 6, 8)),
                ("D", 4, 2, (0, 2, 4, 6, 8)),
                ("D", 4, 3, (0, 2, 3, 4, 5, 6))):
            borel = classical_borel(family, n, k)
            stop = k if family == "B" else k - 1
            generators = {"q%d" % i: borel.tail_polynomial(k, i)
                          for i in range(1, stop + 1)}
            if family == "D":
                generators["p"] = borel.euler_polynomial(k)
            algebra = GeneratorAlgebra(borel, generators)
            for degree in degrees:
                result = orthogonal_smith(family, n, k, degree)
                matrix = algebra.integral_degree(degree)
                self.assertEqual(result.smith_factors(), matrix.smith_factors)
                self.assertEqual(result.saturation_index(), matrix.saturation_index)
                COUNTS["fresh_BGG_integer_image_checks"] += 1
            scaled = {name: (2 ** polynomial.degree()) * polynomial
                      for name, polynomial in generators.items()}
            character = GeneratorAlgebra(borel, scaled).integral_degree(4)
            result = orthogonal_smith(family, n, k, 4, normalization="character")
            self.assertEqual(result.smith_factors(), character.smith_factors)
            COUNTS["fresh_scaled_character_matrix_checks"] += 1

    def test_summary_and_gluing_example(self):
        result = orthogonal_smith("d", ZZ(4), ZZ(2), ZZ(4))
        self.assertEqual(result.multiplicities, {1: 1, 2: 2})
        self.assertEqual(result.exponent_multiplicities, ((0, 1), (1, 2)))
        self.assertEqual(result.rank, 3)
        self.assertEqual(result.cohomological_degree, 8)
        self.assertEqual(result.saturation_exponent, 2)
        self.assertEqual(result.saturation_index(), 4)
        self.assertFalse(result.is_saturated)
        self.assertEqual(result.unit_gluing_dimension, 1)
        self.assertEqual(result.smith_factors(), (1, 2, 2))
        with self.assertRaises(FrozenInstanceError):
            result.degree = 5
        result.multiplicities.clear()
        self.assertEqual(result.rank, 3)
        encoded = result.as_dict()
        self.assertEqual(encoded["cohomological_degree"], 8)
        self.assertEqual(encoded["normalized_unit_gluing_dimension"], 1)
        self.assertNotIn("smith_factors", encoded)
        json.dumps(encoded)
        character = orthogonal_smith("D", 4, 2, 4, normalization="character")
        self.assertEqual(character.multiplicities, {16: 1, 32: 2})
        self.assertEqual(character.unit_gluing_dimension, 1)

    def test_odd_tail_and_one_head_euler_boundary(self):
        for n, k in ((4, 3), (6, 3), (6, 5)):
            for degree in range(2 * k * (n - k) + 2):
                full = orthogonal_smith("D", n, k, degree)
                sector = "square" if degree % 2 == 0 else "euler"
                self.assertEqual(full.exponent_multiplicities,
                                 orthogonal_smith("D", n, k, degree,
                                                   sector=sector).exponent_multiplicities)
                self.assertEqual(full.unit_gluing_dimension, 0)
        self.assertEqual(orthogonal_smith("D", 5, 4, 4, sector="euler").smith_factors(), (1,))
        self.assertEqual(orthogonal_smith("D", 5, 4, 4).smith_factors(), (1, 2))
        self.assertEqual(orthogonal_smith("D", 5, 4, 6, sector="euler").rank, 0)

    def test_empty_and_constant_images(self):
        for family, n, k in (("B", 2, 1), ("B", 8, 3), ("D", 4, 2)):
            for degree in (2 * k * (n - k) + 1, 10 ** 8):
                result = orthogonal_smith(family, n, k, degree)
                self.assertEqual(result.rank, 0)
                self.assertEqual(result.smith_factors(), ())
                self.assertEqual(result.saturation_index(), 1)
                self.assertTrue(result.is_saturated)
            result = orthogonal_smith(family, n, k, 0, normalization="character")
            self.assertEqual(result.smith_factors(), (1,))
        self.assertEqual(orthogonal_smith("B", 8, 3, 3).rank, 0)

    def test_large_rank_without_weyl_groups_and_output_guards(self):
        with patch("sage_borel.core.WeylGroup", side_effect=AssertionError("Weyl construction")):
            result = orthogonal_smith("D", 200, 80, 1000)
        self.assertGreater(result.rank, 100000)
        self.assertGreater(result.saturation_exponent, 1000000)
        with self.assertRaises(ValueError):
            result.smith_factors()
        with self.assertRaises(ValueError):
            result.saturation_index()
        json.dumps(result.as_dict())
        small = orthogonal_smith("D", 4, 2, 4)
        with self.assertRaises(ValueError):
            small.smith_factors(max_rank=2)
        self.assertEqual(small.smith_factors(max_rank=3), (1, 2, 2))
        with self.assertRaises(ValueError):
            small.saturation_index(max_bits=2)
        self.assertEqual(small.saturation_index(max_bits=3), 4)

    def test_invalid_inputs_are_not_silently_coerced(self):
        for family in ("A", "C", "F4", "B4", None, 1):
            with self.assertRaises(ValueError):
                orthogonal_smith(family, 4, 2, 4)
        for position in (1, 2, 3):
            for bad in (True, 2.5, QQ(3) / 2, "4"):
                args = ["D", 4, 2, 4]
                args[position] = bad
                with self.assertRaises(TypeError):
                    orthogonal_smith(*args)
        for args in (("B", 1, 1, 0), ("D", 3, 2, 0), ("B", 4, 0, 0),
                     ("D", 4, 1, 0), ("B", 4, 4, 0), ("D", 4, 2, -1)):
            with self.assertRaises(ValueError):
                orthogonal_smith(*args)
        for kwargs in ({"sector": "mixed"}, {"normalization": "Spin"}):
            with self.assertRaises(ValueError):
                orthogonal_smith("D", 4, 2, 4, **kwargs)
        with self.assertRaises(ValueError):
            orthogonal_smith("B", 4, 2, 4, sector="euler")
        result = orthogonal_smith("D", 4, 2, 4)
        with self.assertRaises(TypeError):
            result.smith_factors(max_rank=True)
        with self.assertRaises(ValueError):
            result.saturation_index(max_bits=0)


if __name__ == "__main__":
    unittest.main()
