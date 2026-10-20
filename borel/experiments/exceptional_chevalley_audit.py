"""Reconstruct every saved G2/F4 matrix by root/coroot Chevalley products.

No library divided differences, expansion, multiplication or Smith matrices
are used to obtain expected columns. Catalogue polynomials are the shared
inputs, while their independent reconstruction is covered by the catalogue
test suite. Run after exceptional_integral_images.py.
"""
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

from sage.all import QQ, PolynomialRing
from sage.env import SAGE_VERSION
from sage_borel import exceptional_borel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.test_exceptional import RootChevalleyOracle


class ChevalleyAction(RootChevalleyOracle):
    @lru_cache(None)
    def monomial_on(self, powers, start):
        if not any(powers):
            return ((start, QQ.one()),)
        position = next(i for i, exponent in enumerate(powers) if exponent)
        lower = list(powers)
        lower[position] -= 1
        answer = defaultdict(QQ)
        for w, coefficient in self.monomial_on(tuple(lower), start):
            for target, scalars in self.covers(w):
                if scalars[position]:
                    answer[target] += coefficient * scalars[position]
        return tuple((w, c) for w, c in answer.items() if c)

    @lru_cache(None)
    def polynomial_on(self, polynomial, start):
        answer = defaultdict(QQ)
        for powers, coefficient in polynomial.dict().items():
            for w, scalar in self.monomial_on(tuple(powers), start):
                answer[w] += coefficient * scalar
        return tuple((w, c) for w, c in answer.items() if c)

    def multiply_polynomial(self, polynomial, expansion):
        answer = defaultdict(QQ)
        for start, coefficient in expansion.items():
            for w, scalar in self.polynomial_on(polynomial, start):
                answer[w] += coefficient * scalar
        return {w: c for w, c in answer.items() if c}


def main():
    started = perf_counter()
    source = ROOT / "experiments/exceptional-integral-results.json"
    tables = json.loads(source.read_text())
    assert tables["result"] == "passed"
    columns, rows = 0, 0
    for kind in ("G2", "F4"):
        oracle = ChevalleyAction(kind)
        for case in tables["cases"]:
            if case["kind"] != kind:
                continue
            b = exceptional_borel(kind, case["parabolic"])
            algebra = b.generator_algebra(scales={n: QQ(v) for n, v in case["multipliers"].items()})
            original = dict(zip(algebra.generator_names, algebra.generator_polynomials))
            inputs = dict(original)
            for addition in case["added_image_generators"]:
                formal = algebra.domain(addition["expression_in_scaled_catalogue"])
                inputs[addition["name"]] = b.ring(formal(*algebra.generator_polynomials))
            names = case["retained_generators"]
            assert [str(inputs[n]) for n in names] == case["retained_polynomials"]
            domain = PolynomialRing(QQ, len(names), names=names)
            @lru_cache(None)
            def column(powers):
                if not any(powers):
                    return {oracle.W.one(): QQ.one()}
                i = next(i for i, power in enumerate(powers) if power)
                lower = list(powers)
                lower[i] -= 1
                return oracle.multiply_polynomial(
                    inputs[names[i]], column(tuple(lower)))

            for degree in case["degrees"]:
                labels = [oracle.W.from_reduced_word(word) for word in degree["schubert_words"]]
                for j, monomial in enumerate(degree["monomials"]):
                    powers = tuple(domain(monomial).exponents()[0])
                    expected = column(powers)
                    actual = {w: QQ(degree["matrix"][i][j])
                              for i, w in enumerate(labels) if degree["matrix"][i][j]}
                    assert expected == actual, (kind, case["parabolic"], case["form"], degree["degree"], monomial)
                    columns += 1
                rows += 1
            print(kind, case["parabolic"], case["form"], "all matrices agree", flush=True)
            column.cache_clear()
        oracle.polynomial_on.cache_clear()
        oracle.monomial_on.cache_clear()
        oracle.covers.cache_clear()
    result = {"result": "passed", "sage_version": SAGE_VERSION,
              "degree_matrices": rows, "independent_columns": columns,
              "seconds": round(perf_counter()-started, 3),
              "scope": "Every stored source column independently reconstructed by positive-root Chevalley multiplication, including top degree.",
              "source_hash_semantics": "SHA256 of UTF-8 file bytes with CRLF normalized to LF, matching Git text blobs.",
              "source_sha256": {str(p.relative_to(ROOT)): sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                                for p in [source, Path(__file__).resolve(), ROOT / "tests/test_exceptional.py",
                                          ROOT / "src/sage_borel/exceptional.py",
                                          ROOT / "src/sage_borel/data/exceptional_generators.json"]}}
    (ROOT / "experiments/exceptional-chevalley-audit.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
