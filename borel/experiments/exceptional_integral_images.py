"""Complete G2/F4 maximal-parabolic image tables for two explicit Z-forms.

Run: PYTHONPATH=src sage -python experiments/exceptional_integral_images.py
The second form makes individual nonzero generator images primitive; it does
not saturate their products automatically or identify a character lattice.
"""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

from sage.all import QQ, gcd, lcm
from sage.env import SAGE_VERSION
from sage_borel import GeneratorAlgebra, exceptional_borel

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.test_exceptional import RootChevalleyOracle


def reduce_generators(algebra):
    """Discard a generator only with a verified integral polynomial preimage."""
    chosen, relations = {}, []
    entries = sorted(zip(algebra.generator_names, algebra.generator_polynomials),
                     key=lambda item: item[1].total_degree())
    for name, polynomial in entries:
        image = algebra.borel_map.expand(polynomial)
        if chosen:
            smaller = GeneratorAlgebra(algebra.borel_map, chosen)
            degree = smaller.integral_degree(int(polynomial.total_degree()))
            if degree.contains(image):
                preimage = degree.preimage(image)
                assert smaller.evaluate(preimage) == image
                relations.append({"name": name, "integer_image_expression": str(preimage)})
                continue
        chosen[name] = polynomial
    # In the two middle F4 cases neither quadratic can simply be discarded.
    # Their sum, however, is an integral replacement for both. Certify this
    # rather than substituting a rationally equivalent set of generators.
    added = []
    quadratics = [name for name, p in chosen.items() if p.total_degree() == 2]
    if len(quadratics) == 2:
        first, second = quadratics
        candidate = {name: p for name, p in chosen.items() if name not in quadratics}
        candidate["quadratic_sum"] = chosen[first] + chosen[second]
        smaller = GeneratorAlgebra(algebra.borel_map, candidate)
        degree = smaller.integral_degree(2)
        images = [algebra.borel_map.expand(chosen[name]) for name in quadratics]
        if all(degree.contains(image) for image in images):
            for name, image in zip(quadratics, images):
                preimage = degree.preimage(image)
                assert smaller.evaluate(preimage) == image
                relations.append({"name": name, "integer_image_expression": str(preimage)})
            added.append({"name": "quadratic_sum", "expression_in_scaled_catalogue": first + "+" + second})
            chosen = candidate
    return GeneratorAlgebra(algebra.borel_map, chosen), relations, added


def compare_column(oracle, algebra, report, column):
    powers = report.monomials[column].exponents()[0]
    polynomial = algebra.borel_map.ring.one()
    for generator, exponent in zip(algebra.generator_polynomials, powers):
        polynomial *= generator ** exponent
    expected = oracle.expand(polynomial)
    expected = {algebra.borel_map.W.from_reduced_word(w.reduced_word()): c
                for w, c in expected.items() if c}
    actual = {w: report.matrix[i, column]
              for i, w in enumerate(report.schubert_indices) if report.matrix[i, column]}
    assert actual == expected


def main():
    started = perf_counter()
    output = ROOT / "experiments/exceptional-integral-results.json"
    result = {"sage_version": SAGE_VERSION,
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "All degrees of all G2/F4 maximal parabolics, for explicitly selected generator forms; not a representation-ring intersection or novelty claim.",
              "cases": []}
    checked = 0
    for kind, nodes in (("G2", range(1, 3)), ("F4", range(1, 5))):
        oracle = RootChevalleyOracle(kind)
        for node in nodes:
            borel = exceptional_borel(kind, node)
            audit, cleared, primitive = [], {}, {}
            for info in borel.generators():
                image = borel.generator_image(info.name)
                denominator = lcm([c.denominator() for c in image.values()])
                numerator = gcd([c * denominator for c in image.values()])
                cleared[info.name] = denominator
                primitive[info.name] = denominator / numerator if numerator else QQ.one()
                # Independent algorithm for every original catalogue generator.
                expected = oracle.expand(borel.generator_polynomial(info.name))
                assert {borel.W.from_reduced_word(w.reduced_word()): c
                        for w, c in expected.items()} == dict(image)
                checked += 1
                audit.append({"name": info.name, "degree": info.polynomial_degree,
                              "image": image.to_terms(),
                              "least_integer_image_multiplier": str(denominator),
                              "primitive_image_multiplier": str(primitive[info.name]),
                              "zero_image": not bool(image)})
            top = max(d for d in range(int(borel.W.long_element().length()) + 1)
                      if borel.basis(d))
            for form, scales in (("image_denominators_cleared", cleared),
                                 ("individually_primitive_images", primitive)):
                full = borel.generator_algebra(scales=scales)
                algebra, eliminated, added = reduce_generators(full)
                case = {"kind": kind, "parabolic": node, "form": form,
                        "complex_dimension": top, "generators": audit,
                        "multipliers": {name: str(value) for name, value in scales.items()},
                        "retained_generators": list(algebra.generator_names),
                        "retained_polynomials": [str(p) for p in algebra.generator_polynomials],
                        "added_image_generators": added,
                        "eliminated_image_generators": eliminated, "degrees": []}
                for d in range(top + 1):
                    report = algebra.integral_degree(d, include_ambient=True)
                    assert report.rank == report.ambient_dimension
                    # Complete low-degree matrices and all G2 columns, not
                    # merely equality of Smith factors between two algorithms.
                    if kind == "G2" or d <= 6:
                        for column in range(len(report.monomials)):
                            compare_column(oracle, algebra, report, column)
                            checked += 1
                    case["degrees"].append({
                        "degree": d, "rank": report.rank,
                        "ambient_dimension": report.ambient_dimension,
                        "smith_factors": [int(v) for v in report.smith_factors],
                        "saturation_index": int(report.saturation_index),
                        "matrix": [[int(v) for v in row] for row in report.matrix.rows()],
                        "monomials": [str(m) for m in report.monomials],
                        "schubert_words": [list(w.reduced_word()) for w in report.schubert_indices]})
                    print(kind, node, form, d, report.smith_factors, flush=True)
                assert not borel.basis(top + 1)
                result["cases"].append(case)
                output.write_text(json.dumps(result, indent=2) + "\n")
        oracle.monomial.cache_clear()
        oracle.covers.cache_clear()
    result["independent_chevalley_columns"] = checked
    result["degree_rows"] = sum(len(c["degrees"]) for c in result["cases"])
    result["seconds"] = round(perf_counter() - started, 3)
    files = [Path(__file__).resolve(), ROOT / "tests/test_exceptional.py"]
    files += sorted((ROOT / "src/sage_borel").glob("*.py"))
    files += [ROOT / "src/sage_borel/data/exceptional_generators.json"]
    result["source_hash_semantics"] = "SHA256 of UTF-8 file bytes with CRLF normalized to LF, matching Git text blobs."
    result["source_sha256"] = {str(p.relative_to(ROOT)): sha256(p.read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                               for p in files}
    result["result"] = "passed"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASSED", result["degree_rows"], "degree rows;", checked,
          "independent Chevalley columns;", result["seconds"], "seconds", flush=True)


if __name__ == "__main__":
    main()
