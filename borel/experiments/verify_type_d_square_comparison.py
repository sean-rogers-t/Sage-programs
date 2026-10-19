"""Finite exact checks for the type-D square-sector comparison note.

Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_square_comparison.py
The general equality of D and B square Smith lists is not assumed or proved
by this finite experiment.  The false Gysin shortcut is tested independently
using Cartesian type C and exact B/C row rescaling as well as the library BGG.
"""

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import PolynomialRing, QQ, ZZ, matrix, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, classical_borel
from sage_borel.operations import elementary_polynomial
from type_b_cartesian_frontier import CartesianCFrontier


def rectangle_partitions(size, width, height):
    def visit(left, maximum, slots, prefix):
        if not left:
            yield prefix
        elif slots:
            for first in range(min(left, maximum), 0, -1):
                yield from visit(left-first, first, slots-1, prefix+(first,))
    return tuple(visit(size, width, height, ()))


def schur_polynomials(a, k, m):
    ring = PolynomialRing(QQ, k, names=["q%d" % i for i in range(1, k+1)])
    q = ring.gens()
    complete = [ring.one()]
    for degree in range(1, a+k):
        complete.append(ring(sum((-1)**(i-1)*q[i-1]*complete[degree-i]
                                 for i in range(1, min(k, degree)+1))))
    partitions = rectangle_partitions(m, a, k)
    polynomials = []
    for partition in partitions:
        length = len(partition)
        polynomials.append(matrix(ring, length, length,
            lambda i, j: complete[partition[i]-i+j]
            if partition[i]-i+j >= 0 else ring.zero()).det())
    return ring, partitions, polynomials


def coefficient_matrix(borel, expansions):
    words = tuple(sorted({word for expansion in expansions for word in expansion},
                         key=borel._word))
    return matrix(ZZ, len(words), len(expansions),
                  lambda i, j: expansions[j].get(words[i], 0))


def smith_factors(coefficients):
    return tuple(abs(ZZ(value)) for value in coefficients.elementary_divisors() if value)


def exponent(factor):
    power = factor.valuation(2)
    assert factor == 2**power, factor
    return int(power)


def square_checks():
    checks, stored = [], {}
    # Complete square ranges in the first five cases; a separate low-degree
    # sample in rank six.  Rank-four k=1 uses the two omitted fork nodes.
    parameters = [(1, 3, 3), (3, 1, 3), (2, 2, 4),
                  (2, 3, 6), (3, 2, 6), (3, 3, 4)]
    for a, k, maximum in parameters:
        n = a+k
        b = classical_borel("B", n, k)
        d = (BorelMap(["D", n], omitted_nodes=[n-1, n]) if k == 1
             else classical_borel("D", n, k))
        for m in range(maximum+1):
            started = perf_counter()
            ring, partitions, formal = schur_polynomials(a, k, m)
            b_subst = ring.hom([b.tail_polynomial(k, i) for i in range(1, k+1)], b.ring)
            d_subst = ring.hom([d.tail_polynomial(k, i) for i in range(1, k+1)], d.ring)
            b_matrix = coefficient_matrix(b, [b.expand(b_subst(f)) for f in formal])
            d_matrix = coefficient_matrix(d, [d.expand(d_subst(f)) for f in formal])
            b_factors, d_factors = smith_factors(b_matrix), smith_factors(d_matrix)
            b_rows, d_rows = b_matrix.row_module(), d_matrix.row_module()
            assert len(b_factors) == len(d_factors) == len(partitions)
            assert all(row in b_rows for row in d_rows.basis())
            b_exponents, d_exponents = tuple(map(exponent, b_factors)), tuple(map(exponent, d_factors))
            assert all(x <= y for x, y in zip(b_exponents, d_exponents))
            if 2*m <= max(a, k):
                assert all(value == 1 for value in d_factors)
            # This finite observed equality is recorded separately from the
            # assertions that corroborate the proved one-sided comparison.
            agrees = b_factors == d_factors
            stored[a, k, m] = d_factors
            row = {"a": a, "k": k, "n": n,
                   "polynomial_degree": 2*m, "source_partitions": [list(p) for p in partitions],
                   "rank": len(partitions), "B_smith_factors": list(map(str, b_factors)),
                   "D_square_smith_factors": list(map(str, d_factors)),
                   "rectangle_rank_matches": True, "D_rows_contained_in_B_rows": True,
                   "sorted_exponent_lower_bound_holds": True,
                   "within_proved_primitive_range": 2*m <= max(a, k),
                   "finite_square_B_D_equality": agrees,
                   "seconds": round(perf_counter()-started, 6)}
            checks.append(row)
            print("a=%d k=%d d=%d rank=%d B=%s D=%s row inclusion verified" %
                  (a, k, 2*m, len(partitions), b_factors, d_factors), flush=True)
    dual_checks = []
    for (a, k, m), factors in sorted(stored.items()):
        if a < k and (k, a, m) in stored:
            assert factors == stored[k, a, m]
            dual_checks.append({"a": a, "k": k, "polynomial_degree": 2*m,
                                "square_block_duality_matches": True})
    return checks, dual_checks


def gysin_counterexample():
    a = k = 2
    n, m = a+k, 2
    b, c = classical_borel("B", n, k), classical_borel("C", n, k)
    formal_ring, partitions, formal = schur_polynomials(a, k, m)
    b_subst = formal_ring.hom([b.tail_polynomial(k, i) for i in range(1, k+1)], b.ring)
    b_head = prod(b.classical_coordinates()[:a])
    direct = [[b.expand(b_subst(f)) for f in formal],
              [b.expand(b_head*b_subst(f)) for f in formal]]
    cartesian_ring = PolynomialRing(QQ, n, names=["x%d" % i for i in range(1, n+1)])
    x = cartesian_ring.gens()
    cartesian_subst = formal_ring.hom([
        elementary_polynomial(cartesian_ring, [z*z for z in x[-k:]], i)
        for i in range(1, k+1)], cartesian_ring)
    frontier = CartesianCFrontier(c)
    independent = []
    for head_factor in [cartesian_ring.one(), prod(x[:a])]:
        columns = []
        for polynomial in formal:
            expansion = frontier.expand(head_factor*cartesian_subst(polynomial))
            scaled = {b.W.from_reduced_word(word.reduced_word()):
                      2**sum(entry < 0 for entry in word.matrix().list())*value
                      for word, value in expansion.items()}
            columns.append(b._expansion(scaled))
        independent.append(columns)
    for actual_columns, expected_columns in zip(direct, independent):
        assert all(actual == expected for actual, expected in zip(actual_columns, expected_columns))
    factors = [smith_factors(coefficient_matrix(b, columns)) for columns in direct]
    assert factors == [(ZZ(1), ZZ(2)), (ZZ(2), ZZ(2))], factors
    print("False Gysin shortcut reproduced independently: %s -> %s" % tuple(factors), flush=True)
    return {"a": a, "k": k, "input_polynomial_degree": 2*m,
            "multiplied_polynomial_degree": 2*m+a,
            "source_partitions": [list(p) for p in partitions],
            "square_smith_factors": list(map(str, factors[0])),
            "head_top_Chern_times_square_smith_factors": list(map(str, factors[1])),
            "Cartesian_C_rescaled_B_matches_fundamental_weight_B_in_both_matrices": True,
            "primitive_multiplication_assertion_is_false": True}


def main():
    started = perf_counter()
    checks, dual_checks = square_checks()
    counterexample = gysin_counterexample()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, here.with_name("type_b_cartesian_frontier.py"),
              borel/"TYPE_D_SQUARE_COMPARISON.md",
              borel/"src"/"sage_borel"/"core.py",
              borel/"src"/"sage_borel"/"operations.py"]
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION,
              "scope": "Finite exact corroboration of proved partial comparisons; general square equality remains unresolved.",
              "normalization": "unscaled normalized squared-tail inputs; polynomial degree is half cohomological degree",
              "source_sha256": {str(path.relative_to(borel)): sha256(path.read_bytes()).hexdigest()
                                 for path in inputs},
              "degree_checks": checks, "block_duality_checks": dual_checks,
              "failed_Gysin_shortcut": counterexample,
              "summary": {"degree_components_checked": len(checks),
                          "block_duality_components_checked": len(dual_checks),
                          "primitive_range_components_checked": sum(row["within_proved_primitive_range"] for row in checks),
                          "square_B_D_mismatches_observed": sum(not row["finite_square_B_D_equality"] for row in checks),
                          "all_proved_partial_comparison_checks_passed": True,
                          "false_primitive_multiplication_counterexample_confirmed": True},
              "seconds": round(perf_counter()-started, 6)}
    output = here.with_name("type_d_square_comparison_verification.json")
    output.write_text(json.dumps(report, indent=2)+"\n")
    print("Saved %s" % output, flush=True)
    print(report["summary"], flush=True)


if __name__ == "__main__":
    main()
