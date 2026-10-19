"""Check the derivation/core-one upper bound for type-D Smith units.

Run from borel/: sage -python experiments/verify_type_d_unit_upper_bound.py

In R'=H*(Gr(2h,a+2h-1);GF(2)), independently compute the kernel of
multiplication by e_1 on the Frobenius image.  Compare its graded ranks
with q^(h*floor(a/2))*[h+ceil(a/2)-1 choose h]_q.  Also compare e_1 F
with its core-one two-quotient LR matrix.  Small exact ZZ matrices check
the signs, integral row modules, and primitivity.  No Borel engine or
Weyl-group calculation is imported.  These checks concern the proved
unit-count argument, not the still-open full type-D Smith formula.
"""

import argparse
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import GF, Partition, Partitions, PolynomialRing, SymmetricFunctions, ZZ, matrix, q_binomial
from sage.env import SAGE_VERSION


F2 = GF(2)
SF2 = SymmetricFunctions(F2).schur()
SZ = SymmetricFunctions(ZZ).schur()
Q = PolynomialRing(ZZ, "q")
q = Q.gen()


@lru_cache(None)
def square(mu):
    return SF2(Partition(mu)) ** 2


@lru_cache(None)
def e1_square(mu):
    return SF2([1]) * square(mu)


@lru_cache(None)
def e1_psi_integer(mu):
    return SZ([1]) * SZ(Partition(mu)).frobenius(2)


@lru_cache(None)
def box(width, height, degree):
    if width < 0 or height < 0 or degree < 0:
        return ()
    return tuple(tuple(p) for p in Partitions(degree, max_part=width, max_length=height))


@lru_cache(None)
def product(alpha, beta):
    return SZ(Partition(alpha)) * SZ(Partition(beta))


def unbalanced_matrix(a, h, degree, source):
    u, v = (a + 1) // 2, a // 2
    factors = [(alpha, beta) for first_degree in range(degree + 1)
               for alpha in box(v - 1, h + 1, first_degree)
               for beta in box(u, h - 1, degree - first_degree)]
    coefficients = matrix(ZZ, len(factors), len(source),
                          lambda i, j: product(*factors[i])[source[j]])
    return coefficients, factors


def kernel_check(a, h, degree, exact_integer=False):
    width, k = a - 1, 2 * h
    source = tuple(Partitions(degree, max_part=width, max_length=k))
    target_even = tuple(Partitions(2 * degree, max_part=width, max_length=k))
    target_odd = tuple(Partitions(2 * degree + 1, max_part=width, max_length=k))
    frobenius = matrix(F2, len(target_even), len(source),
                       lambda i, j: square(tuple(source[j]))[target_even[i]])
    multiplied = matrix(F2, len(target_odd), len(source),
                        lambda i, j: e1_square(tuple(source[j]))[target_odd[i]])
    lr, factorizations = unbalanced_matrix(a, h, degree, source)
    rf, re, rl = frobenius.rank(), multiplied.rank(), lr.change_ring(F2).rank()
    u, v = (a + 1) // 2, a // 2
    candidate = Q(q ** (h * v) * q_binomial(h + u - 1, h, q))
    expected = int(candidate[degree])
    assert re == rl, (a, h, degree, re, rl)
    assert rf - re == expected, (a, h, degree, rf, re, expected)
    row = {"a": a, "k": k, "source_half_degree": degree,
           "source_rank": len(source), "frobenius_rank": int(rf),
           "e1_frobenius_rank": int(re), "core_one_LR_rank": int(rl),
           "kernel_rank": int(rf - re), "gaussian_kernel_rank": expected,
           "core_one_rectangles": [[v - 1, h + 1], [u, h - 1]],
           "kernel_and_core_one_ranks_match": True}
    if exact_integer:
        actual = matrix(ZZ, len(target_odd), len(source),
                        lambda i, j: e1_psi_integer(tuple(source[j]))[target_odd[i]])
        # A core-one partition has h+1 even beads and h-1 odd beads.
        # Its ordinary SXP sign is a single row sign, which does not
        # affect the exact integer row module.
        expected_rows, actual_rows = {}, {}
        for i, pair in enumerate(factorizations):
            expected_rows[pair] = tuple(lr.row(i))
        for i, lam in enumerate(target_odd):
            values = tuple(actual.row(i))
            if lam.core(2) != Partition([1]):
                assert not any(values), (a, h, degree, lam, values)
                continue
            beads = [int(lam[j]) + k - 1 - j if j < len(lam) else k - 1 - j for j in range(k)]
            odd = sorted(((b - 1) // 2 for b in beads if b % 2), reverse=True)
            even = sorted((b // 2 for b in beads if b % 2 == 0), reverse=True)
            assert len(even) == h + 1 and len(odd) == h - 1, (a, h, degree, lam, beads)
            alpha = tuple(value - (len(even) - 1 - j) for j, value in enumerate(even))
            beta = tuple(value - (len(odd) - 1 - j) for j, value in enumerate(odd))
            pair = (tuple(x for x in alpha if x), tuple(x for x in beta if x))
            assert pair in expected_rows, (a, h, degree, lam, pair)
            expected_row = expected_rows[pair]
            assert values == expected_row or values == tuple(-x for x in expected_row), (a, h, degree, lam, values, expected_row)
            actual_rows[pair] = values
        assert set(actual_rows) == set(expected_rows), (a, h, degree)
        equal = actual.row_module() == lr.row_module()
        assert equal, (a, h, degree)
        factors = (() if not all(lr.dimensions()) else
                   tuple(abs(ZZ(d)) for d in lr.smith_form()[0].diagonal() if d))
        assert all(d == 1 for d in factors), (a, h, degree, factors)
        row.update(exact_integer_row_modules_equal=True,
                   core_one_sign_rows_checked=len(actual_rows),
                   integer_LR_smith_factors=list(map(int, factors)),
                   integer_LR_row_lattice_primitive=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-a", type=int, default=7)
    parser.add_argument("--max-h", type=int, default=3)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("type_d_unit_upper_bound_verification.json"))
    args = parser.parse_args()
    started = perf_counter()
    proof_path = Path(__file__).resolve().parent.parent / "TYPE_D_FROBENIUS.md"
    proof_hash = sha256(proof_path.read_bytes()).hexdigest()
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "proof_document": {"file": proof_path.name, "sha256": proof_hash},
              "Borel_or_Weyl_engine_imported": False,
              "scope": "Finite exact verification of the type-D unit upper-bound/core-one argument; no full type-D Smith classification follows from these tests.",
              "arguments": vars(args) | {"output": str(args.output)}, "checks": []}
    for a in range(1, args.max_a + 1):
        for h in range(1, args.max_h + 1):
            for degree in range((a - 1) * h + 1):
                exact = h <= 2 and a <= 5
                report["checks"].append(kernel_check(a, h, degree, exact_integer=exact))
            print("a={} k={} kernel/core-one checks passed".format(a, 2 * h), flush=True)
            report["seconds"] = round(perf_counter() - started, 6)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    report["summary"] = {
        "kernel_degree_components": len(report["checks"]),
        "core_one_rank_comparisons": len(report["checks"]),
        "integer_row_module_and_primitivity_checks": sum("exact_integer_row_modules_equal" in row for row in report["checks"]),
        "individual_core_one_sign_rows": sum(row.get("core_one_sign_rows_checked", 0) for row in report["checks"]),
        "all_assertions_passed": True}
    assert sha256(proof_path.read_bytes()).hexdigest() == proof_hash, "Proof document changed during verification; rerun after editing."
    report["seconds"] = round(perf_counter() - started, 6)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"summary": report["summary"], "seconds": report["seconds"]}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
