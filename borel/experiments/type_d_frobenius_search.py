"""Exact finite checks for the type-D mod-2 Springer image.

Run from borel/: PYTHONPATH=src sage -python experiments/type_d_frobenius_search.py

The independent part uses only ordinary Grassmannian Schur multiplication
over GF(2).  In R=H*(Gr(k,a+k);GF(2)), compute the sum of the Frobenius
image and e_k times the Frobenius image for Gr(k,a+k-1).  For even k,
compare their intersection with a proposed Gaussian formula and explicit
Schur-square basis.  These finite tests do not prove that the integral
gluing graph consists only of unit-unit pairs or prove the conjectural
type-D square-sector Smith comparison.
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
SCHUR = SymmetricFunctions(F2).schur()
Q = PolynomialRing(ZZ, "q")
q = Q.gen()


@lru_cache(maxsize=None)
def square(partition):
    return SCHUR(Partition(partition)) ** 2


def intersection_candidate(a, k):
    if k % 2:
        return Q.zero()
    h, u, v = k // 2, (a + 1) // 2, a // 2
    return Q(q ** (h * (v + 1)) * q_binomial(h + u - 1, u - 1, q))


def grassmannian_data(a, k, degree, check_basis=False):
    """Return exact mod-2 ranks; degree is polynomial degree in x."""
    target = tuple(Partitions(degree, max_part=a, max_length=k))
    plus_source = tuple(Partitions(degree // 2, max_part=a, max_length=k)) if degree % 2 == 0 else ()
    minus_source = (tuple(Partitions((degree - k) // 2, max_part=a - 1, max_length=k))
                    if degree >= k and (degree - k) % 2 == 0 else ())
    plus = matrix(F2, len(target), len(plus_source),
                  lambda i, j: square(tuple(plus_source[j]))[target[i]])
    minus = matrix(F2, len(target), len(minus_source),
                   lambda i, j: square(tuple(minus_source[j]))[
                       Partition(x - 1 for x in target[i] if x > 1)]
                   if len(target[i]) == k else F2.zero())
    rp, rm, total = plus.rank(), minus.rank(), plus.augment(minus).rank()
    result = {"a": a, "k": k, "polynomial_degree": degree,
              "square_rank": int(rp), "euler_rank": int(rm),
              "full_unit_rank": int(total), "intersection_rank": int(rp + rm - total)}
    if k % 2 == 0 and degree % 2 == 0:
        expected = int(intersection_candidate(a, k)[degree // 2])
        result.update(gaussian_intersection_rank=expected,
                      gaussian_matches=expected == rp + rm - total)
        assert result["gaussian_matches"], result
    if check_basis:
        h, v = k // 2, a // 2
        candidates = [mu for mu in Partitions(degree // 2, max_part=a, length=h)
                      if mu[-1] >= v + 1]
        columns, identities = [], []
        for mu in candidates:
            nu = Partition(x - 1 for x in mu if x > 1)
            left = {tuple(lam): int(c) for lam, c in square(tuple(mu))
                    if len(lam) <= k and (not lam or lam[0] <= a)}
            right = {tuple(x + 1 for x in list(lam) + [0] * (k - len(lam))): int(c)
                     for lam, c in square(tuple(nu))
                     if len(lam) <= k and (not lam or lam[0] <= a - 1)}
            assert left == right and left, (a, k, mu, left, right)
            columns.append([left.get(tuple(lam), 0) for lam in target])
            identities.append({"mu": list(map(int, mu)), "euler_preimage": list(map(int, nu)),
                               "exact_identity_matches": True})
        basis = matrix(F2, len(target), len(columns), lambda i, j: columns[j][i])
        assert basis.rank() == len(columns) == rp + rm - total, result
        result.update(explicit_intersection_basis=identities,
                      explicit_basis_rank=int(basis.rank()))
    return result


def source_hash(path):
    return {"file": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-a", type=int, default=7)
    parser.add_argument("--max-h", type=int, default=3)
    parser.add_argument("--seconds-per-degree", type=float, default=45)
    parser.add_argument("--skip-fresh-bgg", action="store_true")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("type_d_frobenius_search_results.json"))
    args = parser.parse_args()
    started = perf_counter()
    base = Path(__file__).parent
    recorded_file = base / "orthogonal_sector_search_results.json"
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "degree_convention": "polynomial degree d is cohomological degree 2d",
              "scope": "Finite exact corroboration. Gaussian intersections and unit-only integral gluing are not inferred as general theorems from this experiment.",
              "arguments": vars(args) | {"output": str(args.output)},
              "source_hashes": [source_hash(Path(__file__)), source_hash(recorded_file)],
              "grassmannian_checks": [], "recorded_D_checks": [], "fresh_D_checks": []}

    def persist():
        report["seconds"] = round(perf_counter() - started, 6)
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    for a in range(1, args.max_a + 1):
        for h in range(1, args.max_h + 1):
            k = 2 * h
            for degree in range(0, a * k + 1, 2):
                report["grassmannian_checks"].append(grassmannian_data(a, k, degree, check_basis=True))
            print("Grassmannian a={} k={} Gaussian {} passed".format(a, k, intersection_candidate(a, k)), flush=True)
            persist()

    recorded = json.loads(recorded_file.read_text())
    for case in recorded["cases"]:
        a, k = case["a"], case["k"]
        for row in case["degrees"]:
            if row["status"] != "computed":
                continue
            expected = grassmannian_data(a, k, row["polynomial_degree"])
            actual = row["combined"]["power_of_two_exponents"].count(0)
            expected.update(n=case["n"], actual_full_unit_rank=actual,
                            actual_square_unit_rank=row["square"]["power_of_two_exponents"].count(0),
                            actual_euler_unit_rank=row["euler"]["power_of_two_exponents"].count(0))
            assert expected["square_rank"] == expected["actual_square_unit_rank"], expected
            assert expected["euler_rank"] == expected["actual_euler_unit_rank"], expected
            assert expected["full_unit_rank"] == actual, expected
            expected["all_unit_ranks_match"] = True
            report["recorded_D_checks"].append(expected)
    persist()

    if not args.skip_fresh_bgg:
        from sage_borel import GeneratorAlgebra, classical_borel
        from orthogonal_gluing_flags import analyze
        from orthogonal_sector_search import image_matrix
        from orthogonal_smith_search import DegreeTimeout, deadline
        for name in ("orthogonal_gluing_flags.py", "orthogonal_sector_search.py", "orthogonal_smith_search.py"):
            report["source_hashes"].append(source_hash(base / name))
        parameters = [(7, 2, [6, 8, 10, 12, 14, 16]), (7, 4, [6, 8, 10, 12]),
                      (8, 4, [8, 10, 12]), (9, 4, [12, 14, 16])]
        for n, k, degrees in parameters:
            borel = classical_borel("D", n, k)
            algebra = GeneratorAlgebra(borel, {"q{}".format(i): borel.tail_polynomial(k, i)
                                              for i in range(1, k + 1)})
            p = borel.euler_polynomial(k)
            for degree in degrees:
                case_started = perf_counter()
                expected = grassmannian_data(n - k, k, degree)
                try:
                    with deadline(args.seconds_per_degree):
                        plus = tuple(algebra.evaluate(m) for m in algebra.monomials(degree))
                        minus = tuple(borel.expand(p * algebra.substitute(m)) for m in algebra.monomials(degree - k))
                        indices = tuple(sorted({w for value in plus + minus for w in value}, key=borel._word))
                        actual = analyze(image_matrix(borel, plus, indices)[1], image_matrix(borel, minus, indices)[1])
                    expected.update(status="computed", n=n, **actual)
                    expected["all_unit_ranks_match"] = (
                        expected["square_rank"] == actual["plus_exponents"].count(0)
                        and expected["euler_rank"] == actual["minus_exponents"].count(0)
                        and expected["full_unit_rank"] == actual["actual_combined_exponents"].count(0))
                    assert expected["all_unit_ranks_match"], expected
                    expected["unit_only_gluing_observed"] = all(
                        pair["plus_exponent"] == pair["minus_exponent"] == 0
                        for pair in actual["mixed_flag_pairs"])
                    expected["gaussian_matches_full_gluing_dimension"] = expected["intersection_rank"] == actual["gluing_dimension"]
                    print("Fresh D{} k={} d={} pairs {}".format(n, k, degree, actual["mixed_flag_pairs"]), flush=True)
                except DegreeTimeout:
                    expected.update(n=n, status="timeout", seconds_per_degree=args.seconds_per_degree)
                    print("Fresh D{} k={} d={} timed out".format(n, k, degree), flush=True)
                expected["seconds"] = round(perf_counter() - case_started, 6)
                report["fresh_D_checks"].append(expected)
                persist()

    report["summary"] = {
        "grassmannian_degree_checks": len(report["grassmannian_checks"]),
        "explicit_schur_square_identities": sum(len(row["explicit_intersection_basis"]) for row in report["grassmannian_checks"]),
        "recorded_D_unit_rank_comparisons": len(report["recorded_D_checks"]),
        "fresh_D_computed": sum(row["status"] == "computed" for row in report["fresh_D_checks"]),
        "fresh_D_timeouts": sum(row["status"] == "timeout" for row in report["fresh_D_checks"]),
        "fresh_nonunit_gluing_observed": sum(not row["unit_only_gluing_observed"] for row in report["fresh_D_checks"] if row["status"] == "computed"),
        "all_assertions_passed": True}
    persist()
    print(json.dumps(report["summary"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
