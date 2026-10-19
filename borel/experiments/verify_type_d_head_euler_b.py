"""Independent finite checks for the relocated B head-Euler argument.

Checks actual H>K prefix lattices using only Schur products, and signed
level-zero identities using ordinary symmetric-function plethysm.
Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_head_euler_b.py
"""
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import GF, Partition, Partitions, QQ, SymmetricFunctions, ZZ, matrix
from sage.env import SAGE_VERSION
from verify_type_b_filtration import schur_product_matrix, factors
from verify_type_d_candidate import type_b_counts, unit_intersection

SF = SymmetricFunctions(ZZ)
S, COMPLETE = SF.schur(), SF.homogeneous()
QSF = SymmetricFunctions(QQ)
QS, POWER = QSF.schur(), QSF.power()


def partitions(width, height, degree):
    return tuple(tuple(p) for p in Partitions(degree, max_part=width, max_length=height))


def corner(H, K, t, gamma):
    n, b = H+K, H-t
    if not b:
        return 0
    r, c, eps = n//2, (n+1)//2, n % 2
    for j in range(max(0, b-r), min(b, c, (b+eps)//2)+1):
        bounds = (c-j,)*j+(r-b+j,)*(b-j)
        if all(value <= bound for value, bound in zip(gamma[t:], bounds)):
            return j
    return None


def bead_level(H, K, gamma):
    n = H+K
    beads = [gamma[i]+H-i-1 for i in range(H)]
    d = sum(value < n//2 for value in beads)
    return sum(value >= (n+1)//2+d for value in beads)


def actual_head_check(H, K, degree):
    source = partitions(K, H, degree)
    padded = [p+(0,)*(H-len(p)) for p in source]
    index = {p: i for i, p in enumerate(padded)}
    result = []
    for t in range(H+1):
        chosen = [(p, corner(H, K, t, p)) for p in padded
                  if corner(H, K, t, p) is not None]
        assert all((corner(H, K, t, p) is not None) == (bead_level(H, K, p) <= t)
                   for p in padded)
        chosen.sort(key=lambda item: (sum(value*value for value in item[0]), item[0]))
        schur_rows, h_rows = [], []
        for gamma, j in chosen:
            mu = tuple(value for value in gamma[:t] if value)
            alpha = tuple(value for value in gamma[t:t+j] if value)
            beta = tuple(value for value in gamma[t+j:] if value)
            value = S(COMPLETE(mu))*S(alpha)*S(beta)
            schur_rows.append([ZZ(value[Partition(p)]) for p in source])
            coordinates = [ZZ(0)]*len(source)
            for p, coefficient in COMPLETE(value):
                if p and p[0] > K:
                    continue
                pp = tuple(p)+(0,)*(H-len(p))
                coordinates[index[pp]] += coefficient
                if pp != gamma and coefficient:
                    assert sum(x*x for x in pp) > sum(x*x for x in gamma)
            assert coordinates[index[gamma]] == 1
            h_rows.append(coordinates)
        constructed = matrix(ZZ, len(chosen), len(source),
                             [x for row in schur_rows for x in row])
        h_matrix = matrix(ZZ, len(chosen), len(source), [x for row in h_rows for x in row])
        minor = h_matrix.matrix_from_columns([index[p] for p, _j in chosen])
        assert minor.is_square() and minor.det() == 1
        balanced, _metadata = schur_product_matrix(K, H, degree, t, source)
        assert constructed.row_module() == balanced.row_module()
        assert all(value == 1 for value in factors(balanced))
        expected = sum(count for exponent, count in type_b_counts(H, K, 2*degree).items()
                       if exponent <= t)
        assert len(chosen) == expected
        result.append({"threshold": t, "rank": len(chosen),
                       "feasible_corner_bead_criterion": True,
                       "unit_triangular_minor": True,
                       "exact_integer_balanced_module_equality": True,
                       "prefix_primitive_over_ZZ": True})
    return {"actual_head": H, "tail": K, "half_polynomial_degree": degree,
            "prefixes": result}


@lru_cache(None)
def adams2(p):
    value = QS(QS(p).plethysm(POWER([2])))
    assert all(coefficient.denominator() == 1 for _q, coefficient in value)
    return S(value)


def level_zero_check(H, K, half_degree):
    h = H//2
    target = partitions(K, H, 2*half_degree)
    plus = partitions(K, H, half_degree)
    minus = partitions(K-1, H, half_degree-h) if half_degree >= h else ()
    def plus_coefficient(p, gamma):
        return ZZ(adams2(p)[Partition(gamma)])
    def minus_coefficient(p, gamma):
        padded = list(gamma)+[0]*(H-len(gamma))
        if min(padded) == 0:
            return ZZ(0)
        previous = tuple(value-1 for value in padded if value > 1)
        return ZZ(adams2(p)[Partition(previous)])
    A = matrix(ZZ, len(target), len(plus), lambda i, j: plus_coefficient(plus[j], target[i]))
    B = matrix(ZZ, len(target), len(minus), lambda i, j: minus_coefficient(minus[j], target[i]))
    M = A.augment(B)
    assert all(value == 1 for value in factors(A))
    assert all(value == 1 for value in factors(B))
    expected = unit_intersection(K, H, 2*half_degree)
    rational_overlap = A.rank()+B.rank()-M.rank()
    mod2_overlap = A.change_ring(GF(2)).rank()+B.change_ring(GF(2)).rank()-M.change_ring(GF(2)).rank()
    assert rational_overlap == mod2_overlap == expected
    assert M.rank() == M.change_ring(GF(2)).rank()
    v = K//2
    explicit = [theta for theta in partitions(K, h, half_degree)
                if len(theta) == h and theta[-1] >= v+1]
    assert len(explicit) == expected
    for theta in explicit:
        earlier = tuple(value-1 for value in theta if value > 1)
        assert theta in plus and earlier in minus
        assert A.column(plus.index(theta)) == (-1)**h*B.column(minus.index(earlier))
    return {"actual_head": H, "tail": K, "half_polynomial_degree": half_degree,
            "rational_overlap": int(rational_overlap), "mod2_overlap": int(mod2_overlap),
            "Gaussian_overlap": expected, "exact_signed_overlap_identities": len(explicit),
            "combined_image_primitive_at_two": True,
            "combined_integer_nonunit_factors": [int(value) for value in factors(M) if value != 1]}


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    dependencies = [here, here.with_name("verify_type_b_filtration.py"),
                    here.with_name("verify_type_d_candidate.py")]
    proof = here.parent.parent / "TYPE_D_HEAD_EULER_B_COMPARISON.md"
    source_hashes = {path.name: sha256(path.read_bytes()).hexdigest()
                     for path in dependencies}
    proof_hashes = {proof.name: sha256(proof.read_bytes()).hexdigest()}
    report = {"generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "source_sha256": source_hashes,
              "proof_sha256": proof_hashes,
              "scope": "Finite independent Schur/LR and ordinary plethysm checks; no all-rank theorem inferred.",
              "actual_head_checks": [], "level_zero_checks": []}
    cases = [(2,1), (3,1), (3,2), (4,1), (4,2), (4,3),
             (5,1), (5,2), (5,3), (6,2), (6,3)]
    for H, K in cases:
        for degree in range(H*K+1):
            report["actual_head_checks"].append(actual_head_check(H, K, degree))
        print("H=%d K=%d: all actual-head prefix modules passed" % (H, K), flush=True)
    for H, K in [(2,K) for K in range(1,7)]+[(4,K) for K in range(1,6)]+[(6,2), (6,3), (6,4)]:
        for degree in range(H*K//2+1):
            report["level_zero_checks"].append(level_zero_check(H, K, degree))
        print("H=%d K=%d: all level-zero signed/2-local checks passed" % (H, K), flush=True)
    report["status"] = "passed"
    report["summary"] = {"actual_head_degree_components": len(report["actual_head_checks"]),
                         "actual_head_prefix_module_comparisons": sum(len(row["prefixes"]) for row in report["actual_head_checks"]),
                         "level_zero_degree_components": len(report["level_zero_checks"]),
                         "exact_signed_overlap_identities": sum(row["exact_signed_overlap_identities"] for row in report["level_zero_checks"]),
                         "counterexamples": 0}
    assert source_hashes == {path.name: sha256(path.read_bytes()).hexdigest()
                             for path in dependencies}
    assert proof_hashes == {proof.name: sha256(proof.read_bytes()).hexdigest()}
    report["inputs_stable_during_run"] = True
    report["seconds"] = round(perf_counter()-started, 6)
    here.with_name("type_d_head_euler_b_verification.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(report["summary"]), flush=True)


if __name__ == "__main__":
    main()
