"""Finite exact corroboration of the general type-B integral squeeze.

Run: sage -python experiments/verify_type_b_general_smith.py
No Borel engine, Weyl group, or divided difference is imported.
"""

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import Partition, SymmetricFunctions, ZZ, matrix
from sage.env import SAGE_VERSION
from verify_type_b_high_layers import add, candidate, padded_partitions, same
from verify_type_b_filtration import schur_product_matrix


def corner(A, ell, t, gamma):
    b = ell - t
    if not b:
        return 0
    n = A + ell
    r, c, epsilon = n//2, (n+1)//2, n % 2
    for j in range((b + epsilon)//2 + 1):
        bounds = (c-j,)*j + (r-b+j,)*(b-j)
        if all(x <= y for x, y in zip(gamma[t:], bounds)):
            return j
    return None


def bead_statistic(A, ell, gamma):
    n = A + ell
    r, c = n//2, (n+1)//2
    beads = tuple(gamma[i] + ell-i-1 for i in range(ell))
    d = sum(x < r for x in beads)
    return sum(x >= c+d for x in beads)


def bead_checks():
    checks, comparisons = [], 0
    for n in range(2, 19):
        for ell in range(1, min(6, n//2) + 1):
            A = n - ell
            counts = [[0]*(A*ell+1) for _ in range(ell+1)]
            for gamma in padded_partitions(A, ell):
                e = bead_statistic(A, ell, gamma)
                counts[e][sum(gamma)] += 1
                for t in range(ell+1):
                    assert (corner(A, ell, t, gamma) is not None) == (e <= t), (
                        n, ell, t, gamma, e)
                    comparisons += 1
            for e, polynomial in enumerate(counts):
                assert same(polynomial, candidate(A, ell, e)), (n, ell, e)
            checks.append({"n": n, "ell": ell, "partitions": sum(map(sum, counts)),
                           "all_memberships_and_weight_counts_match": True})
    return checks, comparisons


def integer_module_checks():
    sf = SymmetricFunctions(ZZ)
    schur, complete = sf.schur(), sf.homogeneous()
    checks = []
    parameters = [
        (5, 5, 1, [6, 12, 16]),
        (6, 5, 1, [6, 15, 20]),
        (6, 6, 1, [7, 18, 24]),
        (7, 6, 1, [7, 21, 26]),
        (7, 7, 1, [8, 24, 32]),
        (8, 7, 2, [8, 18, 26]),
        (8, 8, 1, [9, 14, 18]),
    ]
    for A, ell, t, degrees in parameters:
        for degree in degrees:
            # Partitions are the actual ambient complete-coordinate basis.
            padded = [p for p in padded_partitions(A, ell) if sum(p) == degree]
            source = tuple(tuple(x for x in p if x) for p in padded)
            source_index = {p: i for i, p in enumerate(padded)}
            selected = [(p, corner(A, ell, t, p)) for p in padded
                        if corner(A, ell, t, p) is not None]
            selected.sort(key=lambda row: (sum(x*x for x in row[0]), row[0]))
            selected_columns = [source_index[p] for p, _ in selected]
            h_rows, s_rows = [], []
            used_corners = set()
            for gamma, j in selected:
                used_corners.add(j)
                mu = tuple(x for x in gamma[:t] if x)
                alpha = tuple(x for x in gamma[t:t+j] if x)
                beta = tuple(x for x in gamma[t+j:] if x)
                value = schur(complete(mu)) * schur(alpha) * schur(beta)
                s_rows.append([ZZ(value[Partition(p)]) for p in source])
                h_value = complete(value)
                h_row = [ZZ(0)]*len(source)
                energy = sum(x*x for x in gamma)
                for partition, coefficient in h_value:
                    if partition and partition[0] > A:
                        continue
                    assert len(partition) <= ell
                    coordinate = tuple(partition) + (0,)*(ell-len(partition))
                    h_row[source_index[coordinate]] += coefficient
                    if coordinate != gamma and coefficient:
                        assert sum(x*x for x in coordinate) > energy
                assert h_row[source_index[gamma]] == 1
                h_rows.append(h_row)
            dim, size = len(source), len(selected)
            h_matrix = matrix(ZZ, size, dim, [v for row in h_rows for v in row])
            minor = h_matrix.matrix_from_columns(selected_columns)
            assert all(minor[i, j] == 0 for i in range(size) for j in range(i))
            assert minor.det() == 1
            constructed = matrix(ZZ, size, dim, [v for row in s_rows for v in row])
            balanced, metadata = schur_product_matrix(A, ell, degree, t, source)
            equal = constructed.row_module() == balanced.row_module()
            assert equal, (A, ell, t, degree)
            predicted_rank = sum(int(candidate(A, ell, e)[degree])
                                 if degree < len(candidate(A, ell, e)) else 0
                                 for e in range(t+1))
            assert predicted_rank == size
            check = {"A": A, "ell": ell, "t": t, "m": degree,
                     "source_rank": dim, "prefix_rank": size,
                     "used_unbalanced_corners": sorted(used_corners),
                     "unit_triangular_minor": True,
                     "integer_balanced_and_constructed_modules_equal": equal,
                     "candidate_cumulative_rank": predicted_rank,
                     "balanced_nonzero_rows": metadata["nonzero_product_rows"]}
            checks.append(check)
            print(json.dumps(check), flush=True)
    return checks


def main():
    started = perf_counter()
    here = Path(__file__).resolve().parent
    report = {"status": "running", "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION,
              "scope": "finite corroboration; general proof is TYPE_B_GENERAL_SMITH.md",
              "Borel_or_Weyl_engine_imported": False,
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    dependencies = [here/"verify_type_b_high_layers.py",
                    here/"verify_type_b_filtration.py",
                    here.parent/"TYPE_B_GENERAL_SMITH.md",
                    here.parent/"TYPE_B_DEFECT_COMPARISON.md"]
    report["dependencies"] = [
        {"file": str(p.relative_to(here.parent)), "sha256": sha256(p.read_bytes()).hexdigest()}
        for p in dependencies]
    report["bead_checks"], count = bead_checks()
    report["integer_module_checks"] = integer_module_checks()
    report["summary"] = {
        "bead_parameter_cases": len(report["bead_checks"]),
        "corner_membership_comparisons": count,
        "exact_integer_LR_module_comparisons": len(report["integer_module_checks"]),
        "unit_triangular_minors": len(report["integer_module_checks"]),
        "counterexamples": 0}
    report["status"] = "passed"
    report["seconds"] = round(perf_counter()-started, 6)
    destination = here/"type_b_general_smith_verification.json"
    destination.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"status": report["status"], **report["summary"],
                      "seconds": report["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
