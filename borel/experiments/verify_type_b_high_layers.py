"""Exact finite checks for TYPE_B_HIGH_LAYERS.md.

Run with Sage Python. Polynomial and graph checks use integer arithmetic
without the Borel engine; a separate oracle compares exact ZZ LR modules.
Finite checks corroborate the displayed proofs, not arbitrary-rank claims.
"""

from collections import defaultdict
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
from time import perf_counter


@lru_cache(None)
def gaussian(n, k):
    if not 0 <= k <= n:
        return (0,)
    if not k or k == n:
        return (1,)
    return tuple(add(gaussian(n - 1, k - 1), gaussian(n - 1, k), k))


def add(a, b, shift=0, sign=1):
    result = list(a) + [0] * max(0, len(b) + shift - len(a))
    for i, value in enumerate(b):
        result[i + shift] += sign * value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def multiply(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def same(a, b):
    return add(a, b, sign=-1) == [0]


def units(n, ell):
    r, c = n // 2, (n + 1) // 2
    result = [0]
    for j in range(ell // 2 + 1):
        result = add(result, multiply(gaussian(r, j),
                                     gaussian(c - j, ell - 2*j)),
                     j * (c - ell + j))
    return result


def candidate(A, ell, e):
    n = A + ell
    return add([0], multiply(gaussian(n // 2, e),
                             units(n - 2*e, ell - e)),
               e * ((n + 1) // 2))


def cumulative(A, ell, t):
    result = [0]
    for e in range(t + 1):
        result = add(result, candidate(A, ell, e))
    return result


def prefix_three_formula(A, ell):
    n = A + ell
    r, c = n // 2, (n + 1) // 2
    result = [0]
    for j in range(ell - 2):
        s = ell - j
        small = list(gaussian(c + s - 1, s))
        if n % 2:
            correction = add([1] * (c - 1), [1] * c, 1)
            small = add(small, correction, (s - 1)*(c - 1) - 1, -1)
        else:
            small = add(small, [1] * c, (s - 1)*(c - 1), -1)
        result = add(result, multiply(gaussian(r - ell + j, j), small),
                     j*c)
    return result


def padded_partitions(A, ell):
    return [tuple(reversed(p))
            for p in combinations_with_replacement(range(A + 1), ell)]


def graph_components(A, ell):
    """All rows (5), including isolated source coordinates."""
    n = A + ell
    r, c = n // 2, (n + 1) // 2
    parts = padded_partitions(A, ell)
    index = {p: i for i, p in enumerate(parts)}
    parent = list(range(len(parts)))
    anchors = set()

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for mu in padded_partitions(A, ell - 3):
        for u in range(c - 1):
            for v in range(u + 1):
                for w in range(r):
                    start = tuple(sorted(mu + (u, v, w), reverse=True))
                    i = index[start]
                    if not v:
                        anchors.add(i)
                    else:
                        end = tuple(sorted(mu + (u + 1, v - 1, w),
                                           reverse=True))
                        left, right = find(i), find(index[end])
                        if left != right:
                            parent[right] = left
    anchored = {find(i) for i in anchors}
    groups = defaultdict(list)
    for i, p in enumerate(parts):
        groups[find(i)].append(p)
    return [(values, root in anchored) for root, values in groups.items()]


def graph_rank(A, ell, groups):
    result = [0] * (A * ell + 1)
    for values, anchored in groups:
        degree = sum(values[0])
        assert all(sum(v) == degree for v in values)
        result[degree] += len(values) - (not anchored)
    return result


def lr_module_checks():
    """Compare with independent LR products, not their graph expansion."""
    from sage.all import Partition, SymmetricFunctions, ZZ, matrix
    from verify_type_b_filtration import schur_product_matrix

    sf = SymmetricFunctions(ZZ)
    schur, complete = sf.schur(), sf.homogeneous()
    checks = []
    for A, ell in [(4, 4), (5, 5), (6, 5), (6, 6)]:
        groups = graph_components(A, ell)
        degrees = sorted({0, A, A+1, A*ell//2,
                          A*ell//2+1, A*ell-1, A*ell})
        all_parts = padded_partitions(A, ell)
        for degree in degrees:
            source = tuple(tuple(x for x in p if x)
                           for p in all_parts if sum(p) == degree)
            dim = len(source)
            change_rows = []
            for p in source:
                value = schur(complete(p))
                change_rows.append([ZZ(value[Partition(x)]) for x in source])
            change = matrix(ZZ, dim, dim,
                            [x for row in change_rows for x in row])
            assert abs(change.det()) == 1
            padded_index = {p: i for i, p in
                            enumerate(p for p in all_parts if sum(p) == degree)}
            graph_rows = []
            for values, anchored in groups:
                if sum(values[0]) != degree:
                    continue
                root = padded_index[values[0]]
                if anchored:
                    row = [0] * dim
                    row[root] = 1
                    graph_rows.append(row)
                for value in values[1:]:
                    row = [0] * dim
                    row[padded_index[value]], row[root] = 1, -1
                    graph_rows.append(row)
            graph_matrix = matrix(ZZ, len(graph_rows), dim,
                                  [x for row in graph_rows for x in row]) * change
            lr, _ = schur_product_matrix(A, ell, degree, ell - 3, source)
            equal_three = graph_matrix.row_module() == lr.row_module()
            r, c = (A+ell)//2, (A+ell+1)//2
            coordinate_rows = [change[i, :] for p, i in padded_index.items()
                               if p[-2] <= c-1 and p[-1] <= r-1]
            coordinates = matrix(ZZ, len(coordinate_rows), dim,
                                 [x for row in coordinate_rows for x in row])
            lr_two, _ = schur_product_matrix(A, ell, degree, ell - 2, source)
            equal_two = coordinates.row_module() == lr_two.row_module()
            assert equal_three and equal_two, (A, ell, degree)
            checks.append({"A": A, "ell": ell, "m": degree,
                           "unimodular_complete_change": True,
                           "b3_graph_equals_integer_LR_module": equal_three,
                           "b2_coordinates_equal_integer_LR_module": equal_two})
    return checks


def main():
    started = perf_counter()
    report = {"status": "running", "generated_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "finite corroboration of proof, no Borel or Weyl engine",
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "polynomial_checks": [], "graph_checks": []}
    for ell in range(3, 11):
        for A in range(ell, ell + 16):
            formula = prefix_three_formula(A, ell)
            assert same(formula, cumulative(A, ell, ell - 3)), (A, ell)
            r, c = (A+ell)//2, (A+ell+1)//2
            # Independently count the b=2 coordinate shape.
            if A <= ell+2 and ell <= 7:
                coordinates = [0] * (A*ell+1)
                for p in padded_partitions(A, ell):
                    if p[-2] <= c-1 and p[-1] <= r-1:
                        coordinates[sum(p)] += 1
                assert same(coordinates, cumulative(A, ell, ell - 2))
            report["polynomial_checks"].append({"A": A, "ell": ell,
                                               "prefix_matches": True})
    for ell in range(3, 8):
        for A in range(ell, ell + 3):
            groups = graph_components(A, ell)
            ranks = graph_rank(A, ell, groups)
            assert same(ranks, prefix_three_formula(A, ell)), (A, ell)
            report["graph_checks"].append({"A": A, "ell": ell,
                                          "vertices": sum(len(v) for v, _ in groups),
                                          "degree_components": A*ell+1,
                                          "matches": True})
    report["integer_LR_checks"] = lr_module_checks()
    report["summary"] = {
        "polynomial_parameter_cases": len(report["polynomial_checks"]),
        "full_graph_parameter_cases": len(report["graph_checks"]),
        "graph_degree_components": sum(x["degree_components"] for x in report["graph_checks"]),
        "integer_LR_row_module_equalities": 2*len(report["integer_LR_checks"]),
        "counterexamples": 0}
    report["status"] = "passed"
    report["seconds"] = round(perf_counter() - started, 6)
    destination = Path(__file__).with_name("type_b_high_layers_verification.json")
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], **report["summary"],
                      "seconds": report["seconds"]}), flush=True)


if __name__ == "__main__":
    main()
