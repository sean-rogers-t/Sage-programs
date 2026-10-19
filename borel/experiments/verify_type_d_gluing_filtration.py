"""Exact finite tests of the conditional type-D filtered row criterion.

Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_gluing_filtration.py
The all-rank strictness hypotheses remain unproved. These calculations
test actual ZZ module intersections, not just matching Smith lists.
"""
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import GF, ZZ, matrix
from sage.env import SAGE_VERSION
from verify_type_d_even_gluing import three_head_typed_matrices
from verify_type_d_three_row_smith import ThreeRootChernModel, smith_factors
from verify_type_d_two_row_smith import TwoRootChernModel
from verify_type_d_candidate import unit_intersection


def integer_matrix(rows, columns):
    return matrix(ZZ, len(rows), columns, [entry for row in rows for entry in row])


def exponents(coefficients):
    answer = []
    for factor in smith_factors(coefficients):
        exponent = int(factor.valuation(2))
        assert factor == 2**exponent
        answer.append(exponent)
    return answer


def check_component(K, degree, models):
    data = three_head_typed_matrices(K, degree, models)
    actual = data["full"]
    weights = [sum(part > K for part in shape)
               for shape, _kind in data["typed_rows"]]
    assert all(all(entry % 2**level == 0 for entry in actual.row(i))
               for i, level in enumerate(weights))
    unweighted = matrix(ZZ, actual.nrows(), actual.ncols(),
                        lambda i, j: actual[i, j]//2**weights[i])
    plus_rank, minus_rank = data["plus"].ncols(), data["minus"].ncols()
    total_rank = plus_rank+minus_rank
    ambient = ZZ**total_rank
    modules, sectors, rows = [], [], []
    for threshold in range(4):
        indices = [i for i, level in enumerate(weights) if level <= threshold]
        prefix = unweighted.matrix_from_rows(indices)
        plus = prefix.matrix_from_columns(range(plus_rank))
        minus = prefix.matrix_from_columns(range(plus_rank, total_rank))
        assert all(value == 1 for value in smith_factors(plus))
        assert all(value == 1 for value in smith_factors(minus))
        sector_rows = [list(row)+[ZZ(0)]*minus_rank
                       for row in plus.row_module().basis_matrix().rows()]
        sector_rows += [[ZZ(0)]*plus_rank+list(row)
                        for row in minus.row_module().basis_matrix().rows()]
        sector = integer_matrix(sector_rows, total_rank).row_module()
        module = prefix.row_module()
        assert module.saturation() == sector
        assert all(2*row in module for row in sector.basis_matrix().rows())
        modules.append(module)
        sectors.append(sector)
        factors = smith_factors(prefix)
        assert all(value in (1, 2) for value in factors)
        rows.append({"threshold": threshold, "rank": int(module.rank()),
                     "sector_projection_prefixes_primitive": True,
                     "saturation_is_sector_direct_sum": True,
                     "nonunit_factors": [int(value) for value in factors if value != 1]})
    assert sectors[-1] == ambient
    original_overlap = unit_intersection(3, K, degree)
    quotient_dimension = sum(exponents(unweighted))
    assert quotient_dimension == original_overlap
    for threshold in range(4):
        strict = modules[threshold] == sectors[threshold].intersection(modules[-1])
        assert strict
        assert len(rows[threshold]["nonunit_factors"]) == quotient_dimension
        rows[threshold]["exact_module_intersection_strictness"] = True
    detected = sectors[0]+modules[-1] == ambient
    assert detected
    explicit_lift = None
    if original_overlap:
        assert original_overlap == 1
        h, j = K//2, degree//2-K
        gamma = tuple(value for value in (h, h, j) if value)
        delta = tuple(value for value in (h, j) if value)
        plus_index = data["square_partitions"].index(gamma)
        minus_index = plus_rank+data["Euler_partitions"].index(delta)
        relation = unweighted.column(plus_index)+unweighted.column(minus_index)
        assert all(value % 2 == 0 for value in relation)
        assert any((row[plus_index]+row[minus_index]) % 2
                   for row in sectors[0].basis_matrix().rows())
        explicit_lift = {"square_partition": list(gamma), "Euler_partition": list(delta),
                         "unit_relation_already_even_in_unweighted_full_matrix": True,
                         "detected_nontrivially_by_level_zero_sector_prefix": True}
    # The exact weighted-prefix module must be the actual Schubert row module.
    weighted_rows = [list(2**threshold*row) for threshold, module in enumerate(modules)
                     for row in module.basis_matrix().rows()]
    weighted = integer_matrix(weighted_rows, total_rank)
    assert weighted.row_module() == actual.row_module()
    expected, previous = [], 0
    for threshold, sector in enumerate(sectors):
        expected += [threshold]*(sector.rank()-previous)
        previous = sector.rank()
    for _ in range(quotient_dimension):
        expected.remove(0)
        expected.append(1)
    expected.sort()
    assert exponents(weighted) == expected
    return {"a": 3, "k": K, "polynomial_degree": degree,
            "source_sector_ranks": [plus_rank, minus_rank],
            "unweighted_full_quotient_dimension": quotient_dimension,
            "Gaussian_original_unit_overlap": original_overlap,
            "all_quotient_characters_detected_at_level_zero": detected,
            "weighted_prefix_module_equals_actual_full_module": True,
            "conditional_formula_matches_actual_exponents": True,
            "actual_full_exponents": expected, "explicit_unit_lift": explicit_lift,
            "prefixes": rows}


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    dependencies = [here, here.with_name("verify_type_d_even_gluing.py"),
                    here.with_name("verify_type_d_three_row_smith.py"),
                    here.with_name("verify_type_d_two_row_smith.py"),
                    here.with_name("verify_type_d_candidate.py")]
    report = {"generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION,
              "scope": "Finite exact ZZ tests of conditional filtered-row strictness; no all-rank assertion.",
              "source_sha256": {path.name: sha256(path.read_bytes()).hexdigest()
                                 for path in dependencies}, "checks": []}
    for K in (2, 4, 6, 8, 10, 12):
        models = ThreeRootChernModel(K), TwoRootChernModel(K)
        for degree in range(0, 6*K+1, 2):
            report["checks"].append(check_component(K, degree, models))
        print("k=%d: all %d degree components passed exact module tests" %
              (K, 3*K+1), flush=True)
    report["status"] = "passed"
    report["summary"] = {
        "head_three_tail_sizes": 6, "degree_components": len(report["checks"]),
        "exact_prefix_module_intersections": 4*len(report["checks"]),
        "level_zero_quotient_detection_tests": len(report["checks"]),
        "weighted_module_and_conditional_formula_checks": len(report["checks"]),
        "counterexamples": 0}
    report["seconds"] = round(perf_counter()-started, 6)
    destination = here.with_name("type_d_gluing_filtration_verification.json")
    destination.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report["summary"]), flush=True)


if __name__ == "__main__":
    main()
