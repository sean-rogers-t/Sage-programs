"""Exact higher-head full-D gluing and unweighted-prefix finite checks.

Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_higher_gluing_filtration.py
All modules/intersections/saturations are over ZZ. In particular, equal
rational row spaces or equal prefix Smith lists do not count as the exact
strictness assertion M_t=Sat(M_t) intersection M_a. None of the finite
checks proves that assertion outside the recorded parameter range.
"""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import signal
from time import perf_counter

from sage.all import GF, Partition, Partitions, ZZ, gcd, matrix, vector
from sage.env import SAGE_VERSION

from type_d_chern_frontier import BinaryChernModel, typed_matrices
from verify_type_d_even_gluing import (analyze_matrices, cartesian_matrices,
                                      three_head_typed_matrices)


def exact_prefix_checks(a, k, degree, data):
    """Compute every full unweighted typed-row module and its exact saturation."""
    actual = data["full"]
    assert actual.base_ring() is ZZ
    weights = [sum(part > k for part in partition)
               for partition, _kind in data["typed_rows"]]
    assert all(actual[i, j] % 2**weights[i] == 0
               for i in range(actual.nrows()) for j in range(actual.ncols()))
    unweighted = matrix(ZZ, actual.nrows(), actual.ncols(),
                        lambda i, j: actual[i, j]//2**weights[i])
    modules = [unweighted.matrix_from_rows(
        [i for i, weight in enumerate(weights) if weight <= t]).row_module()
        for t in range(a+1)]
    full = modules[-1]
    ambient = ZZ**actual.ncols()
    checks = []
    for t, module in enumerate(modules):
        saturation = module.saturation()
        intersection = saturation.intersection(full)
        assert all(value.base_ring() is ZZ for value in (module, saturation, full, intersection))
        selected = [i for i, weight in enumerate(weights) if weight <= t]
        projection_primitive = []
        for first, last in ((0, data["plus"].ncols()),
                            (data["plus"].ncols(), actual.ncols())):
            projected = unweighted.matrix_from_rows_and_columns(
                selected, list(range(first, last))).row_module()
            assert projected.base_ring() is ZZ
            projection_primitive.append(projected == projected.saturation())
        factors = [int(abs(value)) for value in module.basis_matrix().elementary_divisors() if value]
        checks.append({
            "t": t, "rank": int(module.rank()), "nonunit_factors": [value for value in factors if value != 1],
            "exact_strictness_Mt_equals_SatMt_intersection_Ma": module == intersection,
            "two_times_saturation_contained_in_Mt": all(2*value in module for value in saturation.basis()),
            "square_and_Euler_projection_modules_primitive": projection_primitive})
    # Known original unit-intersection classes. Check their divisibility even
    # after removing ALL Schubert powers of two, a stronger finite assertion.
    h, v = k//2, a//2
    relations, representatives = [], []
    for mu in Partitions(degree//2, max_part=a, max_length=h):
        if len(mu) != h or mu[-1] < v+1:
            continue
        gamma = tuple(Partition(mu).conjugate())
        delta = tuple(Partition([part-1 for part in mu]).conjugate())
        first = data["square_partitions"].index(gamma)
        second = data["plus"].ncols()+data["Euler_partitions"].index(delta)
        relation = vector(GF(2), actual.ncols())
        relation[first] = relation[second] = 1
        relations.append(relation)
        primitive = gcd(list(actual.column(first))) == gcd(list(actual.column(second))) == 1
        even_actual = all(value % 2 == 0 for value in actual.column(first)+actual.column(second))
        even_unweighted = all(value % 2 == 0 for value in unweighted.column(first)+unweighted.column(second))
        representatives.append({"mu": list(map(int, mu)), "square_head_partition": list(map(int, gamma)),
                                "Euler_head_partition": list(map(int, delta)),
                                "both_actual_sector_classes_primitive": primitive,
                                "actual_half_sum_integral": even_actual,
                                "unweighted_half_sum_integral": even_unweighted})
    explicit_space = (matrix(GF(2), relations).row_space() if relations else
                      matrix(GF(2), 0, actual.ncols()).row_space())
    return {"prefixes": checks,
            "level_zero_detects_entire_quotient": modules[0].saturation()+full == ambient,
            "explicit_unit_relations": representatives,
            "unweighted_mod2_kernel_exactly_explicit_unit_relation_space":
                unweighted.change_ring(GF(2)).right_kernel() == explicit_space,
            "base_ring": "ZZ"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("type_d_higher_gluing_filtration_verification.json"))
    parser.add_argument("--component-seconds", type=int, default=90)
    parser.add_argument("--total-seconds", type=int, default=600)
    args = parser.parse_args()
    started = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    inputs = [here, here.with_name("type_d_chern_frontier.py"),
              here.with_name("verify_type_d_even_gluing.py"),
              here.with_name("verify_type_d_three_row_smith.py"),
              here.with_name("verify_type_d_two_row_smith.py"),
              here.with_name("verify_type_d_candidate.py"),
              here.with_name("type_d_cartesian_frontier.py"),
              here.with_name("verify_type_d_square_comparison.py"),
              borel/"src"/"sage_borel"/"core.py", borel/"src"/"sage_borel"/"operations.py"]
    hashes = {str(path.relative_to(borel)): sha256(path.read_bytes()).hexdigest() for path in inputs}
    report = {"schema_version": 1, "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "source_sha256": hashes,
              "scope": "Finite higher-head evidence, including exact ZZ prefix intersections. The general even-tail theorem is proved separately in TYPE_D_EVEN_GLUING_THEOREM.md. The stronger original-head prefix strictness checked here is not inferred from that theorem.",
              "normalization": "q_i=e_i(x_tail^2), p=prod(x_tail); polynomial degree d means cohomological degree 2d",
              "resource_limits": {"component_seconds": args.component_seconds, "total_seconds": args.total_seconds},
              "generalized_model_crosschecks_with_existing_three_head_model": [],
              "completed_components": [], "excluded_components": []}

    def save():
        report["seconds"] = round(perf_counter()-started, 6)
        payload = json.dumps(report, indent=2)+"\n"
        try:
            args.output.write_text(payload)
        except OSError:
            fallback = Path("/tmp/type_d_higher_gluing_filtration_verification.json")
            fallback.write_text(payload)
            print("Report checkpoint fallback: %s" % fallback, flush=True)

    # Compare all small a=3 source columns in the exact same head alphabet.
    for k in (2, 4):
        models = BinaryChernModel(3, k), BinaryChernModel(2, k)
        for degree in range(0, 6*k+1, 2):
            general, previous = typed_matrices(3, k, degree, models), three_head_typed_matrices(k, degree)
            assert general["square_partitions"] == previous["square_partitions"]
            assert general["Euler_partitions"] == previous["Euler_partitions"]
            assert general["full"] == previous["full"]
            report["generalized_model_crosschecks_with_existing_three_head_model"].append(
                {"a": 3, "k": k, "polynomial_degree": degree, "exact_integer_matrix_equality": True})
    save()

    def timeout(_signal, _frame):
        raise TimeoutError("component wall-clock bound reached")

    signal.signal(signal.SIGALRM, timeout)
    cartesian_cache = {}
    for a, k in ((4, 2), (4, 4), (4, 6), (5, 2), (5, 4)):
        models = BinaryChernModel(a, k), BinaryChernModel(a-1, k)
        for degree in range(0, 2*a*k+1, 2):
            stamp = perf_counter()
            remaining = int(args.total_seconds-(stamp-started))
            if remaining <= 0:
                report["excluded_components"].append({"a": a, "k": k, "polynomial_degree": degree,
                                                       "reason": "overall wall-clock bound"})
                continue
            signal.alarm(min(args.component_seconds, remaining))
            try:
                data = typed_matrices(a, k, degree, models)
                row = analyze_matrices(a, k, degree, data)
                row.update(exact_prefix_checks(a, k, degree, data))
                I = row["Gaussian_unit_overlap"]
                assert len(row["explicit_unit_relations"]) == I
                row["all_prefix_nonunit_lists_equal_Gaussian_I_copies_of_two"] = all(
                    check["nonunit_factors"] == [2]*I for check in row["prefixes"])
                cross = (a, k, degree) in ((4, 2, 6), (4, 2, 8), (4, 4, 18), (5, 4, 14))
                if cross:
                    direct = cartesian_matrices(a, k, degree, cartesian_cache, force_head=True)
                    assert direct["full"].row_module() == data["full"].row_module(), (a, k, degree)
                row["independent_FULL_Cartesian_ZZ_row_module_crosscheck"] = cross
                row["seconds"] = round(perf_counter()-stamp, 6)
                report["completed_components"].append(row)
                if not row["FULL_gluing_exhausted_by_unit_overlap"] or not all(
                        check["exact_strictness_Mt_equals_SatMt_intersection_Ma"] for check in row["prefixes"]):
                    print("COUNTEREXAMPLE: %s" % row, flush=True)
                print("a=%d,k=%d,d=%d: W=%d,I=%d, exact ZZ prefixes checked" %
                      (a, k, degree, row["FULL_saturated_gluing_dimension"], I), flush=True)
            except TimeoutError as error:
                report["excluded_components"].append({"a": a, "k": k, "polynomial_degree": degree,
                                                       "reason": str(error), "seconds": round(perf_counter()-stamp, 6)})
            finally:
                signal.alarm(0)
            save()
    assert all(sha256(path.read_bytes()).hexdigest() == hashes[str(path.relative_to(borel))] for path in inputs)
    report["proof_snapshot_sha256_at_finalization"] = {
        name: sha256((borel/name).read_bytes()).hexdigest()
        for name in ("TYPE_D_GLUING_FILTRATION.md", "TYPE_D_CANDIDATE.md", "TYPE_D_SQUARE_THEOREM.md",
                     "TYPE_D_EVEN_GLUING_THEOREM.md", "TYPE_D_RELATIVE_EULER_EXTENSION.md",
                     "TYPE_D_HEAD_EULER_B_COMPARISON.md")}
    rows = report["completed_components"]
    prefixes = [check for row in rows for check in row["prefixes"]]
    report["summary"] = {
        "full_components": len(rows), "exact_ZZ_prefix_modules": len(prefixes),
        "independent_FULL_Cartesian_crosschecks": sum(row["independent_FULL_Cartesian_ZZ_row_module_crosscheck"] for row in rows),
        "generalized_model_exact_matrix_crosschecks": len(report["generalized_model_crosschecks_with_existing_three_head_model"]),
        "FULL_gluing_counterexamples": sum(not row["FULL_gluing_exhausted_by_unit_overlap"] for row in rows),
        "full_candidate_mismatches": sum(not row["full_candidate_matches"] for row in rows),
        "exact_strictness_failures": sum(not check["exact_strictness_Mt_equals_SatMt_intersection_Ma"] for check in prefixes),
        "level_zero_detection_failures": sum(not row["level_zero_detects_entire_quotient"] for row in rows),
        "nonprimitive_sector_prefix_projections": sum(not all(check["square_and_Euler_projection_modules_primitive"]) for check in prefixes),
        "prefix_nonunit_count_failures": sum(not row["all_prefix_nonunit_lists_equal_Gaussian_I_copies_of_two"] for row in rows),
        "explicit_unit_relations": sum(len(row["explicit_unit_relations"]) for row in rows),
        "explicit_unweighted_unit_relation_failures": sum(not relation["unweighted_half_sum_integral"]
            for row in rows for relation in row["explicit_unit_relations"]),
        "unweighted_mod2_kernel_description_failures": sum(
            not row["unweighted_mod2_kernel_exactly_explicit_unit_relation_space"] for row in rows),
        "excluded_components": len(report["excluded_components"])}
    save()
    print(report["summary"], flush=True)


if __name__ == "__main__":
    main()
