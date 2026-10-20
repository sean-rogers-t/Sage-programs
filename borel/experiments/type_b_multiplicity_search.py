#!/usr/bin/env -S sage -python
"""Higher-degree type-B Smith multiplicities from rectangular Schur columns.

For a=n-k, the integral symmetric-function quotient by h_(a+1),...,h_n is
free on s_lambda with lambda inside (a^k).  Its rational characteristic-map
kernel is that ideal, so this integral source basis spans exactly the same
image lattice as all elementary-generator monomials.  This avoids computing
many columns that differ by source relations.

For the same squared-coordinate input, C and B Schubert coefficient matrices
satisfy B_w=2^nu(w)*C_w, where nu is the number of negative signed-coordinate
entries.  This follows from the simple-root denominators in divided
differences and is directly checked in the smaller calculations here.

The experiment also tests (without assuming) the candidate multiplicities
rank_F2(C_{nu<=t})-rank_F2(C_{nu<t}).  Integer Smith factors of these filtered
C matrices are recorded to distinguish mod-2 evidence from primitivity.
"""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from time import perf_counter

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix
from sage.env import SAGE_VERSION

from sage_borel import classical_borel
from sage_borel.operations import elementary_polynomial
from orthogonal_sector_search import encoded_matrix, smith_data
from orthogonal_smith_search import DegreeTimeout, deadline, factor_data
from type_b_cartesian_frontier import CartesianCFrontier


def rectangle_partitions(size, width, length):
    def visit(remaining, maximum, slots, prefix):
        if remaining == 0:
            yield prefix
        elif slots:
            for first in range(min(remaining, maximum), 0, -1):
                yield from visit(remaining - first, first, slots - 1, prefix + (first,))
    return tuple(visit(size, width, length, ()))


def schur_formal(partition, ring, complete):
    if not partition:
        return ring.one()
    length = len(partition)
    return matrix(ring, length, length,
                  lambda i, j: complete[partition[i] - i + j]
                  if partition[i] - i + j >= 0 else ring.zero()).det()


def highest_exponent_candidate(a, k, m, exponents):
    short, long = min(a, k), max(a, k)
    threshold = (long + short - 1) // 2
    expected = sum(len(partition) == short and partition[-1] > threshold
                   for partition in rectangle_partitions(m, long, short))
    actual = list(exponents).count(short)
    return {"highest_exponent": short, "last_part_threshold": threshold,
            "predicted_multiplicity": int(expected), "actual_multiplicity": actual,
            "matches": expected == actual}


def run_case(n, k, args, report, reference, persist):
    case_started = perf_counter()
    a, maximum = n - k, k * (n - k)
    c, b = classical_borel("C", n, k), classical_borel("B", n, k)
    ring = PolynomialRing(QQ, k, names=["q{}".format(i) for i in range(1, k + 1)])
    q = ring.gens()
    complete = [ring.one()]
    for degree in range(1, n):
        complete.append(ring(sum((-1) ** (i - 1) * q[i - 1] * complete[degree - i]
                                 for i in range(1, min(k, degree) + 1))))
    c_map = ring.hom([c.tail_polynomial(k, i) for i in range(1, k + 1)], c.ring)
    b_map = ring.hom([b.tail_polynomial(k, i) for i in range(1, k + 1)], b.ring)
    cartesian_ring = PolynomialRing(QQ, n, names=["x{}".format(i) for i in range(1, n + 1)])
    squared_tail = tuple(x ** 2 for x in cartesian_ring.gens()[-k:])
    cartesian_map = ring.hom([elementary_polynomial(cartesian_ring, squared_tail, i)
                              for i in range(1, k + 1)], cartesian_ring)
    cartesian = CartesianCFrontier(c)
    # Larger rectangles first receive the familiar low degrees and their top
    # class; remaining degrees are attempted while the case budget permits.
    degrees = list(range(maximum + 1)) if n <= 6 else (
        list(range(min(6, maximum) + 1)) + [maximum]
        + [m for m in range(7, maximum) if m != maximum])
    existing = next((case for case in report["cases"] if case["n"] == n and case["k"] == k), None)
    if existing is None:
        case = {"n": n, "k": k, "a": a, "maximum_half_polynomial_degree": maximum,
                "requested_half_polynomial_degrees": degrees, "degrees": [], "complete": False,
                "source_basis": "s_lambda(tail^2), lambda contained in (a^k), |lambda|=m"}
        report["cases"].append(case)
    else:
        case = existing
        case.setdefault("previous_attempts", []).extend(row for row in case["degrees"] if row["status"] != "computed")
        case["degrees"] = [row for row in case["degrees"] if row["status"] == "computed"]
        case.pop("uncomputed_half_polynomial_degrees", None)
        case.pop("stopped_for_budget", None)
    previous_seconds = case.get("total_seconds", 0)
    already_computed = {row["half_polynomial_degree"] for row in case["degrees"]}
    persist()
    timeouts = 0
    for m in degrees:
        if m in already_computed:
            continue
        if perf_counter() - case_started > args.seconds_per_case or timeouts >= 2:
            case["uncomputed_half_polynomial_degrees"] = [d for d in degrees
                if d not in {row["half_polynomial_degree"] for row in case["degrees"]}]
            case["stopped_for_budget"] = True
            break
        row_started = perf_counter()
        partitions = rectangle_partitions(m, a, k)
        try:
            with deadline(args.seconds_per_degree):
                formal = tuple(schur_formal(partition, ring, complete) for partition in partitions)
                expansions = tuple((cartesian.expand(cartesian_map(polynomial)) if args.backend == "cartesian"
                                    else c.expand(c_map(polynomial))) for polynomial in formal)
                indices = tuple(sorted({w for expansion in expansions for w in expansion}, key=c._word))
                weights = [sum(entry < 0 for entry in word.matrix().list()) for word in indices]
                if any(weight != list(word.reduced_word()).count(n) for weight, word in zip(weights, indices)):
                    raise ArithmeticError("signed-coordinate count differs from end-node word count")
                c_matrix = matrix(ZZ, len(indices), len(partitions),
                                  {(i, j): expansion.get(word, 0) for i, word in enumerate(indices)
                                   for j, expansion in enumerate(expansions)})
                b_matrix = matrix(ZZ, len(indices), len(partitions),
                                  lambda i, j: (2 ** weights[i]) * c_matrix[i, j])
                c_factors, b_factors = smith_data(c_matrix), smith_data(b_matrix)
                direct_columns = []
                if n <= args.direct_b_through_rank or 2 * m <= args.direct_b_max_degree:
                    direct_columns = list(range(len(partitions)))
                for column in direct_columns:
                    actual = b.expand(b_map(formal[column]))
                    expected = {b.W.from_reduced_word(word.reduced_word()): b_matrix[row, column]
                                for row, word in enumerate(indices) if b_matrix[row, column]}
                    if actual.as_dict() != expected:
                        raise ArithmeticError("B/C row-scaling identity failed")
                c_direct_columns = (list(range(len(partitions))) if 2 * m <= 4 else
                                    [0] if n <= 6 and partitions else [])
                for column in c_direct_columns:
                    if c.expand(c_map(formal[column])) != expansions[column]:
                        raise ArithmeticError("Cartesian and fundamental-weight C frontiers disagree")
                ranks, predicted, previous = [], [], 0
                for threshold in range(max(weights, default=0) + 1):
                    filtered = c_matrix.matrix_from_rows([row for row, weight in enumerate(weights)
                                                         if weight <= threshold])
                    rank = filtered.change_ring(GF(2)).rank()
                    filtered_factors = smith_data(filtered)
                    ranks.append({"threshold": threshold, "rank_mod_two": int(rank),
                                  "integer_smith_factors": [str(factor) for factor in filtered_factors],
                                  "integer_image_is_primitive": all(factor == 1 for factor in filtered_factors)})
                    predicted.extend([threshold] * (rank - previous))
                    previous = rank
                actual_data = factor_data(b_factors)
            row = {"half_polynomial_degree": m, "polynomial_degree": 2 * m,
                   "cohomological_degree": 4 * m, "status": "computed",
                   "partitions": [list(partition) for partition in partitions],
                   "source_schur_polynomials": [str(polynomial) for polynomial in formal],
                   "rank": len(b_factors), "schubert_support_size": len(indices),
                   "C_smith_factors": [str(factor) for factor in c_factors],
                   "C_image_is_primitive": all(factor == 1 for factor in c_factors),
                   "maximum_absolute_C_entry": max((abs(int(value)) for value in c_matrix.list()), default=0),
                   "C_entries_in_zero_plus_minus_one": all(abs(value) <= 1 for value in c_matrix.list()),
                   "row_negative_counts": weights,
                   "schubert_indices": [list(map(int, word.reduced_word())) for word in indices],
                   "C_matrix": encoded_matrix(c_matrix), "B_matrix": encoded_matrix(b_matrix),
                   "direct_B_columns_checked": len(direct_columns),
                   "direct_C_fundamental_columns_checked": len(c_direct_columns),
                   "coefficient_backend": args.backend,
                   "Cartesian_divided_differences_cumulative": cartesian.divided_differences,
                   "filtered_C_ranks": ranks,
                   "candidate_exponents_from_filtered_mod_two_ranks": predicted,
                   "candidate_matches_actual_smith": predicted == actual_data["power_of_two_exponents"],
                   "seconds": round(perf_counter() - row_started, 6), **actual_data}
            row["highest_exponent_partition_candidate"] = highest_exponent_candidate(
                a, k, m, row["power_of_two_exponents"])
            old = reference.get((n, k, 2 * m))
            if old is not None:
                row["matches_original_elementary_monomial_search"] = row["smith_factors"] == old["smith_factors"]
            print("B{} k={} a={} m={}: exponents{} filtered{} matches{} ({:.3f}s)".format(
                n, k, a, m, row["power_of_two_exponents"], predicted,
                row["candidate_matches_actual_smith"], row["seconds"]), flush=True)
            if not row["candidate_matches_actual_smith"]:
                report["counterexamples"].append({"n": n, "k": k, "a": a, "m": m,
                                                  "actual": row["power_of_two_exponents"], "predicted": predicted})
        except DegreeTimeout as error:
            row = {"half_polynomial_degree": m, "polynomial_degree": 2 * m,
                   "status": "timeout", "seconds": round(perf_counter() - row_started, 6),
                   "message": str(error)}
            timeouts += 1
            print("B{} k={} m={} TIMEOUT".format(n, k, m), flush=True)
        case["degrees"].append(row)
        case["total_seconds"] = round(previous_seconds + perf_counter() - case_started, 6)
        persist()
    case["complete"] = (len(case["degrees"]) == len(degrees)
                        and all(row["status"] == "computed" for row in case["degrees"]))
    case["total_seconds"] = round(previous_seconds + perf_counter() - case_started, 6)
    persist()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds-per-degree", type=float, default=30)
    parser.add_argument("--seconds-per-case", type=float, default=120)
    parser.add_argument("--direct-b-through-rank", type=int, default=6)
    parser.add_argument("--direct-b-max-degree", type=int, default=4)
    parser.add_argument("--backend", choices=("cartesian", "fundamental"), default="cartesian")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--case", action="append", dest="selected_cases",
                        help="specific n:k pair; repeat for a custom grid")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("type_b_multiplicity_search_results.json"))
    args = parser.parse_args(argv)
    started = perf_counter()
    reference_file = Path(__file__).with_name("orthogonal_smith_search_results.json")
    reference = {(case["n"], case["k"], row["polynomial_degree"]): row
                 for case in json.loads(reference_file.read_text())["cases"] if case["family"] == "B"
                 for row in case["degrees"] if row["status"] == "computed"}
    if args.resume and args.output.exists():
        report = json.loads(args.output.read_text())
        for case in report["cases"]:
            for row in case["degrees"]:
                if row["status"] == "computed":
                    row.setdefault("coefficient_backend", "fundamental")
                    row.setdefault("direct_C_fundamental_columns_checked", 0)
                    entries = [abs(int(entry["value"])) for entry in row["C_matrix"]["entries"]]
                    row.setdefault("maximum_absolute_C_entry", max(entries, default=0))
                    row.setdefault("C_entries_in_zero_plus_minus_one", all(value <= 1 for value in entries))
                    row.setdefault("highest_exponent_partition_candidate", highest_exponent_candidate(
                        case["a"], case["k"], row["half_polynomial_degree"], row["power_of_two_exponents"]))
    else:
        report = {"schema_version": 1, "sage_version": SAGE_VERSION,
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Finite higher-degree Smith and filtered-rank experiments, not a closed multiplicity proof.",
              "normalization": "unscaled elementary squared tail inputs; source changed by integral Schur relations only",
              "ambient_basis_enumerated": False,
              "budgets": {"seconds_per_degree": args.seconds_per_degree, "seconds_per_case": args.seconds_per_case},
                  "cases": [], "counterexamples": []}
    report.setdefault("runs", []).append({"backend": args.backend,
                                           "started_utc": datetime.now(timezone.utc).isoformat()})
    def persist():
        report["last_run_seconds"] = round(perf_counter() - started, 6)
        report["total_seconds"] = round(sum(case.get("total_seconds", 0) for case in report["cases"]), 6)
        temporary = args.output.with_suffix(args.output.suffix + ".tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n")
        temporary.replace(args.output)
    selected = ([tuple(map(int, value.split(":"))) for value in args.selected_cases]
                if args.selected_cases else ((5, 2), (6, 2), (6, 3), (7, 2), (7, 3), (8, 2), (8, 3), (8, 4)))
    for n, k in selected:
        run_case(n, k, args, report, reference, persist)
    rows = [row for case in report["cases"] for row in case["degrees"] if row["status"] == "computed"]
    report["summary"] = {"cases": len(report["cases"]), "complete_cases": sum(case["complete"] for case in report["cases"]),
                         "degrees_computed": len(rows), "candidate_counterexamples": len(report["counterexamples"]),
                         "direct_B_columns_checked": sum(row["direct_B_columns_checked"] for row in rows),
                         "direct_C_fundamental_columns_checked": sum(row["direct_C_fundamental_columns_checked"] for row in rows),
                         "Cartesian_degree_reports": sum(row["coefficient_backend"] == "cartesian" for row in rows),
                         "original_monomial_search_mismatches": sum(row.get("matches_original_elementary_monomial_search") is False for row in rows),
                         "nonprimitive_C_images": sum(not row["C_image_is_primitive"] for row in rows),
                         "maximum_absolute_C_entry": max((row["maximum_absolute_C_entry"] for row in rows), default=0),
                         "C_matrices_with_entry_beyond_zero_plus_minus_one": sum(not row["C_entries_in_zero_plus_minus_one"] for row in rows),
                         "nonprimitive_filtered_C_images": sum(not filtered["integer_image_is_primitive"] for row in rows for filtered in row["filtered_C_ranks"]),
                         "highest_exponent_partition_candidate_mismatches": sum(not row["highest_exponent_partition_candidate"]["matches"] for row in rows),
                         "timeouts": sum(row["status"] == "timeout" for case in report["cases"] for row in case["degrees"]),
                         "odd_B_factors": any(row["odd_factor_present"] for row in rows)}
    persist()
    print("Saved {}".format(report["summary"]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
