#!/usr/bin/env -S sage -python
"""Compare square and Euler sectors of orthogonal integral Borel images.

The type-D tail algebra is Z[q_1,...,q_k,p]/(p^2-q_k), with p the Euler
product and deg(p)=k.  Its degree-d lattice is consequently the sum of the
square-only columns q^alpha and the columns p*q^alpha.  This script computes
both sectors in the SAME integral Schubert basis, compares their Smith
factors with B_n and B_(n-1), and measures the index caused by gluing their
rational eigenspaces.  It does not assume a general Smith formula.
"""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from time import perf_counter

from sage.all import ZZ, matrix, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, GeneratorAlgebra, classical_borel
from orthogonal_smith_search import DegreeTimeout, deadline, factor_data, save


def smith_data(coefficients):
    if not all(coefficients.dimensions()):
        return ()
    diagonal, left, right = coefficients.smith_form()
    if left * coefficients * right != diagonal:
        raise ArithmeticError("Smith transformation failed its exact check")
    return tuple(abs(ZZ(diagonal[i, i])) for i in range(min(diagonal.dimensions()))
                 if diagonal[i, i])


def image_matrix(borel, expansions, indices=None):
    if indices is None:
        indices = tuple(sorted({w for expansion in expansions for w in expansion}, key=borel._word))
    coefficients = matrix(ZZ, len(indices), len(expansions),
                          {(row, column): expansion.get(word, 0)
                           for row, word in enumerate(indices)
                           for column, expansion in enumerate(expansions)})
    return indices, coefficients


def sector_summary(factors):
    return {"rank": len(factors), "saturation_index": str(prod(factors)),
            **factor_data(factors)}


def type_b_algebra(n, k):
    borel = (classical_borel("B", n, k) if n > k else
             BorelMap(["B", n], omitted_nodes=[]))
    return GeneratorAlgebra(borel, {"q{}".format(i): borel.tail_polynomial(k, i)
                                  for i in range(1, k + 1)})


def encoded_matrix(coefficients):
    return {"rows": coefficients.nrows(), "columns": coefficients.ncols(),
            "entries": [{"row": int(row), "column": int(column), "value": str(value)}
                        for (row, column), value in sorted(coefficients.dict().items()) if value]}


def run_case(n, k, args, report, reference, started):
    case_started = perf_counter()
    borel = classical_borel("D", n, k)
    square = GeneratorAlgebra(borel, {"q{}".format(i): borel.tail_polynomial(k, i)
                                    for i in range(1, k + 1)})
    euler = borel.euler_polynomial(k)
    big_b, small_b = type_b_algebra(n, k), type_b_algebra(n - 1, k)
    top = int(borel.longest_representative().length())
    bound = top if n <= args.complete_through_rank else min(top, args.max_degree)
    degrees = [degree for degree in range(bound + 1)
               if square.monomials(degree) or degree >= k and square.monomials(degree - k)]
    case = {"n": n, "k": k, "a": n - k, "maximum_polynomial_degree": bound,
            "euler_parity": "odd" if k % 2 else "even", "degrees": [], "complete": False}
    report["cases"].append(case)
    save(report, args.output, started)
    for degree in degrees:
        degree_started = perf_counter()
        try:
            with deadline(args.seconds_per_degree):
                square_report = square.integral_degree(degree)
                square_expansions = tuple(square.evaluate(monomial) for monomial in square.monomials(degree))
                euler_expansions = tuple(
                    borel.expand(euler * square.substitute(monomial))
                    for monomial in (square.monomials(degree - k) if degree >= k else ()))
                indices = tuple(sorted({word for expansion in square_expansions + euler_expansions
                                        for word in expansion}, key=borel._word))
                _, q_matrix = image_matrix(borel, square_expansions, indices)
                _, p_matrix = image_matrix(borel, euler_expansions, indices)
                combined_matrix = q_matrix.augment(p_matrix)
                q_factors = square_report.smith_factors
                p_factors = smith_data(p_matrix)
                combined_factors = smith_data(combined_matrix)
                b_factors = big_b.integral_degree(degree).smith_factors
                smaller_factors = (small_b.integral_degree(degree - k).smith_factors
                                   if degree >= k else ())
            independent_ranks = len(combined_factors) == len(q_factors) + len(p_factors)
            product_index = ZZ(prod(q_factors) * prod(p_factors))
            index = ZZ(prod(combined_factors))
            if independent_ranks and index % product_index:
                raise ArithmeticError("combined index is not divisible by the sector indices")
            row = {"polynomial_degree": degree, "status": "computed",
                   "square": sector_summary(q_factors), "euler": sector_summary(p_factors),
                   "combined": sector_summary(combined_factors),
                   "B_n_reference": sector_summary(b_factors),
                   "B_n_minus_one_shifted_reference": sector_summary(smaller_factors),
                   "square_matches_B_n": q_factors == b_factors,
                   "euler_matches_B_n_minus_one_shifted": p_factors == smaller_factors,
                   "sector_ranks_add": independent_ranks,
                   "integral_eigenspace_gluing_index": str(index // product_index) if independent_ranks else None,
                   "schubert_support_size": len(indices),
                   "square_columns": q_matrix.ncols(), "euler_columns": p_matrix.ncols(),
                   "seconds": round(perf_counter() - degree_started, 6)}
            recorded = reference.get((n, k, degree))
            if recorded is not None:
                row["agrees_with_full_generator_algebra_search"] = (
                    row["combined"]["smith_factors"] == recorded["smith_factors"])
            if args.include_matrices or n == 4 and k == 2 and degree == 4:
                row["schubert_indices"] = [list(map(int, word.reduced_word())) for word in indices]
                row["square_matrix"] = encoded_matrix(q_matrix)
                row["euler_matrix"] = encoded_matrix(p_matrix)
            print("D{} k={} d={}: Q{} E{} all{} glue{} B={} smallerB={} ({:.3f}s)".format(
                n, k, degree, row["square"]["power_of_two_multiplicities"],
                row["euler"]["power_of_two_multiplicities"],
                row["combined"]["power_of_two_multiplicities"],
                row["integral_eigenspace_gluing_index"], row["square_matches_B_n"],
                row["euler_matches_B_n_minus_one_shifted"], row["seconds"]), flush=True)
        except DegreeTimeout as error:
            row = {"polynomial_degree": degree, "status": "timeout", "message": str(error),
                   "seconds": round(perf_counter() - degree_started, 6)}
            case["degrees"].append(row)
            case["total_seconds"] = round(perf_counter() - case_started, 6)
            save(report, args.output, started)
            print("D{} k={} d={} TIMEOUT".format(n, k, degree), flush=True)
            return
        case["degrees"].append(row)
        case["total_seconds"] = round(perf_counter() - case_started, 6)
        save(report, args.output, started)
    case["complete"] = True
    save(report, args.output, started)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-rank", type=int, default=4)
    parser.add_argument("--max-rank", type=int, default=6)
    parser.add_argument("--max-degree", type=int, default=12)
    parser.add_argument("--complete-through-rank", type=int, default=4)
    parser.add_argument("--seconds-per-degree", type=float, default=45)
    parser.add_argument("--include-matrices", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("orthogonal_sector_search_results.json"))
    args = parser.parse_args(argv)
    if args.min_rank < 4 or args.max_rank < args.min_rank or args.max_degree < 0:
        parser.error("invalid rank or degree bound")
    started = perf_counter()
    original = Path(__file__).with_name("orthogonal_smith_search_results.json")
    reference = {}
    if original.exists():
        reference = {(case["n"], case["k"], row["polynomial_degree"]): row
                     for case in json.loads(original.read_text())["cases"] if case["family"] == "D"
                     for row in case["degrees"] if row["status"] == "computed"}
    if args.resume and args.output.exists():
        report = json.loads(args.output.read_text())
    else:
        report = {"schema_version": 1, "sage_version": SAGE_VERSION,
                  "generated_utc": datetime.now(timezone.utc).isoformat(),
                  "scope": "Finite sector and integral gluing experiment; no general Smith theorem assumed.",
                  "degree_convention": "polynomial degree d is cohomological degree 2d",
                  "normalization": "unscaled elementary squared tail generators and Euler product",
                  "ambient_basis_enumerated": False, "cases": [], "runs": []}
    report["runs"].append({"arguments": vars(args) | {"output": str(args.output)},
                           "started_utc": datetime.now(timezone.utc).isoformat()})
    existing = {(case["n"], case["k"]) for case in report["cases"]}
    for n in range(args.min_rank, args.max_rank + 1):
        for k in range(2, n):
            if (n, k) not in existing:
                run_case(n, k, args, report, reference, started)
    rows = [row for case in report["cases"] for row in case["degrees"] if row["status"] == "computed"]
    report["summary"] = {
        "cases_recorded": len(report["cases"]), "degrees_computed": len(rows),
        "square_B_mismatches": sum(not row["square_matches_B_n"] for row in rows),
        "euler_smaller_B_mismatches": sum(not row["euler_matches_B_n_minus_one_shifted"] for row in rows),
        "sector_rank_additivity_failures": sum(not row["sector_ranks_add"] for row in rows),
        "nontrivial_gluing_degrees": sum(row["integral_eigenspace_gluing_index"] != "1" for row in rows),
        "full_original_search_mismatches": sum(row.get("agrees_with_full_generator_algebra_search") is False for row in rows),
        "timeouts": sum(row["status"] == "timeout" for case in report["cases"] for row in case["degrees"]),
        "odd_factors_observed": any(row[sector]["odd_factor_present"] for row in rows
                                    for sector in ("square", "euler", "combined")),
    }
    save(report, args.output, started)
    print("Saved {}".format(report["summary"]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
