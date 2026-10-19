#!/usr/bin/env -S sage -python
"""Exact finite Smith searches for the unscaled orthogonal tail generators.

From borel/, run ``PYTHONPATH=src sage -python
experiments/orthogonal_smith_search.py``.  Type B uses all elementary symmetric
polynomials in the last k squared coordinates; type D additionally includes
the tail Euler product.  No Cartan denominators are cleared and no ambient
Schubert basis is requested.  Results are persisted after every degree.

This is an experiment, not a proof of a general Smith formula.  A per-degree
POSIX timer bounds expensive BGG calculations.  Degree-zero validation includes
every chosen generator, even those above the requested calculation bound.
"""

import argparse
from collections import Counter
from contextlib import contextmanager
from datetime import datetime, timezone
import json
from pathlib import Path
import signal
from time import perf_counter

from sage.env import SAGE_VERSION

from sage_borel import GeneratorAlgebra, classical_borel


class DegreeTimeout(TimeoutError):
    pass


@contextmanager
def deadline(seconds):
    def expired(signum, frame):
        raise DegreeTimeout("per-degree calculation exceeded {} seconds".format(seconds))
    previous = signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


def factor_data(factors):
    valuations = []
    odd_parts = []
    for factor in factors:
        factor = int(factor)
        exponent = 0
        while factor % 2 == 0:
            factor //= 2
            exponent += 1
        valuations.append(exponent)
        odd_parts.append(factor)
    return {
        "smith_factors": [str(int(factor)) for factor in factors],
        "power_of_two_exponents": valuations,
        "power_of_two_multiplicities": {
            str(exponent): multiplicity for exponent, multiplicity in
            sorted(Counter(valuations).items())},
        "odd_factor_present": any(part != 1 for part in odd_parts),
        "odd_parts": odd_parts,
    }


def one_head_prediction(family, n, k, degree):
    """Elementary odd/even quadric prediction, independent of Smith calls."""
    if n - k != 1:
        return None
    if family == "B":
        if degree % 2 or degree > 2 * n - 1:
            return []
        return [str(2 if degree >= n else 1)]
    if degree == k:
        return ["1", "2"] if k % 2 == 0 else ["1"]
    if degree % 2 or degree > 2 * k:
        return []
    return [str(1 if degree < k else 2)]


def save(report, output, started):
    report["last_run_seconds"] = round(perf_counter() - started, 6)
    report["total_seconds"] = round(sum(case.get("total_seconds", 0)
                                        for case in report["cases"]), 6)
    report["total_seconds_basis"] = "sum of case elapsed times across resumed runs; excludes process imports"
    report["updated_utc"] = datetime.now(timezone.utc).isoformat()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(output)


def run_case(family, n, k, args, report, output, started):
    case_started = perf_counter()
    borel = classical_borel(family, n, k)
    generators = {"q{}".format(i): borel.tail_polynomial(k, i)
                  for i in range(1, k + 1)}
    if family == "D":
        generators["euler"] = borel.euler_polynomial(k)
    algebra = GeneratorAlgebra(borel, generators)
    top = int(borel.longest_representative().length())
    if family == "B":
        selected_top = 2 * k * (n - k)
        maximum = selected_top if n <= args.complete_through_rank else min(
            selected_top, args.max_degree)
        degrees = list(range(0, maximum + 1, 2))
    else:
        maximum = top if n <= args.complete_through_rank else min(top, args.max_degree)
        degrees = [degree for degree in range(maximum + 1) if algebra.monomials(degree)]
    case = {
        "family": family, "n": n, "k": k, "a": n - k,
        "omitted_node": n - k,
        "ambient_top_polynomial_degree": top,
        "maximum_polynomial_degree": maximum,
        "requested_degrees": degrees,
        "generator_names": list(algebra.generator_names),
        "generator_degrees": list(algebra.generator_degrees),
        "normalization": "unscaled e_i(tail^2), with Euler(tail) additionally in D",
        "ambient_basis_enumerated": False,
        "degrees": [], "complete": False,
    }
    report["cases"].append(case)
    save(report, output, started)
    for degree in degrees:
        degree_started = perf_counter()
        try:
            with deadline(args.seconds_per_degree):
                analysis = algebra.integral_degree(degree)
            row = {
                "polynomial_degree": degree,
                "cohomological_degree": 2 * degree,
                "status": "computed",
                "domain_dimension": int(analysis.domain_dimension),
                "rank": int(analysis.rank),
                "schubert_support_size": len(analysis.schubert_indices),
                "saturation_index": str(analysis.saturation_index),
                "seconds": round(perf_counter() - degree_started, 6),
            }
            row.update(factor_data(analysis.smith_factors))
            prediction = one_head_prediction(family, n, k, degree)
            if prediction is not None:
                row["one_head_quadric_prediction"] = prediction
                row["one_head_quadric_prediction_matches"] = row["smith_factors"] == prediction
            print("{}{} k={} a={} d={}: rank{} 2-exponents{} ({:.3f}s)".format(
                family, n, k, n - k, degree, row["rank"],
                row["power_of_two_multiplicities"], row["seconds"]), flush=True)
        except DegreeTimeout as error:
            row = {"polynomial_degree": degree, "status": "timeout",
                   "seconds": round(perf_counter() - degree_started, 6), "message": str(error)}
            case["degrees"].append(row)
            print("{}{} k={} d={} TIMEOUT after {:.3f}s".format(
                family, n, k, degree, row["seconds"]), flush=True)
            case["total_seconds"] = round(perf_counter() - case_started, 6)
            save(report, output, started)
            return
        except Exception as error:
            row = {"polynomial_degree": degree, "status": "error",
                   "seconds": round(perf_counter() - degree_started, 6),
                   "exception": type(error).__name__, "message": str(error)}
            case["degrees"].append(row)
            case["total_seconds"] = round(perf_counter() - case_started, 6)
            save(report, output, started)
            raise
        case["degrees"].append(row)
        case["total_seconds"] = round(perf_counter() - case_started, 6)
        save(report, output, started)
    case["complete"] = True
    case["total_seconds"] = round(perf_counter() - case_started, 6)
    save(report, output, started)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--families", nargs="+", choices=("B", "D"), default=["B", "D"])
    parser.add_argument("--min-rank", type=int, default=2)
    parser.add_argument("--max-rank", type=int, default=6)
    parser.add_argument("--max-degree", type=int, default=12)
    parser.add_argument("--complete-through-rank", type=int, default=4)
    parser.add_argument("--seconds-per-degree", type=float, default=45)
    parser.add_argument("--resume", action="store_true", help="retain all recorded cases and fill missing cases")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("orthogonal_smith_search_results.json"))
    args = parser.parse_args(argv)
    if (args.min_rank < 2 or args.max_rank < args.min_rank or args.max_degree < 0
            or args.seconds_per_degree <= 0):
        parser.error("invalid rank, degree or deadline")
    started = perf_counter()
    if args.resume and args.output.exists():
        report = json.loads(args.output.read_text(encoding="utf-8"))
        report["runs"].append({"started_utc": datetime.now(timezone.utc).isoformat(),
                               "arguments": vars(args) | {"output": str(args.output)}})
    else:
        report = {
            "schema_version": 1,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "sage_version": SAGE_VERSION,
            "arithmetic": "ZZ matrices in the integral Schubert basis; exact Smith forms",
            "degree_convention": "polynomial degree d is cohomological degree 2d",
            "scope": "Finite calculations; no general Smith formula or novelty claim follows.",
            "ambient_basis_enumerated": False,
            "runs": [{"started_utc": datetime.now(timezone.utc).isoformat(),
                      "arguments": vars(args) | {"output": str(args.output)}}],
            "cases": [],
        }
    existing = {(case["family"], case["n"], case["k"]) for case in report["cases"]}
    for n in range(args.min_rank, args.max_rank + 1):
        for family in args.families:
            if family == "D" and n < 4:
                continue
            for k in range(1 if family == "B" else 2, n):
                if (family, n, k) not in existing:
                    run_case(family, n, k, args, report, args.output, started)
    rows = [row for case in report["cases"] for row in case["degrees"]]
    report["summary"] = {
        "cases_recorded": len(report["cases"]),
        "cases_complete": sum(case["complete"] for case in report["cases"]),
        "degrees_computed": sum(row["status"] == "computed" for row in rows),
        "timeouts": sum(row["status"] == "timeout" for row in rows),
        "odd_factors_observed": any(row.get("odd_factor_present", False) for row in rows),
        "one_head_prediction_mismatches": sum(
            row.get("one_head_quadric_prediction_matches") is False for row in rows),
    }
    save(report, args.output, started)
    print("Saved {}".format(report["summary"]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
