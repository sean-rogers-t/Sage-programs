#!/usr/bin/env -S sage -python
"""Compute a bounded grid of joint Levi-factor images, using exact arithmetic.

Run from the ``borel`` directory with ``PYTHONPATH=src sage -python
experiments/two_step_symplectic.py``.  This computes a fresh report and writes
it after every completed case.  It does not use previously recorded ranks or
relations as expected answers.  No stability theorem or novelty claim follows
from this finite experiment.
"""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from time import perf_counter

from sage.env import SAGE_VERSION

from sage_borel.families import two_step_symplectic


def _relation_data(algebra, relation):
    return {
        "polynomial": str(relation),
        "terms": [
            {"exponents": list(map(int, exponents)), "coefficient": str(coefficient)}
            for exponents, coefficient in sorted(relation.dict().items(), reverse=True)
        ],
    }


def run_case(n, b, k, max_degree, include_ambient=False):
    started = perf_counter()
    family = two_step_symplectic(n, b, k)
    algebra = family.joint_algebra()
    predicted_hilbert = family.predicted_hilbert_prefix(max_degree)
    predicted_relations = family.predicted_relations()
    setup_seconds = perf_counter() - started
    rows = []
    for degree in range(max_degree + 1):
        degree_started = perf_counter()
        analysis = algebra.degree(degree, include_ambient=include_ambient)
        if analysis.image_dimension != predicted_hilbert[degree]:
            raise ArithmeticError("BGG image rank disagrees with presentation Hilbert coefficient")
        row = {
            "polynomial_degree": degree,
            "cohomological_degree": 2 * degree,
            "domain_dimension": int(analysis.domain_dimension),
            "image_dimension": int(analysis.image_dimension),
            "kernel_dimension": int(analysis.kernel_dimension),
            "schubert_support_size": len(analysis.schubert_indices),
            "predicted_image_dimension": predicted_hilbert[degree],
            "relations": [_relation_data(algebra, relation)
                          for relation in analysis.relations],
            "seconds": round(perf_counter() - degree_started, 6),
        }
        if include_ambient:
            row["ambient_dimension"] = int(analysis.ambient_dimension)
        rows.append(row)
        print("C{} b={} k={} degree {}: domain={}, image={}, kernel={} ({:.3f}s)".format(
            n, b, k, degree, row["domain_dimension"], row["image_dimension"],
            row["kernel_dimension"], row["seconds"]), flush=True)
    relation_degrees = [row["polynomial_degree"] for row in rows
                        if row["kernel_dimension"]]
    presentation = []
    for r, relation in enumerate(predicted_relations, family.a + 1):
        data = _relation_data(algebra, relation)
        data["polynomial_degree"] = 2 * r
        data["bgg_evaluated_within_bound"] = 2 * r <= max_degree
        if data["bgg_evaluated_within_bound"]:
            if algebra.evaluate(relation):
                raise ArithmeticError("presentation relation has nonzero BGG image")
            data["bgg_image_zero"] = True
        presentation.append(data)
    return {
        "parameters": family.family_parameters(),
        "generator_names": list(algebra.generator_names),
        "generator_degrees": list(map(int, algebra.generator_degrees)),
        "input_polynomials": {name: str(polynomial)
                              for name, polynomial in family.joint_generators().items()},
        "maximum_polynomial_degree": max_degree,
        "first_observed_relation_degree": min(relation_degrees) if relation_degrees else None,
        "domain_hilbert_prefix": [row["domain_dimension"] for row in rows],
        "image_hilbert_prefix": [row["image_dimension"] for row in rows],
        "kernel_hilbert_prefix": [row["kernel_dimension"] for row in rows],
        "predicted_image_hilbert_prefix": list(predicted_hilbert),
        "presentation_relations": presentation,
        "predicted_first_relation_degree": 2 * (family.a + 1),
        "all_bgg_ranks_match_presentation_through_bound": True,
        "setup_seconds": round(setup_seconds, 6),
        "total_seconds": round(perf_counter() - started, 6),
        "degrees": rows,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=8,
                        help="largest polynomial degree (cohomological degree is twice this)")
    parser.add_argument("--n-count", type=int, default=5,
                        help="number of n values starting at b+k+1, for each b,k")
    parser.add_argument("--b", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--k", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--include-ambient", action="store_true",
                        help="also enumerate the complete Schubert basis in each degree")
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("two_step_symplectic_results.json"))
    args = parser.parse_args(argv)
    if args.max_degree < 0 or args.n_count < 1 or min(args.b + args.k) < 1:
        parser.error("degree must be nonnegative; n-count, b and k must be positive")
    report = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sage_version": SAGE_VERSION,
        "arithmetic": "QQ (exact rational arithmetic)",
        "degree_convention": "polynomial degree d is cohomological degree 2d",
        "scope": "Finite experiments only; null first relation means none through the reported bound.",
        "grid": {"b_values": args.b, "k_values": args.k,
                 "n_rule": "b+k+1 through b+k+n_count, inclusive",
                 "n_count": args.n_count, "maximum_polynomial_degree": args.max_degree},
        "ambient_basis_enumerated": args.include_ambient,
        "cases": [],
        "complete": False,
    }
    started = perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for b in args.b:
        for k in args.k:
            for n in range(b + k + 1, b + k + 1 + args.n_count):
                report["cases"].append(run_case(n, b, k, args.max_degree, args.include_ambient))
                report["total_seconds"] = round(perf_counter() - started, 6)
                args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    report["complete"] = True
    report["total_seconds"] = round(perf_counter() - started, 6)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("Saved {} exact cases to {} ({:.3f}s)".format(
        len(report["cases"]), args.output, report["total_seconds"]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
