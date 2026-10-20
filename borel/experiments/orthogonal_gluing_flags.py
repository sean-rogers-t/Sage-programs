#!/usr/bin/env -S sage -python
"""Test a two-flag formula for Smith factors of integral eigenspace gluing.

For sector lattices L+ and L- with saturated bases S+ and S-, extract the
elementary-2 graph W describing the saturation of their sum.  In sector Smith
bases, pull back the flags spanned by exponents <=i.  Mixed intersection
differences give pairs (i,j).  Replace each paired sector exponents i,j by
min(i,j),max(i,j)+1, leaving unpaired exponents unchanged.

The script recomputes every nontrivial gluing example in the recorded D grid,
checks trivial examples directly from their recorded Smith lists, tests
random higher-dimensional elementary-2 gluing models, and exhaustively checks
small mod-2 flag stabilizers lift to integral lattice automorphisms.  These
finite checks support a proof but do not replace it.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from random import Random
from time import perf_counter

from sage.all import GF, QQ, ZZ, diagonal_matrix, identity_matrix, matrix
from sage.env import SAGE_VERSION

from sage_borel import GeneratorAlgebra, classical_borel
from orthogonal_sector_search import image_matrix
from orthogonal_smith_search import deadline


F2 = GF(2)


def exponents_of_two(factors):
    answer = []
    for factor in factors:
        factor = ZZ(abs(factor))
        exponent = int(factor.valuation(2))
        if factor != 2 ** exponent:
            raise ArithmeticError("this formula test requires 2-primary sector lattices")
        answer.append(exponent)
    return answer


def decomposition(coefficients):
    rows, columns = coefficients.dimensions()
    if not rows or not columns:
        return (), matrix(ZZ, rows, 0)
    diagonal, left, right = coefficients.smith_form()
    if left * coefficients * right != diagonal:
        raise ArithmeticError("Smith transformations failed their exact check")
    factors = tuple(abs(ZZ(diagonal[i, i])) for i in range(min(rows, columns)) if diagonal[i, i])
    basis = matrix(ZZ, left.inverse()).matrix_from_columns(range(len(factors)))
    return factors, basis


def analyze(plus, minus, include_graph=False):
    plus_factors, plus_saturation = decomposition(plus)
    minus_factors, minus_saturation = decomposition(minus)
    combined_factors, combined_saturation = decomposition(plus.augment(minus))
    rp, rm = len(plus_factors), len(minus_factors)
    if len(combined_factors) != rp + rm:
        raise ArithmeticError("sector rational spans are not independent")
    sector_coordinates = matrix(ZZ, combined_saturation.change_ring(QQ).solve_right(
        plus_saturation.augment(minus_saturation)))
    rank = rp + rm
    if rank:
        twice_coordinates = matrix(ZZ, 2 * sector_coordinates.change_ring(QQ).inverse())
    else:
        twice_coordinates = matrix(ZZ, 0, 0)
    graph = matrix(F2, twice_coordinates).column_space().basis_matrix().transpose()
    graph_rank = graph.ncols()
    plus_graph, minus_graph = graph[:rp, :], graph[rp:, :]
    if plus_graph.rank() != graph_rank or minus_graph.rank() != graph_rank:
        raise ArithmeticError("the gluing graph does not inject into both sectors")
    if abs(sector_coordinates.det()) != 2 ** graph_rank:
        raise ArithmeticError("gluing quotient is not elementary 2-primary")
    ep, em = exponents_of_two(plus_factors), exponents_of_two(minus_factors)
    actual = exponents_of_two(combined_factors)
    dimensions = {}
    for i in range(max(ep, default=0) + 1):
        high_plus = plus_graph.matrix_from_rows([row for row, exponent in enumerate(ep) if exponent > i])
        for j in range(max(em, default=0) + 1):
            high_minus = minus_graph.matrix_from_rows([row for row, exponent in enumerate(em) if exponent > j])
            dimensions[i, j] = graph_rank - high_plus.stack(high_minus).rank()
    plus_remaining, minus_remaining = Counter(ep), Counter(em)
    pairs, replacements = [], []
    for (i, j), dimension in sorted(dimensions.items()):
        multiplicity = (dimension - dimensions.get((i - 1, j), 0)
                        - dimensions.get((i, j - 1), 0) + dimensions.get((i - 1, j - 1), 0))
        if multiplicity < 0:
            raise ArithmeticError("negative mixed flag multiplicity")
        if multiplicity:
            pairs.append({"plus_exponent": i, "minus_exponent": j,
                          "multiplicity": int(multiplicity)})
            plus_remaining[i] -= multiplicity
            minus_remaining[j] -= multiplicity
            replacements.extend([min(i, j), max(i, j) + 1] * multiplicity)
    if any(count < 0 for count in list(plus_remaining.values()) + list(minus_remaining.values())):
        raise ArithmeticError("flag pairing consumes too many sector exponents")
    if sum(pair["multiplicity"] for pair in pairs) != graph_rank:
        raise ArithmeticError("mixed flag multiplicities do not exhaust the gluing graph")
    predicted = sorted(replacements + [exponent for exponent, count in plus_remaining.items()
                                      for _ in range(count)]
                       + [exponent for exponent, count in minus_remaining.items() for _ in range(count)])
    result = {
        "plus_exponents": ep, "minus_exponents": em,
        "actual_combined_exponents": actual, "predicted_combined_exponents": predicted,
        "matches": predicted == actual,
        "gluing_dimension": graph_rank, "gluing_index": str(2 ** graph_rank),
        "mixed_flag_pairs": pairs,
        "flag_intersection_dimensions": [{"i": i, "j": j, "dimension": int(value)}
                                          for (i, j), value in sorted(dimensions.items())],
    }
    if include_graph:
        result["graph_in_sector_smith_coordinates"] = [list(map(int, row)) for row in graph.rows()]
    return result


def random_unimodular(size, random):
    result = identity_matrix(ZZ, size)
    for _ in range(3 * size):
        if size < 2:
            break
        first, second = random.sample(range(size), 2)
        result.add_multiple_of_row(first, second, random.choice((-1, 1)))
    return result


def synthetic_model(random):
    rp, rm = random.randint(1, 5), random.randint(1, 5)
    graph_rank = random.randint(1, min(rp, rm))
    plus_saturation, minus_saturation = matrix(ZZ, rp + rm, rp), matrix(ZZ, rp + rm, rm)
    for i in range(rp):
        plus_saturation[i, i] = 2 if i < graph_rank else 1
        if i < graph_rank:
            plus_saturation[rp + i, i] = -1
    for i in range(rm):
        minus_saturation[rp + i, i] = 1
    plus_exponents = [random.randint(0, 4) for _ in range(rp)]
    minus_exponents = [random.randint(0, 4) for _ in range(rm)]
    plus = (plus_saturation * random_unimodular(rp, random)
            * diagonal_matrix(ZZ, [2 ** exponent for exponent in plus_exponents]))
    minus = (minus_saturation * random_unimodular(rm, random)
             * diagonal_matrix(ZZ, [2 ** exponent for exponent in minus_exponents]))
    ambient_change = random_unimodular(rp + rm, random)
    return ambient_change * plus, ambient_change * minus, graph_rank


def flag_stabilizer_lifts(exponents):
    size = len(exponents)
    generators = []
    for row in range(size):
        for column in range(size):
            if row != column and exponents[row] <= exponents[column]:
                generator = identity_matrix(ZZ, size)
                generator[row, column] = 1
                generators.append(generator)
    identity = identity_matrix(ZZ, size)
    lifts = {tuple(identity.change_ring(F2).list()): identity}
    pending = list(lifts)
    while pending:
        key = pending.pop()
        lift = lifts[key]
        for generator in generators:
            candidate = lift * generator
            reduced = tuple(candidate.change_ring(F2).list())
            if reduced not in lifts:
                lifts[reduced] = candidate
                pending.append(reduced)
    expected = set()
    for entries in product((0, 1), repeat=size * size):
        candidate = matrix(F2, size, size, entries)
        if candidate.rank() == size and all(not candidate[row, column]
                for row in range(size) for column in range(size)
                if exponents[row] > exponents[column]):
            expected.add(tuple(candidate.list()))
    diagonal = diagonal_matrix(QQ, [2 ** exponent for exponent in exponents])
    for lift in lifts.values():
        preserved = diagonal.inverse() * lift * diagonal
        if abs(lift.det()) != 1 or any(value.denominator() != 1 for value in preserved.list()):
            raise ArithmeticError("an enumerated flag stabilizer did not lift integrally")
    return {"sector_exponents": list(exponents), "mod_two_stabilizer_order": len(expected),
            "integral_lifts_constructed": len(lifts), "surjective": set(lifts) == expected}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic-trials", type=int, default=200)
    parser.add_argument("--seed", type=int, default=271828)
    parser.add_argument("--seconds-per-degree", type=float, default=45)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("orthogonal_gluing_flags_results.json"))
    args = parser.parse_args(argv)
    if args.synthetic_trials < 0 or args.seconds_per_degree <= 0:
        parser.error("trial count must be nonnegative and deadline positive")
    started = perf_counter()
    source = Path(__file__).with_name("orthogonal_sector_search_results.json")
    recorded = json.loads(source.read_text())
    report = {"schema_version": 1, "sage_version": SAGE_VERSION,
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Exact finite tests of the abstract two-flag Smith formula, not a proof.",
              "recorded_sector_source": source.name,
              "recorded_sector_sha256": sha256(source.read_bytes()).hexdigest(),
              "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "random_seed": args.seed,
              "synthetic_seed": args.seed,
              "synthetic_parameters": {"sector_rank_range": [1, 5], "sector_exponent_range": [0, 4],
                                       "gluing_dimension_rule": "uniform from 1 through min(sector ranks)",
                                       "basis_changes": "independent random integral elementary shears in both sectors and the ambient lattice"},
              "actual_D_cases": [], "recorded_trivial_gluing_checks": [],
              "synthetic_models": [], "flag_stabilizer_lifts": []}
    def persist():
        report["total_seconds"] = round(perf_counter() - started, 6)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    for case in recorded["cases"]:
        n, k = case["n"], case["k"]
        nontrivial = [row for row in case["degrees"] if row["status"] == "computed"
                      and row["integral_eigenspace_gluing_index"] != "1"]
        for row in case["degrees"]:
            if row["status"] == "computed" and row["integral_eigenspace_gluing_index"] == "1":
                predicted = sorted(row["square"]["power_of_two_exponents"]
                                   + row["euler"]["power_of_two_exponents"])
                report["recorded_trivial_gluing_checks"].append(
                    {"n": n, "k": k, "polynomial_degree": row["polynomial_degree"],
                     "matches": predicted == row["combined"]["power_of_two_exponents"]})
        if not nontrivial:
            continue
        borel = classical_borel("D", n, k)
        square = GeneratorAlgebra(borel, {"q{}".format(i): borel.tail_polynomial(k, i)
                                        for i in range(1, k + 1)})
        p = borel.euler_polynomial(k)
        for row in nontrivial:
            degree = row["polynomial_degree"]
            with deadline(args.seconds_per_degree):
                plus = tuple(square.evaluate(monomial) for monomial in square.monomials(degree))
                minus = tuple(borel.expand(p * square.substitute(monomial))
                              for monomial in square.monomials(degree - k))
                indices = tuple(sorted({word for expansion in plus + minus for word in expansion}, key=borel._word))
                result = analyze(image_matrix(borel, plus, indices)[1],
                                 image_matrix(borel, minus, indices)[1], include_graph=True)
            result.update({"n": n, "k": k, "polynomial_degree": degree})
            report["actual_D_cases"].append(result)
            print("D{} k={} d={} pairs{} predicted{} actual{} matches{}".format(
                n, k, degree, result["mixed_flag_pairs"], result["predicted_combined_exponents"],
                result["actual_combined_exponents"], result["matches"]), flush=True)
            persist()
    random = Random(args.seed)
    for trial in range(args.synthetic_trials):
        plus, minus, expected_graph_rank = synthetic_model(random)
        result = analyze(plus, minus)
        if result["gluing_dimension"] != expected_graph_rank:
            raise ArithmeticError("synthetic graph dimension changed unexpectedly")
        result["trial"] = trial
        report["synthetic_models"].append(result)
        if not result["matches"]:
            persist()
            raise ArithmeticError("two-flag formula failed on a synthetic model")
    for exponents in ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 2)):
        report["flag_stabilizer_lifts"].append(flag_stabilizer_lifts(exponents))
    report["summary"] = {
        "actual_D_nontrivial_gluing_cases": len(report["actual_D_cases"]),
        "actual_D_formula_mismatches": sum(not row["matches"] for row in report["actual_D_cases"]),
        "recorded_trivial_gluing_cases": len(report["recorded_trivial_gluing_checks"]),
        "recorded_trivial_gluing_mismatches": sum(not row["matches"] for row in report["recorded_trivial_gluing_checks"]),
        "synthetic_models": len(report["synthetic_models"]),
        "synthetic_formula_mismatches": sum(not row["matches"] for row in report["synthetic_models"]),
        "synthetic_maximum_gluing_dimension": max((row["gluing_dimension"] for row in report["synthetic_models"]), default=0),
        "synthetic_models_with_gluing_dimension_at_least_two": sum(row["gluing_dimension"] >= 2 for row in report["synthetic_models"]),
        "synthetic_models_with_both_paired_exponents_positive": sum(
            any(pair["plus_exponent"] > 0 and pair["minus_exponent"] > 0 for pair in row["mixed_flag_pairs"])
            for row in report["synthetic_models"]),
        "synthetic_models_with_multiple_pair_types": sum(len(row["mixed_flag_pairs"]) > 1 for row in report["synthetic_models"]),
        "flag_stabilizer_patterns_checked": len(report["flag_stabilizer_lifts"]),
        "flag_stabilizer_integral_lifts_checked": sum(row["integral_lifts_constructed"] for row in report["flag_stabilizer_lifts"]),
        "flag_stabilizer_surjectivity_failures": sum(not row["surjective"] for row in report["flag_stabilizer_lifts"]),
    }
    persist()
    failures = ("actual_D_formula_mismatches", "recorded_trivial_gluing_mismatches",
                "synthetic_formula_mismatches", "flag_stabilizer_surjectivity_failures")
    if any(report["summary"][name] for name in failures):
        raise ArithmeticError("gluing verification failed; inspect the saved report")
    print("Saved {} in {:.3f}s".format(report["summary"], report["total_seconds"]), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
