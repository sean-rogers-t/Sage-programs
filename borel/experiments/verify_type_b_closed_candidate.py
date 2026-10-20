"""Finite exact checks of the all-rank type-B multiplicity formula.

The unit-factor formula is an audited repository deduction from published
inputs (novelty unestablished) and is evaluated by a positive Gaussian
sum, independently checked against its alternating expression.  The full
Gaussian expression is now derived in TYPE_B_GENERAL_SMITH.md. These
finite checks are corroboration, not its proof. No Borel engine or Weyl
group is imported.
"""

from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import PolynomialRing, ZZ

from verify_type_b_three_row_units import gaussian, q, R


@lru_cache(None)
def alternating_units(n, short):
    if short < 0 or n < 0 or short > n:
        return R.zero()
    if short == 0:
        return R.one()
    first, second = n // 2, (n + 1) // 2
    return sum(((-1) ** j * q ** (j * second - j * (j - 1) // 2)
                * gaussian(first, j) * gaussian(n - 2 * j, short - j)
                for j in range(min(first, short) + 1)), R.zero())


@lru_cache(None)
def positive_units(n, short):
    if short < 0 or n < 0 or short > n:
        return R.zero()
    if short == 0:
        return R.one()
    first, second = n // 2, (n + 1) // 2
    return sum((q ** (j * (second - short + j)) * gaussian(first, j)
                * gaussian(second - j, short - 2 * j)
                for j in range(min(first, short // 2) + 1)), R.zero())


def candidate_multiplicity(n, short, exponent):
    if not 0 <= exponent <= short:
        return R.zero()
    return (q ** (exponent * ((n + 1) // 2)) * gaussian(n // 2, exponent)
            * positive_units(n - 2 * exponent, short - exponent))


def main():
    started = perf_counter()
    here = Path(__file__).resolve().parent
    report = {"suite_result": "running",
              "status_of_formula": {"unit_factors_e_zero": "audited deduction in TYPE_B_UNIT_SMITH.md from published inputs; novelty not established",
                                    "counting_polynomial": "classical Birkhoff/Hall enumeration; see PROOF_STATUS.md",
                                    "general_intermediate_factors": "complete repository derivation in TYPE_B_GENERAL_SMITH.md; finite verification is not a proof"},
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "Borel_or_Weyl_engine_imported": False,
              "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "inputs": [], "Smith_checks": [], "prefix_checks": [], "counterexamples": []}
    used_unit_parameters = set()
    def note_parameters(n, short):
        used_unit_parameters.update((n - 2 * exponent, short - exponent)
                                    for exponent in range(short + 1))
    dependency = here / "verify_type_b_three_row_units.py"
    report["source_dependency"] = {"filename": dependency.name,
                                   "sha256": sha256(dependency.read_bytes()).hexdigest()}
    for filename in ("orthogonal_smith_search_results.json", "type_b_multiplicity_search_results.json"):
        path = here / filename
        raw = path.read_bytes()
        report["inputs"].append({"filename": filename, "sha256": sha256(raw).hexdigest(), "bytes": len(raw)})
        for case in json.loads(raw)["cases"]:
            if case.get("family", "B") != "B":
                continue
            n, short = case["n"], min(case["a"], case["k"])
            note_parameters(n, short)
            polynomials = [candidate_multiplicity(n, short, exponent) for exponent in range(short + 1)]
            for row in case["degrees"]:
                if row["status"] != "computed":
                    continue
                degree = row["polynomial_degree"] // 2
                counts = [int(polynomial[degree]) for polynomial in polynomials]
                predicted = [2 ** exponent for exponent, count in enumerate(counts) for _ in range(count)]
                actual = [int(value) for value in row["smith_factors"]]
                check = {"filename": filename, "n": n, "ell": short, "m": degree,
                         "matches": predicted == actual and all(count >= 0 for count in counts),
                         "candidate_counts": counts, "actual_Smith": actual}
                report["Smith_checks"].append(check)
                if not check["matches"]:
                    report["counterexamples"].append(check)
                for prefix in row.get("filtered_C_ranks", []):
                    threshold = prefix["threshold"]
                    expected = sum(counts[:threshold + 1])
                    prefix_check = {"filename": filename, "n": n, "ell": short, "m": degree,
                                    "threshold": threshold, "candidate_rank": expected,
                                    "observed_rank": prefix["rank_mod_two"],
                                    "matches": expected == prefix["rank_mod_two"]}
                    report["prefix_checks"].append(prefix_check)
                    if not prefix_check["matches"]:
                        report["counterexamples"].append(prefix_check)
    larger = here / "type_b_schur_prefix_search_results.json"
    if larger.exists():
        raw = larger.read_bytes()
        report["inputs"].append({"filename": larger.name, "sha256": sha256(raw).hexdigest(), "bytes": len(raw)})
        for case in json.loads(raw)["cases"]:
            n, short = case["n"], min(case["head"], case["tail"])
            note_parameters(n, short)
            polynomials = [candidate_multiplicity(n, short, exponent) for exponent in range(short + 1)]
            for row in case["degrees"]:
                degree = row["half_polynomial_degree"]
                for prefix in row["prefixes"]:
                    threshold = prefix["threshold"]
                    expected = sum(int(polynomial[degree]) for polynomial in polynomials[:threshold + 1])
                    check = {"filename": larger.name, "n": n, "ell": short, "m": degree,
                             "threshold": threshold, "candidate_rank": expected,
                             "observed_rank": prefix["rank"], "matches": expected == prefix["rank"]}
                    report["prefix_checks"].append(check)
                    if not check["matches"]:
                        report["counterexamples"].append(check)
    consistency = []
    ring = PolynomialRing(ZZ, names=("x", "t"))
    x, t = ring.gens()
    for n in range(2, 25):
        at_one = (1 + x) ** (n % 2) * (1 + (1 + t) * x + x ** 2) ** (n // 2)
        for short in range(1, min(6, n // 2) + 1):
            note_parameters(n, short)
            polynomials = [candidate_multiplicity(n, short, exponent) for exponent in range(short + 1)]
            nonnegative = all(coefficient >= 0 for polynomial in polynomials for coefficient in polynomial.list())
            total_rank = sum(polynomials, R.zero()) == gaussian(n, short)
            at_one_matches = all(int(polynomial(1)) == at_one.monomial_coefficient(x ** short * t ** exponent)
                                 for exponent, polynomial in enumerate(polynomials))
            check = {"n": n, "ell": short, "nonnegative": nonnegative,
                     "sum_is_Gaussian_rank": total_rank, "q_one_bivariate_identity": at_one_matches}
            consistency.append(check)
            if not (nonnegative and total_rank and at_one_matches):
                report["counterexamples"].append(check)
    report["polynomial_consistency_checks"] = consistency
    report["positive_vs_alternating_unit_identity_checks"] = []
    for n, short in sorted(used_unit_parameters):
        positive, alternating = positive_units(n, short), alternating_units(n, short)
        check = {"n": n, "ell": short, "matches": positive == alternating,
                 "positive_generating_polynomial": str(positive)}
        report["positive_vs_alternating_unit_identity_checks"].append(check)
        if not check["matches"]:
            report["counterexamples"].append(check)
    report["suite_result"] = "passed" if not report["counterexamples"] else "failed"
    report["summary"] = {"Smith_lists_checked": len(report["Smith_checks"]),
                         "filtered_prefixes_checked": len(report["prefix_checks"]),
                         "larger_LR_prefixes_checked": sum(row["filename"] == larger.name for row in report["prefix_checks"]),
                         "parameter_cases_for_polynomial_consistency": len(consistency),
                         "positive_vs_alternating_unit_identities": len(used_unit_parameters),
                         "counterexamples": len(report["counterexamples"])}
    report["seconds"] = round(perf_counter() - started, 6)
    output = here / "type_b_closed_candidate_verification.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"suite_result": report["suite_result"], "summary": report["summary"],
                      "counterexamples": report["counterexamples"][:3], "seconds": report["seconds"]}))
    if report["counterexamples"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
