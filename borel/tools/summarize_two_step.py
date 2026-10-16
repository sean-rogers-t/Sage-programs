"""Render the completed two-step experiment as a compact Markdown summary."""

import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    report = json.loads((root / "experiments/two_step_symplectic_results.json").read_text(encoding="utf-8"))
    if not report["complete"]:
        raise RuntimeError("experiment has not completed")
    cases = report["cases"]
    if not all(case["all_bgg_ranks_match_presentation_through_bound"] for case in cases):
        raise RuntimeError("a computed rank does not match the presentation")
    comparisons = sum(len(case["degrees"]) for case in cases)
    relation_checks = sum(sum(r["bgg_evaluated_within_bound"] for r in case["presentation_relations"])
                          for case in cases)
    lines = [
        "# Two-step symplectic experiment", "",
        "The completed grid uses `b,k` in `{1,2}`, `a=1,...,5`, `n=a+b+k`,",
        "and polynomial degrees `0,...,8` (cohomological degrees `0,...,16`).",
        "The selected image uses the joint `GL(b) x Sp(2k)` generators.", "",
        "All **{} cases and {} homogeneous ranks** agree with the independent".format(len(cases), comparisons),
        "formal presentation, and all **{} presentation relations** within the".format(relation_checks),
        "degree bound have zero BGG image. The separate test suite also compares",
        "whole kernel spaces with coordinate-coinvariant Groebner reductions.", "",
        "| n | b | k | a | First relation found through degree 8 | Degree-8 domain | Degree-8 image |",
        "| ---: | ---: | ---: | ---: | --- | ---: | ---: |",
    ]
    for case in cases:
        parameters = case["parameters"]
        degree8 = next(row for row in case["degrees"] if row["polynomial_degree"] == 8)
        first = case["first_observed_relation_degree"]
        lines.append("| {} | {} | {} | {} | {} | {} | {} |".format(
            parameters["n"], parameters["b"], parameters["k"], parameters["a"],
            first if first is not None else "None through 8",
            degree8["domain_dimension"], degree8["image_dimension"]))
    lines.extend([
        "", "The first computed relation degrees are 4, 6 and 8 for a=1,2,3.",
        "For a=4,5 there is no relation through degree 8. Finite computations",
        "alone establish only these bounded facts. The",
        "[invariant-theory derivation](TWO_STEP_DERIVATION.md) separately proves the",
        "full ideal `(h_(a+1)(Z),...,h_n(Z))` over QQ, and hence the sharp first",
        "relation degree `2(a+1)`. No novelty claim is made.", "",
        "For n=3,b=k=1, a mixed relation is `c1^4+c1^2*e1+e1^2=0`, and",
        "`e1^3=0`. The image Hilbert polynomial has coefficients",
        "`(1,1,2,2,2,2,1,1)`.", "",
        "The [complete exact report](two_step_symplectic_results.json) contains input",
        "polynomials, generator degrees, every computed kernel basis, structured",
        "rational coefficients and timing information. No complete ambient Weyl",
        "group or ambient Schubert basis was enumerated for this grid.", "",
        "Run from `borel/`:", "",
        "```bash", "PYTHONPATH=src sage -python experiments/two_step_symplectic.py",
        "python tools/summarize_two_step.py", "```", "",
        "The recorded SageMath {} run took {:.2f} seconds; timing is machine dependent.".format(
            report["sage_version"], report["total_seconds"]), "",
        "The script accepts degree/rank-grid bounds; edit or use its `--help` options",
        "for a new experiment. This summary requires the recorded degree-8 grid.", "",
    ])
    (root / "experiments/README.md").write_text("\n".join(lines), encoding="utf-8")
    print("Summarized {} cases and {} degrees.".format(len(cases), comparisons))


if __name__ == "__main__":
    main()
