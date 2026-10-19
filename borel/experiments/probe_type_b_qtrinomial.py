"""A finite q-trinomial probe; the general proof is TYPE_B_GENERAL_SMITH.md."""

import json
from hashlib import sha256
from pathlib import Path
import sys
from time import perf_counter

from verify_type_b_three_row_units import R, gaussian, q


def units(n, ell):
    if ell < 0 or ell > n:
        return R.zero()
    c, r = (n + 1) // 2, n // 2
    return sum(((-1)**j*q**(j*c-j*(j-1)//2)*gaussian(r, j)
                *gaussian(n-2*j, ell-j)
                for j in range(min(ell, r) + 1)), R.zero())


def layer(n, ell, exponent):
    return (q**(exponent*((n+1)//2))*gaussian(n//2, exponent)
            *units(n-2*exponent, ell-exponent))


def main():
    if "--unit-frontier" in sys.argv:
        from verify_type_b_filtration import factors, schur_product_matrix
        reports = []
        started = perf_counter()
        for width, ell in ((6, 6), (8, 4)):
            n = width + ell
            predicted = units(n, ell)
            checked = 0
            degrees = []
            for m in range(width*ell+1):
                coefficients, _ = schur_product_matrix(width, ell, m, 0)
                smith = factors(coefficients)
                assert all(value == 1 for value in smith), (width, ell, m, smith)
                assert len(smith) == predicted[m], (width, ell, m, len(smith), predicted[m])
                checked += 1
                degrees.append({"half_polynomial_degree": m,
                                "rank": len(smith),
                                "smith_factors": [str(value) for value in smith]})
            report = {"A": width, "ell": ell, "degrees": checked,
                      "unit_polynomial": str(predicted), "all_images_primitive": True,
                      "degree_reports": degrees}
            reports.append(report)
            print(json.dumps({key: value for key, value in report.items()
                              if key != "degree_reports"}), flush=True)
        here = Path(__file__).resolve().parent
        result = {"unit_frontier_result": "passed", "cases": reports,
                  "scope": "Independent exact ZZ Littlewood--Richardson P0 Smith forms in all degrees; no Borel or Weyl engine. The general multiplicity formula is derived separately in TYPE_B_GENERAL_SMITH.md; this computation is finite corroboration.",
                  "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
                  "helper_sha256": sha256(here.joinpath("verify_type_b_filtration.py").read_bytes()).hexdigest(),
                  "seconds": perf_counter()-started}
        output = here / "type_b_unit_frontier_verification.json"
        output.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps({"unit_frontier_result": "passed", "cases": len(reports),
                          "degrees": sum(row["degrees"] for row in reports),
                          "report": str(output), "seconds": result["seconds"]}), flush=True)
        return
    data = json.loads(Path(__file__).with_name(
        "type_b_multiplicity_search_results.json").read_text())
    checked, bad = 0, []
    for case in data["cases"]:
        n, ell = int(case["n"]), min(int(case["a"]), int(case["k"]))
        for row in case["degrees"]:
            if row["status"] != "computed":
                continue
            m = int(row["half_polynomial_degree"])
            exponents = [int(factor).bit_length()-1
                         for factor in row["smith_factors"]]
            predicted = [int(layer(n, ell, e)[m]) for e in range(ell+1)]
            actual = [exponents.count(e) for e in range(ell+1)]
            checked += 1
            if predicted != actual:
                bad.append({"n": n, "ell": ell, "m": m,
                            "predicted": predicted, "actual": actual})
    print(json.dumps({"degrees_checked": checked, "mismatches": bad}))
    for e in range(5):
        print("n8,ell4,e%d: %s" % (e, layer(8, 4, e)))


if __name__ == "__main__":
    main()
