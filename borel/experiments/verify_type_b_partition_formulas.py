"""Independent checks of closed type-B partition formulas.

No sage_borel, Weyl group, or divided-difference code is imported.  The
three-row multiplicities are compared with recorded Smith data; its middle
prefix is independently computed by integral Schur-LR products.  The highest
factor is checked against Gaussian coefficients, explicit partition counts,
and primitive coordinate row modules from the one-line Pieri formula.
"""

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter

from sage.all import Partition, ZZ, matrix
from sage.env import SAGE_VERSION

from verify_type_b_filtration import box_by_degree, box_degree, factors, schur_product
from verify_type_b_three_row_units import gaussian, q, unit_polynomial


def rectangle(width, height, degree):
    return box_degree(box_by_degree(width, height), degree)


def highest_count(width, height, degree):
    threshold = (width + height - 1) // 2
    return sum(len(partition) == height and partition[-1] > threshold
               for partition in rectangle(width, height, degree))


def three_row_counts(width, degree):
    partitions = rectangle(width, 3, degree)
    upper, lower = (width + 2) // 2, (width + 1) // 2
    r0 = int(unit_polynomial(width)[degree])
    r1 = 0
    for partition in partitions:
        padded = partition + (0,) * (3 - len(partition))
        r1 += padded[1] <= upper and padded[2] <= lower
    r3 = len(partitions)
    r2 = r3 - highest_count(width, 3, degree)
    multiplicities = (r0, r1 - r0, r2 - r1, r3 - r2)
    if any(value < 0 for value in multiplicities):
        raise ArithmeticError("negative three-row multiplicity")
    return (r0, r1, r2, r3), tuple(map(int, multiplicities))


def three_row_middle_lr(width, degree):
    upper, lower = (width + 2) // 2, (width + 1) // 2
    source = rectangle(width, 3, degree)
    factorizations, raw_columns = set(), 0
    for first in range(width + 1):
        for second in range(upper + 1):
            third = degree - first - second
            if 0 <= third <= lower:
                raw_columns += 1
                factorizations.add(tuple(sorted((value,) for value in (first, second, third) if value)))
    rows = []
    for factorization in sorted(factorizations):
        product = schur_product(factorization)
        rows.append(tuple(ZZ(product[Partition(partition)]) for partition in source))
    coefficients = matrix(ZZ, len(rows), len(source), [value for row in rows for value in row])
    smith = factors(coefficients)
    expected = three_row_counts(width, degree)[0][1]
    assert len(smith) == expected, ("middle LR rank", width, degree, expected, smith)
    assert all(value == 1 for value in smith), ("middle LR primitivity", width, degree, smith)
    return {"rank": len(smith), "integer_smith_factors": [str(value) for value in smith],
            "raw_triples": raw_columns, "distinct_products": len(factorizations)}


def one_line_pieri_coordinate_oracle(width, height, degree):
    """Horizontal-strip rows, independently specified by interlacing parts."""
    threshold = (width + height - 1) // 2
    source = rectangle(width, height, degree)
    rows = []
    for strip in range(min(threshold, degree) + 1):
        for mu in rectangle(width, height - 1, degree - strip):
            mu_padded = mu + (0,) * (height - 1 - len(mu))
            row = []
            for partition in source:
                padded = partition + (0,) * (height - len(partition))
                row.append(int(all(padded[i] >= mu_padded[i] >= padded[i + 1]
                                   for i in range(height - 1))))
            rows.append(tuple(row))
    coefficients = matrix(ZZ, len(rows), len(source), [value for row in rows for value in row])
    allowed = [column for column, partition in enumerate(source)
               if len(partition) < height or partition[-1] <= threshold]
    coordinates = matrix(ZZ, len(allowed), len(source),
                         {(row, column): 1 for row, column in enumerate(allowed)})
    assert coefficients.row_module() == coordinates.row_module(), (
        "Pieri coordinate row module", width, height, degree)
    smith = factors(coefficients)
    assert len(smith) == len(allowed) and all(value == 1 for value in smith)
    return {"rank": len(smith), "source_dimension": len(source),
            "Pieri_rows": len(rows), "all_integer_smith_factors_one": True,
            "exact_coordinate_row_module_matches": True}


def main():
    started = perf_counter()
    here = Path(__file__).resolve().parent
    report = {"suite_result": "running", "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "Borel_or_Weyl_engine_imported": False,
              "scope": "Closed normalized B multiplicities for min-block three and highest-factor partition/Pieri checks.",
              "checks": {}, "three_row_cases": [], "highest_factor_cases": [],
              "Pieri_coordinate_cases": [], "recorded_comparisons": [], "input_datasets": [],
              "source_files": []}
    for filename in (Path(__file__).name, "verify_type_b_filtration.py", "verify_type_b_three_row_units.py"):
        path = here / filename
        raw = path.read_bytes()
        report["source_files"].append({"filename": filename, "sha256": sha256(raw).hexdigest(), "bytes": len(raw)})
    counts = {"three_row_degrees_with_nonnegative_closed_multiplicities": 0,
              "three_row_middle_LR_rank_and_primitivity": 0,
              "highest_partition_vs_Gaussian_coefficients": 0,
              "one_line_Pieri_exact_coordinate_row_modules": 0,
              "closed_three_row_vs_recorded_full_Smith_lists": 0,
              "closed_three_row_middle_vs_recorded_C_prefix": 0,
              "highest_factor_vs_recorded_Smith_lists": 0}
    for width in range(3, 13):
        rows = []
        for degree in range(3 * width + 1):
            prefix_ranks, multiplicities = three_row_counts(width, degree)
            assert sum(multiplicities) == int(gaussian(width + 3, 3)[degree])
            counts["three_row_degrees_with_nonnegative_closed_multiplicities"] += 1
            row = {"half_polynomial_degree": degree, "prefix_ranks": list(prefix_ranks),
                   "multiplicities_of_1_2_4_8": list(multiplicities)}
            if width <= 8:
                row["middle_LR_oracle"] = three_row_middle_lr(width, degree)
                counts["three_row_middle_LR_rank_and_primitivity"] += 1
            rows.append(row)
        generating = [sum((row["multiplicities_of_1_2_4_8"][exponent]
                           * q ** row["half_polynomial_degree"] for row in rows), q.parent().zero())
                      for exponent in range(4)]
        report["three_row_cases"].append({"A": width,
                                          "generating_polynomials_for_1_2_4_8": list(map(str, generating)),
                                          "degrees": rows})
    for height in range(1, 7):
        for width in range(height, 11):
            threshold = (width + height - 1) // 2
            predicted = q ** (height * (threshold + 1)) * gaussian(width - threshold - 1 + height, height)
            rows = []
            for degree in range(width * height + 1):
                counted = highest_count(width, height, degree)
                assert counted == int(predicted[degree]), ("highest Gaussian", width, height, degree)
                counts["highest_partition_vs_Gaussian_coefficients"] += 1
                rows.append({"half_polynomial_degree": degree, "highest_factor_multiplicity": counted})
            report["highest_factor_cases"].append({"A": width, "ell": height,
                                                   "threshold": threshold,
                                                   "generating_polynomial": str(predicted), "degrees": rows})
            if width <= 6:
                oracle_rows = []
                for degree in range(width * height + 1):
                    oracle_rows.append({"half_polynomial_degree": degree,
                                        **one_line_pieri_coordinate_oracle(width, height, degree)})
                    counts["one_line_Pieri_exact_coordinate_row_modules"] += 1
                report["Pieri_coordinate_cases"].append({"A": width, "ell": height, "degrees": oracle_rows})
    for filename in ("orthogonal_smith_search_results.json", "type_b_multiplicity_search_results.json"):
        path = here / filename
        raw = path.read_bytes()
        data = json.loads(raw)
        report["input_datasets"].append({"filename": filename, "sha256": sha256(raw).hexdigest(), "bytes": len(raw)})
        for case in data["cases"]:
            if case.get("family", "B") != "B":
                continue
            height, width = min(case["a"], case["k"]), max(case["a"], case["k"])
            for row in case["degrees"]:
                if row["status"] != "computed":
                    continue
                d = int(row["polynomial_degree"])
                assert d % 2 == 0
                degree = d // 2
                actual = [int(value) for value in row["smith_factors"]]
                highest = sum(value == 2 ** height for value in actual)
                assert highest == highest_count(width, height, degree)
                counts["highest_factor_vs_recorded_Smith_lists"] += 1
                comparison = {"input_file": filename, "A": width, "ell": height,
                              "half_polynomial_degree": degree, "highest_multiplicity": highest}
                if height == 3:
                    prefix, multiplicities = three_row_counts(width, degree)
                    expected = [2 ** exponent for exponent, count in enumerate(multiplicities)
                                for _ in range(count)]
                    assert actual == expected, ("full three-row Smith", filename, width, degree, expected, actual)
                    counts["closed_three_row_vs_recorded_full_Smith_lists"] += 1
                    comparison["full_Smith_list_matches"] = True
                    if "filtered_C_ranks" in row:
                        available = [entry for entry in row["filtered_C_ranks"] if entry["threshold"] <= 1]
                        observed_prefix = available[-1]["rank_mod_two"] if available else 0
                        assert observed_prefix == prefix[1]
                        counts["closed_three_row_middle_vs_recorded_C_prefix"] += 1
                report["recorded_comparisons"].append(comparison)
    assert not any(name == "sage_borel" or name.startswith("sage_borel.") for name in sys.modules)
    report["suite_result"] = "passed"
    report["checks"] = counts
    report["seconds"] = round(perf_counter() - started, 6)
    output = here / "type_b_partition_formula_verification.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"suite_result": "passed", "checks": counts, "seconds": report["seconds"]}))


if __name__ == "__main__":
    main()
