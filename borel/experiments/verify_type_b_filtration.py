"""Independent Schur-LR verification of the proposed type-B filtration.

No Borel calculator, Weyl group, or divided difference is imported.  Products
on the right side of the incidence/2-quotient formula are expanded over ZZ in
the ordinary Grassmannian Schur basis.  Their exact integer row modules are
compared, by Hermite/module equality, with the recorded symplectic prefixes.
"""

from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import Partition, Partitions, SymmetricFunctions, ZZ, matrix
from sage.env import SAGE_VERSION


SCHUR = SymmetricFunctions(ZZ).schur()


@lru_cache(None)
def box_by_degree(width, height):
    if not width or not height:
        return (((),),)
    return tuple(tuple(tuple(partition) for partition in Partitions(degree, max_part=width,
                                                                   max_length=height))
                 for degree in range(width * height + 1))


def box_degree(box, degree):
    return box[degree] if 0 <= degree < len(box) else ()


@lru_cache(None)
def schur_product(factors):
    product = SCHUR.one()
    for factor in factors:
        if factor:
            product *= SCHUR(factor)
    return product


def schur_product_matrix(width, short, degree, threshold, source=None):
    source = (tuple(tuple(partition) for partition in Partitions(degree, max_part=width,
                                                                 max_length=short))
              if source is None else tuple(source))
    residual = short - threshold
    first, second = (residual + 1) // 2, residual // 2
    first_width = (width + short + 1) // 2 - first
    second_width = (width + short) // 2 - second
    mu_box = box_by_degree(width, threshold)
    alpha_box = box_by_degree(first_width, first)
    beta_box = box_by_degree(second_width, second)
    factors, raw_count = set(), 0
    for mu_degree in range(degree + 1):
        for alpha_degree in range(degree - mu_degree + 1):
            beta_degree = degree - mu_degree - alpha_degree
            for mu in box_degree(mu_box, mu_degree):
                for alpha in box_degree(alpha_box, alpha_degree):
                    for beta in box_degree(beta_box, beta_degree):
                        raw_count += 1
                        factors.add(tuple(sorted(factor for factor in (mu, alpha, beta) if factor)))
    rows = []
    for factorization in sorted(factors):
        value = schur_product(factorization)
        row = tuple(ZZ(value[Partition(partition)]) for partition in source)
        if any(row):
            rows.append(row)
    coefficients = matrix(ZZ, len(rows), len(source),
                          [value for row in rows for value in row])
    return coefficients, {"source_dimension": len(source), "raw_product_count": raw_count,
                          "distinct_factorizations": len(factors),
                          "nonzero_product_rows": len(rows),
                          "mu_box": [width, threshold],
                          "alpha_box": [first_width, first],
                          "beta_box": [second_width, second]}


def factors(coefficients):
    if not all(coefficients.dimensions()):
        return ()
    diagonal = coefficients.smith_form()[0]
    return tuple(abs(ZZ(diagonal[i, i])) for i in range(min(diagonal.dimensions())) if diagonal[i, i])


def decode_matrix(encoded):
    return matrix(ZZ, encoded["rows"], encoded["columns"],
                  {(entry["row"], entry["column"]): ZZ(entry["value"]) for entry in encoded["entries"]})


def main():
    started = perf_counter()
    here = Path(__file__).resolve().parent
    path = here / "type_b_multiplicity_search_results.json"
    raw = path.read_bytes()
    dataset = json.loads(raw)
    report = {"suite_result": "running", "generated_utc": datetime.now(timezone.utc).isoformat(),
              "sage_version": SAGE_VERSION, "Borel_or_Weyl_engine_imported": False,
              "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "input_dataset": {"filename": path.name, "sha256": sha256(raw).hexdigest(), "bytes": len(raw)},
              "scope": "Exact ZZ Schur-product row-module equality against stored C prefixes; no fresh BGG calculation.",
              "orientation": "Use the actual head a and tail k. Replace each stored tail Schur index lambda by its conjugate head index; the common degree sign does not change a row module.",
              "checks": [], "counterexamples": []}
    output = here / "type_b_filtration_verification.json"
    for case in dataset["cases"]:
        head, tail = int(case["a"]), int(case["k"])
        for row in case["degrees"]:
            if row["status"] != "computed":
                continue
            degree = row["half_polynomial_degree"]
            source = tuple(tuple(Partition(partition).conjugate()) for partition in row["partitions"])
            actual = decode_matrix(row["C_matrix"])
            for threshold in range(head + 1):
                product_rows, metadata = schur_product_matrix(tail, head, degree, threshold, source)
                prefix = actual.matrix_from_rows([i for i, weight in enumerate(row["row_negative_counts"])
                                                  if weight <= threshold])
                product_module, actual_module = product_rows.row_module(), prefix.row_module()
                equal = product_module == actual_module
                smith = factors(product_rows)
                check = {"head_size": head, "tail_size": tail, "half_polynomial_degree": degree,
                         "threshold": threshold, "row_modules_equal": equal,
                         "LR_rank": int(product_module.rank()), "BGG_prefix_rank": int(actual_module.rank()),
                         "LR_smith_factors": [str(value) for value in smith],
                         "LR_image_is_primitive": all(value == 1 for value in smith), **metadata}
                report["checks"].append(check)
                if not equal:
                    check["LR_Hermite_basis"] = [list(map(int, values)) for values in product_module.basis_matrix().rows()]
                    check["BGG_Hermite_basis"] = [list(map(int, values)) for values in actual_module.basis_matrix().rows()]
                    report["counterexamples"].append(check)
                    output.write_text(json.dumps(report, indent=2) + "\n")
                    if len(report["counterexamples"]) <= 3:
                        print("MISMATCH {}".format(check), flush=True)
            print("head={} tail={} m={} allprefixeschecked".format(head, tail, degree), flush=True)
    report["suite_result"] = "passed" if not report["counterexamples"] else "failed"
    report["summary"] = {"prefix_row_module_comparisons": len(report["checks"]),
                         "row_module_mismatches": len(report["counterexamples"]),
                         "nonprimitive_LR_images": sum(not row["LR_image_is_primitive"] for row in report["checks"]),
                         "Schur_product_cache_size": schur_product.cache_info().currsize}
    report["seconds"] = round(perf_counter() - started, 6)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"suite_result": report["suite_result"], "summary": report["summary"], "seconds": report["seconds"]}))
    if report["counterexamples"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
