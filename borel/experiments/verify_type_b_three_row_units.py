"""Verify the closed three-row type-B unit formula without the Borel engine.

From borel/: sage -python experiments/verify_type_b_three_row_units.py

The independent oracle expands Schur products s_(r,s)*h_t by Sage's
Littlewood--Richardson arithmetic over GF(2).  A second oracle directly
squares Schur classes in finite ordinary Grassmannian cohomology.  Neither
oracle imports sage_borel, Weyl groups, or divided differences.  Available
recorded type-B Smith data are read for additional corroboration only.
"""

from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import GF, Partitions, PolynomialRing, SymmetricFunctions, ZZ, matrix
from sage.env import SAGE_VERSION


R = PolynomialRing(ZZ, "q")
q = R.gen()


@lru_cache(None)
def gaussian(n, k):
    """Gaussian polynomial via its usual Pascal recurrence."""
    if k < 0 or k > n:
        return R.zero()
    if not k or k == n:
        return R.one()
    return gaussian(n - 1, k - 1) + q**k * gaussian(n - 1, k)


def unit_polynomial(width):
    b, parity = divmod(width, 2)
    if parity:
        return gaussian(b + 3, 3) + q**(b + 1)*gaussian(b + 2, 2)
    return (gaussian(b + 3, 3)
            + q**(b + 1)*(gaussian(b + 2, 2) - q**(2*b))
            - q**(2*b + 1)*sum((q**j for j in range(b)), R.zero()))


def lr_product_rank(width, degree, schur):
    """Rank of the two-quotient Schur-product matrix, using LR arithmetic."""
    b = width // 2
    c = width - b
    products = []
    for r in range(b + 1):
        for s in range(r + 1):
            t = degree - r - s
            if 0 <= t <= c:
                alpha = tuple(part for part in (r, s) if part)
                beta = (t,) if t else ()
                products.append(schur(alpha)*schur(beta))
    indices = tuple(Partitions(degree, max_part=width, max_length=3))
    allowed = set(tuple(partition) for partition in indices)
    assert all(tuple(partition) in allowed
               for product in products for partition in product.support())
    coefficients = matrix(GF(2), len(indices), len(products),
                          lambda row, column: products[column][indices[row]])
    return int(coefficients.rank()), len(products)


def direct_frobenius_rank(width, degree, schur):
    """Independent squaring in the finite ordinary Grassmannian Schur basis."""
    source = tuple(Partitions(degree, max_part=width, max_length=3))
    target = tuple(Partitions(2*degree, max_part=width, max_length=3))
    if not target:
        return 0
    products = tuple(schur(partition)**2 for partition in source)
    coefficients = matrix(GF(2), len(target), len(source),
                          lambda row, column: products[column][target[row]])
    return int(coefficients.rank())


def main():
    started = perf_counter()
    schur = SymmetricFunctions(GF(2)).schur()
    counts = {"closed_formula_vs_LR_product_ranks": 0,
              "closed_formula_vs_direct_frobenius": 0,
              "closed_formula_vs_recorded_B_Smith_lists": 0}
    cases = []
    oracle = {}
    for width in range(3, 13):
        predicted = unit_polynomial(width)
        assert all(value >= 0 for value in predicted.list())
        rows = []
        for m in range(3*width + 1):
            expected = int(predicted[m])
            rank, columns = lr_product_rank(width, m, schur)
            assert expected == rank, ("LR", width, m, expected, rank)
            counts["closed_formula_vs_LR_product_ranks"] += 1
            if width <= 6:
                direct = direct_frobenius_rank(width, m, schur)
                assert direct == expected, ("Frobenius", width, m, expected, direct)
                counts["closed_formula_vs_direct_frobenius"] += 1
            oracle[width, m] = expected
            rows.append({"half_polynomial_degree": m,
                         "polynomial_degree": 2*m,
                         "unit_count": expected,
                         "LR_product_columns": columns})
        cases.append({"block_sizes": [width, 3],
                      "unit_generating_polynomial": str(predicted),
                      "degrees": rows})

    inputs = []
    recorded = []
    here = Path(__file__).resolve().parent
    for filename in ("orthogonal_smith_search_results.json",
                     "type_b_multiplicity_search_results.json"):
        path = here / filename
        if not path.exists():
            continue
        raw = path.read_bytes()
        data = json.loads(raw)
        inputs.append({"filename": filename, "sha256": sha256(raw).hexdigest(),
                       "bytes": len(raw)})
        for case in data["cases"]:
            if case.get("family", "B") != "B":
                continue
            a, k = int(case["a"]), int(case["k"])
            if 3 not in (a, k):
                continue
            width = k if a == 3 else a
            predicted = unit_polynomial(width)
            for row in case["degrees"]:
                if row["status"] != "computed":
                    continue
                d = int(row["polynomial_degree"])
                assert d % 2 == 0
                observed = sum(int(factor) == 1 for factor in row["smith_factors"])
                expected = int(predicted[d // 2])
                assert expected == observed, (filename, a, k, d, expected, observed)
                counts["closed_formula_vs_recorded_B_Smith_lists"] += 1
                recorded.append({"input_file": filename, "block_sizes": [a, k],
                                 "polynomial_degree": d, "unit_count": observed})

    report = {
        "suite_result": "passed",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sage_version": SAGE_VERSION,
        "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "oracle_scope": {
            "other_block_widths": list(range(3, 13)),
            "polynomial_degrees": "all even degrees 0 through 6*A for each width A",
            "direct_Frobenius_widths": [3, 4, 5, 6],
            "expected_values": "Schur Littlewood--Richardson arithmetic over GF(2)",
            "Borel_or_Weyl_engine_imported": False,
        },
        "checks": counts,
        "input_datasets": inputs,
        "cases": cases,
        "recorded_Smith_comparisons": recorded,
        "scope": "Closed counts of unit factors only. Higher 2-primary multiplicities are not classified by this verification.",
        "seconds": round(perf_counter() - started, 6),
    }
    output = here / "type_b_three_row_unit_verification.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"suite_result": "passed", "checks": counts,
                      "seconds": report["seconds"], "report": str(output)}))


if __name__ == "__main__":
    main()
