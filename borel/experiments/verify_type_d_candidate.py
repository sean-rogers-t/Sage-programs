"""Compare the proved full type-D formula with exact recorded Smith lists.

Run from the repository root with ordinary Python 3:
    python borel/experiments/verify_type_d_candidate.py

No Sage, Borel engine, or earlier verifier is imported. The arbitrary-rank proof is in
TYPE_D_EVEN_GLUING_THEOREM.md; finite agreement is not its proof.
"""

from collections import Counter
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import platform
from time import perf_counter


@lru_cache(maxsize=None)
def gaussian(n, k):
    if k < 0 or k > n:
        return ()
    if k == 0 or k == n:
        return (1,)
    left, right = gaussian(n - 1, k), gaussian(n - 1, k - 1)
    result = [0] * max(len(left), len(right) + n - k)
    for i, value in enumerate(left):
        result[i] += value
    for i, value in enumerate(right):
        result[i + n - k] += value
    return tuple(result)


def multiply(left, right):
    if not left or not right:
        return ()
    result = [0] * (len(left) + len(right) - 1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return tuple(result)


@lru_cache(maxsize=None)
def unit_polynomial(n, length):
    assert 0 <= 2 * length <= n
    r, c = n // 2, (n + 1) // 2
    result = [0] * (length * (n - length) + 1)
    for j in range(length // 2 + 1):
        product = multiply(gaussian(r, j), gaussian(c - j, length - 2 * j))
        shift = j * (c - length + j)
        for i, value in enumerate(product):
            result[i + shift] += value
    return tuple(result)


def type_b_counts(a, k, degree):
    """The proved B formula; keys are 2-adic exponents, not factors."""
    if degree < 0 or degree % 2:
        return Counter()
    if a == 0:
        return Counter({0: 1}) if degree == 0 else Counter()
    n, length, m = a + k, min(a, k), degree // 2
    r, c = n // 2, (n + 1) // 2
    result = Counter()
    for exponent in range(length + 1):
        product = multiply(gaussian(r, exponent),
                           unit_polynomial(n - 2 * exponent, length - exponent))
        index = m - exponent * c
        if 0 <= index < len(product) and product[index]:
            result[exponent] = product[index]
    return result


def unit_intersection(a, k, degree):
    if k % 2 or degree % 2:
        return 0
    h, u, v = k // 2, (a + 1) // 2, a // 2
    product = gaussian(h + u - 1, h)
    index = degree // 2 - h * (v + 1)
    return product[index] if 0 <= index < len(product) else 0


def full_candidate(a, k, degree):
    result = type_b_counts(a, k, degree) + type_b_counts(a - 1, k, degree - k)
    intersection = unit_intersection(a, k, degree)
    result[0] -= intersection
    result[1] += intersection
    assert all(value >= 0 for value in result.values()), (a, k, degree, result)
    return +result


def expanded(counts):
    return sorted(exponent for exponent, count in counts.items() for _ in range(count))


def main():
    started = perf_counter()
    base = Path(__file__).resolve().parent
    paths = (base / "orthogonal_sector_search_results.json",
             base / "type_d_frobenius_search_results.json")
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "scope": "Finite verification of the full type-D formula proved in TYPE_D_EVEN_GLUING_THEOREM.md; no theorem inferred from finite data. Historical input timeouts remain excluded here and are superseded separately by the new even-gluing report.",
        "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_hashes": [{"file": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}
                         for path in paths],
        "degree_convention": "polynomial degree d is cohomological degree 2d",
        "checks": [], "skipped_timeouts": [],
    }

    def compare(a, k, degree, actual, source, plus=None, minus=None, gluing=None):
        predicted = expanded(full_candidate(a, k, degree))
        assert predicted == actual, (a, k, degree, predicted, actual)
        row = {"a": a, "k": k, "n": a + k, "polynomial_degree": degree,
               "source": source, "actual_exponents": actual,
               "candidate_exponents": predicted, "all_full_factors_match": True,
               "proved_unit_intersection": unit_intersection(a, k, degree)}
        if plus is not None:
            row["square_matches_full_B_formula"] = plus == expanded(type_b_counts(a, k, degree))
            assert row["square_matches_full_B_formula"], row
        if minus is not None:
            row["euler_matches_shifted_B_formula"] = minus == expanded(type_b_counts(a - 1, k, degree - k))
            assert row["euler_matches_shifted_B_formula"], row
        if gluing is not None:
            row["full_gluing_is_unit_gluing_in_this_case"] = gluing == row["proved_unit_intersection"]
            assert row["full_gluing_is_unit_gluing_in_this_case"], row
        report["checks"].append(row)

    historical = json.loads(paths[0].read_text())
    for case in historical["cases"]:
        for row in case["degrees"]:
            if row["status"] != "computed":
                continue
            gluing_index = int(row["integral_eigenspace_gluing_index"])
            assert gluing_index > 0 and gluing_index & (gluing_index - 1) == 0
            compare(case["a"], case["k"], row["polynomial_degree"],
                    row["combined"]["power_of_two_exponents"], paths[0].name,
                    row["square"]["power_of_two_exponents"],
                    row["euler"]["power_of_two_exponents"], gluing_index.bit_length() - 1)
    fresh = json.loads(paths[1].read_text())
    for row in fresh["fresh_D_checks"]:
        if row["status"] != "computed":
            report["skipped_timeouts"].append({key: row[key] for key in
                                               ("n", "a", "k", "polynomial_degree", "status")})
            continue
        compare(row["a"], row["k"], row["polynomial_degree"],
                row["actual_combined_exponents"], paths[1].name,
                row["plus_exponents"], row["minus_exponents"], row["gluing_dimension"])
    report["summary"] = {"full_smith_lists_compared": len(report["checks"]),
                         "timeouts_excluded": len(report["skipped_timeouts"]),
                         "counterexamples": 0}
    report["status"] = "passed"
    report["seconds"] = round(perf_counter() - started, 6)
    destination = base / "type_d_candidate_verification.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], **report["summary"],
                      "seconds": report["seconds"]}))


if __name__ == "__main__":
    main()
