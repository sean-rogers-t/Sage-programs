"""Independent pure-Python corroboration of the general type-B count.

Run from any directory with Python 3:
    python borel/experiments/verify_type_b_general_audit.py

This deliberately imports no Sage, library implementation, or earlier
verification code. It checks corner membership against the beta-number
statistic and compares the full weighted distributions with the positive
Gaussian expression. These are finite checks, not a proof of all ranks.
"""

from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import platform
from time import perf_counter


PARAMETERS = ((14, 7), (15, 7), (16, 8), (17, 8), (18, 9),
              (19, 9), (20, 10), (21, 10), (22, 11))


@lru_cache(maxsize=None)
def gaussian(n, k):
    """Exact coefficients of [n choose k]_q, in increasing degree."""
    if k < 0 or k > n:
        return ()
    if k == 0 or k == n:
        return (1,)
    left = gaussian(n - 1, k)
    right = gaussian(n - 1, k - 1)
    result = [0] * max(len(left), len(right) + n - k)
    for degree, coefficient in enumerate(left):
        result[degree] += coefficient
    for degree, coefficient in enumerate(right):
        result[degree + n - k] += coefficient
    return tuple(result)


def multiply(left, right):
    """Convolve integer coefficient lists without a polynomial package."""
    result = [0] * (len(left) + len(right) - 1) if left and right else []
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i + j] += first * second
    return result


def unit_polynomial(n, length):
    """Positive formula U_(n,length); all used inputs obey 2*length<=n."""
    assert 0 <= 2 * length <= n
    r, c = n // 2, (n + 1) // 2
    result = [0] * (length * (n - length) + 1)
    for j in range(length // 2 + 1):
        product = multiply(gaussian(r, j), gaussian(c - j, length - 2 * j))
        shift = j * (c - length + j)
        for degree, coefficient in enumerate(product):
            result[degree + shift] += coefficient
    return result


def in_corner(gamma, t, r, c, epsilon):
    """Directly test the union of allowed partition corners."""
    b = len(gamma) - t
    return b == 0 or any(
        all(gamma[t + i] <= (c - j if i < j else r - b + j)
            for i in range(b))
        for j in range((b + epsilon) // 2 + 1)
    )


def check_case(n, ell):
    A = n - ell
    r, c, epsilon = n // 2, (n + 1) // 2, n % 2
    counts = [[0] * (A * ell + 1) for _ in range(ell + 1)]
    partition_count = membership_count = 0
    for increasing_beads in combinations(range(n), ell):
        beads = increasing_beads[::-1]
        gamma = tuple(beads[i] - ell + i + 1 for i in range(ell))
        d = sum(value < r for value in beads)
        exponent = sum(value >= c + d for value in beads)
        counts[exponent][sum(gamma)] += 1
        partition_count += 1
        for t in range(ell + 1):
            assert in_corner(gamma, t, r, c, epsilon) == (exponent <= t), (
                n, ell, t, gamma, exponent)
            membership_count += 1

    for exponent in range(ell + 1):
        predicted = [0] * (A * ell + 1)
        product = multiply(gaussian(r, exponent),
                           unit_polynomial(n - 2 * exponent, ell - exponent))
        shift = exponent * c
        for degree, coefficient in enumerate(product):
            predicted[degree + shift] = coefficient
        assert counts[exponent] == predicted, (n, ell, exponent)

    # The total is also checked against the ordinary rectangle Gaussian.
    assert [sum(row[m] for row in counts) for m in range(A * ell + 1)] == list(
        gaussian(n, ell))
    return {
        "n": n, "ell": ell, "A": A,
        "partitions": partition_count,
        "corner_membership_comparisons": membership_count,
        "exact_weighted_polynomial_comparisons": ell + 1,
        "all_exponent_weight_distributions": counts,
        "corner_criterion_passed": True,
        "all_weighted_polynomials_passed": True,
        "total_rectangle_gaussian_passed": True,
    }


def main():
    started = perf_counter()
    script = Path(__file__).resolve()
    proof = script.parent.parent / "TYPE_B_GENERAL_SMITH.md"
    proof_start_hash = sha256(proof.read_bytes()).hexdigest()
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "scope": "finite independent corroboration; no general-proof claim",
        "sage_imported": False,
        "library_or_prior_verifier_imported": False,
        "source_file": script.name,
        "source_sha256": sha256(script.read_bytes()).hexdigest(),
        "audited_proof_file": "TYPE_B_GENERAL_SMITH.md",
        "audited_proof_sha256_at_start": proof_start_hash,
        "cases": [],
    }
    for n, ell in PARAMETERS:
        case = check_case(n, ell)
        report["cases"].append(case)
        print(json.dumps({key: value for key, value in case.items()
                          if key != "all_exponent_weight_distributions"}), flush=True)
    report["audited_proof_sha256_at_end"] = sha256(proof.read_bytes()).hexdigest()
    report["proof_unchanged_during_run"] = (
        proof_start_hash == report["audited_proof_sha256_at_end"])
    report["summary"] = {
        "parameter_cases": len(report["cases"]),
        "partitions": sum(case["partitions"] for case in report["cases"]),
        "corner_membership_comparisons": sum(
            case["corner_membership_comparisons"] for case in report["cases"]),
        "exact_weighted_polynomial_comparisons": sum(
            case["exact_weighted_polynomial_comparisons"] for case in report["cases"]),
        "total_rectangle_gaussian_comparisons": len(report["cases"]),
        "counterexamples": 0,
    }
    report["seconds"] = round(perf_counter() - started, 6)
    report["status"] = "passed"
    destination = script.with_name("type_b_general_audit_verification.json")
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], **report["summary"],
                      "seconds": report["seconds"],
                      "report": str(destination)}), flush=True)


if __name__ == "__main__":
    main()
