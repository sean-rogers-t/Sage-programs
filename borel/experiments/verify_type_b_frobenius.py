"""Independently verify type-B Smith consequences with Schur LR arithmetic.

Run from borel/:
  sage -python experiments/verify_type_b_frobenius.py

Expected unit-factor counts use ordinary Grassmannian Frobenius over GF(2),
computed by Sage's Schur Littlewood--Richardson multiplication.  No BorelMap,
Weyl group, divided differences or library image calculations are imported.
The remaining assertions concern ranks, bounds and symmetry of the recorded
integer data; the mathematical proofs are in TYPE_B_SMITH.md.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import GF, Partitions, SymmetricFunctions, matrix
from sage.env import SAGE_VERSION


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=Path(__file__).with_name("orthogonal_smith_search_results.json"))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("type_b_frobenius_verification.json"))
    args = parser.parse_args()
    started = perf_counter()
    source = args.input.read_bytes()
    data = json.loads(source)
    schur = SymmetricFunctions(GF(2)).schur()
    counts = {"rectangle_rank_checks": 0, "frobenius_unit_count_checks": 0,
              "power_of_two_and_exponent_bound_checks": 0,
              "stable_range_checks": 0, "block_duality_checks": 0}
    records = {}
    case_count = 0
    for case in data["cases"]:
        if case["family"] != "B":
            continue
        case_count += 1
        a, k, n = (int(case[key]) for key in ("a", "k", "n"))
        assert a + k == n and min(a, k) >= 1
        for result in case["degrees"]:
            if result["status"] != "computed":
                continue
            degree = int(result["polynomial_degree"])
            assert degree % 2 == 0
            m = degree // 2
            partitions = tuple(Partitions(m, max_part=a, max_length=k))
            target = tuple(Partitions(2 * m, max_part=a, max_length=k))
            factors = tuple(int(value) for value in result["smith_factors"])
            assert len(partitions) == int(result["rank"]) == len(factors), (n, k, degree)
            counts["rectangle_rank_checks"] += 1

            squares = tuple(schur(partition) ** 2 for partition in partitions)
            frobenius = matrix(GF(2), len(target), len(partitions),
                               lambda row, column: squares[column][target[row]])
            assert frobenius.rank() == factors.count(1), (n, k, degree)
            counts["frobenius_unit_count_checks"] += 1

            assert all(value > 0 and value & (value - 1) == 0
                       and value <= 2 ** min(a, k) for value in factors), (n, k, degree)
            counts["power_of_two_and_exponent_bound_checks"] += 1
            if degree <= max(a, k):
                assert all(value == 1 for value in factors), (n, k, degree)
                counts["stable_range_checks"] += 1
            records[n, a, k, degree] = factors

    for (n, a, k, degree), factors in records.items():
        # Check each pair once, including self-dual blocks.
        if a <= k and (n, k, a, degree) in records:
            assert factors == records[n, k, a, degree], (n, a, k, degree)
            counts["block_duality_checks"] += 1

    report = {
        "suite_result": "passed",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sage_version": SAGE_VERSION,
        "input_file": args.input.name,
        "input_sha256": sha256(source).hexdigest(),
        "input_bytes": len(source),
        "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "type_B_cases": case_count,
        "maximum_rank": max(key[0] for key in records),
        "maximum_polynomial_degree": max(key[3] for key in records),
        "checks": counts,
        "independent_expected_values": (
            "Schur Littlewood--Richardson products over GF(2), restricted to "
            "the ordinary complex Grassmannian rectangle; no Borel engine used"),
        "scope": (
            "Finite verification of proved consequences; does not classify "
            "all higher two-adic multiplicities or establish novelty"),
        "seconds": round(perf_counter() - started, 6),
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
