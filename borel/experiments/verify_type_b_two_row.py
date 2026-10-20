"""Verify the closed type-B two-row Smith theorem.

Run from borel/:
  sage -python experiments/verify_type_b_two_row.py

This script does not import sage_borel. It checks exact identities in two
independent polynomial variables, then compares the closed partition statistic
with the Smith factors recorded by the Borel calculations. The general proof
is in TYPE_B_TWO_ROW_SMITH.md; finite verification is not its substitute.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import PolynomialRing, ZZ, matrix
from sage.env import SAGE_VERSION


def exponent(A, r, s):
    u, v = (A + 1) // 2, A // 2
    if r <= u and s <= v:
        return 0
    return 2 if s > u else 1


def expected(A, m):
    return tuple(sorted(2 ** exponent(A, m - s, s)
                        for s in range(max(0, m - A), m // 2 + 1)))


def nonzero_smith(mat):
    diagonal = mat.smith_form(transformation=False).diagonal()
    return tuple(abs(int(value)) for value in diagonal if value)


def verify_polynomials(maximum):
    ring = PolynomialRing(ZZ, names=("y1", "y2"))
    y1, y2 = ring.gens()
    limit = 4 * maximum + 2
    h = [ring(sum(y1 ** (d - i) * y2 ** i for i in range(d + 1)))
         for d in range(limit + 1)]
    hsquare = [ring(sum(y1 ** (2 * (d - i)) * y2 ** (2 * i)
                       for i in range(d + 1))) for d in range(2 * maximum + 1)]
    counts = {"two_variable_plethysm_identities": 0,
              "pfaffian_telescoping_identities": 0,
              "high_degree_single_theta_identities_mod_relation": 0,
              "closed_factors_vs_two_row_smith_matrices": 0,
              "finite_high_complete_relation_identities": 0}
    def complete(d):
        return h[d] if d >= 0 else ring.zero()

    for A in range(2, maximum + 1):
        # The tail-rank relation h_j(y^2)=0 for j>A follows from these
        # first two relations by the exact two-variable recurrence.
        relations = {A + 1: (ring.one(), ring.zero()),
                     A + 2: (ring.zero(), ring.one())}
        for degree in range(A + 3, 2 * A + 1):
            previous, before = relations[degree - 1], relations[degree - 2]
            relations[degree] = tuple((y1 ** 2 + y2 ** 2) * previous[i]
                                     - (y1 * y2) ** 2 * before[i] for i in (0, 1))
        for degree, (left, right) in relations.items():
            if degree <= 2 * A:
                assert hsquare[degree] == left * hsquare[A + 1] + right * hsquare[A + 2]
                counts["finite_high_complete_relation_identities"] += 1

        for m in range(2 * A + 1):
            source = tuple(range(max(0, m - A), m // 2 + 1))
            schurs = [complete(2 * m - j) * complete(j)
                      - complete(2 * m - j + 1) * complete(j - 1)
                      for j in range(m + 1)]
            pfaffians = [sum((-1) ** (j - i) * schurs[i] for i in range(j + 1))
                         for j in range(m + 1)]
            assert pfaffians[m] == (-1) ** m * hsquare[m]
            for j in range(m + 1):
                direct = complete(2 * m - j) * complete(j) + 2 * sum(
                    (-1) ** shift * complete(2 * m - j + shift) * complete(j - shift)
                    for shift in range(1, j + 1))
                assert pfaffians[j] == direct
                counts["pfaffian_telescoping_identities"] += 1

            for s in source:
                r = m - s
                polynomial = (y1 * y2) ** (2 * s) * hsquare[r - s]
                assert polynomial == sum((-1) ** j * schurs[j]
                                         for j in range(2 * s, m + 1))
                counts["two_variable_plethysm_identities"] += 1
                if m > A:
                    assert polynomial == pfaffians[2 * s - 1] + hsquare[m]
                    counts["high_degree_single_theta_identities_mod_relation"] += 1

            if m <= A:
                # Build the original overlapping Schur columns, not their
                # disjoint-support diagonalization in the written proof.
                mat = matrix(ZZ, m + 1, len(source), lambda j, column:
                    (-1) ** j * 2 ** int(2 * m - j > A)
                    if j >= 2 * source[column] else 0)
            else:
                mat = matrix.diagonal(ZZ, [2 ** (1 + int(2 * s - 1 > A)) for s in source])
            assert nonzero_smith(mat) == expected(A, m), (A, m)
            counts["closed_factors_vs_two_row_smith_matrices"] += 1
    return counts


def verify_recorded(directory):
    checks, inputs = 0, []
    per_input = {}
    for name in ("orthogonal_smith_search_results.json", "type_b_multiplicity_search_results.json"):
        path = directory / name
        if not path.exists():
            continue
        raw = path.read_bytes()
        data = json.loads(raw)
        checked = 0
        for case in data["cases"]:
            if case.get("family", "B") != "B" or min(case["a"], case["k"]) != 2:
                continue
            A = max(case["a"], case["k"])
            for row in case["degrees"]:
                if row["status"] != "computed":
                    continue
                degree = int(row["polynomial_degree"])
                assert degree % 2 == 0
                actual = tuple(int(value) for value in row["smith_factors"])
                assert actual == expected(A, degree // 2), (name, case["n"], case["k"], degree)
                checked += 1
        checks += checked
        per_input[name] = checked
        inputs.append({"name": name, "sha256": sha256(raw).hexdigest(), "bytes": len(raw)})
    return checks, per_input, inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-long-block", type=int, default=12)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("type_b_two_row_verification.json"))
    args = parser.parse_args()
    started = perf_counter()
    counts = verify_polynomials(args.maximum_long_block)
    checks, per_input, inputs = verify_recorded(Path(__file__).parent)
    report = {"suite_result": "passed", "sage_version": SAGE_VERSION,
              "generated_utc": datetime.now(timezone.utc).isoformat(),
              "maximum_long_block_symbolic_verification": args.maximum_long_block,
              "symbolic_checks": counts, "recorded_Borel_Smith_degrees_checked": checks,
              "recorded_checks_by_input": per_input, "input_files": inputs,
              "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "independence": "No sage_borel imports; polynomial and integer Smith arithmetic only",
              "seconds": round(perf_counter() - started, 6)}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
