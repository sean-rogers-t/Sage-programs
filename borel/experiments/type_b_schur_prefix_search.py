"""Independent ordinary-Schur verification of type-B coefficient prefixes.

Run from borel/:
  sage -python experiments/type_b_schur_prefix_search.py

No BorelMap, Weyl group or divided-difference calculation is imported. The
integral prefix description proved in TYPE_B_FILTRATION.md is evaluated with
ordinary Littlewood--Richardson multiplication to search additional larger
rectangles for nonprimitive prefixes. Recorded BGG oracle comparisons live
in the separate verify_type_b_filtration.py suite. A passed finite search
is not a proof of general prefix primitivity.
"""

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from sage.all import Partitions, SymmetricFunctions, ZZ, matrix
from sage.env import SAGE_VERSION


class SchurPrefixes:
    def __init__(self):
        self.schur = SymmetricFunctions(ZZ).schur()
        self.rectangles = {}
        self.products = {}

    def partitions(self, size, width, height):
        key = (size, width, height)
        if key not in self.rectangles:
            self.rectangles[key] = tuple(tuple(map(int, p)) for p in
                Partitions(size, max_part=width, max_length=height)) if height else (
                    ((),) if size == 0 else ())
        return self.rectangles[key]

    def product(self, partitions):
        key = tuple(sorted(p for p in partitions if p))
        if key not in self.products:
            value = self.schur.one()
            for partition in key:
                value *= self.schur(partition)
            self.products[key] = {tuple(map(int, p)): int(c) for p, c in value}
        return self.products[key]

    def prefix(self, A, head, m, threshold):
        n, b = A + head, head - threshold
        b0, b1 = (b + 1) // 2, b // 2
        A0, A1 = (n + 1) // 2 - b0, n // 2 - b1
        target = self.partitions(m, A, head)
        rows = set()
        for mu_degree in range(m + 1):
            for alpha_degree in range(m - mu_degree + 1):
                beta_degree = m - mu_degree - alpha_degree
                for mu in self.partitions(mu_degree, A, threshold):
                    for alpha in self.partitions(alpha_degree, A0, b0):
                        for beta in self.partitions(beta_degree, A1, b1):
                            coefficients = self.product((mu, alpha, beta))
                            row = tuple(coefficients.get(p, 0) for p in target)
                            if any(row):
                                rows.add(row)
        mat = matrix(ZZ, len(rows), len(target), list(rows))
        diagonal = mat.smith_form(transformation=False).diagonal()
        factors = tuple(abs(int(value)) for value in diagonal if value)
        return {"threshold": threshold, "rank": len(factors),
                "source_rank": len(target), "distinct_nonzero_LR_rows": len(rows),
                "smith_factors": list(factors), "primitive": all(value == 1 for value in factors)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("type_b_schur_prefix_search_results.json"))
    args = parser.parse_args()
    started = perf_counter()
    engine = SchurPrefixes()
    reports, counterexamples = [], []
    larger_checks = 0
    for A, head, maximum in ((5, 5, 15), (6, 5, 12), (6, 6, 12)):
        case_report = {"family": "B", "source": "larger_LR_only_search",
                       "tail": A, "head": head, "n": A + head,
                       "maximum_half_polynomial_degree": maximum, "degrees": []}
        for m in range(maximum + 1):
            prefixes = [engine.prefix(A, head, m, threshold) for threshold in range(head - 1)]
            # The omitted last two levels are already proved primitive
            # (highest-factor note and full quaternionic source lattice).
            for prefix in prefixes:
                if not prefix["primitive"]:
                    counterexamples.append({"tail": A, "head": head, "m": m, **prefix})
                larger_checks += 1
            case_report["degrees"].append({"half_polynomial_degree": m, "prefixes": prefixes})
        reports.append(case_report)

    report = {"suite_result": "passed" if not counterexamples else "counterexamples_found",
              "generated_utc": datetime.now(timezone.utc).isoformat(), "sage_version": SAGE_VERSION,
              "scope": "Finite exact integral Schur-product verification; general primitivity is not inferred",
              "independence": "Ordinary LR multiplication and integer Smith forms; no sage_borel/Weyl/BGG imports",
              "source_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "checks": {"larger_LR_only_prefix_Smith_checks": larger_checks},
              "cases": reports, "counterexamples": counterexamples,
              "seconds": round(perf_counter() - started, 6)}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "cases"}, sort_keys=True))


if __name__ == "__main__":
    main()
