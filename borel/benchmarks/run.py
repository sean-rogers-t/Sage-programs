"""Reproducible, bounded benchmarks; run from borel/ with SageMath.

PYTHONPATH=src sage -python benchmarks/run.py --output benchmarks/results.json
Use --skip-legacy for just the new library.  Historical comparisons are
deliberately restricted to ranks 3 and 4: list(W) at rank 8 has 10,321,920
elements, before the repeated products are even performed.
"""

import argparse
from contextlib import redirect_stdout
import io
from itertools import combinations
import json
from pathlib import Path
import platform
from time import perf_counter

import sage.all as sage
from sage.env import SAGE_VERSION
from sage.repl.load import load

from sage_borel import BorelMap, __version__


HERE = Path(__file__).resolve().parent


def timed(function):
    start = perf_counter()
    result = function()
    return result, perf_counter() - start


def new_case(kind, k, degree):
    rank = int(kind[1:])
    omitted = rank + 1 - k if kind.startswith("A") else rank - k
    b, setup = timed(lambda: BorelMap(kind, omitted_nodes=omitted))
    answer, elapsed = timed(lambda: b.tail_elementary(k, degree))
    # Independent literal calculation in orthonormal coordinates, where the
    # tail polynomial has binomial(k,degree) monomials.  Expanding it first in
    # fundamental weights would defeat the very optimization being measured.
    coordinate_ring = sage.PolynomialRing(sage.QQ, rank + 1 if kind.startswith("A") else rank,
                                         names="x")
    x = coordinate_ring.gens()
    power = 1 if kind.startswith("A") else 2
    monomials = {}
    for subset in combinations(range(len(x) - k, len(x)), degree):
        exponent = [0] * len(x)
        for j in subset:
            exponent[j] = power
        monomials[tuple(exponent)] = sage.QQ.one()
    literal = coordinate_ring(monomials)
    classes = b.basis(degree if kind.startswith("A") else 2 * degree)
    for w in classes:
        coefficient = answer[w]
        value = literal
        for node in reversed(w.reduced_word()):
            images = list(x)
            if node < rank or kind.startswith("A"):
                images[node - 1], images[node] = images[node], images[node - 1]
                root = x[node - 1] - x[node]
            elif kind.startswith("C"):
                images[-1] = -images[-1]
                root = 2 * x[-1]
            elif kind.startswith("D"):
                images[-2], images[-1] = -images[-1], -images[-2]
                root = x[-2] + x[-1]
            else:
                raise ValueError("benchmark oracle implements A/C/D only")
            value, remainder = (value - value(*images)).quo_rem(root)
            assert not remainder
        assert value == coefficient
    return {
        "cartan_type": kind, "tail_size": k, "elementary_index": degree,
        "polynomial_degree": degree if kind.startswith("A") else 2 * degree,
        "setup_seconds": round(setup, 6), "calculation_seconds": round(elapsed, 6),
        "weyl_group_order": str(b.W.cardinality()),
        "stats": dict(b.last_stats), "terms": answer.to_terms(),
        "literal_coefficient_checks": len(classes),
        "checked_complete_degree_basis": True,
    }


def legacy_comparison(n, k, degree):
    namespace = dict(vars(sage))
    repo = HERE.parents[1]
    load(str(repo / "chevalley.sage"), namespace)
    load(str(repo / "elem.sage"), namespace)
    with redirect_stdout(io.StringIO()):
        old, old_seconds = timed(lambda: namespace["cohom_0"](n, k, degree))
    b, setup = timed(lambda: BorelMap(["C", n], omitted_nodes=n - k))
    new, new_seconds = timed(lambda: b.tail_elementary(k, degree))
    expected = {}
    for coefficient, w in old:
        local = b.W.from_reduced_word(w.reduced_word())
        expected[local] = expected.get(local, 0) + coefficient
    expected = {w: coefficient for w, coefficient in expected.items() if coefficient}
    assert dict(new) == expected, "legacy comparison disagrees"
    return {"cartan_type": "C" + str(n), "tail_size": k,
            "elementary_index": degree, "legacy_seconds": round(old_seconds, 6),
            "new_setup_seconds": round(setup, 6),
            "new_calculation_seconds": round(new_seconds, 6),
            "calculation_speedup": round(old_seconds / new_seconds, 2),
            "equal_results": True, "output_terms": len(new)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-legacy", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    report = {"sage_version": SAGE_VERSION, "library_version": __version__,
              "python_version": platform.python_version(), "platform": platform.platform(),
              "timing_method": "one run in one Sage process; setup separated; seconds are machine dependent",
              "new_cases": [], "legacy_comparisons": []}

    def save():
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")

    for kind, k, degree in (("C8", 6, 2), ("C12", 10, 5), ("C20", 18, 5),
                            ("D12", 10, 5), ("A15", 13, 5)):
        case = new_case(kind, k, degree)
        report["new_cases"].append(case)
        save()
        print("new {} k={} i={}: {:.6f}s, {} terms, frontier {}".format(
            kind, k, degree, case["calculation_seconds"], case["stats"]["output_terms"],
            case["stats"]["max_frontier"]), flush=True)
    if not args.skip_legacy:
        for n, k, degree in ((3, 2, 1), (4, 3, 1)):
            case = legacy_comparison(n, k, degree)
            report["legacy_comparisons"].append(case)
            save()
            print("legacy C{}: {:.6f}s -> {:.6f}s ({:.2f}x), exact agreement".format(
                n, case["legacy_seconds"], case["new_calculation_seconds"],
                case["calculation_speedup"]), flush=True)
    save()


if __name__ == "__main__":
    main()
