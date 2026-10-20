"""Independent Cartesian corroboration of selected higher D/B square layers.

Run from borel/: PYTHONPATH=src sage -python experiments/verify_type_d_higher_square.py
This experiment proves no general equality. It compares exact coefficient-row
modules on rectangle Schur input bases, independently of the library's
fundamental-weight polynomial divided differences. Small cases cross-check
those two coordinate implementations. A full-flag all-even counterexample
guards against discarding the two-block symmetry in a proposed proof.
"""

from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import signal
from time import perf_counter

from sage.all import PolynomialRing, QQ, ZZ
from sage.env import SAGE_VERSION
from sage_borel import BorelMap, classical_borel
from sage_borel.operations import elementary_polynomial
from type_b_cartesian_frontier import CartesianCFrontier
from type_d_cartesian_frontier import CartesianDFrontier
from verify_type_d_square_comparison import (
    coefficient_matrix, schur_polynomials, smith_factors,
)


def rescaled_B(c, frontier, polynomial):
    return c._expansion({
        word: 2**sum(entry < 0 for entry in word.matrix().list())*value
        for word, value in frontier.expand(polynomial).items()
    })


def timeout_handler(_signal, _frame):
    raise TimeoutError("component computation exceeded its wall-clock budget")


def main():
    started = perf_counter()
    here = Path(__file__).resolve()
    borel = here.parent.parent
    sources = [here, here.with_name("type_d_cartesian_frontier.py"),
               here.with_name("type_b_cartesian_frontier.py"),
               here.with_name("verify_type_d_square_comparison.py"),
               borel/"src"/"sage_borel"/"core.py",
               borel/"src"/"sage_borel"/"operations.py"]
    report = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sage_version": SAGE_VERSION,
        "scope": "Finite exact evidence only; general higher square comparison is not inferred.",
        "normalization": "q_i=e_i(x_tail^2), polynomial degree is half cohomological degree",
        "independence": "Cartesian finite-geometric divided differences; B uses exact C/B row scaling.",
        "source_sha256": {str(p.relative_to(borel)): sha256(p.read_bytes()).hexdigest()
                          for p in sources},
        "resource_policy": {"per_component_wall_seconds": 120,
                            "overall_wall_seconds": 240,
                            "excluded_by_design": [
                                "An exploratory a=k=5, degree-24 run is not counted unless completed here.",
                                "These checks concern square sectors, not the full D Euler-sector gluing."
                            ]},
        "selected_components": [], "excluded_components": [],
    }
    output = here.with_name("type_d_higher_square_verification.json")

    def save():
        report["seconds"] = round(perf_counter()-started, 6)
        output.write_text(json.dumps(report, indent=2)+"\n")

    signal.signal(signal.SIGALRM, timeout_handler)
    configurations = [(2, 2, tuple(range(5))), (4, 4, (8, 12, 16)),
                      (4, 5, (8,)), (5, 5, (8,))]
    for a, k, degrees in configurations:
        n = a+k
        d, c = classical_borel("D", n, k), classical_borel("C", n, k)
        fd, fc = CartesianDFrontier(d), CartesianCFrontier(c)
        ring = PolynomialRing(QQ, n, names=["x%d" % i for i in range(1, n+1)])
        x = ring.gens()
        for m in degrees:
            row = {"a": a, "k": k, "polynomial_degree": 2*m}
            remaining = int(240-(perf_counter()-started))
            if remaining <= 0:
                row["reason"] = "overall computation budget exhausted"
                report["excluded_components"].append(row)
                save()
                continue
            component_started = perf_counter()
            signal.alarm(min(120, remaining))
            try:
                formal_ring, partitions, formal = schur_polynomials(a, k, m)
                substitution = formal_ring.hom([
                    elementary_polynomial(ring, [value*value for value in x[a:]], i)
                    for i in range(1, k+1)], ring)
                d_columns, b_columns = [], []
                for f in formal:
                    polynomial = substitution(f)
                    actual_d = fd.expand(polynomial)
                    if n == 4:
                        weight_substitution = formal_ring.hom([
                            d.tail_polynomial(k, i) for i in range(1, k+1)], d.ring)
                        assert actual_d == d.expand(weight_substitution(f))
                    d_columns.append(actual_d)
                    b_columns.append(rescaled_B(c, fc, polynomial))
                md = coefficient_matrix(d, d_columns)
                mb = coefficient_matrix(c, b_columns)
                assert md.rank() == mb.rank() == len(partitions)
                assert md.row_module() == mb.row_module()
                factors = smith_factors(md)
                assert factors == smith_factors(mb)
                assert all(factor == 2**factor.valuation(2) and
                           factor.valuation(2) <= min(a, k) for factor in factors)
                row.update({
                    "source_partitions": [list(p) for p in partitions],
                    "rank": len(partitions),
                    "D_smith_factors": list(map(str, factors)),
                    "B_smith_factors": list(map(str, factors)),
                    "D_nonzero_coefficient_rows": md.nrows(),
                    "B_nonzero_coefficient_rows": mb.nrows(),
                    "exact_integer_row_modules_equal": True,
                    "every_D_exponent_at_most_min_a_k": True,
                    "fundamental_weight_D_crosscheck": n == 4,
                    "seconds": round(perf_counter()-component_started, 6),
                })
                report["selected_components"].append(row)
                print("Cartesian D/B a=%d k=%d d=%d rank=%d factors=%s" %
                      (a, k, 2*m, len(partitions), factors), flush=True)
            except TimeoutError as error:
                row["reason"] = str(error)
                row["seconds"] = round(perf_counter()-component_started, 6)
                report["excluded_components"].append(row)
                print("Excluded timed-out component %s" % row, flush=True)
            finally:
                signal.alarm(0)
            save()

    # This is deliberately outside the selected two-block-symmetric problem.
    c, d = BorelMap(["C", 3]), BorelMap(["D", 3])
    ring = PolynomialRing(QQ, 3, names=["x1", "x2", "x3"])
    x = ring.gens()
    polynomial = x[0]**4*x[1]**2
    de = CartesianDFrontier(d).expand(polynomial)
    be = rescaled_B(c, CartesianCFrontier(c), polynomial)
    coordinates = d.classical_coordinates()
    assert de == d.expand(coordinates[0]**4*coordinates[1]**2)
    ds = smith_factors(coefficient_matrix(d, [de]))
    bs = smith_factors(coefficient_matrix(c, [be]))
    assert ds == (ZZ(4),) and bs == (ZZ(2),)
    report["full_flag_counterexample"] = {
        "n": 3, "input": "x1^4*x2^2", "polynomial_degree": 6,
        "B_smith_factors": list(map(str, bs)),
        "D_smith_factors": list(map(str, ds)),
        "fundamental_weight_D_crosscheck": True,
        "conclusion": "A claim for arbitrary all-even full-flag inputs is false; block symmetry must be retained.",
    }
    report["summary"] = {
        "selected_components_completed": len(report["selected_components"]),
        "selected_components_excluded": len(report["excluded_components"]),
        "selected_row_module_mismatches": 0,
        "selected_components_satisfying_exponent_bound":
            sum(row["every_D_exponent_at_most_min_a_k"]
                for row in report["selected_components"]),
        "fundamental_weight_D_components_crosschecked":
            sum(row["fundamental_weight_D_crosscheck"] for row in report["selected_components"])+1,
        "all_completed_checks_passed": True,
    }
    save()
    print(report["summary"], flush=True)


if __name__ == "__main__":
    main()
