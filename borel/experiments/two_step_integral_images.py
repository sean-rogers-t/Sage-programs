"""Reproducible bounded joint-image experiment; no general conjecture assumed.

Run from borel/: PYTHONPATH=src sage -python experiments/two_step_integral_images.py
Outputs exact degree tables plus independently checked Cartesian BGG matrices.
"""

from itertools import combinations
import hashlib
import json
from pathlib import Path
import sys
import time

from sage.all import ZZ, matrix, prod
from sage.env import SAGE_VERSION

from sage_borel import BorelMap, GeneratorAlgebra, two_step_symplectic

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.oracle import OrthonormalBGG
from tests.test_integral import determinantal_smith_factors


def elementary(values, j):
    return sum(prod(c) for c in combinations(values, j))


def algebras(kind, n, b, k):
    a = n-b-k
    if kind == "C":
        flag = two_step_symplectic(n, b, k)
        return flag, flag.middle_algebra(), flag.tail_algebra(), flag.joint_algebra()
    flag = BorelMap("%s%d" % (kind, n), omitted_nodes=(a, a+b))
    x = flag.classical_coordinates()
    middle = {"c%d" % j: elementary(x[a:a+b], j) for j in range(1, b+1)}
    tail = {"e%d" % j: elementary([u*u for u in x[-k:]], j)
            for j in range(1, k+1 if kind == "B" else k)}
    if kind == "D":
        tail["q"] = prod(x[-k:])
    return (flag, GeneratorAlgebra(flag, middle), GeneratorAlgebra(flag, tail),
            GeneratorAlgebra(flag, dict(middle, **tail)))


def summary(lattice):
    return {"rank": int(lattice.rank),
            "smith_factors": [int(x) for x in lattice.smith_factors],
            "saturation_index": int(lattice.saturation_index)}


def independently_check(algebra, report, oracle):
    inputs = [oracle.from_fundamental(f) for f in algebra.generator_polynomials]
    columns = []
    for monomial in report.monomials:
        powers = monomial.exponents()[0]
        polynomial = prod(f**e for f, e in zip(inputs, powers))
        raw = oracle.expansion(polynomial, report.degree)
        columns.append({algebra.borel_map.W.from_reduced_word(w.reduced_word()): ZZ(c)
                        for w, c in raw.items()})
    assert set(report.schubert_indices) == set().union(*(set(c) for c in columns))
    expected = matrix(ZZ, len(report.schubert_indices), len(columns),
                      {(i, j): col.get(w, 0) for i, w in enumerate(report.schubert_indices)
                       for j, col in enumerate(columns)})
    assert expected == report.matrix
    # Small source spaces permit the independent determinantal-divisor check.
    assert determinantal_smith_factors(expected) == report.smith_factors
    return len(columns)


def main():
    start = time.time()
    cases = [("C", 3, 1, 1, 9), ("C", 4, 1, 2, 12),
             ("C", 4, 2, 1, 14), ("C", 4, 1, 1, 12), ("C", 5, 2, 2, 8),
             ("B", 3, 1, 1, 9), ("B", 4, 1, 2, 12), ("D", 4, 1, 2, 12)]
    output = []
    independent_matrices = independent_columns = 0
    for kind, n, b, k, max_degree in cases:
        flag, middle, tail, joint = algebras(kind, n, b, k)
        # Keep the independent Weyl enumeration bounded; D4 has order 192.
        oracle = OrthonormalBGG("%s%d" % (kind, n)) if n <= 4 else None
        predicted = flag.predicted_hilbert_prefix(max_degree) if kind == "C" else None
        rows = []
        for degree in range(max_degree+1):
            m, t, j = (algebra.integral_degree(degree) for algebra in (middle, tail, joint))
            comparison = m.compare(t)
            joint_comparison = comparison.sum.compare(j)
            assert joint_comparison.left_in_right
            assert comparison.sum.rank + comparison.intersection.rank == m.rank+t.rank
            assert comparison.sum.quotient_by(m).as_dict() == comparison.right.quotient_by(comparison.intersection).as_dict()
            assert comparison.sum.quotient_by(t).as_dict() == comparison.left.quotient_by(comparison.intersection).as_dict()
            if predicted is not None:
                assert j.rank == predicted[degree]
            extra_check = (kind == "B" and n == 3 and degree == 7) or (kind == "D" and degree == 8)
            if oracle is not None and (degree <= 4 or extra_check):
                for algebra, report in ((middle, m), (tail, t), (joint, j)):
                    independent_columns += independently_check(algebra, report, oracle)
                    independent_matrices += 1
            rows.append({"degree": degree, "middle": summary(m), "tail": summary(t),
                         "sum": summary(comparison.sum),
                         "intersection": summary(comparison.intersection), "joint": summary(j),
                         "middle_in_tail": comparison.left_in_right,
                         "tail_in_middle": comparison.right_in_left,
                         "joint_mod_sum": j.lattice().quotient_by(comparison.sum).as_dict(),
                         "sum_mod_middle": comparison.sum.quotient_by(m).as_dict(),
                         "sum_mod_tail": comparison.sum.quotient_by(t).as_dict()})
        case = {"type": kind, "n": n, "a": n-b-k, "b": b, "k": k,
                "omitted_nodes": [int(i) for i in flag.omitted_nodes],
                "generator_names": list(joint.generator_names),
                "generator_degrees": [int(i) for i in joint.generator_degrees],
                "degrees": rows}
        output.append(case)
        print(json.dumps({key: case[key] for key in ("type", "n", "a", "b", "k")}) +
              " joint indices=" + str([r["joint"]["saturation_index"] for r in rows]) +
              " quotient torsion=" + str([r["joint_mod_sum"]["torsion_factors"] for r in rows]), flush=True)
    report = {"result": "passed", "sage_version": SAGE_VERSION,
              "seconds": round(time.time()-start, 3),
              "command": "PYTHONPATH=src sage -python experiments/two_step_integral_images.py",
              "scope": "Finite exact computations, not a general classification or novelty claim.",
              "normalization": "Middle elementary functions in Cartesian x; tail elementary functions in x^2; type D replaces the last elementary function by the Euler product q.",
              "independent_cartesian_bgg_matrices": independent_matrices,
              "independent_cartesian_bgg_columns": independent_columns,
              "independent_smith_check": "Gcds of all minors on every independently checked matrix.",
              "degree_rows": sum(len(c["degrees"]) for c in output), "cases": output,
              "source_sha256": {str(path.relative_to(ROOT)).replace("\\", "/"):
                                hashlib.sha256(path.read_bytes()).hexdigest()
                                for path in sorted(list((ROOT / "src/sage_borel").glob("*.py")) +
                                                   [Path(__file__).resolve(), ROOT / "tests/oracle.py",
                                                    ROOT / "tests/test_integral.py"])}}
    (ROOT / "two-step-integral-results.json").write_text(json.dumps(report, indent=2)+"\n")
    print("Verified %d degree rows; %d independent matrices in %.3fs" %
          (report["degree_rows"], independent_matrices, report["seconds"]), flush=True)


if __name__ == "__main__":
    main()
