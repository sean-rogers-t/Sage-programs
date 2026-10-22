"""Check selected numerical examples in the exposition with SageMath.

From borel/: PYTHONPATH=src sage -python ../articles/classical-equivariant-exposition/verify_examples.py
These checks corroborate the displayed examples, not the general proofs.
"""
from hashlib import sha256
import json
from pathlib import Path

from sage.all import ZZ, PolynomialRing, matrix
from sage_borel import BorelMap, GeneratorAlgebra, classical_borel, orthogonal_smith

ROOT = Path(__file__).resolve().parent


def main():
    cases = []
    b = classical_borel("B", 3, 1)
    a = GeneratorAlgebra(b, {"q": b.tail_polynomial(1, 1)})
    for degree, expected in ((0, (1,)), (2, (1,)), (4, (2,))):
        actual = tuple(int(s) for s in a.integral_degree(degree).smith_factors)
        formula = tuple(int(s) for s in orthogonal_smith("B", 3, 1, degree).smith_factors())
        assert actual == formula == expected
        cases.append(dict(family="B", n=3, k=1, polynomial_degree=degree, factors=actual))

    d = classical_borel("D", 4, 2)
    a = GeneratorAlgebra(d, {"q": d.tail_polynomial(2, 1), "p": d.euler_polynomial(2)})
    report = a.integral_degree(4, include_ambient=True)
    assert tuple(report.smith_factors) == (1, 2, 2)
    assert report.rank == 3 and report.ambient_dimension == 4
    assert tuple(orthogonal_smith("D", 4, 2, 4).smith_factors()) == (1, 2, 2)
    cases.append(dict(family="D", n=4, k=2, polynomial_degree=4,
                      factors=[1, 2, 2], ambient_rank=4, free_cokernel_rank=1))

    b = BorelMap("A3", omitted_nodes=1)
    w1, w2, w3 = b.weights
    eps = [w1, w2 - w1, w3 - w2, -w3]
    roots = [eps[0] - x for x in eps[1:]]
    f1 = sum(roots)
    f2 = sum(roots[i] * roots[j] for i in range(3) for j in range(i + 1, 3))
    f3 = roots[0] * roots[1] * roots[2]
    a = GeneratorAlgebra(b, {"E1": f1, "E2": f2, "E3": f3})
    actual = [tuple(int(s) for s in a.integral_degree(d).smith_factors) for d in range(4)]
    assert actual == [(1,), (4,), (2,), (4,)]
    cases.append(dict(family="A", group="SL4", flag="P1", factors_by_degree=actual))

    ring = PolynomialRing(ZZ, "v")
    v = ring.gen()
    b = 2*v - v**2 + v**3
    cols = [[int(((b**j) % v**4)[i]) for i in range(4)] for j in range(4)]
    inclusion = matrix(ZZ, cols).transpose()
    expected = matrix(ZZ, [[1, 0, 0, 0], [0, 2, 0, 0], [0, -1, 4, 0], [0, 1, -4, 8]])
    assert inclusion == expected and abs(inclusion.det()) == 64
    assert tuple(inclusion.elementary_divisors()) == (1, 1, 8, 8)
    assert 1 * 2 * 4 * 8 == 64
    cases.append(dict(family="C", group="Sp4", flag="P1",
                      K_matrix=[list(map(int, row)) for row in inclusion.rows()],
                      K_factors=[1, 1, 8, 8], H_factors=[1, 2, 4, 8], index=64))

    result = dict(result="passed", cases=cases,
                  methods="BGG Schubert image matrices versus closed orthogonal formulas; root-difference source on SL4; direct truncated projective K multiplication",
                  scope="Selected displayed examples only; not an independent audit of all general proofs.",
                  source_sha256=sha256(Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
                  tex_sha256=sha256((ROOT/"main.tex").read_bytes().replace(b"\r\n", b"\n")).hexdigest())
    (ROOT/"examples-verification.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
