"""Independent SymPy/stdlib corroboration; imports no Sage or sage_borel.

Run with a Python installation containing SymPy. This independently builds
the degree-five integer lattices, not just the displayed witness identity.
"""
from itertools import combinations
from math import comb, gcd
from pathlib import Path
import json

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form


ROOT = Path(__file__).resolve().parents[1]


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first, length-1):
                yield (first,)+rest


def main():
    t1, t2, t3, h = sp.symbols("t1 t2 t3 h")
    ts = (t1, t2, t3)
    weights = ts+(-sum(ts),)
    a, b, c = [sp.expand(sum(sp.prod(v) for v in combinations(weights, j)))
               for j in (2, 3, 4)]
    relation = h**4+a*h**2-b*h+c

    def reduce(f):
        return sp.rem(sp.Poly(f, h), sp.Poly(relation, h)).as_expr().expand()

    generators = (4*h, 2*h**2, 4*h**3+2*a*h)
    u, v, w = generators
    nu = 2*a*a*h
    identities = (
        u*u-8*v, u*v-2*w+a*u,
        u*w+4*a*v-4*b*u+16*c,
        v*v+2*a*v-b*u+4*c,
        nu-v*w-a*w+4*b*v-2*c*u,
        2*nu-a*a*u,
    )
    assert all(reduce(f) == 0 for f in identities)
    assert sp.Poly(a*a, *ts).coeff_monomial(t1**4) == 1
    fs = (4*h, 6*h*h+a, 4*h**3+2*a*h-b)
    for j, f in enumerate(fs, 1):
        for i, value in enumerate(weights):
            direct = sum(sp.prod(value-weights[k] for k in chosen)
                         for chosen in combinations([k for k in range(4) if k != i], j))
            assert sp.expand(f.subs(h, value)-direct) == 0

    def image_generators(degree):
        result = []
        for i in range(degree+1):
            for j in range(degree//2+1):
                for k in range(degree//3+1):
                    used = i+2*j+3*k
                    if used <= degree:
                        reduced = reduce(u**i*v**j*w**k)
                        for exponents in compositions(degree-used, 3):
                            result.append(sp.expand(reduced*sp.prod(z**e for z, e in zip(ts, exponents))))
        return result

    def column(f, degree):
        polynomial = sp.Poly(reduce(f), h, *ts)
        return sp.Matrix([polynomial.coeff_monomial(
            h**i*sp.prod(z**e for z, e in zip(ts, exponents)))
            for i in range(min(3, degree)+1)
            for exponents in compositions(degree-i, 3)])

    degree = 5
    Jcols = sp.Matrix.hstack(*(column(f, degree) for f in image_generators(degree)))
    IJcols = sp.Matrix.hstack(*(column(z*f, degree) for z in ts for f in image_generators(degree-1)))
    J = hermite_normal_form(Jcols)
    IJ = hermite_normal_form(IJcols)
    assert J.shape == IJ.shape == (52, 52)
    transition = J.inv()*IJ
    assert all(value.q == 1 for value in transition)
    diagonal = smith_normal_form(transition, domain=sp.ZZ)
    factors = [abs(int(diagonal[i, i])) for i in range(52)]
    assert factors == [1]*51+[2]
    nuJ = J.inv()*column(nu, degree)
    nuIJ = IJ.inv()*column(nu, degree)
    assert all(value.q == 1 for value in nuJ)
    assert any(value.q == 2 for value in nuIJ)
    assert all((2*value).q == 1 for value in nuIJ)

    arithmetic = 0
    for n in range(2, 81):
        indices = [1]
        for d in range(1, n):
            value = 0
            for j in range(1, d+1):
                value = gcd(value, comb(n, j)*indices[d-j])
            assert value == n//gcd(n, d)
            indices.append(value)
            arithmetic += 1
    report = dict(result="passed", engine="SymPy "+sp.__version__+" and Python stdlib",
                  sage_imported=False, runtime_imported=False,
                  independently_rebuilt_degree5_ambient_rank=52,
                  degree5_quotient_nonunit_factors=[2],
                  witness_in_image=True, witness_not_in_parameter_times_image=True,
                  twice_witness_in_parameter_times_image=True,
                  fixed_point_generator_checks=12, polynomial_identities=6,
                  ordinary_family_positive_degree_rows=arithmetic,
                  source_scope="Symbolic quotient and integer-lattice checks; separate written arguments prove the general statements and identify the literal source.")
    (ROOT/"experiments/local-image-independent.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
