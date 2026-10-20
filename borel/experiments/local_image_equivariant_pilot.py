"""Exact checks for LOCAL_IMAGE_BOUND.md and EQUIVARIANT_IMAGE_PILOT.md.

Run from borel/: PYTHONPATH=src sage -python experiments/local_image_equivariant_pilot.py
No new runtime API or universal equivariant backend is asserted.
"""
from itertools import combinations
from pathlib import Path
import json

from sage.all import (PolynomialRing, QQ, ZZ, binomial, gcd, matrix, prod,
                      vector, xgcd)
from sage_borel import BorelMap, irreducible_character_lattice


ROOT = Path(__file__).resolve().parents[1]


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first,) + rest


def dd(borel, f, word):
    for node in reversed(tuple(word)):
        f = borel.divided_difference(f, node)
    return f


def top_certificate(kind):
    b = BorelMap(kind)
    w0 = b.W.long_element()
    N = int(w0.length())
    t, F, tried = ZZ(0), b.ring.zero(), 0
    for exponents in compositions(N, b.rank):
        monomial = prod(x**e for x, e in zip(b.weights, exponents))
        value = ZZ(dd(b, monomial, w0.reduced_word()))
        common, old, new = xgcd(t, value)
        F = old*F + new*monomial
        t = common
        tried += 1
        if t == 1:
            break
    assert t > 0 and dd(b, F, w0.reduced_word()) == t
    assert all(c.denominator() == 1 for c in F.coefficients())
    checked = 0
    for count in range(b.rank + 1):
        for omitted in combinations(b.nodes, count):
            retained = tuple(i for i in b.nodes if i not in omitted)
            for w in b.W:
                if any(w.has_descent(i, side="right") for i in retained):
                    continue
                fw = dd(b, F, (w.inverse()*w0).reduced_word())
                assert all(b.reflect(fw, i) == fw for i in retained)
                assert b.expand(fw).as_dict() == {w: QQ(t)}
                checked += 1
    return dict(cartan_type=kind, torsion_index=int(t), top_polynomial=str(F),
                top_monomials_checked=tried, parabolic_schubert_lifts=checked)


def projective_checks():
    cases = []
    for n in range(2, 9):
        b = BorelMap("A%d" % (n-1), omitted_nodes=1)
        w = b.weights
        epsilon = (w[0],) + tuple(w[i]-w[i-1] for i in range(1, n-1)) + (-w[-1],)
        differences = tuple(epsilon[0]-a for a in epsilon[1:])
        indices = []
        for j in range(1, n):
            f = sum(prod(c) for c in combinations(differences, j))
            image = b.expand(f)
            assert image == binomial(n, j)*b.expand(w[0]**j)
            indices.append(int(binomial(n, j)))
        lattice = irreducible_character_lattice("A%d" % (n-1), (1,)+(0,)*(n-2))
        assert abs(lattice.linear_basis_matrix.det()) == n
        assert int(lattice.linear_basis_matrix.elementary_divisors()[-1]) == n
        factors = [ZZ(1)]
        for d in range(1, n):
            factors.append(gcd([ZZ(indices[j-1])*factors[d-j] for j in range(1, d+1)]))
            assert factors[-1] == n//gcd(n, d)
        cases.append(dict(n=n, generator_coefficients=indices,
                          image_factors=[int(a) for a in factors]))
    # Independent arithmetic range: dynamic gcd over all source monomials.
    arithmetic_rows = 0
    for n in range(2, 81):
        factors = [ZZ(1)]
        for d in range(1, n):
            value = gcd([binomial(n, j)*factors[d-j] for j in range(1, d+1)])
            assert value == n//gcd(n, d)
            factors.append(value)
            arithmetic_rows += 1
    return dict(bgg_cases=cases, arithmetic_n_max=80,
                arithmetic_positive_degree_rows=arithmetic_rows)


def equivariant_sl4():
    S = PolynomialRing(ZZ, ("t1", "t2", "t3"))
    t = tuple(S.gens()) + (-sum(S.gens()),)
    a, b, c = [sum(prod(cs) for cs in combinations(t, j)) for j in (2, 3, 4)]
    P = PolynomialRing(S, "x")
    x = P.gen()
    H = P.quotient(x**4+a*x**2-b*x+c, "h")
    h = H.gen()
    v1, v2, v3 = 4*h, 2*h**2, 4*h**3+2*a*h
    f1, f2, f3 = 4*h, 6*h**2+a, v3-b
    assert 3*f2-f1**2-3*a == v2
    assert f2 == 3*v2+a and f3 == v3-b
    # Fixed-point reconstruction of the full Levi-invariant generators.
    for j, f in enumerate((f1, f2, f3), 1):
        for i in range(4):
            direct = sum(prod(t[i]-t[k] for k in chosen)
                         for chosen in combinations([k for k in range(4) if k != i], j))
            assert f.lift()(t[i]) == direct
    identities = [
        (v1**2, 8*v2),
        (v1*v2, 2*v3-a*v1),
        (v1*v3, -4*a*v2+4*b*v1-16*c),
        (v2**2, -2*a*v2+b*v1-4*c),
    ]
    assert all(left == right for left, right in identities)
    nu = 2*a**2*h
    assert nu == v2*v3+a*v3-4*b*v2+2*c*v1
    assert 2*nu == a**2*v1
    assert any(ZZ(coefficient) % 2 for coefficient in (a**2).coefficients())

    def coordinates(element, degree):
        basis = [(i, exponents) for i in range(min(3, degree)+1)
                 for exponents in compositions(degree-i, 3)]
        f = element.lift()
        coefficients = [{tuple(e): v for e, v in f[i].dict().items()}
                        for i in range(min(3, degree)+1)]
        return vector(ZZ, [coefficients[i].get(exponents, 0) for i, exponents in basis])

    def generators(degree):
        result = []
        for e1 in range(degree+1):
            for e2 in range(degree//2+1):
                for e3 in range(degree//3+1):
                    used = e1+2*e2+3*e3
                    if used > degree:
                        continue
                    for exponents in compositions(degree-used, 3):
                        result.append(prod(z**e for z, e in zip(S.gens(), exponents))
                                      *v1**e1*v2**e2*v3**e3)
        return result

    degree_reports = []
    for degree in range(6):
        cols = generators(degree)
        J = matrix(ZZ, [coordinates(f, degree) for f in cols]).row_module()
        if degree:
            Icols = [z*f for z in S.gens() for f in generators(degree-1)]
            IJ = matrix(ZZ, [coordinates(f, degree) for f in Icols]).row_module()
            # Integral matrix of IJ in a Z-basis of J, not rational membership.
            relations = matrix(ZZ, [J.coordinate_vector(v) for v in IJ.basis()])
            nonzero = [int(abs(v)) for v in relations.elementary_divisors() if abs(v) > 1]
            free_rank = int(J.rank()-IJ.rank())
        else:
            nonzero, free_rank = [], 1
        assert free_rank == (1 if degree <= 3 else 0), (degree, J.rank(), IJ.rank() if degree else 0, free_rank)
        assert nonzero == ([2] if degree == 5 else [])
        if degree == 5:
            nv = coordinates(nu, degree)
            assert nv in J and nv not in IJ and 2*nv in IJ
        degree_reports.append(dict(polynomial_degree=degree, ambient_rank=len(coordinates(H.zero(), degree)),
                                   image_rank=int(J.rank()), specialized_image_free_rank=free_rank,
                                   specialized_image_torsion=nonzero))
    return dict(group="SL4", representation="defining", space="P3", source="Sym_Z(Q)^(S3)",
                torus_parameter_ring="Z[t1,t2,t3]; t4=-t1-t2-t3",
                a=str(a), b=str(b), c=str(c),
                quotient_relation="h^4+a*h^2-b*h+c",
                image_generators=["4*h", "2*h^2", "4*h^3+2*a*h"],
                witness="nu=2*a^2*h", witness_degree=5, kernel_degree5="Z/2",
                fixed_point_checks=12, multiplication_identities=len(identities)+2,
                graded_integer_lattice_checks=degree_reports)


def main():
    report = dict(result="passed",
                  scope="Finite exact corroboration of the separate proofs; no universal equivariant backend or novelty claim.",
                  torsion_index_lifts=[top_certificate(k) for k in ("A1", "A2", "A3", "B2", "B3", "C3", "G2")],
                  projective_family=projective_checks(),
                  equivariant_pilot=equivariant_sl4())
    (ROOT/"experiments/local-image-equivariant-pilot.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
