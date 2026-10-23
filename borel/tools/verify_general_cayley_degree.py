"""Independent E6 support audit and exact generic normal-layer examples.

No packaged exceptional characters, degree constants or old facet tables
are imported. Generic fibers below are over QQ(b1,b2), not finite fields.
"""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from sage.all import (QQ, ZZ, RootSystem, Polyhedron, matrix, vector,
                      factorial, PolynomialRing, prod, gcd)
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def e6_minimal():
    P = RootSystem(['E', 6]).weight_lattice()
    weights = [tuple(w.to_vector()) for w in P.fundamental_weights()[1].orbit()]
    assert len(weights) == 27
    assert [d for d in matrix(ZZ, weights).elementary_divisors() if d] == [1]*6
    moment = sum((vector(QQ,w).column()*vector(QQ,w).row() for w in weights),
                 matrix(QQ,6,6))
    assert moment.det() != 0
    poly = Polyhedron(vertices=weights)
    assert poly.interior_contains(vector(QQ,[0]*6))
    counts = Counter(); cones = Counter(); total = 0
    for face in poly.facets():
        V = [vector(QQ,v) for v in face.vertices()]
        hs = face.ambient_Hrepresentation()
        assert len(hs) == 1
        raw = list(hs[0]); scale = ZZ.one()
        for a in raw: scale = scale.lcm(a.denominator())
        raw = [ZZ(a*scale) for a in raw]
        common = gcd(raw[1:]); height = raw[0]/common
        if len(V) == 6:
            assert matrix([v-V[0] for v in V[1:]]).rank() == 5
            volume = abs(matrix(V).det())
            kind = 'simplex'
        else:
            assert len(V) == 10
            center = sum(V)/10
            remaining = {tuple(v) for v in V}; pairs = []
            while remaining:
                a = vector(QQ,min(remaining)); b = 2*center-a
                assert tuple(b) in remaining and a != b
                remaining.remove(tuple(a)); remaining.remove(tuple(b)); pairs.append((a,b))
            assert len(pairs) == 5
            assert matrix([a-center for a,b in pairs]).rank() == 5
            # Triangulate the facet along the first opposite pair.
            from itertools import product
            volume = sum(abs(matrix(list(pairs[0])+list(choice)).det())
                         for choice in product(*pairs[1:]))
            kind = 'crosspolytope'
        counts[(kind,int(height))] += 1
        cones[(kind,int(volume))] += 1
        total += volume
    assert counts == {('simplex',1):72, ('crosspolytope',2):27}
    assert cones == {('simplex',3):72, ('crosspolytope',32):27}
    assert total == factorial(6)*poly.volume() == 1080
    # Crosspolytope character on a finite torus cover is sum(z_i+1/z_i).
    # Its critical points have signs +/-1; five signs cannot sum to zero.
    from itertools import product
    critical_values = sorted({2*sum(s) for s in product([-1,1],repeat=5)})
    assert 0 not in critical_values
    return dict(weights=27,lattice='P(E6)',weight_lattice_index=1,
                facets=[dict(kind=k,height=h,count=n) for (k,h),n in sorted(counts.items())],
                normalized_volume=int(total),boundary_correction=0,degree=1080,
                cross_cover_critical_values=critical_values,
                independent_constructions='Sage root-system orbit; actual facet enumeration; every cone triangulated; separate polytope volume')


def generic_examples():
    B = PolynomialRing(QQ,['b1','b2']); K = B.fraction_field()
    b1,b2 = K.gens()
    R = PolynomialRing(K,['x','y','z'],order='degrevlex')
    x,y,z = R.gens(); rows=[]
    # chi=(x^m+x^-m)(y+y^-1+c)+a(x+x^-1).
    for m,c,a,expected in [(1,1,0,8),(1,2,0,6),(1,3,0,8),(2,2,0,12),(2,2,1,14)]:
        numerator=(x**(2*m)+1)*(y*y+c*y+1)+a*(x**(m+1)+x**(m-1))*y
        denominator=x**m*y
        I=R.ideal([x*numerator.derivative(x)-m*numerator-b1*denominator,
                   y*numerator.derivative(y)-numerator-b2*denominator,z*x*y-1])
        length=I.vector_space_dimension()
        assert length == expected
        assert len(I.normal_basis()) == length
        rows.append(dict(m=m,c=c,a=a,base_field='QQ(b1,b2)',generic_degree=int(length),
                         normalized_volume=8*m,boundary_correction=8*m-int(length)))
    # The reduced normal equations at y=-1: m=1 gives -b1*s;
    # m=2 has -a*s-b1*s^2 before terms of higher order.
    S=PolynomialRing(QQ,['s','u','a']); s,u,a=S.gens()
    for m in [1,2]:
        # Multiply also by y=-1+u to stay polynomial; a monomial unit.
        F=(1+s**(2*m))*u*u+a*(s**(m-1)+s**(m+1))*(-1+u)
        if m==1: F=F.subs(a=0)
        N=s*F.derivative(s)-m*F
        assert N.subs(u=0) == (0 if m==1 else a*s-a*s**3)
    return rows


def main():
    result=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
                E6_minimal=e6_minimal(),generic_character_examples=generic_examples(),
                finite_field_inference=False,novelty_claim=False,external_review=False,
                scope='Independent finite E6 geometry and characteristic-zero generic rank-two fibers; theorem proofs are in the note, not certified by this script.')
    paths=['GENERAL_CAYLEY_DEGREE.md','tools/verify_general_cayley_degree.py']
    result['source_sha256']={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths}
    (ROOT/'general-cayley-degree-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__': main()
