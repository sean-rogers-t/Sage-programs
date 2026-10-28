"""Construct universal integral lifts by difference trace, not module HNF.

Research certificate builder only; not a public API. The all-rank proof is
TYPE_A_K_INTEGRAL_RECURRENCE.md. Run with sage -python.
"""
from pathlib import Path
from math import comb, gcd
from hashlib import sha256
import gzip
import json
from sage.all import ZZ, PolynomialRing, xgcd

ROOT = Path(__file__).resolve().parents[1]


def coefficient(poly, index, degree):
    ring = poly.parent()
    out = ring.zero()
    for exponent, value in poly.dict().items():
        if exponent[index] == degree:
            e = list(exponent)
            e[index] = 0
            out += value * ring.monomial(*e)
    return out


def power_sums(elementary, bound, m):
    ring = elementary[0].parent()
    result = [ring(m)]
    for r in range(1, bound + 1):
        value = sum((-1)**(k-1)*elementary[k]*result[r-k]
                    for k in range(1, min(r-1, m)+1))
        if r <= m:
            value += (-1)**(r-1)*r*elementary[r]
        result.append(ring(value))
    return result


def sparse(poly, indices):
    other = set(range(poly.parent().ngens())) - set(indices)
    assert all(all(e[i] == 0 for i in other) for e in poly.dict())
    return [[list(e[i] for i in indices), int(value)]
            for e, value in sorted(poly.dict().items())]


def build_case(m, kind, indices, terms=None):
    ring = PolynomialRing(ZZ, names=(['c%d'%i for i in range(1,m+1)]
        + ['z%d'%i for i in range(1,m+1)] + ['Y','t']))
    variables = list(ring.gens())
    c = [ring.one()] + variables[:m]
    z = [ring.one()] + variables[m:2*m]
    Y,t = variables[-2:]
    b = [comb(m,j) for j in range(m+1)]
    U = [ring.one()]
    contents = [1]
    for d in range(1,m):
        current = 0
        polynomial = ring.zero()
        for j in range(1,d+1):
            candidate = b[j]*contents[d-j]
            g,a,beta = xgcd(ZZ(current),ZZ(candidate))
            polynomial = a*polynomial + beta*z[j]*U[d-j]
            current = int(g)
        assert current == m//gcd(m,d)
        U.append(polynomial)
        contents.append(current)

    degree = sum(indices)
    if kind == 'relation':
        assert terms
        weights={sum((j+1)*e for j,e in enumerate(exponents)) for exponents,value in terms}
        assert len(weights)==1
        degree=weights.pop()
        R=sum(value*ring.prod(z[j+1]**e for j,e in enumerate(exponents))
              for exponents,value in terms)
    elif kind == 'replacement':
        assert len(indices)==1
        j=indices[0]
        R=z[j]-(b[j]//contents[j])*U[j]
    elif kind in ('basis_lower','basis_overflow'):
        assert len(indices)==2
        i,j=indices
        R=U[i]*U[j]
        if kind=='basis_lower':
            assert degree<m and (contents[i]*contents[j])%contents[degree]==0
            R-=(contents[i]*contents[j]//contents[degree])*U[degree]
        else:
            assert degree>=m
    else:
        R=ring.one()
        for j in indices:
            assert 1<=j<m
            R*=z[j]
    if kind == 'lower':
        assert degree < m
        leading=1
        for j in indices:leading*=b[j]
        common = gcd(leading,b[degree])
        R = (b[degree]//common)*R - (leading//common)*z[degree]
    elif kind=='overflow':
        assert degree>=m
    assert R
    at_binomial=variables[:]
    at_binomial[m:2*m]=[ring(x) for x in b[1:]]
    leading=int(ring.hom(at_binomial,ring)(R))
    assert degree>=m or leading==0

    G = ring.zero()
    if leading:
        q,r = divmod(degree-m,m)
        content = m//gcd(m,degree)
        assert leading % content == 0
        G = (leading//content)*z[m]**q*U[r]
    V = z[m]*G-R

    jets = [ring.one()] + [sum((-1)**k*comb(m-k,j-k)*c[k]*Y**(j-k)
                             for k in range(j+1)) for j in range(1,m+1)]
    substitution = variables[:]
    substitution[m:2*m] = jets[1:]
    evaluate_jets = ring.hom(substitution,ring)
    E = evaluate_jets(R-z[m]*G)
    at_origin = variables[:]
    at_origin[:m] = [ring.zero()]*m
    assert ring.hom(at_origin,ring)(E) == 0

    jet_at_root = z[:m]+[ring.zero()]
    translate = variables[:]
    translate[m:2*m] = [sum(comb(m-j+r,r)*jet_at_root[j-r]*t**r
                               for r in range(j+1)) for j in range(1,m+1)]
    shifted_V = ring.hom(translate,ring)(V)
    source_moments = power_sums(jet_at_root,degree,m)
    coefficient_moments = power_sums(c,degree,m)
    W = sum(coefficient(E,2*m,r)*source_moments[r]
            for r in range(degree+1))
    W += sum((-1)**r*coefficient_moments[r]*coefficient(shifted_V,2*m+1,r)
             for r in range(1,degree+1))
    assert all(sum(e[:m]) > 0 for e in W.dict())

    # A builder check, separate from the exported independent replay.
    L = Y
    eliminate = variables[:]
    eliminate[m-1] = (-1)**(m+1)*sum((-1)**k*c[k]*L**(m-k) for k in range(m))
    eliminate[m:2*m] = jets[1:m]+[ring.zero()]
    test = ring.hom(eliminate,ring)(m*R-W)
    assert test == 0
    used = list(range(2*m-1))
    return dict(m=m, kind=kind, indices=list(indices), degree=degree,
                relation=sparse(R,used), witness=sparse(W,used),
                term_count=len(W.dict()),
                source='Integral Euclid, Hasse translation, and Newton trace recurrence; no HNF or ideal membership solver.')


def main():
    # These are implementation controls for the universal recurrence, not
    # rank-by-rank mathematical coverage. Mixed ranks and long words matter.
    requests = [(3,'lower',(1,1)), (3,'overflow',(1,2)),
                (3,'overflow',(2,2)), (4,'overflow',(2,3)),
                (6,'lower',(2,2)), (6,'overflow',(2,4)),
                (6,'overflow',(5,5)), (7,'overflow',(6,6)),
                (8,'overflow',(3,5)), (9,'overflow',(4,5)),
                (10,'overflow',(4,6)), (3,'overflow',(2,2,2,2)),
                (4,'overflow',(3,3,3)), (6,'lower',(1,1,1)),
                (6,'replacement',(2,)),(6,'basis_lower',(2,3)),
                (6,'basis_overflow',(2,5)),
                (4,'relation',(), [([3,0,0],1),([1,1,0],-2),([0,0,1],-4)])]
    cases=[]
    for request in requests:
        case=build_case(*request); cases.append(case)
        print(json.dumps({k:case[k] for k in ('m','kind','indices','degree','term_count')}),flush=True)
    data=dict(result='passed',cases=cases,
        method='Explicit integral difference-trace recurrence; finite implementation controls only.',
        source_sha256={'experiments/type_a_k_integral_recurrence.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    path=ROOT/'experiments/type-a-k-integral-recurrence.json.gz'
    with gzip.GzipFile(filename=str(path),mode='wb',mtime=0) as stream:
        stream.write(json.dumps(data,separators=(',',':')).encode())


if __name__ == '__main__':
    main()
