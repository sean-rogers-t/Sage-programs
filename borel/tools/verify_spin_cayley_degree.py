"""Exact weight and original-torus checks for Spin tensor-power degrees."""
from collections import defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

from sage.all import QQ, PolynomialRing, matrix, prod, vector
from sage.env import SAGE_VERSION
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(None)
def weights(r,m):
    """Doubled spin weights; multiplicities from tensor convolution."""
    cube=list(product([-1,1],repeat=r))
    out={(0,)*r:1}
    for _ in range(m):
        nxt=defaultdict(int)
        for k,a in out.items():
            for w in cube: nxt[tuple(x+y for x,y in zip(k,w))]+=a
        out=dict(nxt)
    return out


def gram(r,W):
    return matrix(QQ,r,r,lambda i,j:sum(a*QQ(k[i]*k[j])/4 for k,a in W.items()))


def moments(W,zs,b,m):
    r=len(zs); p=m%2
    powers=[{e:z**e for e in {(k[i]-p)//2 for k in W}} for i,z in enumerate(zs)]
    values=[b.parent().zero() for _ in zs]
    for k,a in W.items():
        term=a*b**p*prod(powers[i][(k[i]-p)//2] for i in range(r))
        for i in range(r): values[i]+=QQ(k[i])/2*term
    return vector(b.parent(),values)


def scalar_polynomial(r,m,ys):
    R=PolynomialRing(QQ,'c'); c=R.gen()
    return prod(4*c**(2*m)-y*y for y in ys)-4**r*c**(2*r*m-2)


def scalar_checks():
    rows=[]
    for r in range(1,9):
        for m in range(1,7):
            ys=[QQ(i+2) for i in range(r)]
            P=scalar_polynomial(r,m,ys)
            assert P.degree()==2*r*m and P.leading_coefficient()==4**r
            assert P[0]!=0 and P.gcd(P.derivative()).degree()==0
            assert P[0]==(prod(-y*y for y in ys)-(4 if r*m==1 else 0))
            c=P.parent().gen()
            assert P(-c)==P
            rows.append(dict(rank=r,tensor_power=m,degree=int(P.degree()),
                             squarefree=True,zero_root_excluded=True))
    return rows


def weight_checks():
    rows=[]
    for r in range(1,5):
        for m in range(1,5):
            W=weights(r,m); G=gram(r,W)
            assert sum(W.values())==2**(r*m)
            assert G==matrix.identity(QQ,r)*(QQ(m)*2**(r*m)/4)
            # Actual torus point: z_i=a_i^2, b=product(a_i).
            aa=[QQ(i+2) for i in range(r)]; zz=[a*a for a in aa]; b=prod(aa)
            c=prod(1+z for z in zz)/(2**r*b)
            direct=G.inverse()*moments(W,zz,b,m)
            assert direct==vector(QQ,[2*c**m*(z-1)/(z+1) for z in zz])
            assert G.inverse()*moments(W,zz,-b,m)==(-1)**m*direct
            # All roots at once in a finite algebra, no numerical root solving.
            ys=[QQ(i+2) for i in range(r)]; P=scalar_polynomial(r,m,ys)
            A=P.parent().quotient(P,'cbar'); cb=A.gen()
            xs=[y/(2*cb**m) for y in ys]
            zbar=[(1+x)/(1-x) for x in xs]
            bbar=prod(1+z for z in zbar)/(2**r*cb)
            assert bbar*bbar==prod(zbar)
            assert cb**2*prod(1-x*x for x in xs)==1
            assert prod(1+z for z in zbar)/(2**r*bbar)==cb
            image=G.change_ring(A).inverse()*moments(W,zbar,bbar,m)
            assert image==vector(A,ys)
            rows.append(dict(rank=r,tensor_power=m,dimension=2**(r*m),
                             weight_count=len(W),fiber_algebra_dimension=int(P.degree()),
                             trace_and_cover_verified=True,central_parity_verified=True))
    return rows


def original_torus_fiber(r,m,descended=False):
    """Use original Laurent weight moments, never the Cayley fiber equation."""
    W=weights(r,m); G=gram(r,W); p=m%2
    if descended: assert p==0
    names=['z'+str(i+1) for i in range(r)]+([] if descended else ['b'])+['u']
    R=PolynomialRing(QQ,names,order='degrevlex')
    zs=R.gens()[:r]; b=R.one() if descended else R.gen(r); u=R.gens()[-1]
    shift=(m+1)//2; denominator=prod(z**shift for z in zs)
    equations=[]
    for i in range(r):
        moment=sum(a*QQ(k[i])/2*b**p*
                   prod(z**((ki-p)//2+shift) for z,ki in zip(zs,k))
                   for k,a in W.items())
        equations.append(moment-G[i,i]*(i+2)*denominator)
    if not descended: equations.append(b*b-prod(zs))
    equations.append(u*b*prod(zs)-1)
    gb=R.ideal(equations).groebner_basis()
    length=staircase_length([tuple(f.lm().exponents()[0]) for f in gb])
    expected=r*m if descended else 2*r*m
    assert length==expected,(r,m,descended,length,expected)
    print('Original Laurent fiber:',r,m,'SO' if descended else 'Spin',length,flush=True)
    return dict(rank=r,tensor_power=m,group='SO' if descended else 'Spin',
                torus_fiber_length=length,Groebner_basis_size=len(gb),
                saturation='torus coordinates only',coefficient_field='QQ')


def main():
    scalar=scalar_checks()
    trace=weight_checks()
    cases=[(r,m) for r in [1,2] for m in range(1,5)]+[(3,1),(3,2),(4,1)]
    original=[original_torus_fiber(r,m) for r,m in cases]
    descended=[original_torus_fiber(r,m,True) for r,m in cases if m%2==0]
    paths=['SPIN_CAYLEY_DEGREES.md','tools/verify_spin_cayley_degree.py',
           'tools/monomial_staircase.py']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        theorem='deg(theta_(S tensor m))=2*floor(N/2)*m for Spin_N(C), N>=3',
        even_dimension_representation='full spin S+ direct-sum S-, not a single half-spin',
        odd_dimension_representation='irreducible spin; KM uses two equivalent copies',
        finite_flat_over_explicit_open=True,even_power_SO_degree='floor(N/2)*m',
        scalar_fiber_cases=len(scalar),independent_weight_fiber_cases=len(trace),
        original_Spin_Laurent_Groebner_fibers=len(original),
        original_SO_Laurent_Groebner_fibers=len(descended),
        scalar_fibers=scalar,independent_weight_fibers=trace,
        original_Spin_fibers=original,original_SO_fibers=descended,
        novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
        scope='Finite exact corroboration of a separate uniform proof from KM formulas. No single-half-spin, minimum-degree or global-finiteness claim.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'spin-cayley-degrees-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('source_sha256','scalar_fibers','independent_weight_fibers','original_Spin_fibers',
         'original_SO_fibers')},indent=2))


if __name__=='__main__': main()
