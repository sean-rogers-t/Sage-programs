"""Polynomial spin characters: independent weight and Laurent fiber checks."""
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
import json

from sage.all import QQ, ZZ, PolynomialRing, matrix, prod, vector, QuadraticField
from sage.env import SAGE_VERSION
from verify_spin_cayley_degree import weights, gram
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]
Z=PolynomialRing(ZZ,'u'); u=Z.gen()


def examples():
    return [('mixed_1_2',u+u**2),('odd_1_3',u+u**3),
            ('even_2_4',u**2+u**4),('shifted_cube',(u+1)**3),
            ('removed_trivial',u**2-1),('mixed_removed',(u**2-1)*(u+1)),
            ('removed_square',(u**2-1)**2),('odd_removed',u**3-u)]


def mixed_weights(r,F):
    out=defaultdict(int)
    for m,a in enumerate(F.list()):
        for k,v in weights(r,m).items(): out[k]+=int(a)*v
    assert all(v>=0 for v in out.values())
    return {k:v for k,v in out.items() if v}


def parity(k):
    p=k[0]%2
    assert all(ki%2==p for ki in k)
    return p


def moment(W,zs,b):
    r=len(zs); A=b.parent()
    powers=[{e:z**e for e in {(k[i]-parity(k))//2 for k in W}}
            for i,z in enumerate(zs)]
    out=[A.zero() for _ in zs]
    for k,a in W.items():
        p=parity(k)
        value=a*b**p*prod(powers[i][(k[i]-p)//2] for i in range(r))
        for i in range(r): out[i]+=QQ(k[i])/2*value
    return vector(A,out)


def polynomial(r,F,ys):
    R=PolynomialRing(QQ,'c'); c=R.gen(); K=2**r
    h=c*F.derivative()(K*c)/F.derivative()(K)
    h_over_c=h.quo_rem(c)[0]
    P=prod(4*h*h-y*y for y in ys)-4**r*c**(2*r-2)*h_over_c**(2*r)
    return P,h


def scalar_checks():
    rows=[]
    for r in range(1,7):
        for label,F in examples():
            ys=[QQ(i+2) for i in range(r)]; P,h=polynomial(r,F,ys)
            assert P.degree()==2*r*F.degree()
            assert P.leading_coefficient()==4**r*h.leading_coefficient()**(2*r)
            assert P[0]==prod(-y*y for y in ys)-(4*h[1]**2 if r==1 else 0)
            assert P[0]!=0 and P.gcd(h).degree()==0
            assert P.gcd(P.derivative()).degree()==0
            rows.append(dict(rank=r,character=label,polynomial=str(F),degree=int(P.degree()),
                             squarefree=True,no_c_or_h_denominator_roots=True))
    return rows


def weight_checks():
    rows=[]
    for r in range(1,4):
        K=2**r
        for label,F in examples():
            W=mixed_weights(r,F); G=gram(r,W)
            assert sum(W.values())==F(K)
            assert G==matrix.identity(QQ,r)*(QQ(K)/4)*F.derivative()(K)
            aa=[QQ(i+2) for i in range(r)]; zs=[a*a for a in aa]; b=prod(aa)
            c=prod(1+z for z in zs)/(K*b)
            h=c*F.derivative()(K*c)/F.derivative()(K)
            image=G.inverse()*moment(W,zs,b)
            assert image==vector(QQ,[2*h*(z-1)/(z+1) for z in zs])
            deck=G.inverse()*moment(W,zs,-b)
            even=F(-u)==F(u); odd_nonconstant=F(-u)+F(u)==2*F[0]
            assert (deck==image)==even
            assert (deck==-image)==odd_nonconstant
            ys=[QQ(i+2) for i in range(r)]; P,hp=polynomial(r,F,ys)
            A=P.parent().quotient(P,'cb'); cb=A.gen(); hb=hp(cb)
            xs=[y/(2*hb) for y in ys]; zz=[(1+x)/(1-x) for x in xs]
            bb=prod(1+z for z in zz)/(K*cb)
            assert bb*bb==prod(zz)
            assert cb*cb*prod(1-x*x for x in xs)==1
            assert G.change_ring(A).inverse()*moment(W,zz,bb)==vector(A,ys)
            rows.append(dict(rank=r,character=label,dimension=int(F(K)),
                             fiber_algebra_dimension=int(P.degree()),
                             descends_to_SO=bool(even),central_parity_checked=True))
    return rows


def original_fiber(r,F,label,descended=False):
    W=mixed_weights(r,F); G=gram(r,W)
    if descended: assert F(-u)==F(u)
    names=['z'+str(i+1) for i in range(r)]+([] if descended else ['b'])+['v']
    R=PolynomialRing(QQ,names,order='degrevlex')
    zs=R.gens()[:r]; b=R.one() if descended else R.gen(r); v=R.gens()[-1]
    shift=(int(F.degree())+1)//2; denominator=prod(z**shift for z in zs)
    eqs=[]
    for i in range(r):
        numerator=sum(a*QQ(k[i])/2*b**parity(k)*
            prod(z**((ki-parity(k))//2+shift) for z,ki in zip(zs,k))
            for k,a in W.items())
        eqs.append(numerator-G[i,i]*(i+2)*denominator)
    if not descended: eqs.append(b*b-prod(zs))
    eqs.append(v*b*prod(zs)-1)
    gb=R.ideal(eqs).groebner_basis()
    length=staircase_length([tuple(f.lm().exponents()[0]) for f in gb])
    expected=r*int(F.degree())*(1 if descended else 2)
    assert length==expected,(r,label,descended,length)
    print('Original polynomial-spin fiber:',r,label,'SO' if descended else 'Spin',length,flush=True)
    return dict(rank=r,character=label,group='SO' if descended else 'Spin',
                degree=length,Groebner_basis_size=len(gb),
                saturation='torus coordinates only; no Cayley substitution',coefficient_field='QQ')


def collapsed_level_checks():
    rows=[]
    F=u+u**2
    for r in range(1,7):
        K=2**r; c=QQ(-1)/(2*K)
        # t+1/t=2c, all other torus parameters identity.
        A=QuadraticField(c*c-1,'v'); t=A(c)+A.gen()
        zs=[t*t]+[A.one()]*(r-1); b=t
        W=mixed_weights(r,F); G=gram(r,W)
        chi=prod(1+z for z in zs)/b
        assert chi==QQ(-1)/2
        assert F.derivative()(chi)==0
        assert G.change_ring(A).inverse()*moment(W,zs,b)==vector(A,[0]*r)
        assert c!=0 and t!=1 and t!=-1
        rows.append(dict(rank=r,character='u+u^2',spin_character_value='-1/2',
                         direct_weight_image_zero=True,noncentral_point=True))
    return rows


def main():
    scalar=scalar_checks(); trace=weight_checks(); ex=examples()
    cases=[(1,i) for i in range(8)]+[(2,i) for i in [0,1,4,5,6]]+[(3,0)]
    original=[original_fiber(r,ex[i][1],ex[i][0]) for r,i in cases]
    so=[original_fiber(r,ex[i][1],ex[i][0],True) for r,i in cases
        if ex[i][1](-u)==ex[i][1](u)]
    collapsed=collapsed_level_checks()
    paths=['SPIN_POLYNOMIAL_CAYLEY_DEGREES.md','tools/verify_spin_polynomial_degree.py',
           'SPIN_CAYLEY_DEGREES.md','tools/verify_spin_cayley_degree.py','tools/monomial_staircase.py']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        theorem='deg(theta_F)=2*floor(N/2)*deg(F) for actual polynomial spin characters',
        finite_flat_over_explicit_open=True,SO_descent_iff_F_even=True,
        operations='direct sum: maximum degree; tensor product: sum of degrees, within this family',
        scalar_fiber_cases=len(scalar),independent_weight_fiber_cases=len(trace),
        original_Spin_Laurent_Groebner_fibers=len(original),original_SO_Laurent_Groebner_fibers=len(so),
        collapsed_character_level_checks=len(collapsed),
        scalar_fibers=scalar,weight_fibers=trace,original_Spin_fibers=original,
        original_SO_fibers=so,collapsed_levels=collapsed,
        negative_coefficient_examples_actual_by='S self-dual: remove one trivial summand from S tensor S, then take actual sums and tensor products',
        novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
        scope='Finite exact corroboration of a uniform deduction from KM. No all-representation degree classification or global finiteness assertion.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'spin-polynomial-cayley-degrees-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('source_sha256','scalar_fibers','weight_fibers','original_Spin_fibers',
         'original_SO_fibers','collapsed_levels')},indent=2))


if __name__=='__main__': main()
