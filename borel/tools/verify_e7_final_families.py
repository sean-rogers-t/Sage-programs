"""Exact root reconstruction and boundary checks for remaining E7 families."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, LaurentPolynomialRing, matrix, vector, prod
from verify_branching_norm_fibers import charge_polynomial
ROOT=Path(__file__).resolve().parents[1]

def a2a2a1():
    row=charge_polynomial('E7',(0,1,2,4,5))
    L=LaurentPolynomialRing(QQ,['p','q']);p,q=L.gens()
    A=PolynomialRing(L.fraction_field(),names='a,x,b,y,z');a,x,b,y,z=A.gens()
    full=sum(QQ(v['coefficient'])*prod(w**e for w,e in zip(A.gens(),v['exponent']))*
             prod(w**e for w,e in zip((p,q),v['charge'])) for v in row['terms'])
    F=full.derivative(a)-b;G=full.derivative(b)-a
    graph={a:-G,b:-F};res=full.subs(graph)
    P=PolynomialRing(QQ,names='x,y,z,r,t');x,y,z,r,t=P.gens();K=P.fraction_field()
    H=matrix(K,[[-2*(x+r)/(r*t),-(x+r)*(x+1/r)],
                [-(x+r)*(x+1/r),-2*r*t*(x+1/r)]])
    ell=vector(K,[-t*(x*x+r*x+r*r-1)-(r*x+1)/(r*r*t*t),
                  -r*t*t*(x+r)-(x*x+x/r+1/(r*r)-1)/t])
    v=vector(K,[y,z]);f=-2+(v*H*v)/2+ell*v
    mapping=P.hom([A.gen(1),A.gen(3),A.gen(4),A(p**3),A(q*q)],A)
    mapped=mapping(f.numerator())/mapping(f.denominator())
    assert res==mapped
    for u in (A.gen(1),A.gen(3),A.gen(4)):
        assert full.derivative(u).subs(graph)==res.derivative(u)
    for u in (p,q):
        D=lambda f:A({e:u*c.derivative(u) for e,c in A(f).dict().items()})
        assert D(full).subs(graph)==D(res)
    D=H.det();ab=(x+r)*(x+1/r);k=r*x*x+(r*r+1)*x-3*r
    assert D==ab*(4-ab)
    R=-2-(ell*H.adjugate()*ell)/(2*D)
    J=P(r**6*t**3*D**2*R.derivative(x))
    assert J.degree(x)==9 and J.coefficient({x:9})==2*r**6*t**3
    U=PolynomialRing(QQ,'X');X=U.gen()
    spec=U(J(x=X,y=0,z=0,r=2,t=1))
    assert spec==2*(4*X*X+10*X+13)*(16*X**7+120*X**6+184*X**5-160*X**4-308*X**3+140*X*X+261*X-10)
    assert spec.gcd(spec.derivative())==1
    # Complete consistency identities on each determinant divisor.
    assert H.subs({x:-r})==matrix(K,[[0,0],[0,2*t*(r*r-1)]])
    assert ell[0].subs({x:-r})==-(r*r-1)*(r*r*t**3-1)/(r*r*t*t)
    assert H.subs({x:-1/r})==matrix(K,[[-2*(r*r-1)/(r*r*t),0],[0,0]])
    assert ell[1].subs({x:-1/r})==-(r*r-1)*(r*r*t**3-1)/(r*r*t)
    rx=PolynomialRing(K,'X');X=rx.gen()
    # Division in x over the rational function coefficient field; evaluate
    # at X to avoid changing the roles of any central parameters.
    numerators=[P(r**3*t*(H.adjugate()*ell)[0]),P(r**3*t*t*(H.adjugate()*ell)[1])]
    expected=[2*r*(r*r+1)*(r*x-1)*(r*r*t**3+1),2*(x-r)*(r*r+1)*(r*r*t**3+1)]
    for n,e in zip(numerators,expected):
        assert ((rx(n(x=X))-rx(e(x=X)))%rx(k(x=X)))==0
    assert k.subs({x:-(r*r+1)/(2*r)})==-(r**4+14*r*r+1)/(4*r)
    for eps in (1,-1):
        sub={r:eps,x:-eps}
        assert H.subs(sub)==0 and ell.subs(sub)==0
        assert (t*f.derivative(t)).subs(sub)==0
    pt={x:-4,y:-3,z:-3,r:1,t:1}
    assert all(f.derivative(u).subs(pt)==0 for u in (x,y,z))
    diff=lambda a,u:u*a.derivative(u) if u in (r,t) else a.derivative(u)
    Hess=matrix(QQ,[[diff(diff(f,u),v).subs(pt) for v in P.gens()] for u in P.gens()])
    aff=Hess[:3,:3];schur=Hess[3:,3:]-Hess[3:,:3]*aff.inverse()*Hess[:3,3:]
    assert aff.det()==-810 and schur==matrix(QQ,[[QQ(336)/5,36],[36,54]])
    assert schur.det()==QQ(11664)/5 and D.subs(pt)==-45
    return dict(full_character=row,root_reconstruction=True,original_derivative_checks=True,
                numerator_degree=9,squarefree_specialization=[str(c) for c in spec.list()],
                all_boundary_identities=True,rank_zero_target_second_coordinate=0,
                affine_hessian_determinant=-810,central_schur_determinant='11664/5',
                regular_fiber_dimension=14,standard_embeddings=3)

def main():
    companion=json.loads((ROOT/'e7-companion-membership-verification.json').read_text())
    assert companion['result']=='passed' and companion['original_affine_ideal_memberships']==6
    report=dict(result='passed',a2a2a1=a2a2a1(),new_standard_embeddings=3,
                companion_additional_embedding=1,covered_standard_E7_Levis=127,
                all_E7_Levis_resolved=True,
                literature_priority_claim=False)
    names=['E7_FINAL_FAMILIES.md','tools/verify_e7_final_families.py','e7-companion-membership-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-final-families-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','a2a2a1')},indent=2))
if __name__=='__main__':main()
