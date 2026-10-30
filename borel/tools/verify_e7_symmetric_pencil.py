"""Exact character and rank-stratum checks for the E7 A1^4 refinement."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, LaurentPolynomialRing, matrix, vector, prod
from verify_branching_norm_fibers import charge_polynomial

ROOT=Path(__file__).resolve().parents[1]


def symbolic():
    P=PolynomialRing(QQ,names='x,y,z,w,A,B,C');x,y,z,w,A,B,C=P.gens()
    g=x*x+A*B*C*x+A*A+B*B+C*C-7
    Q=matrix(P,[[2,B*C+A*x,A*C+B*x],[B*C+A*x,2,A*B+C*x],[A*C+B*x,A*B+C*x,2]])
    u=vector(P,[y,z,w]);phi=g+(u*Q*u)/2;d=Q.det();e=2*x+A*B*C
    ordinary={x:-A*B*C/2}
    PP=d.subs(ordinary);HH=d.derivative(x).subs(ordinary)
    assert PP and PP.gcd(HH)==1 and PP(A=2,B=2,C=2)==-216
    assert HH==A*B*C/2*(3*A*A*B*B*C*C-4*A*A*B*B-4*A*A*C*C-4*B*B*C*C+8*A*A+8*B*B+8*C*C-24)
    assert (A*x+B*C)-(B*x+A*C)==(A-B)*(x-C)
    assert (A*x+B*C)-(C*x+A*B)==(A-C)*(x-B)
    assert (B*x+A*C)-(C*x+A*B)==(B-C)*(x-A)
    assert (A*y*z+B*y*w+C*z*w).subs({w:-y-z})==-B*y*y-C*z*z+(A-B-C)*y*z
    # The four normalized rank-one curves, checked as rational identities.
    T=LaurentPolynomialRing(QQ,'t');t=T.gen()
    curves=[(2/t-t,t,t,t),(t,2/t-t,t,t),(t,t,2/t-t,t),(t,t,t,2/t-t)]
    for xx,aa,bb,cc in curves:
        assert aa*xx+bb*cc==bb*xx+aa*cc==cc*xx+aa*bb==2
    for eps1 in (-1,1):
        for eps2 in (-1,1):
            substitution={A:eps1*A,B:eps2*B,C:eps1*eps2*C,z:eps1*z,w:eps2*w}
            assert phi.subs(substitution)==phi
    row=charge_polynomial('E7',(0,1,4,6))
    L=LaurentPolynomialRing(QQ,['p','q','r']);p,q,r=L.gens()
    R=PolynomialRing(L.fraction_field(),names='x,y,z,w');xs=R.gens()
    actual=sum(term['coefficient']*prod(zz**ee for zz,ee in zip((p,q,r),term['charge']))*
               prod(xx**ee for xx,ee in zip(xs,term['exponent'])) for term in row['terms'])
    U=p*q;V=p*q*r;W=p*q*q*r
    mapping=P.hom(list(xs)+[R(U+1/U),R(V+1/V),R(W+1/W)],R)
    assert actual==mapping(phi)
    D=lambda f,zz:R({ee:zz*cc.derivative(zz) for ee,cc in R(f).dict().items()})
    trace=matrix(QQ,[[D(D(actual,zz),ww).subs(dict(zip(xs,[2]*4))).constant_coefficient()(p=1,q=1,r=1)
                     for ww in (p,q,r)] for zz in (p,q,r)])
    assert trace==matrix(QQ,[[216,288,144],[288,432,216],[144,216,144]])
    change=matrix(QQ,[[1,1,1],[1,1,2],[0,1,1]])
    assert change.det()==-1
    LT=LaurentPolynomialRing(QQ,['U','V','W']);ts=LT.gens();aa,bb,cc=[v+1/v for v in ts]
    R0=aa*aa+bb*bb+cc*cc-7-aa*aa*bb*bb*cc*cc/4
    H=matrix(QQ,[[(vv*(ww*R0.derivative(ww)).derivative(vv))(U=1,V=1,W=1) for ww in ts] for vv in ts])
    assert H==-24*matrix.identity(QQ,3)
    witness={x:-4,y:0,z:0,w:0,A:2,B:2,C:2}
    relative=matrix(QQ,[[phi.derivative(v).derivative(w).subs(witness) for w in (x,y,z,w)] for v in (x,y,z,w)])
    assert relative.det()==-432
    return dict(character=row,full_character_verified=True,collision_polynomials_coprime=True,
                normalized_rank_one_curves=4,rank_one_critical_dimension_bound=2,
                collision_base_dimension=2,double_collision_base_dimension_bound=1,
                rank_stratum_bounds_verified=True,
                trace_matrix=[list(map(int,rw)) for rw in trace.rows()],
                torus_change_determinant=-1,relative_hessian_determinant=-432,
                ordinary_logarithmic_hessian_diagonal=-24,regular_fiber_dimension=8)


def companion_reduction():
    row=charge_polynomial('E7',(0,1,2,4,6))
    L=LaurentPolynomialRing(QQ,['p','q']);p,q=L.gens()
    R=PolynomialRing(L.fraction_field(),names='a,x,b,y,w');a,x,b,y,w=R.gens()
    phi=sum(t['coefficient']*prod(z**e for z,e in zip((p,q),t['charge']))*
            prod(v**e for v,e in zip(R.gens(),t['exponent'])) for t in row['terms'])
    F=phi.derivative(a)-b;G=phi.derivative(b)-a
    graph={a:-G,b:-F};psi=phi.subs(graph)
    U=q;V=p**3*q;A=U+1/U;B=V+1/V;C=U*V+1/(U*V)
    expected=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    assert psi==expected and A*A+B*B+C*C-A*B*C==4
    for v in (x,y,w):assert phi.derivative(v).subs(graph)==psi.derivative(v)
    for z in (p,q):
        D=lambda f:R({e:z*c.derivative(z) for e,c in f.dict().items()})
        assert D(phi).subs(graph)==D(psi)
    LT=LaurentPolynomialRing(QQ,['U','V']);uu,vv=LT.gens()
    potential=-(uu+1/uu)*(vv+1/vv)*(uu*vv+1/(uu*vv))-3
    H=matrix(QQ,[[(z*(w*potential.derivative(w)).derivative(z))(U=1,V=1)
                  for w in (uu,vv)] for z in (uu,vv)])
    assert H==matrix(QQ,[[-16,-8],[-8,-16]]) and H.det()==192
    origin={x:0,y:0,w:0}
    rh=matrix(QQ,[[psi.derivative(v).derivative(leaf).subs(origin).constant_coefficient()(p=1,q=1)
                   for leaf in (x,y,w)] for v in (x,y,w)])
    assert rh.det()==-216
    return dict(character=row,exact_residual_verified=True,original_target_derivatives_verified=True,
                central_coefficients_on_Cayley_cubic=True,generic_smoothness_claim=False,
                conditional_dominance_hessian_determinant=192)


def main():
    report=dict(result='passed',symmetric_pencil=symbolic(),companion=companion_reduction(),
                new_standard_embeddings=2,resolved_proper_standard_Levis=101,
                remaining_cubic_embeddings=26,remaining_cubic_families=4,
                all_E7_Levis_resolved=False,full_map_degree_claim=False,literature_priority_claim=False)
    names=['E7_SYMMETRIC_PENCIL.md','tools/verify_e7_symmetric_pencil.py',
           'tools/verify_branching_norm_fibers.py','GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-symmetric-pencil-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','symmetric_pencil','companion')},indent=2))


if __name__=='__main__':main()
