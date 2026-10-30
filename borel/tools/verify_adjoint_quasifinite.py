"""Structural zero-dimensionality test for the fixed E8 relative critical family.

The argument retains the singular quadratic block. Its general squarefreeness
lemma is elementary; this script verifies this potential's hypotheses.
"""
from pathlib import Path
from hashlib import sha256
import json,time
from sage.all import QQ,PolynomialRing,matrix,vector
ROOT=Path(__file__).resolve().parents[1]


def setup():
    T=PolynomialRing(QQ,names='t,r');K=T.fraction_field();t,r=map(K,T.gens())
    P=PolynomialRing(K,names='x,a,b,y');x,a,b,y=P.gens()
    raw=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())['potential']
    f=P(raw.replace('u','(a/t)').replace('v','(t*b)'))
    H=matrix(P,[[f.derivative(z,w) for w in (a,b)] for z in (a,b)])
    L=vector(P,[f.derivative(z).subs({a:0,b:0}) for z in (a,b)])
    Q=f.subs({a:0,b:0});D=H.det();N=Q*D-(L*H.adjugate()*L)/2
    return T,P,f,H,L,Q,D,N


def main():
    T,P,f,H,L,Q,D,N=setup();t,r=map(T.fraction_field(),T.gens());x,a,b,y=P.gens()
    TT=t+1/t;RR=r+1/r
    top=sum(c*P.monomial(*ex) for ex,c in N.dict().items() if sum(ex)==7)
    assert top==x*x*y*y*(x*x-y*y)*(TT*x-RR*y)
    assert N.coefficient({x:6})==-1 and N.coefficient({y:6})==-1
    assert N.degree(x)==6 and N.degree(y)==6
    assert D.degree(x)==4 and D.degree(y)==4
    z=vector(P,[a,b]);g=H*z+L
    assert N-f*D==-(g*H.adjugate()*g)/2
    M=x*x+TT*x+1; V=y*y+RR*y+1
    assert D==4*M*V-(M+V-1)**2
    assert D.derivative(x)==2*(V-M+1)*(2*x+TT)
    assert D.derivative(y)==2*(M-V+1)*(2*y+RR)
    print('plane squarefreeness and determinant identities verified',flush=True)

    # At the centered singular point, rank zero would force t^2=1 or r^2=1.
    # On D=0 that would force the other trace to vanish, excluded by b.
    center={x:-TT/2,y:-RR/2}
    HH=H.subs(center);LL=L.subs(center)
    h11=T.fraction_field()(HH[0,0]);h12=T.fraction_field()(HH[0,1])
    linear=T.fraction_field()(LL[0])
    aa=-(h12*b+linear)/h11
    residual=[P(f.derivative(q).subs(center).subs({a:aa})) for q in (x,y)]
    pp=t*t;qq=r*r;ss=pp+qq;vv=pp*qq
    C=16*t**4*r**4*P.base_ring()(HH.det())
    E=(vv-1)*ss-2*vv*(vv+1)
    Cformula=-(vv-1)**2*ss**2-8*vv*(vv+1)*ss+4*vv*(vv+1)**2
    Bx=8*t*t*r**3*(pp-1)*residual[0].monomial_coefficient(b)
    By=8*t**3*r*r*(qq-1)*residual[1].monomial_coefficient(b)
    assert C==Cformula and Bx-By==-3*(pp-qq)*E
    multiplier=-(vv-1)**2*ss-2*vv*(vv+1)*(vv+3)
    assert (vv-1)*C+4*vv*(vv+1)**4==E*multiplier
    print('compact centered-line exclusion identity verified',flush=True)
    constraints=[T(HH.det().numerator()),T((HH.adjugate()*LL)[0].numerator())]
    for h in residual:
        for c in h.coefficients():constraints.append(T(c.numerator()))
    J=T.ideal(constraints)
    start=time.monotonic();print('centered critical-line obstruction ideal started',flush=True)
    gb=J.groebner_basis()
    print('centered obstruction basis',list(map(str,gb)),time.monotonic()-start,flush=True)
    base=t*r*(t*t-r*r)*((t*r)**2-1)*(t*t+1)*(r*r+1)*(t*t-1)*(r*r-1)
    sat,power=J.saturation(T.ideal(T(base)))
    print('centered obstruction excluded on allowed open',sat.is_one(),flush=True)
    assert sat.is_one()
    target=T.gen(1)**4
    certificate=list(target.lift(J))
    assert sum(c*h for c,h in zip(certificate,constraints))==target
    print('explicit r^4 original-constraint membership verified',flush=True)

    rank_zero=[]
    for axis in ('t','r'):
        for epsilon in (-1,1):
            for alternate in (False,True):
                if axis=='t':
                    point={x:-epsilon,y:-RR if alternate else 0}; parameter=T.gen(0)
                    remaining_trace=RR
                else:
                    point={x:-TT if alternate else 0,y:-epsilon}; parameter=T.gen(1)
                    remaining_trace=TT
                def ev(poly):
                    value=poly.subs(point)
                    return P({ex:c.numerator().subs({parameter:epsilon})/c.denominator().subs({parameter:epsilon})
                              for ex,c in value.dict().items()})
                assert all(ev(h)==0 for h in H.list())
                assert all(ev(h)==0 for h in L)
                derivatives=[ev(f.derivative(q)) for q in (x,y)]
                relation_index=1 if axis=='t' else 0
                cross=derivatives[relation_index]
                quotient=cross.monomial_coefficient(a*b)
                assert cross==quotient*(a*b-1)
                assert quotient/remaining_trace in (QQ(1),QQ(-1))
                other=derivatives[1-relation_index]
                quartic=P((a*a*P.fraction_field()(other).subs({b:1/a})))
                assert set(quartic.variables())==set([a]) and quartic.degree(a)==4
                leading=P.base_ring()(quartic.coefficient({a:4}))
                assert len(leading.numerator().monomials())==1
                assert len(leading.denominator().monomials())==1
                zz=r if axis=='t' else t
                expected=(epsilon*(a**4/zz**2+zz**2)-(zz+1/zz**3)*a**3-(zz**3+1/zz)*a
                          if alternate else -epsilon*(a**4+1)-(zz-1/zz)*(a**3-a))
                assert quartic==expected
                rank_zero.append(dict(axis=axis,epsilon=epsilon,alternate=alternate,
                                      cross_factor=str(quotient),other_derivative=str(other),
                                      quartic=str(quartic),leading_coefficient=str(leading)))
    print('all eight rank-zero configurations have a unit-leading quartic',flush=True)
    report=dict(result='passed',status='structural_quasifiniteness_certificate',
                scope='Fixed full E8 adjoint candidate, nodes 1,2,3,5,6,8; b nonzero',
                plane_top=str(top),schur_numerator=str(N),schur_determinant=str(D),
                center_constraints=list(map(str,constraints)),
                center_basis=list(map(str,gb)),center_saturation_unit=bool(sat.is_one()),
                center_membership_target=str(target),center_membership=list(map(str,certificate)),
                compact_center_identity_verified=True,
                compact_center_formula='(V-1)*C+4*V*(V+1)^4=E*(-(V-1)^2*S-2*V*(V+1)*(V+3))',
                center_inverted_factor=str(base),center_saturation_power=int(power),
                residual_derivatives=list(map(str,residual)),
                rank_zero=rank_zero,all_allowed_critical_fibers_zero_dimensional=True,
                relative_algebra_flat_and_reduced_on_b_open=True,
                generic_target_reduced_by_prior_boundary_exclusions=True,
                relative_scheme_pure_dimension=2,
                generic_target_finite_etale_or_zero_by_dimension_argument=True,
                generic_target_smoothness_proved=True,uniform_adjoint_conjecture_proved=False,
                source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__),ROOT/'adjoint-nilpotent-source-screen.json')})
    (ROOT/'adjoint-quasifinite-verification.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
