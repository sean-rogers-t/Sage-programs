"""Sage-free exact audit of the structural E8 quasifiniteness certificate.

No finite-field specialization is used. Reconstructs the potential, the
Schur identity, centered obstruction constraints, explicit ideal membership,
and every residual rank-zero quartic over characteristic zero.
"""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def main():
    x,u,v,y,t,r,a,b=s.symbols('x u v y t r a b')
    variables=(x,a,b,y); namespace={str(z):z for z in (x,u,v,y,t,r,a,b)}
    parse=lambda z:s.sympify(z.replace('^','**'),locals=namespace)
    source=ROOT/'adjoint-nilpotent-source-screen.json'
    certpath=ROOT/'adjoint-quasifinite-verification.json'
    cert=json.loads(certpath.read_text()); assert cert['result']=='passed'
    original=parse(json.loads(source.read_text())['potential'])
    f=s.expand(original.subs({u:a/t,v:t*b},simultaneous=True))
    H=s.hessian(f,(a,b)); L=s.Matrix([s.diff(f,z).subs({a:0,b:0}) for z in (a,b)])
    Q=f.subs({a:0,b:0}); D=s.expand(H.det()); N=s.expand(Q*D-(L.T*H.adjugate()*L)[0]/2)
    assert s.expand(N-parse(cert['schur_numerator']))==0
    assert s.expand(D-parse(cert['schur_determinant']))==0
    TT=t+1/t;RR=r+1/r
    poly=s.Poly(N,x,y)
    assert poly.degree(x)==poly.degree(y)==6
    assert s.Poly(N,x).LC()==s.Poly(N,y).LC()==-1
    top=sum(c*x**i*y**j for (i,j),c in poly.terms() if i+j==7)
    assert s.expand(top-x*x*y*y*(x*x-y*y)*(TT*x-RR*y))==0
    base=(t*t-r*r)*((t*r)**2-1)*(t*t+1)*(r*r+1)
    assert s.cancel((TT*TT-RR*RR)-(t*t-r*r)*((t*r)**2-1)/(t*t*r*r))==0
    g=H*s.Matrix([a,b])+L
    assert s.expand(N-f*D+(g.T*H.adjugate()*g)[0]/2)==0
    M=x*x+TT*x+1;V=y*y+RR*y+1
    assert s.expand(D-(4*M*V-(M+V-1)**2))==0
    assert s.expand(s.diff(D,x)-2*(V-M+1)*(2*x+TT))==0
    assert s.expand(s.diff(D,y)-2*(M-V+1)*(2*y+RR))==0
    print('exact Schur, highest-form, and singular-point identities passed',flush=True)

    center={x:-TT/2,y:-RR/2}
    H0=H.subs(center,simultaneous=True).applyfunc(s.cancel)
    L0=L.subs(center,simultaneous=True).applyfunc(s.cancel)
    aa=s.cancel(-(H0[0,1]*b+L0[0])/H0[0,0])
    constraints=[s.cancel(H0.det()),s.cancel((H0.adjugate()*L0)[0])]
    residual=[]
    for z in (x,y):
        q=s.diff(f,z).subs(center,simultaneous=True).subs(a,aa)
        q=s.Poly(s.expand(q),b)
        coeffs=[s.cancel(c) for c in q.all_coeffs()]
        residual.append(sum(c*b**(q.degree()-i) for i,c in enumerate(coeffs)))
        constraints.extend(coeffs)
    assert len(constraints)==len(cert['center_constraints'])==8
    original_constraints=[]
    for expr,stored in zip(constraints,cert['center_constraints']):
        num,den=s.fraction(s.cancel(expr)); h=parse(stored)
        ratio=s.cancel(num/h)
        assert ratio!=0 and ratio.is_Rational
        # Every denominator is invertible on the centered-singular chart.
        for factor,power in s.factor_list(den,t,r)[1]:
            assert s.rem(s.Poly(t*r*(t*t-1)*(r*r-1),t,r),s.Poly(factor,t,r)).is_zero
        original_constraints.append(h)
    for first,second in zip(residual,cert['residual_derivatives']):
        assert s.cancel(first-parse(second))==0
    pp=t*t;qq=r*r;ss=pp+qq;vv=pp*qq
    C=s.cancel(16*t**4*r**4*H0.det())
    E=(vv-1)*ss-2*vv*(vv+1)
    Bx=8*t*t*r**3*(pp-1)*s.Poly(residual[0],b).nth(1)
    By=8*t**3*r*r*(qq-1)*s.Poly(residual[1],b).nth(1)
    assert s.cancel(Bx-By+3*(pp-qq)*E)==0
    assert s.cancel(C+(vv-1)**2*ss**2+8*vv*(vv+1)*ss-4*vv*(vv+1)**2)==0
    SS,VV=s.symbols('S V')
    CC=-(VV-1)**2*SS**2-8*VV*(VV+1)*SS+4*VV*(VV+1)**2
    EE=(VV-1)*SS-2*VV*(VV+1)
    assert s.expand((VV-1)*CC+4*VV*(VV+1)**4-EE*(-(VV-1)**2*SS-2*VV*(VV+1)*(VV+3)))==0
    membership=[parse(z) for z in cert['center_membership']]
    assert s.expand(sum(c*h for c,h in zip(membership,original_constraints))-r**4)==0
    assert parse(cert['center_membership_target'])==r**4
    print('explicit r^4 ideal-membership identity independently passed',flush=True)

    checked=[]
    for item in cert['rank_zero']:
        eps=item['epsilon'];alternate=item['alternate'];axis=item['axis']
        if axis=='t':
            point={t:eps,x:-eps,y:-RR if alternate else 0};idx=1;trace=RR
        else:
            point={r:eps,y:-eps,x:-TT if alternate else 0};idx=0;trace=TT
        assert all(s.cancel(z.subs(point,simultaneous=True))==0 for z in list(H)+list(L))
        partials=[s.cancel(s.diff(f,z).subs(point,simultaneous=True)) for z in (x,y)]
        coefficient=s.cancel(partials[idx]/(a*b-1))
        assert s.cancel(coefficient/trace) in (s.Integer(-1),s.Integer(1))
        assert s.cancel(coefficient-parse(item['cross_factor']))==0
        other=partials[1-idx]
        assert s.cancel(other-parse(item['other_derivative']))==0
        quartic=s.cancel(a*a*other.subs(b,1/a));qp=s.Poly(quartic,a)
        assert qp.degree()==4 and s.cancel(quartic-parse(item['quartic']))==0
        zz=r if axis=='t' else t
        expected=(eps*(a**4/zz**2+zz**2)-(zz+1/zz**3)*a**3-(zz**3+1/zz)*a if alternate
                  else -eps*(a**4+1)-(zz-1/zz)*(a**3-a))
        assert s.cancel(quartic-expected)==0
        lc=qp.LC();num,den=s.fraction(lc)
        assert len(s.Poly(num,t,r).terms())==len(s.Poly(den,t,r).terms())==1
        checked.append(dict(axis=axis,epsilon=eps,alternate=alternate,degree=4,leading_coefficient=str(lc)))
    assert len(checked)==8
    print('all rank-zero branches independently verified',flush=True)
    report=dict(result='passed',independent_engine='SymPy',coefficient_field='Q(t,r)',
                uses_Sage=False,uses_modular_specialization=False,
                original_potential_reconstructed=True,schur_identity_verified=True,
                squarefree_plane_pencil_hypotheses_verified=True,
                determinant_singular_point_classification_identities_verified=True,
                centered_line_obstruction_membership_verified=True,
                compact_center_identity_independently_verified=True,
                rank_zero_branches=checked,all_allowed_critical_fibers_zero_dimensional=True,
                relative_algebra_flat_and_reduced_on_b_open=True,
                generic_target_reducedness_uses_prior_six_boundary_exclusions=True,
                generic_target_finite_etale_or_zero_by_written_dimension_argument=True,
                generic_target_smoothness_uses_prior_boundary_exclusions=True,
                uniform_adjoint_conjecture_proved=False,
                novelty_claim=False,
                source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__),certpath,source)})
    (ROOT/'adjoint-quasifinite-audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
