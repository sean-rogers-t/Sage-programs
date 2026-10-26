"""Exact tests of a uniform Jacobian identity and a failed spectral inference.

The A2 calculation tests a proposed general implication; no Levi census
or new smoothness coverage is claimed.
"""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def checks():
    a,z,x,v=s.symbols('a z x v')
    eigen=[z*a,z/a,z**-2]
    roots=[s.cancel(eigen[i]/eigen[j]) for i in range(3) for j in range(3) if i!=j]
    external=[s.cancel(eigen[i]/eigen[j]) for i,j in ((0,2),(1,2),(2,0),(2,1))]
    Phi=x*x+(z**3+z**-3)*x
    assert s.cancel(2+sum(roots)-Phi.subs(x,a+1/a))==0
    S=z**3+z**-3;R=z**3-z**-3
    Pint=1-v*(x*x-2)+v*v
    Pext=(1-v*z**3*x+v*v*z**6)*(1-v*z**-3*x+v*v*z**-6)
    P=s.expand(Pint*Pext)
    assert s.cancel(P.subs(x,a+1/a)-s.prod(1-v*r for r in roots))==0
    assert s.cancel(Pext.subs(x,a+1/a)-s.prod(1-v*r for r in external))==0
    E=s.factor(Pext.subs(v,1))
    assert s.cancel(E-(x-S)**2)==0
    chars=[z*x+z**-2,x/z+z*z]
    J=s.Matrix([[s.diff(f,x),z*s.diff(f,z)] for f in chars])
    assert s.cancel(J.det()-2*(S-x))==0
    assert s.cancel(J.det()**2-4*E)==0
    assert s.cancel(chars[0]*chars[1]-1-Phi)==0
    # Moment hierarchy from the SAME full root polynomial, using Newton's
    # recurrence rather than evaluating root sums repeatedly.
    coefficients=[s.expand(P).coeff(v,i) for i in range(7)]
    moments={}
    for m in range(1,5):
        moments[m]=s.expand(-m*coefficients[m]-sum(coefficients[j]*moments[m-j] for j in range(1,m)))
        assert s.cancel(moments[m].subs(x,a+1/a)-sum(r**m for r in roots))==0
    assert s.cancel(moments[1]-(Phi-2))==0
    phi2=(x*x-2)**2+(S*S-2)*(x*x-2)
    assert s.cancel(moments[2]+2-phi2)==0
    # Nonzero remainder on the ENTIRE generic critical incidence.
    critical=-S/2
    remainder=s.factor(s.diff(phi2,x).subs(x,critical))
    assert s.cancel(remainder+s.Rational(3,2)*S*R*R)==0
    mu=z*s.diff(Phi,z)
    assert s.cancel(mu.subs(x,critical)+s.Rational(3,2)*S*R)==0
    assert s.cancel(E.subs(x,critical)-s.Rational(9,4)*S*S)==0
    point={z:2,x:-s.Rational(65,16)}
    assert s.diff(Phi,x).subs(point)==0
    assert remainder.subs(z,2)==-s.Rational(773955,1024)
    assert mu.subs(point)==-s.Rational(12285,128)
    assert E.subs(point)==s.Rational(38025,256)
    # Internal discriminant is NOT needed for the quotient Jacobian unit.
    assert J.det().subs({x:2,z:2})==s.Rational(49,4)
    # The raw spectrum itself does not become stationary on the critical locus.
    first=s.cancel(s.diff(P,x).subs(v,0))
    assert first==0
    first_v=s.diff(P,x,v).subs(v,0)
    assert s.cancel(first_v.subs(x,critical))==0
    second_v=s.diff(P,x,v,2).subs(v,0)/2
    assert s.cancel(second_v.subs(x,critical))!=0
    # Retain the previous cubic control: it already has an etale chart with
    # f as one coordinate on an open containing the ACTUAL incidence.
    q,y,k,c=s.symbols('q y k c')
    f=(z+1/z)*(1+y*q*q)-c/k
    chart=s.Matrix([[s.diff(g,t) if t!=z else z*s.diff(g,z)
                    for t in (q,y,z)] for g in (q,y,f)])
    assert s.cancel(chart.det()-(z-1/z)*(1+y*q*q))==0
    return dict(A2_character_from_roots=True,
        A2_full_and_external_spectral_products=True,
        A2_Jacobian_square_equals_four_external_discriminant=True,
        Newton_moments_checked=4,
        A2_second_Adams_derived_gradient_remainder='-3*S*R^2/2',
        A2_exact_specialization=dict(z='2',x='-65/16',mu='-12285/128',
            external_discriminant='38025/256',second_Adams_gradient='-773955/1024'),
        criticality_propagates_to_Adams_characters=False,
        failure_on_generic_critical_incidence=True,
        internal_discriminant_retained=True,
        abstract_control_has_etale_chart_with_potential_coordinate=True)


def main():
    report=dict(result='passed',engine='SymPy root products and Newton recurrence',
        checks=checks(),uniform_Jacobian_formula_is_Steinberg_corollary=True,
        uniform_adjoint_conjecture_proved=False,new_Levi_coverage=False,
        adjoint_counterexample=False,novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-character-discriminant-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
