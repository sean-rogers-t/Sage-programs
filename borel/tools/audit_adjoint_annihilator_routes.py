"""Exact audit of failed radial descent and central-norm annihilator routes.

No exceptional nilradical annihilator or new smoothness result is asserted.
"""
from pathlib import Path
from hashlib import sha256
from collections import Counter
import json
import sympy as s
from audit_affine_branching_certificates import cartan, root_data

ROOT=Path(__file__).resolve().parents[1]


def checks():
    p,q,t,r,w,z=s.symbols('p q t r w z')
    variables=(p,q,t,r,w,z)
    a=p+q+1/(p*q);b=1/p+1/q+p*q
    u=r+w+1/(r*w);v=1/r+1/w+r*w
    x=t+1/t;S=z**3+z**-3;U=z**3-z**-3
    phi=a*b+u*v+x*x-2+S*x+(z**-2+z*x)*a*v+(z*z+x/z)*b*u
    # Reconstruct all 72 roots independently of the polynomial formula.
    C=cartan('E6');positive=root_data(C)[1]
    roots=list(positive)+[tuple(-i for i in alpha) for alpha in positive]
    exponents=[(h[0]+h[2],h[2],h[1],h[4]+h[5],h[5],
                sum(i*j for i,j in zip(h,(2,3,4,6,4,2)))) for h in roots]
    expected=Counter(exponents);expected[(0,)*6]=6
    actual=Counter()
    for term in s.expand(phi).as_ordered_terms():
        powers=term.as_powers_dict();ex=tuple(int(powers.get(a,0)) for a in variables)
        actual[ex]+=term/s.prod(a**i for a,i in zip(variables,ex))
    assert actual==expected
    # At both SL3 trace pairs zero, their eigenvalues are cube roots of 1.
    # Each quotient eigenvalue therefore occurs three times.
    xx,ss,uu=s.symbols('x S U')
    T3=xx**3-3*xx
    E=-(ss-xx)**2*uu**6*(ss-T3)**6
    logx=2/(xx-ss)+18*(1-xx*xx)/(ss-T3)
    logz=6*uu/(ss-xx)+18*ss/uu+18*uu/(ss-T3)
    assert s.cancel(s.diff(E,xx)/E-logx)==0
    assert s.cancel((3*uu*s.diff(E,ss)+3*ss*s.diff(E,uu))/E-logz)==0
    # The cover metric has SL2 norm 2 and central norm 6.
    radial_F=2*xx+(xx*xx-4)*logx+uu*logz/2
    # Replace U^2 by S^2-4 using a polynomial remainder, including denominators.
    def relation_reduce(f):
        n,d=s.fraction(s.cancel(f));rel=uu**2-ss**2+4
        return s.factor(s.rem(n,rel,uu)/s.rem(d,rel,uu))
    defect=relation_reduce(radial_F.subs(xx,-ss/2))
    assert s.cancel(defect-2*(ss*ss+320)/(3*ss))==0
    assert relation_reduce(E.subs(xx,-ss/2))==-s.Rational(9,2**20)*ss**8*(ss-2)**9*(ss+2)**9
    mu=-s.Rational(3,2)*ss*uu
    radial_phi=relation_reduce((logz/6*mu).subs(xx,-ss/2))
    assert s.cancel(radial_phi-(22*(-ss*ss/4-2)+12))==0
    assert defect.subs(ss,s.Rational(65,8))==s.Rational(1647,52)
    # The restricted target map is nonconstant: every nonzero p(mu)
    # stays nonzero in C(z); this is the all-polynomials argument in the note.
    mu_z=-s.Rational(3,2)*S*U
    assert s.factor(z*s.diff(mu_z,z))!=0
    # Smaller actual A2 control, with its full discriminant vector field.
    ax=xx+(xx*xx-4)/(xx-ss);bz=-uu/(xx-ss)
    a2defect=relation_reduce((2*ax+3*uu*bz).subs(xx,-ss/2))
    assert s.cancel(a2defect-2*(ss-2)*(ss+2)/(3*ss))==0
    # A generic finite central algebra is not a finite whole fiber algebra.
    # E6 contact has z^2+1=0: norm of the tempting source localization is zero.
    Z=s.Matrix([[0,-1],[1,0]])
    assert Z*Z+s.eye(2)==s.zeros(2)
    # An etale center does not control nilpotents along excess variables.
    c=s.symbols('c');h=z**12+c*z**6/6-1
    assert s.factor(s.resultant(h,s.diff(h,z),z))!=0
    return dict(E6_full_root_character_reconstructed=True,E6_root_count=72,
        E6_external_discriminant_slice='-R^6*(S-x)^2*(S-x^3+3*x)^6',
        E6_radial_gradient_defect='2*(S^2+320)/(3*S)',
        E6_exact_defect_at_z2='1647/52',
        E6_radial_eigenidentity_retained=True,
        E6_branch_target_dominant=True,
        no_nonzero_target_multiple_repairs_radial_descent=True,
        A2_radial_gradient_defect='2*(S^2-4)/(3*S)',
        E6_contact_source_norm_zero=True,
        altered_Levi_central_polynomial_generically_squarefree=True)


def main():
    report=dict(result='passed',engine='SymPy and independent root exponent reconstruction',
        checks=checks(),exceptional_target_annihilator_derived=False,
        uniform_generic_reducedness_proved=False,adjoint_counterexample=False,
        new_Levi_coverage=False,novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
          for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-annihilator-routes-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['checks'],indent=2))


if __name__=='__main__':main()
