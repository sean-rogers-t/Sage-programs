"""Independent SymPy syzygy certificates and exact target multiplicities.

The arbitrary-pair and arbitrary-central-dimension proof is written in
the note; these finite polynomial checks corroborate it.
"""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s
from audit_affine_branching_certificates import cartan, root_data

ROOT=Path(__file__).resolve().parents[1]


def intersection(first,second,variables):
    t=s.Dummy('intersection')
    G=s.groebner([t*f for f in first]+[(1-t)*g for g in second],t,*variables)
    return [p.as_expr() for p in G.polys if not p.as_expr().has(t)]


def checks():
    x,d,e=s.symbols('x d e');counts=[]
    for m in range(1,6):
        u=s.symbols('u:'+str(m));v=s.symbols('v:'+str(m))
        w=1-sum(a*b for a,b in zip(u,v));P=x*x-d*e
        F0=x*w
        for h in u+v:
            n=d*w*h;Fi=P*h
            assert s.expand(e*n-(x*h*F0-w*Fi))==0
        # The intersection proof uses P*x in J, including all boundaries.
        assert s.expand(P*x-(P*F0+sum(x*b*(P*a) for a,b in zip(u,v))))==0
        counts.append(m)
    x,u,v,z=s.symbols('x u v z');w=1-u*v
    delta=(z-1)**2*(z+1)**3;d=(z-1)*(z+1);e=(z-1)*(z+1)**2
    P=x*x-delta;I=[x*w,P*u,P*v]
    actual=intersection(intersection([x,u,v],[P,w],(x,u,v,z)),[x,delta],(x,u,v,z))
    G=s.groebner(I,x,u,v,z)
    assert G==s.groebner(actual,x,u,v,z)
    mu=z-w*z*s.diff(delta,z);n=d*w*u
    exact=(mu-1)*(mu+1)**2
    assert G.reduce(s.expand(exact*n))[1]==0
    assert G.reduce(s.expand((mu*mu-1)*n))[1]!=0
    assert G.reduce(s.expand((mu*mu-1)**2*n))[1]==0
    assert G.reduce(s.expand(n*n))[1]!=0
    assert G.reduce(s.expand(n**3))[1]==0
    # Full E6 character, independently reconstructed from Cartan roots.
    p,q,t,r,h,z=s.symbols('p q t r h z');vars=(p,q,t,r,h,z)
    a=p+q+1/(p*q);b=1/p+1/q+p*q
    u=r+h+1/(r*h);v=1/r+1/h+r*h;x=t+1/t
    S=z**3+z**-3
    phi=a*b+u*v+x*x-2+S*x+(z**-2+z*x)*a*v+(z*z+x/z)*b*u
    positive=root_data(cartan('E6'))[1]
    roots=list(positive)+[tuple(-i for i in root) for root in positive]
    exponents=[(r[0]+r[2],r[2],r[1],r[4]+r[5],r[5],
                sum(i*j for i,j in zip(r,(2,3,4,6,4,2)))) for r in roots]
    root_character=6+sum(s.prod(v**k for v,k in zip(vars,ex)) for ex in exponents)
    assert s.expand(root_character-phi)==0
    # Global polynomial transformations, with no inversion of S or 1-u*v.
    a,b,x,u,v,z,aa,bb,X=s.symbols('a b x u v z aa bb X')
    S=z**3+z**-3;R=z**3-z**-3;w=1-u*v
    A=z**-2+z*x;B=z*z+x/z
    phi=a*b+u*v+x*x-2+S*x+A*a*v+B*b*u
    reduced=w*(x*x+S*x)-2
    assert s.expand(phi-(a+B*u)*(b+A*v)-reduced)==0
    change={a:aa-B*u,b:bb-A*v}
    changed=s.expand(phi.subs(change,simultaneous=True))
    assert s.expand(changed-aa*bb-reduced)==0
    shifted=s.expand(reduced.subs(x,X-S/2))
    assert s.expand(shifted-(w*(X*X-S*S/4)-2))==0
    # A z-dependent coordinate change changes the original central derivative
    # only by the explicitly displayed original gradient multiples.
    original_mu=z*s.diff(phi,z)
    correction=z*s.diff(phi,a).subs(change,simultaneous=True)*s.diff(change[a],z)
    correction+=z*s.diff(phi,b).subs(change,simultaneous=True)*s.diff(change[b],z)
    assert s.expand(z*s.diff(changed,z)-original_mu.subs(change,simultaneous=True)-correction)==0
    shift_correction=-z*s.diff(S,z)/2*s.diff(reduced,x).subs(x,X-S/2)
    assert s.expand(z*s.diff(shifted,z)-(z*s.diff(reduced,z)).subs(x,X-S/2)-shift_correction)==0
    mu=-s.Rational(3,2)*S*R*w;P=X*X-S*S/4;F0=2*X*w
    for c in (u,v):
        n=S*w*c
        assert s.expand(mu*n-(-3*R*w*X*c*F0+6*R*w*w*(P*c)))==0
    # S has simple zeros on Gm, so delta has exactly multiplicity two.
    assert s.gcd(z**6+1,6*z**5)==1
    return dict(universal_syzygies_pair_counts=counts,
        entire_three_ideal_intersection_verified=True,
        mixed_multiplicity_exact_target_polynomial='(c-1)*(c+1)^2',
        squarefree_support_polynomial_alone_fails=True,
        sufficient_power_and_exact_nilpotence_checked=True,
        E6_full_character_root_reconstruction=True,
        E6_global_polynomial_change_and_original_derivative_chain_rule=True,
        E6_explicit_original_ideal_membership_certificate=True,
        E6_delta_zeros_exactly_double=True)


def main():
    report=dict(result='passed',engine='SymPy explicit syzygies and elimination, separate Cartan roots',
        checks=checks(),aligned_family_exact_annihilator_proved_in_note=True,
        all_exceptional_annihilators_proved=False,new_Levi_coverage=False,
        generic_smoothness_coverage_already_known=True,novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
          for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'aligned-critical-annihilator-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['checks'],indent=2))


if __name__=='__main__':main()
