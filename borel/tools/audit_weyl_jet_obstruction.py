"""Native weight enumeration and exact jets for the SL3 Weyl-symmetry control.

Uses Gelfand--Tsetlin interlacing, not Sage character decompositions.
The tensor-power extension and group-fiber geometry are proved in the note.
"""
from collections import Counter
from pathlib import Path
from hashlib import sha256
import json
import sympy as s
from fractions import Fraction
from audit_affine_branching_certificates import cartan,root_data

ROOT=Path(__file__).resolve().parents[1]
SUMMANDS={(0,0):13,(1,1):4,(3,0):6,(0,3):6,(2,2):3,
          (4,1):3,(1,4):3,(6,0):1,(0,6):1,(3,3):2}


def weights(p,q):
    # GL3 partition (p+q,q,0), restricted to determinant one.
    total=p+2*q;out=Counter()
    for m1 in range(q,p+q+1):
        for m2 in range(q+1):
            for n in range(m2,m1+1):
                out[(n,m1+m2-n,total-m1-m2)]+=1
    assert sum(out.values())==(p+1)*(q+1)*(p+q+2)//2
    return out


def root_projectors():
    records=[]
    for kind in ('A2','G2','F4','E6','E7','E8'):
        C=((2,-1),(-1,2)) if kind=='A2' else cartan(kind)
        Q,positive=root_data(C)
        roots=list(positive)+[tuple(-a for a in b) for b in positive]
        def qvec(a):return tuple(sum(row[j]*a[j] for j in range(len(a))) for row in Q)
        def dot(a,b):return sum(x*y for x,y in zip(a,b))
        lengths={a:dot(a,qvec(a)) for a in roots};maximum=max(lengths.values())
        longs=[a for a in positive if lengths[a]==maximum]
        for alpha in longs:
            coroot=tuple(2*v/maximum for v in qvec(alpha))
            for beta in roots:
                k=dot(beta,coroot);assert k.denominator==1;k=int(k)
                plus=Fraction(k**4-k*k,24)+Fraction(k**3-k,12)
                minus=Fraction(k**4-k*k,24)-Fraction(k**3-k,12)
                assert plus==int(beta==alpha)
                assert minus==int(beta==tuple(-a for a in alpha))
        records.append(dict(type=kind,roots=len(roots),positive_long_roots=len(longs)))
    return records


def checks():
    a,b,x,z,p,c=s.symbols('a b x z p c')
    W=Counter()
    for (i,j),mult in SUMMANDS.items():
        for (w1,w2,w3),m in weights(i,j).items():
            W[(w1-w2,w1+w2-2*w3)]+=mult*m
    assert sum(W.values())==640
    U=a*b-1;V=a**3+b**3-3*a*b+5
    poly=s.expand(V*V-6*U*U)
    A=z*x+z**-2;B=x/z+z*z
    f=s.expand(poly.subs({a:A,b:B},simultaneous=True))
    torus=sum(m*p**i*z**j for (i,j),m in W.items())
    assert s.expand(torus-f.subs(x,p+1/p))==0
    def ch(i,j):
        return sum(m*p**(w1-w2)*z**(w1+w2-2*w3)
                   for (w1,w2,w3),m in weights(i,j).items())
    assert s.expand(ch(3,0)**2-ch(6,0)-ch(4,1)-ch(2,2)-ch(0,3))==0
    assert s.expand(ch(3,0)*ch(0,3)-ch(3,3)-ch(2,2)-ch(1,1)-1)==0
    assert s.expand(ch(3,0)*ch(1,1)-ch(4,1)-ch(2,2)-ch(3,0)-ch(1,1))==0
    assert s.expand(ch(1,1)**2-ch(2,2)-ch(3,0)-ch(0,3)-2*ch(1,1)-1)==0
    # Full S3 symmetry in the defining eigenvalues, using the actual weights.
    t1,t2=s.symbols('t1 t2')
    char=sum(mult*sum(m*t1**(w1-w3)*t2**(w2-w3)
              for (w1,w2,w3),m in weights(i,j).items())
              for (i,j),mult in SUMMANDS.items())
    assert s.expand(char-char.subs({t1:t2,t2:t1},simultaneous=True))==0
    assert s.expand(char-char.subs(t2,1/(t1*t2)))==0
    S=z**3+z**-3;R=z**3-z**-3
    expected=S**4+2*S*(S*S-6)*x**3+3*x**4+6*S*x**5+S*S*x**6
    assert s.expand(f-expected)==0
    Aunit=S*(S*S-6)+2*x+5*S*x*x+S*S*x**3
    assert s.expand(s.diff(f,x)-6*x*x*Aunit)==0
    D=lambda g:s.expand(z*s.diff(g,z))
    M=12*S**3*R
    assert s.expand(D(f).subs(x,0)-M)==0
    assert s.expand(D(M)-144*S*S*(S*S-3))==0
    # Direct traces from diagonal representation weights independently
    # reproduce the vanishing jets and the nonvanishing third jet.
    for degree in (1,2):
        trace=sum(m*(s.I**i)*z**j*i**degree for (i,j),m in W.items())
        assert s.expand(trace)==0
    trace3=sum(m*(s.I**i)*z**j*i**3 for (i,j),m in W.items())
    assert s.expand(trace3+96*s.I*S*(S*S-6))==0
    traceH=sum(m*(s.I**i)*z**j*j for (i,j),m in W.items())
    assert s.expand(traceH-M)==0
    assert sum(m*j*j for (i,j),m in W.items())==10368
    assert sum(m*i*i for (i,j),m in W.items())==3456
    assert sum(m*i*j for (i,j),m in W.items())==0
    # Generic-target polynomial is squarefree and misses all declared bad
    # coefficient loci. This is over QQ(c), not at sample target values.
    field=s.QQ.frac_field(c)
    g=s.Poly(s.expand(z**12*(M-c)),z,domain=field)
    assert g.degree()==24 and g.TC()==-12
    assert s.gcd(g,g.diff()).degree()==0
    bad=s.Poly(s.expand(z**12*S*S*(S*S-6)),z,domain=field)
    assert s.gcd(g,bad).degree()==0
    # Explicit trace normalization agrees with differentiation at identity.
    assert f.subs({x:2,z:1})==640
    assert D(D(f)).subs({x:2,z:1})==10368
    # Tensor-power identities, independently checked for several exponents;
    # their all-n proof uses the product rule and tensor trace formula.
    for n in range(1,5):
        restricted=s.expand(n*S**(4*n-4)*M)
        assert s.expand(restricted-12*n*S**(4*n-1)*R)==0
        assert s.expand(D(restricted)-144*n*S**(4*n-2)*(n*S*S-(4*n-1)))==0
    # The actual adjoint does NOT have the cancellation on this source stratum.
    adj=s.expand(U.subs({a:A,b:B},simultaneous=True))
    assert s.expand(adj-(x*x+S*x))==0
    assert s.expand(s.diff(adj,x).subs(x,0)-S)==0
    return dict(representation_dimension=640,irreducible_summands={str(k):v for k,v in SUMMANDS.items()},
                native_weight_count_with_multiplicity=640,full_Weyl_invariance_verified=True,
                trace_form_HH=10368,derived_trace_form=3456,
                first_two_derived_jets_vanish=True,third_derived_jet_nonzero=True,
                generic_target_polynomial_degree=24,generic_target_polynomial_squarefree=True,
                generic_target_avoids_all_chart_denominators=True,
                tensor_power_identities_checked=[1,2,3,4],
                adjoint_linear_jet_does_not_cancel=True,
                adjoint_long_root_jet_projectors=root_projectors())


def main():
    report=dict(result='passed',engine='Native Gelfand--Tsetlin weights and SymPy exact arithmetic',
                checks=checks(),actual_simple_group_representation=True,
                representation_is_reducible=True,representation_is_adjoint=False,
                genericity_is_within_proper_Levi_central_stratum=True,
                infinite_tensor_power_family_proved_in_note=True,
                uniform_adjoint_conjecture_disproved=False,
                new_adjoint_coverage=False,literature_priority_claim=False,
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for p in (Path(__file__),ROOT/'WEYL_JET_OBSTRUCTION.md',
                                         ROOT/'tools/audit_affine_branching_certificates.py')})
    (ROOT/'weyl-jet-obstruction-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
