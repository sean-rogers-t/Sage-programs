"""Exact algebra supporting the all-Levi type A nilradical theorem.

Finite checks corroborate the written arbitrary-rank proof. This does
not claim new generic-smoothness coverage or an exceptional-group result.
"""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def check():
    count=0
    for m in range(2,8):
        aa=s.symbols('a:'+str(m));cc=s.symbols('c:'+str(m))
        alpha,beta,u=s.symbols('alpha beta u')
        A=sum(aa)+alpha;B=sum(c*a for c,a in zip(cc,aa))+beta
        phi=A*B-1
        F=[s.diff(phi,a) for a in aa]
        assert all(s.expand(F[i]-(B+cc[i]*A))==0 for i in range(m))
        change={aa[-1]:u-alpha-sum(aa[:-1])}
        L=sum((cc[i]-cc[-1])*aa[i] for i in range(m-1))+beta-cc[-1]*alpha
        assert s.expand(F[-1].subs(change)-(2*cc[-1]*u+L))==0
        for i in range(m-1):
            assert s.expand((F[i]-F[-1]).subs(change)-(cc[i]-cc[-1])*u)==0
        star=-L/(2*cc[-1])
        assert all(s.cancel(F[i].subs(change).subs(u,star)+(cc[i]-cc[-1])*L/(2*cc[-1]))==0 for i in range(m-1))
        # The central derivative in original coordinates factors by A on
        # incidence, independently of any chosen central cocharacter.
        dA,dB=s.symbols('dA dB')
        assert s.expand((dA*B+A*dB)-A*(dB-cc[-1]*dA)-dA*F[-1])==0
        count+=1
    # A root sum independently gives trace(V)*trace(V^-1)-1 in every
    # dimension; check the exponent dictionary, avoiding specializations.
    root_counts=[]
    for n in range(2,11):
        e=[tuple(int(i==j) for i in range(n)) for j in range(n)]
        product={}
        for v in e:
            for w in e:
                exponent=tuple(a-b for a,b in zip(v,w))
                product[exponent]=product.get(exponent,0)+1
        zero=(0,)*n;product[zero]-=1
        assert product[zero]==n-1
        assert all(v==1 for k,v in product.items() if k!=zero)
        assert len(product)-1==n*(n-1)
        root_counts.append(n*(n-1))
    x,y,z,w,a=s.symbols('x y z w a')
    S=z*z+z**-2;R=z*z-z**-2
    phi=x*x+y*y+S*x*y-1
    assert s.cancel(s.diff(phi,x).subs(y,-S*x/2)+(R*R*x)/2)==0
    mu=z*s.diff(phi,z)
    assert s.cancel(mu.subs(y,-S*x/2)+R*S*x*x)==0
    base=x*(z**4-1)**2;nil=x*(z**4-1)
    assert s.rem(nil*nil,base,x,z)==0
    assert s.rem(nil,base,x,z)!=0
    mu_num=-(z**4+1)*(z**4-1)*x*x
    assert s.rem(mu_num*nil,base,x,z)==0
    assert s.rem(mu_num,base,x,z)!=0
    assert s.rem(mu_num**2,base,x,z)==0
    # One singleton: retain ALL parameter boundaries. The Laurent ideal
    # is the principal squarefree product below, with z,w already units.
    d=w*w-z*z
    L=w*w*d*a+z**4*w**6-1
    polynomial=d*L
    gcd=polynomial
    for variable in (a,z,w):gcd=s.gcd(gcd,s.diff(polynomial,variable))
    assert gcd==1
    assert s.expand(L.subs(z,w)-(w**10-1))==0
    assert s.expand(L.subs(z,-w)-(w**10-1))==0
    # On every equal-square center component with k>=1 singleton blocks,
    # this determinant-compatible curve shows g=beta-c*alpha is not zero.
    t=s.symbols('t')
    for m in range(2,7):
        for k in range(1,5):
            c=t**(-2*k);alpha=k*t**(-2*m);beta=k*t**(2*m)
            assert s.expand(beta-c*alpha)==k*t**(2*m)-k*t**(-2*k-2*m)
            assert beta-c*alpha!=0
    return dict(full_type_A_root_sum_dimensions=list(range(2,11)),
        root_counts=root_counts,universal_linear_elimination_checks=count,
        all_boundary_coefficients_retained=True,
        SL4_actual_incidence_nonreduced=True,
        SL4_nilradical_generator='x*(z^4-1)',
        SL4_central_target_nonzero_nilpotent=True,
        SL4_target_annihilates_nilradical=True,
        SL5_singleton_incidence_squarefree=True,
        singleton_equal_square_component_nonvanishing_curves=20)


def main():
    report=dict(result='passed',engine='SymPy coefficient identities and root exponent dictionaries',
        checks=check(),all_type_A_nilradical_theorem_written=True,
        generic_type_A_reducedness_already_followed_from_quadratic_theorem=True,
        uniform_all_simple_reducedness_proved=False,new_Levi_coverage=False,
        novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'type-a-incidence-nilradical-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
