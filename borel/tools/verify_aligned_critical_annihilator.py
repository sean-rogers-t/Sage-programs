"""Independent exact radicals, intersections and colons for the written theorem.

Sage reconstructs the full E6 character from ambient roots. Original
gradient ideals retain all Laurent boundaries and nilpotents.
"""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, LaurentPolynomialRing, RootSystem, matrix, vector

ROOT=Path(__file__).resolve().parents[1]


def family_checks():
    records=[]
    for m,kind in [(1,'simple'),(1,'square'),(1,'mixed'),(2,'mixed'),(3,'square')]:
        A=PolynomialRing(QQ,names=['x']+['u'+str(i) for i in range(m)]
                         +['v'+str(i) for i in range(m)]+['s','t'])
        x=A.gen();us=A.gens()[1:1+m];vs=A.gens()[1+m:1+2*m];s,t=A.gens()[-2:]
        w=1-sum(u*v for u,v in zip(us,vs));M=list(us)+list(vs)
        if kind=='simple':delta=s;d=s;maximum=1
        elif kind=='square':delta=s*s;d=s;maximum=2
        else:delta=s*s*t**3;d=s*t;maximum=3
        e=delta//d;P=x*x-delta
        I=A.ideal([x*w]+[P*y for y in M])
        first=A.ideal([x]+M);second=A.ideal(P,w);third=A.ideal(x,delta)
        assert I==first.intersection(second).intersection(third)
        nil=[d*w*y for y in M];radical=I+A.ideal(nil)
        assert I.radical()==radical
        assert I.quotient(radical)==A.ideal(x,e)
        assert all(g**maximum in I for g in nil)
        if maximum>1:assert nil[0]**(maximum-1) not in I
        records.append(dict(pairs=m,delta=str(delta),radical_and_colon_verified=True,
                            exact_nilpotence_index=maximum))
    # A sharp target certificate that remembers multiplicity, not just support.
    A=PolynomialRing(QQ,names='x,u,v,z');x,u,v,z=A.gens()
    w=1-u*v;delta=(z-1)**2*(z+1)**3;d=(z-1)*(z+1)
    I=A.ideal(x*w,(x*x-delta)*u,(x*x-delta)*v)
    n=d*w*u;mu=z-w*z*delta.derivative(z)  # Q=z, D=z*d/dz
    exact=(mu-1)*(mu+1)**2
    assert exact*n in I
    assert (mu*mu-1)*n not in I
    assert (mu*mu-1)**2*n in I
    # The excluded identically-zero delta has a persistent embedded nilpotent.
    bad=A.ideal(x*w,x*x*u,x*x*v)
    assert x not in bad and x*x in bad
    assert bad.quotient(A.ideal(x))==A.ideal(x,w)
    # With Q=z, mu=z. Ann(x) intersects QQ[z] trivially; multiplying by
    # a target polynomial cannot remove the embedded component.
    assert z*x not in bad
    return records


def e6_checks():
    # Ambient Sage roots, independently of the archived character formula.
    L=LaurentPolynomialRing(QQ,names='p,q,t,r,s,z');p,q,t,r,s,z=L.gens()
    lattice=RootSystem(['E',6]).root_lattice()
    C=matrix(QQ,lattice.cartan_type().cartan_matrix())
    columns=[(1,0,1,0,0,0),(0,0,1,0,0,0),(0,1,0,0,0,0),
             (0,0,0,0,1,1),(0,0,0,0,0,1),(2,3,4,6,4,2)]
    M=matrix(QQ,columns).transpose();root_sum=L(6)
    for alpha in lattice.roots():
        coeff=vector(QQ,[alpha.coefficient(i) for i in range(1,7)])
        charge=M.transpose()*C*coeff;mon=L(1)
        for variable,n in zip(L.gens(),charge):mon*=variable**int(n)
        root_sum+=mon
    a=p+q+1/(p*q);b=1/p+1/q+p*q
    u=r+s+1/(r*s);v=1/r+1/s+r*s;x=t+1/t;S=z**3+z**-3
    phi=a*b+u*v+x*x-2+S*x+(z**-2+z*x)*a*v+(z*z+x/z)*b*u
    assert phi==root_sum
    # ORIGINAL affine character variables; inverse z is r with z*r=1.
    A=PolynomialRing(QQ,names='a,b,x,u,v,z,r');a,b,x,u,v,z,r=A.gens()
    S=z**3+r**3;R=z**3-r**3;w=1-u*v
    phi=a*b+u*v+x*x-2+S*x+(r*r+z*x)*a*v+(z*z+r*x)*b*u
    J=A.ideal([phi.derivative(y) for y in (a,b,x,u,v)]+[z*r-1])
    mu=z*phi.derivative(z)-r*phi.derivative(r)
    generators=[S*w*u,S*w*v]
    K=J+A.ideal(generators)
    assert J.radical()==K
    assert J.quotient(K)==J+A.ideal(x,S)
    assert all(g not in J and g*g in J and mu*g in J for g in generators)
    # Check explicit membership certificate, after a global polynomial change.
    a1=a+(z*z+r*x)*u;b1=b+(r*r+z*x)*v;y=x+S/2
    P=y*y-S*S/4
    assert phi-a1*b1-w*P+2 in A.ideal(z*r-1)
    assert mu+(QQ(3)/2)*S*R*w in J
    for h in (u,v):
        assert -(QQ(3)/2)*S*R*w*(S*w*h)-(
            -3*R*w*y*h*(2*y*w)+6*R*w*w*(P*h))==0
    return dict(full_E6_character_from_ambient_roots=True,
        original_Laurent_gradient_radical_verified=True,
        original_Laurent_annihilator_colon_verified=True,
        nonzero_nilradical_square_zero=True,
        central_target_annihilates_entire_nilradical=True,
        global_unit_free_coordinate_change_verified=True)


def main():
    T=PolynomialRing(QQ,'c');K=T.fraction_field()
    A=PolynomialRing(K,names='u,v,z');u,v,z=A.gens()
    contact=A.ideal(z*z+1)
    assert contact.radical()==contact and contact.dimension()==2
    altered=A.ideal(u*u*v,u*v*v)
    assert u*v not in altered and (u*v)**2 in altered
    report=dict(result='passed',engine='Sage ambient roots and Singular original-ideal operations',
        family_checks=family_checks(),E6=e6_checks(),native_audit_imported=False,
        excluded_delta_zero_embedded_control_verified=True,
        mixed_target_support_polynomial_alone_fails=True,
        E6_contact_positive_excess_retained=True,
        altered_Levi_embedded_nilpotent_retained=True,
        aligned_family_exact_annihilator_proved_in_note=True,
        all_exceptional_annihilators_proved=False,new_Levi_coverage=False,
        generic_smoothness_coverage_already_known=True,novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
            for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'aligned-critical-annihilator-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
