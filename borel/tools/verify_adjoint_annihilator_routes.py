"""Independent ambient-root and ideal audit; imports no native audit code."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import (QQ, PolynomialRing, CyclotomicField, RootSystem,
                      matrix, vector)

ROOT=Path(__file__).resolve().parents[1]


def checks():
    base=CyclotomicField(3);omega=base.gen()
    P=PolynomialRing(base,'T');T=P.gen()
    L=base.extension(T*T+QQ(65)/16*T+1,'t');t=L.gen()
    values=[L(omega),L(omega**2),t,L(omega),L(omega**2),L(2)]
    # Independent Sage E6 root lattice, not archived branching terms.
    lattice=RootSystem(['E',6]).root_lattice()
    C=matrix(QQ,lattice.cartan_type().cartan_matrix())
    columns=[(1,0,1,0,0,0),(0,0,1,0,0,0),(0,1,0,0,0,0),
             (0,0,0,0,1,1),(0,0,0,0,0,1),(2,3,4,6,4,2)]
    M=matrix(QQ,columns).transpose();metric=M.transpose()*C*M
    assert metric==matrix(QQ,[[2,1,0,0,0,0],[1,2,0,0,0,0],
                             [0,0,2,0,0,0],[0,0,0,2,1,0],
                             [0,0,0,1,2,0],[0,0,0,0,0,6]])
    inverse=metric.inverse();drift=vector(L,[0]*6)
    total=L(6);external=L(1);all_discriminant=L(1);count=0
    for alpha in lattice.positive_roots():
        coeff=vector(QQ,[alpha.coefficient(i) for i in range(1,7)])
        charge=M.transpose()*C*coeff
        e=L(1)
        for value,n in zip(values,charge):e*=value**int(n)
        assert e!=1
        total+=e+1/e
        pair=(1-e)*(1-1/e);all_discriminant*=pair
        if coeff[3]!=0:external*=pair
        drift+=(e+1)/(e-1)*vector(L,inverse*charge)
        count+=2
    x=t+1/t;z=L(2);S=z**3+z**-3;R=z**3-z**-3
    assert count==72 and x==-S/2
    assert total==-S*S/4-2
    expected=-R**6*(S-x)**2*(S-x**3+3*x)**6
    assert external==expected and external!=0 and all_discriminant!=0
    # At both A2 trace pairs zero, the FULL derivative Phi_x has no
    # first derivatives in those four coordinates. Its torus differential
    # is 2(t-t^-1)dlog(t)+3R dlog(z).
    defect=2*(t-1/t)*drift[2]+3*R*drift[5]
    assert defect==QQ(1647)/52
    mu=3*R*x
    assert mu==-QQ(12285)/128
    assert drift[5]*mu==22*total+12
    # Recheck the entire critical branch from the full character, before
    # restricting any gradient or changing the central derivative.
    P=PolynomialRing(QQ,names='a,b,x,u,v,z');a,b,x,u,v,z=P.gens();F=P.fraction_field()
    S=z**3+F(z)**-3;R=z**3-F(z)**-3
    phi=a*b+u*v+x*x-2+S*x+(F(z)**-2+z*x)*a*v+(z*z+x/F(z))*b*u
    branch={a:F(0),b:F(0),u:F(0),v:F(0),x:-S/2}
    assert all(phi.derivative(vv).subs(branch)==0 for vv in (a,b,x,u,v))
    assert (z*phi.derivative(z)).subs(branch)==-QQ(3)/2*S*R
    assert (-QQ(3)/2*S*R).derivative(z)!=0
    # Preserve the altered-character embedded nilpotent over generic c.
    A=PolynomialRing(QQ,'c');c=A.gen();K=A.fraction_field()
    B=PolynomialRing(K,names='u,v,z');u,v,z=B.gens()
    central=z**12+c*z**6/6-1
    I=B.ideal(u*u*v,u*v*v,central)
    assert u*v not in I and (u*v)**2 in I
    assert I.radical()==B.ideal(u*v,central)
    assert B.ideal(central,central.derivative(z)).is_one()
    assert I.quotient(B.ideal(u*v))==B.ideal(u,v,central)
    # Generic E6 contact excess remains a polynomial algebra over a
    # reduced finite central algebra; z^2+1 cannot be inverted.
    contact=B.ideal(z*z+1)
    assert contact.radical()==contact and contact.dimension()==2
    assert (contact+B.ideal(1-u*(z*z+1))).is_one()
    return dict(independent_Sage_E6_roots=72,
        E6_full_root_trace_and_external_discriminant_verified=True,
        E6_root_metric_verified=True,
        E6_radial_defect_exact='1647/52',
        E6_full_critical_branch_and_dominant_target_verified=True,
        E6_regular_semisimple_source_verified=True,
        E6_radial_eigenidentity_retained=True,
        altered_Levi_etale_center_retains_embedded_nilpotent=True,
        E6_contact_positive_excess_retained=True)


def main():
    report=dict(result='passed',engine='Sage ambient E6 roots over an exact number field and Singular ideals',
        checks=checks(),native_audit_imported=False,
        exceptional_target_annihilator_derived=False,uniform_generic_reducedness_proved=False,
        adjoint_counterexample=False,new_Levi_coverage=False,novelty_claim=False,
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
          for p in (Path(__file__),ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md')})
    (ROOT/'adjoint-annihilator-routes-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['checks'],indent=2))


if __name__=='__main__':main()
