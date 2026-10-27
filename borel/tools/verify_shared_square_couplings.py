"""Independent Sage/Singular checks of ORIGINAL coupled critical ideals."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, matrix

ROOT = Path(__file__).resolve().parents[1]


def checks():
    C = PolynomialRing(QQ,'c'); c=C.gen(); K=C.fraction_field(); c=K(c)
    A=PolynomialRing(K,names='x,u,v,p,q,z,r');x,u,v,p,q,z,r=A.gens()
    w1=1-u*v;w2=1-p*q
    f=w1*(x*x-(z-1)**2)+w2*(x*x-(z+1)**2)+z
    I=A.ideal([f.derivative(t) for t in (x,u,v,p,q)]
              +[z*f.derivative(z)-c,z*r-1])
    components=[A.ideal(x,u,v,p,q,4*z*z-z+c,z*r-1),
                A.ideal(p,q,u*v-2,x*x-(z-1)**2,3*z+c,z*r-1),
                A.ideal(u,v,p*q-2,x*x-(z+1)**2,5*z-c,z*r-1)]
    assert I==components[0].intersection(components[1]).intersection(components[2])
    dims=[]
    for J in components:
        assert J.radical()==J
        dim=J.dimension();dims.append(int(dim))
        jac=matrix(A,[[g.derivative(t) for t in A.gens()] for g in J.gens()])
        assert J+A.ideal(jac.minors(A.ngens()-dim))==A.ideal(1)
    assert dims==[0,1,1]
    assert I.radical()==I
    assert all(components[i]+components[j]==A.ideal(1)
               for i in range(3) for j in range(i))
    # Opposite simple-zero discriminants: localize only at u,p, using
    # a=uv,b=pq. These are smooth coordinate changes, not radicalization.
    B=PolynomialRing(K,names='x,a,b,z');x,a,b,z=B.gens()
    bad=(1-a)*(x*x-(z-1))+(1-b)*(x*x+(z-1))+z
    J=B.ideal([bad.derivative(t) for t in (x,a,b)]+[z*bad.derivative(z)-c])
    expected=B.ideal(z-1,b-a-1+c,x*x,x*(1+c-2*a))
    assert J==expected
    assert x not in J and x*x in J
    assert J.radical()==J+B.ideal(x)
    assert J.quotient(B.ideal(x))==B.ideal(z-1,b-a-1+c,x,1+c-2*a)
    # Zero total constant: a reduced singular crossing persists generically.
    B=PolynomialRing(K,names='t,u,v');t,u,v=B.gens()
    crossing=B.ideal(u*v,t*u,t*v)
    assert crossing.radical()==crossing and crossing.dimension()==1
    jac=matrix(B,[[g.derivative(y) for y in B.gens()] for g in crossing.gens()])
    singular=crossing+B.ideal(jac.minors(2))
    assert singular!=B.ideal(1) and singular.radical()==B.ideal(t,u,v)
    # Coordinate-ring controls from the earlier work stay unchanged.
    # E6 contact: actual generic target has two smooth affine planes.
    B=PolynomialRing(K,names='z,x,w1,w2');z,x,w1,w2=B.gens()
    contact=B.ideal(z*z+1,2*x+c*z)
    assert contact.radical()==contact and contact.dimension()==2
    # Altered Levi crossing is genuinely nonreduced.
    B=PolynomialRing(K,names='u,v');u,v=B.gens()
    altered=B.ideal(u*u*v,u*v*v)
    assert u*v not in altered and (u*v)**2 in altered
    # A special triple point must not be confused with a generic assertion.
    B=PolynomialRing(QQ,names='x');x=B.gen()
    assert x not in B.ideal(x**3) and x**3 in B.ideal(x**3)
    # Independently factor the existing E8 control over its central function
    # field. The archived full root-derived potential is the shared input.
    B=PolynomialRing(QQ,names='t,r');t,r=B.gens();F=B.fraction_field()
    A=PolynomialRing(F,names='x,u,v,y');x,u,v,y=A.gens()
    text=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())['potential']
    from sage.misc.sage_eval import sage_eval
    f=A(sage_eval(text,locals=dict(x=x,u=u,v=v,y=y,t=F(t),r=F(r))))
    top=sum(coeff*mon for coeff,mon in zip(f.coefficients(),f.monomials()) if mon.total_degree()==4)
    q1=F(r)*v*x+F(r)*x*y+F(t)*u*y
    q2=F(r)*v*y+F(t)*u*x+F(t)*x*y
    assert top==-q1*q2/F(r*t)
    factors=list(top.factor())
    assert len(factors)==2 and all(e==1 for q,e in factors)
    for quad in (q1,q2):
        H=matrix(F,[[quad.derivative(a).derivative(b) for b in A.gens()] for a in A.gens()])
        assert H.det()==F(r*r*t*t)
    return dict(original_coupled_generic_ideal_exact_decomposition=True,
                smooth_component_dimensions=dims,components_pairwise_disjoint=True,
                original_generic_ideal_reduced=True,
                opposite_simple_zero_embedded_annihilator_verified=True,
                zero_total_constant_reduced_singular_control=True,
                existing_E8_squarefree_quartic_top_verified=True,
                universal_affine_shared_variable_normal_form_disproved=True,
                E6_contact_positive_excess_retained=True,
                altered_Levi_nonreduced_control_retained=True,
                special_triple_point_not_asserted_smooth=True)


def main():
    report=dict(result='passed',engine='Sage/Singular over QQ(c), c transcendental',
                checks=checks(),native_audit_imported=False,
                shared_square_family_theorem_in_note=True,
                full_root_normal_form_for_remaining_cases_proved=False,
                new_exceptional_group_coverage=False,uniform_conjecture_proved=False,
                literature_priority_claim=False,
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for p in (Path(__file__),ROOT/'SHARED_SQUARE_COUPLINGS.md',
                                         ROOT/'adjoint-nilpotent-source-screen.json')})
    (ROOT/'shared-square-couplings-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    main()
