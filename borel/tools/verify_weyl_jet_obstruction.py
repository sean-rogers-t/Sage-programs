"""Independent Sage representation and Singular generic-target scheme audit."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import QQ, PolynomialRing, LaurentPolynomialRing, WeylCharacterRing,RootSystem

ROOT=Path(__file__).resolve().parents[1]


def checks():
    W=WeylCharacterRing('A2',style='coroots')
    adj=W(1,1);V=W(3,0)+W(0,3)+adj+4*W(0,0)
    rho=V*V-6*adj*adj
    expected=(13*W(0,0)+4*adj+6*W(3,0)+6*W(0,3)+3*W(2,2)
              +3*W(4,1)+3*W(1,4)+W(6,0)+W(0,6)+2*W(3,3))
    assert rho==expected and rho.degree()==640
    assert all(v>0 for v in rho.coefficients())
    assert W(3,0)**2==W(6,0)+W(4,1)+W(2,2)+W(0,3)
    assert W(3,0)*W(0,3)==W(3,3)+W(2,2)+adj+W(0,0)
    assert W(3,0)*adj==W(4,1)+W(2,2)+W(3,0)+adj
    assert adj*adj==W(2,2)+W(3,0)+W(0,3)+2*adj+W(0,0)
    # Independent weight evaluation from Sage's ambient A2 weights.
    L=LaurentPolynomialRing(QQ,names='p,z');p,z=L.gens()
    torus=L(0);norm=0
    for weight,m in rho.weight_multiplicities().items():
        w1,w2,w3=weight.to_vector()
        i=int(w1-w2);j=int(w1+w2-2*w3)
        torus+=m*p**i*z**j;norm+=m*j*j
    x=p+1/p;a=z*x+z**-2;b=x/z+z*z
    assert torus==(a**3+b**3-3*a*b+5)**2-6*(a*b-1)**2
    assert norm==10368
    # Characteristic-zero generic TARGET, with every justified inverse
    # represented explicitly. No generic-source specialization is used.
    C=PolynomialRing(QQ,'c');c=C.gen();K=C.fraction_field();c=K(c)
    A=PolynomialRing(K,names='x,z,r,w');x,z,r,w=A.gens();F=A.fraction_field()
    zz=F(z);xx=F(x)
    a=zz*xx+zz**-2;b=xx/zz+zz**2
    phi=(a**3+b**3-3*a*b+5)**2-6*(a*b-1)**2
    S=zz**3+zz**-3;R=zz**3-zz**-3
    unit=S*(S*S-6)+2*xx+5*S*xx*xx+S*S*xx**3
    def polynomial(value):
        # Fractions may retain cancellable monomial denominators in Sage.
        quotient,remainder=value.numerator().quo_rem(value.denominator())
        assert remainder==0
        return A(quotient)
    U=polynomial(zz**9*unit)
    fx=polynomial(zz**9*phi.derivative(x))
    mu=polynomial(zz**12*(zz*phi.derivative(z)-c))
    target=polynomial(zz**12*(12*S**3*R-c))
    I=A.ideal(fx,mu,z*r-1,U*w-1)
    double=A.ideal(x*x,target,z*r-1,U*w-1)
    assert I==double
    assert x not in I and x*x in I
    reduced=I+A.ideal(x)
    assert I.radical()==reduced
    assert I.quotient(A.ideal(x))==reduced
    assert I.dimension()==0
    length=int(I.vector_space_dimension());reduced_length=int(reduced.vector_space_dimension())
    assert (length,reduced_length)==(48,24)
    # Separability and admissibility checked without relying on a radical
    # algorithm: Euclidean identities over the rational target field.
    Z=PolynomialRing(K,'z');z=Z.gen()
    g=12*(z**24+2*z**18-2*z**6-1)-c*z**12
    assert g.gcd(g.derivative())==1
    bad=(z**6+1)*((z**6+1)**2-6*z**6)
    assert g.gcd(bad)==1 and g(0)!=0
    projector_checks=[]
    for kind in ('A2','G2','F4','E6','E7','E8'):
        ambient=RootSystem(kind).ambient_space()
        roots=[alpha.to_vector() for alpha in ambient.roots()]
        maximum=max(a.dot_product(a) for a in roots)
        longs=[a for a in roots if a.dot_product(a)==maximum]
        for alpha in longs:
            for beta in roots:
                k=2*beta.dot_product(alpha)/maximum
                plus=(k**4-k*k)/24+(k**3-k)/12
                minus=(k**4-k*k)/24-(k**3-k)/12
                assert plus==int(beta==alpha) and minus==int(beta==-alpha)
        projector_checks.append(dict(type=kind,roots=len(roots),long_roots=len(longs)))
    return dict(representation_dimension=int(rho.degree()),positive_irreducible_decomposition=True,
                full_character_from_independent_Sage_weights=True,trace_form_HH=int(norm),
                original_localized_gradient_target_ideal_equal_double_structure=True,
                nilradical_generated_by_x=True,nilradical_square_zero_nonzero=True,
                exact_nilradical_annihilator_equals_reduction_ideal=True,
                generic_covered_chart_length=length,reduced_covered_chart_length=reduced_length,
                target_polynomial_separable=True,all_source_chart_units_justified=True,
                adjoint_long_root_jet_projectors=projector_checks)


def main():
    report=dict(result='passed',engine='Sage Weyl characters and Singular ideals over QQ(c)',
                checks=checks(),native_audit_imported=False,
                actual_simple_group_representation=True,representation_is_reducible=True,
                representation_is_adjoint=False,
                genericity_is_within_proper_Levi_central_stratum=True,
                infinite_tensor_power_family_proved_in_note=True,
                uniform_adjoint_conjecture_disproved=False,
                new_adjoint_coverage=False,literature_priority_claim=False,
                source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for p in (Path(__file__),ROOT/'WEYL_JET_OBSTRUCTION.md')})
    (ROOT/'weyl-jet-obstruction-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
