"""Independent native audit, including universal nonzero-target ideal equality.

No Sage or elimination helper is imported. The character certificate is
rebuilt from Cartan reflections and exterior subset weight models.
"""
from hashlib import sha256
from pathlib import Path
import json
import sympy as s
from audit_affine_branching_certificates import cartan, root_data
from audit_branching_norm_fibers import reconstruct, mod_rank, prod_mod

ROOT=Path(__file__).resolve().parents[1]


def symbolic():
    a,b,x,p,q,r,v,z=s.symbols('a b x p q r v z')
    phi=a*b+x*x+p*v-2+z*a*x*v+z*z*b*r+z**3*x*q+z**4*a*p+z**5*x*b+z**6*v+x*b*p/z+a*q/z**2+x*r/z**3+b*v/z**4+a*x/z**5+p/z**6
    # Solve the original derivative equations without using the helper.
    graph=s.solve([s.diff(phi,q),s.diff(phi,r),s.diff(phi,a),s.diff(phi,b)],(a,b,q,r),dict=True)[0]
    reduced=s.expand(phi.subs(graph,simultaneous=True))
    pair=s.solve([s.diff(reduced,p),s.diff(reduced,v)],(p,v),dict=True)[0]
    graph={key:s.expand(value.subs(pair,simultaneous=True)) for key,value in graph.items()};graph.update(pair)
    psi=s.expand(phi.subs(graph,simultaneous=True))
    expected=-(x*x-1)*(x*x+(z**15+z**-15)*x)-3
    assert s.expand(psi-expected)==0
    for var in (a,b,p,q,r,v):
        assert s.expand(s.diff(phi,var).subs(graph,simultaneous=True))==0
    assert s.expand(s.diff(phi,x).subs(graph,simultaneous=True)-s.diff(psi,x))==0
    assert s.expand((z*s.diff(phi,z)).subs(graph,simultaneous=True)-z*s.diff(psi,z))==0
    weights={a:-10,b:-20,x:-15,p:-24,q:-18,r:-12,v:-6}
    variables=(a,b,x,p,q,r,v,z)
    for expression,weight in [(phi,0)]+[(value,weights[key]) for key,value in graph.items()]:
        # A monomial Laurent shift suffices; no rational-function simplification.
        polynomial=s.Poly(s.expand(expression*z**60),*variables)
        for exp,coeff in polynomial.terms():
            assert (sum(exp[i]*weights[var] for i,var in enumerate(variables[:-1]))+exp[-1]-60-weight)%30==0
    Y,M,N,c,h,k=s.symbols('Y M N c h k')
    B=1-Y;C=1-3*Y
    F=900*Y*B**2*(4*Y-1)*(Y*Y-3*Y+1)-c*c*C*C
    # Invariant presentation after mu15: M=x*z^15, N=x*z^-15, MN=Y.
    # Only the TARGET c is inverted. No source divisor is saturated away.
    I=s.groebner([M*N-Y,C*(M+N)-2*Y*(2*Y-1),15*B*(M-N)-c,c*h-1],M,N,Y,h,c,order='grevlex',domain=s.QQ)
    J=s.groebner([F,30*B*C*M-30*B*Y*(2*Y-1)-c*C,
                  30*B*C*N-30*B*Y*(2*Y-1)+c*C,c*h-1],M,N,Y,h,c,order='grevlex',domain=s.QQ)
    assert I==J
    # All potentially divided source expressions are already units.
    for value in (0,1,s.Rational(1,3)):
        test=s.groebner([F,c*h-1,Y-value],Y,h,c,domain=s.QQ)
        assert list(test)==[1]
    P=900*Y*(1-Y)**2*(4*Y-1)*(Y*Y-3*Y+1)-k*(1-3*Y)**2
    delta=4096*k**5+4485075*k**4+1691730000*k**3+292161500000*k*k+23287500000000*k+683437500000000
    assert s.expand(s.discriminant(P,Y)+30611001600000000*k*delta)==0
    assert s.gcd(delta,s.diff(delta,k))==1 and delta.subs(k,0)!=0
    # Ramification locations of the rational target map, a second route.
    num=P.subs(k,0);den=(1-3*Y)**2
    ram=48*Y**5-165*Y**4+162*Y**3-70*Y**2+14*Y-1
    assert s.expand(s.diff(num,Y)*den-num*s.diff(den,Y)-900*(Y-1)*(3*Y-1)*ram)==0
    assert s.gcd(ram,s.diff(ram,Y))==1
    assert s.rem(s.resultant(P,ram,Y),delta,k)==0
    return dict(original_derivatives_and_reduced_potential=True,central_mu30_action_and_graph=True,
                universal_nonzero_target_ideal_equality=True,only_target_inverted=True,
                no_hidden_roots_at_divided_coordinates=True,discriminant_verified=True,
                independent_rational_ramification_identity=True,
                degree=6,exceptional_nonzero_targets=10,regular_orbit_dimension=28)


def root_audit(row):
    A=cartan('E8');positive=root_data(A)[1]
    roots=list(positive)+[tuple(-t for t in rr) for rr in positive]
    H=[10,15,20,30,24,18,12,6];prime=row['prime'];params=row['torus_coroot_parameters']
    z=row['central_parameter'];x=row['x'];c=row['target_c']
    B=[[sum(rr[i]*rr[j] for rr in roots) for j in range(8)] for i in range(8)]
    assert sum(H[i]*B[i][j]*H[j] for i in range(8) for j in range(8))==1800
    values=[prod_mod([pow(t,e,prime) for t,e in zip(params,rr)],prime) for rr in roots]
    assert 1 not in values
    target=[c*pow(1800,-1,prime)*h%prime for h in H]
    gradient=[sum(v*rr[i] for v,rr in zip(values,roots))%prime for i in range(8)]
    assert gradient==[sum(B[i][j]*target[j] for j in range(8))%prime for i in range(8)]
    Hess=[[sum(v*rr[i]*rr[j] for v,rr in zip(values,roots))%prime for j in range(8)] for i in range(8)]
    assert mod_rank(Hess,prime)==8
    outside=sum(sum(rr[i]*target[i] for i in range(8))%prime!=0 for rr in roots)
    assert outside==212
    Y=x*x%prime
    assert (900*Y*(1-Y)**2*(4*Y-1)*(Y*Y-3*Y+1)-c*c*(1-3*Y)**2)%prime==0
    assert (15*x*(1-x*x)*(pow(z,15,prime)-pow(z,-15,prime))-c)%prime==0
    return dict(prime=prime,all_240_roots_reconstructed=True,trace_norm=1800,
                original_trace_target_verified=True,Cartan_rank=8,root_rank=outside,
                full_rank=220,kernel_dimension=28,no_source_root_value_equals_one=True)


def main():
    source=json.loads((ROOT/'multiplier-branching-verification.json').read_text())
    report=dict(result='passed',uses_Sage=False,
                independent_character=reconstruct(source['charge_record']),
                symbolic=symbolic(),root_point=root_audit(source['root_point']),
                all_exceptional_Levis_resolved=False,full_map_generic_degree_claim=False,
                literature_priority_claim=False)
    names=['MULTIPLIER_BRANCHING_FIBERS.md','tools/audit_multiplier_branching_fibers.py',
           'tools/audit_affine_branching_certificates.py','tools/audit_branching_norm_fibers.py',
           'multiplier-branching-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'multiplier-branching-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
