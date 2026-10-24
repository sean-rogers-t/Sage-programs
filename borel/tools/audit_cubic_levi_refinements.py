"""Native universal polynomial and root audit of the cubic refinements."""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import sympy as s
from audit_affine_branching_certificates import cartan,root_data
from audit_branching_norm_fibers import mod_rank,prod_mod

ROOT=Path(__file__).resolve().parents[1]
H=[[-2,0,-1,0,0,0],[2,3,4,6,4,2],[0,0,0,0,-1,-2]]


def universal_fricke():
    X,Y,Z,A,B,C=s.symbols('X Y Z A B C');variables=(X,Y,Z)
    domain=s.QQ.poly_ring(A,B,C)
    gradients=[2*X+Y*Z+A,2*Y+X*Z+B,2*Z+X*Y+C]
    basis=[X*X-Z*Z+(A*X-C*Z)/2,X*Y+2*Z+C,
           Y*Y-Z*Z+(B*Y-C*Z)/2,X*Z+2*Y+B,Y*Z+2*X+A,
           Z**3+C*Z*Z/2+B*X+A*Y-4*Z-2*C+A*B/2]
    computed=s.groebner(gradients,*variables,order='grevlex',domain=domain)
    assert set(computed.exprs)==set(basis)
    polys=[s.Poly(f,*variables,domain=domain) for f in basis]
    leads=[p.LM(order='grevlex').exponents for p in polys]
    assert all(p.LC(order='grevlex')==1 for p in polys)
    assert set(leads)=={(2,0,0),(1,1,0),(0,2,0),(1,0,1),(0,1,1),(0,0,3)}
    monomial=lambda e:s.prod(v**n for v,n in zip(variables,e))
    for i,j in combinations(range(6),2):
        common=tuple(max(a,b) for a,b in zip(leads[i],leads[j]))
        spoly=monomial(tuple(a-b for a,b in zip(common,leads[i])))*basis[i]
        spoly-=monomial(tuple(a-b for a,b in zip(common,leads[j])))*basis[j]
        assert s.reduced(spoly,basis,*variables,order='grevlex',domain=domain)[1]==0
    return dict(monic_basis_over_polynomial_parameter_ring=True,S_pairs_checked=15,
                free_basis=['1','X','Y','Z','Z^2'],rank=5,
                rank_is_over_central_parameters_not_target_parameters=True)


def exact_product_transfer():
    # All nonzero targets on a polynomial coefficient model, not only a numerical fiber.
    x,y,p,z,r,t,c,d,h,k=s.symbols('x y p z r t c d h k')
    S=z+t;T=p+r;phi=-(x*x+T*x)*(y*y+S*y)
    generators=[s.diff(phi,x),s.diff(phi,y),p*s.diff(phi,p)-r*s.diff(phi,r)-c,
                z*s.diff(phi,z)-t*s.diff(phi,t)-d,p*r-1,z*t-1,c*h-1,d*k-1]
    scalar=-T*T*S*S/16
    expected=[2*x+T,2*y+S,p*s.diff(scalar,p)-r*s.diff(scalar,r)-c,
              z*s.diff(scalar,z)-t*s.diff(scalar,t)-d,p*r-1,z*t-1,c*h-1,d*k-1]
    variables=(x,y,h,k,p,z,r,t,c,d)
    G=s.groebner(generators,*variables,order='grevlex',domain=s.QQ)
    E=s.groebner(expected,*variables,order='grevlex',domain=s.QQ)
    assert G==E
    return dict(product_transfer_on_both_nonzero_targets=True,radical_not_taken=True,
                model='S=z+z^-1, T=p+p^-1; universal powers handled by the written proof')


def master_and_normal_form():
    x,y,w,p,z,q=s.symbols('x y w p z q',nonzero=True)
    xs=(x,y,w);zs=(p,z,q)
    aa=[t*v+t**-2 for t,v in zip(zs,xs)];bb=[t**-1*v+t*t for t,v in zip(zs,xs)]
    chi=sum(a*b for a,b in zip(aa,bb))+s.prod(aa)+s.prod(bb)-3
    product=p*z*q;D=product+1/product;E=1/product-product
    delta=[(t**3-t**-3)/2 for t in zs]
    X=s.symbols('X0:3');U=s.symbols('U0:3')
    Xi=[v+(t**3+t**-3)/2 for v,t in zip(xs,zs)]
    expanded=sum(v*v for v in X)+D*s.prod(X)+E*sum(delta[i]*X[(i+1)%3]*X[(i+2)%3] for i in range(3))
    expanded+=D*sum(X[i]*delta[(i+1)%3]*delta[(i+2)%3] for i in range(3))
    expanded+=E*s.prod(delta)-sum(t*t for t in delta)-3
    assert s.expand(chi-expanded.subs(dict(zip(X,Xi)),simultaneous=True))==0
    # Derivative identity avoids expanding an unnecessary central rational constant.
    replacement={X[i]:(U[i]-E*delta[i])/D for i in range(3)}
    coefficients=[4*delta[(i+1)%3]*delta[(i+2)%3]-2*E*delta[i] for i in range(3)]
    for i in range(3):
        actual=D*s.diff(expanded,X[i]).subs(replacement,simultaneous=True)
        expected=2*U[i]+U[(i+1)%3]*U[(i+2)%3]+coefficients[i]
        assert s.cancel(actual-expected)==0
    # Boundary identities in Q(i)(p,z), both components checked.
    for sign in (-1,1):
        de=[s.cancel(t.subs(q,sign*s.I/(p*z))) for t in delta]
        ee=-2*sign*s.I
        Q=s.Matrix(3,3,lambda i,j:2 if i==j else ee*de[3-i-j])
        assert s.cancel(Q.det())==0
        assert s.cancel(ee*s.prod(de)-sum(t*t for t in de)-3+2)==0
    return dict(master_expansion=True,base_unit_relative_gradient_identity=True,
                both_boundary_components_checked=True)


def root_audit(record):
    A=cartan('E6');positive=root_data(A)[1]
    roots=list(positive)+[tuple(-v for v in rr) for rr in positive]
    n=record['refined_factors'];prime=record['prime'];params=record['torus_coroot_parameters']
    B=[[sum(rr[i]*rr[j] for rr in roots) for j in range(6)] for i in range(6)]
    assert s.Matrix(H)*s.Matrix(B)*s.Matrix(H).T==144*s.eye(3)
    I=[1,2,4,5] if n==2 else [1,2,4]
    cover=s.Matrix(6,6,lambda i,j:int(i==I[j]) if j<len(I) else H[j-len(I)][i])
    assert abs(cover.det())==(12 if n==2 else 24)
    vals=[prod_mod([pow(params[j],rr[j],prime) for j in range(6)],prime) for rr in roots]
    assert 1 not in vals
    target=[sum(record['central_target_coefficients'][k]*H[k][j] for k in range(n))%prime for j in range(6)]
    grad=[sum(v*rr[i] for v,rr in zip(vals,roots))%prime for i in range(6)]
    assert grad==[sum(B[i][j]*target[j] for j in range(6))%prime for i in range(6)]
    zero_roots=[]
    for rr in roots:
        central=sum(x*y for x,y in zip(rr,target))%prime
        vanishes=all(sum(rr[j]*H[k][j] for j in range(6))==0 for k in range(n))
        assert (central==0)==vanishes
        zero_roots.append(vanishes)
    hh=[[sum(v*rr[i]*rr[j] for v,rr in zip(vals,roots))%prime for j in range(6)] for i in range(6)]
    rank=mod_rank(hh,prime)+len(roots)-sum(zero_roots)
    assert rank==record['full_rank']
    return dict(refinements=n,prime=prime,full_rank=rank,kernel_dimension=78-rank,
                target_centralizer_checked=True,central_cover_degree=int(abs(cover.det())))


def main():
    data=json.loads((ROOT/'cubic-levi-refinements-verification.json').read_text())
    for name,digest in data['source_sha256'].items():
        assert sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==digest,name
    report=dict(result='passed',uses_Sage=False,universal_fricke=universal_fricke(),
                product_transfer=exact_product_transfer(),normal_form=master_and_normal_form(),
                native_root_checks=[root_audit(r) for r in data['independent_root_points']],
                all_exceptional_Levis_resolved=False,literature_priority_claim=False)
    names=['CUBIC_LEVI_REFINEMENTS.md','tools/audit_cubic_levi_refinements.py',
           'tools/audit_affine_branching_certificates.py','tools/audit_branching_norm_fibers.py',
           'cubic-levi-refinements-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'cubic-levi-refinements-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
