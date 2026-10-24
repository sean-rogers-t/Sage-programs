"""Native independent root/subset and exact-ideal audit of the pencil-scope result."""
from hashlib import sha256
from pathlib import Path
import json
import sympy as s
from audit_affine_branching_certificates import cartan, root_data, inverse
from audit_branching_norm_fibers import reconstruct, mod_rank, prod_mod

ROOT=Path(__file__).resolve().parents[1]


def diagram_audit(kind):
    A=cartan(kind);n=len(A);neighbors={i:{j for j in range(n) if j!=i and A[i][j]} for i in range(n)}
    found=[]
    for mask in range(1,2**n-1):
        selected={i for i in range(n) if mask&(1<<i)}
        blocks=[];remaining=set(selected)
        while remaining:
            block={min(remaining)}
            while True:
                grown=block|set.union(*(neighbors[i]&selected for i in block))
                if grown==block:break
                block=grown
            remaining-=block;blocks.append(block)
        if len(blocks)<2 or sum(len(b)==1 for b in blocks)!=1:continue
        if any(len(neighbors[i]&b)>2 for b in blocks for i in b):continue
        leaves=[b for b in blocks if len(b)>1]
        if all(sum(bool(neighbors[j]&b) for b in leaves) in (0,2) for j in set(range(n))-selected):
            found.append(sorted(i+1 for i in selected))
    assert found==[[i for i in range(1,n+1) if i!=4]]
    return dict(type=kind,candidates=found)


def symbolic():
    a,b,x,u,v,w,z,c=s.symbols('a b x u v w z c')
    phi=a*b+u*w+x*x-2+x*(z*a*u+b*w/z)+v*(z*z*b+a/z**2)+x*(z**3*w+u/z**3)+z**4*a+b/z**4
    af=-z**4*b;vf=-z*z*(b+x*z*u+z**4)
    reduced=s.expand(phi.subs({a:af,v:vf},simultaneous=True))
    S=z**6-z**-6;y=z*z*b
    up=u+x*(b/z+z**3);wp=w+x*(-z**5*b+z**-3)
    assert s.expand(reduced-up*wp-(x*x-1)*(y*y+S*y)+2)==0
    assert s.expand(s.diff(phi,a).subs(v,vf))==0
    assert s.expand(s.diff(phi,v).subs(a,af))==0
    values={x:0,u:0,w:0,a:z*z*S/2,b:-S/(2*z*z),v:-(z**6+z**-6)/2}
    for t in (a,b,x,u,v,w):
        assert s.expand(s.diff(phi,t).subs(values,simultaneous=True))==0
    assert s.expand((z*s.diff(phi,z)).subs(values,simultaneous=True)-3*(z**12-z**-12))==0
    W=s.symbols('W');polynomial=3*W*W-c*W-3
    assert s.discriminant(polynomial,W)==c*c+36
    # Central-kernel action on the original characters, not only z^12.
    weights={a:-4,b:-8,x:-6,w:-9,v:-6,u:-3}
    polynomial_phi=s.Poly(s.expand(phi*z**4),a,b,x,u,v,w,z)
    for exponents,coeff in polynomial_phi.terms():
        assert (sum(k*weights[t] for t,k in zip((a,b,x,u,v,w),exponents[:6]))+exponents[-1]-4)%12==0
    # Uniform target localization, retaining all nilpotents at exceptional targets.
    X,Y,S,D,C,H=s.symbols('X Y S D C H');f=(X*X-1)*(Y*Y+S*Y)
    variables=(X,Y,S,D,C,H)
    I=s.groebner([s.diff(f,X),s.diff(f,Y),(X*X-1)*Y*D-C,C*H-1],*variables,order='grevlex',domain=s.QQ)
    J=s.groebner([X,2*Y+S,S*D-2*C,C*H-1],*variables,order='grevlex',domain=s.QQ)
    assert I==J
    return dict(multiplier_substitution_exact=True,anchored_product_identity=True,
                all_original_derivatives_verified=True,mu12_action_on_original_character=True,
                uniform_nonzero_target_ideal_equality=True,actual_quotient_discriminant='c^2+36',
                generic_regular_orbits=2,regular_orbit_dimension=20)


def root_audit(row):
    A=cartan('E7');positive=root_data(A)[1];roots=list(positive)+[tuple(-t for t in rr) for rr in positive]
    H=[4,6,8,12,9,6,3];p=row['prime'];params=row['torus_coroot_parameters'];c=row['target_c']
    B=[[sum(rr[i]*rr[j] for rr in roots) for j in range(7)] for i in range(7)]
    assert sum(H[i]*B[i][j]*H[j] for i in range(7) for j in range(7))==432
    vals=[prod_mod([pow(t,e,p) for t,e in zip(params,rr)],p) for rr in roots]
    assert 1 not in vals
    target=[c*pow(432,-1,p)*h%p for h in H]
    grad=[sum(v*rr[i] for v,rr in zip(vals,roots))%p for i in range(7)]
    assert grad==[sum(B[i][j]*target[j] for j in range(7))%p for i in range(7)]
    Hess=[[sum(v*rr[i]*rr[j] for v,rr in zip(vals,roots))%p for j in range(7)] for i in range(7)]
    rank=mod_rank(Hess,p);assert rank==7
    outside=sum(sum(rr[i]*target[i] for i in range(7))%p!=0 for rr in roots)
    assert outside==106
    assert (3*(pow(row['central_parameter'],12,p)-pow(row['central_parameter'],-12,p))-c)%p==0
    return dict(prime=p,no_source_root_value_equals_one=True,target_centralizer_checked=True,
                Cartan_rank=rank,root_rank=outside,full_rank=113,kernel_dimension=20,trace_norm=432)


def highest_audit(rows):
    result=[]
    for rec in rows:
        A=cartan(rec['type']);Ai=inverse(A);roots=root_data(A)[1]
        coefficients=[tuple(sum(Ai[i][j]*rr[j] for j in range(len(A))) for i in range(len(A))) for rr in roots]
        theta=max(coefficients,key=sum);labels=s.Matrix(A)*s.Matrix(theta)
        assert list(map(int,theta))==rec['simple_coefficients']
        assert list(labels).index(1)+1==rec['highest_weight_node']
        result.append(dict(type=rec['type'],highest_root_checked_by_reflections=True))
    return result


def main():
    source=json.loads((ROOT/'adjoint-pencil-scope-verification.json').read_text())
    report=dict(result='passed',uses_Sage=False,
                independently_reconstructed_characters=[reconstruct(row) for row in source['charge_records']],
                diagram_checks=[diagram_audit(k) for k in ('E6','E7','E8')],
                highest_root_checks=highest_audit(source['highest_root_obstructions']),
                symbolic=symbolic(),original_root_point=root_audit(source['independent_original_root_point']),
                arbitrary_affine_changes_classified=False,all_exceptional_Levis_resolved=False,
                literature_priority_claim=False)
    names=['ADJOINT_PENCIL_SCOPE.md','tools/audit_adjoint_pencil_scope.py',
           'tools/audit_affine_branching_certificates.py','tools/audit_branching_norm_fibers.py',
           'adjoint-pencil-scope-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'adjoint-pencil-scope-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
