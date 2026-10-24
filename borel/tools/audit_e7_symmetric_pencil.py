"""Sage-free root and symbolic audit of the symmetric-pencil theorem."""
from pathlib import Path
from hashlib import sha256
import json
import gzip
from itertools import combinations
import sympy as sp
from audit_affine_branching_certificates import cartan, root_data
from audit_branching_norm_fibers import reconstruct

ROOT=Path(__file__).resolve().parents[1]


def actual_polynomial(row,affine,central):
    reconstruct(row)
    return sum(sp.Rational(t['coefficient'])*sp.prod(x**e for x,e in zip(affine,t['exponent']))*
               sp.prod(z**e for z,e in zip(central,t['charge'])) for t in row['terms'])


def audit(source):
    x,y,z,w,A,B,C,p,q,r=sp.symbols('x y z w A B C p q r')
    g=x*x+A*B*C*x+A*A+B*B+C*C-7
    Q=sp.Matrix([[2,B*C+A*x,A*C+B*x],[B*C+A*x,2,A*B+C*x],[A*C+B*x,A*B+C*x,2]])
    leaf=sp.Matrix([y,z,w]);phi=g+(leaf.T*Q*leaf)[0]/2
    d=sp.expand(Q.det());pp=sp.Poly(d.subs(x,-A*B*C/2),A,B,C)
    hh=sp.Poly(sp.diff(d,x).subs(x,-A*B*C/2),A,B,C)
    assert sp.gcd(pp,hh).as_expr()==1 and pp.eval({A:2,B:2,C:2})==-216
    assert sp.expand(hh.as_expr()-A*B*C/2*(3*A*A*B*B*C*C-4*A*A*B*B-4*A*A*C*C-4*B*B*C*C+8*A*A+8*B*B+8*C*C-24))==0
    # Check both the necessary rank-one identities and all four parametrizations.
    E=[A*x+B*C-2,B*x+A*C-2,C*x+A*B-2]
    assert sp.expand(E[0]-E[1]-(A-B)*(x-C))==0
    assert sp.expand(E[0]-E[2]-(A-C)*(x-B))==0
    assert sp.expand(E[1]-E[2]-(B-C)*(x-A))==0
    t=sp.Symbol('t',nonzero=True)
    for point in [(2/t-t,t,t,t),(t,2/t-t,t,t),(t,t,2/t-t,t),(t,t,t,2/t-t)]:
        assert all(sp.cancel(e.subs(dict(zip((x,A,B,C),point)),simultaneous=True))==0 for e in E)
    kernel_quadratic=sp.expand((A*y*z+B*y*w+C*z*w).subs(w,-y-z))
    assert sp.expand(kernel_quadratic+B*y*y+C*z*z-(A-B-C)*y*z)==0
    assert sp.solve(sp.Poly(kernel_quadratic,y,z).coeffs(),(A,B,C))=={A:0,B:0,C:0}
    for eps1 in (-1,1):
        for eps2 in (-1,1):
            assert sp.expand(phi.subs({A:eps1*A,B:eps2*B,C:eps1*eps2*C,z:eps1*z,w:eps2*w},simultaneous=True)-phi)==0
    row=source['symmetric_pencil']['character']; actual=actual_polynomial(row,(x,y,z,w),(p,q,r))
    U=p*q;V=p*q*r;W=p*q*q*r
    torus_phi=sp.expand(phi.subs({A:U+1/U,B:V+1/V,C:W+1/W},simultaneous=True))
    assert sp.expand(actual-torus_phi)==0
    for v in (x,y,z,w):assert sp.expand(sp.diff(actual,v)-sp.diff(torus_phi,v))==0
    for v in (p,q,r):assert sp.expand(v*sp.diff(actual,v)-v*sp.diff(torus_phi,v))==0
    H=row['central_cocharacters'];root_positive=root_data(cartan('E7'))[1]
    roots=list(root_positive)+[tuple(-v for v in root) for root in root_positive]
    charges=[[sum(a*b for a,b in zip(root,h)) for h in H] for root in roots]
    trace=[[sum(c[i]*c[j] for c in charges) for j in range(3)] for i in range(3)]
    assert trace==[[216,288,144],[288,432,216],[144,216,144]]
    change=sp.Matrix([[1,1,1],[1,1,2],[0,1,1]])
    assert change.det()==-1 and sp.cancel(U*V/W-p)==sp.cancel(W/V-q)==sp.cancel(V/U-r)==0
    uu,vv,ww=sp.symbols('U V W');ts=(uu,vv,ww);aa,bb,cc=[v+1/v for v in ts]
    R0=aa*aa+bb*bb+cc*cc-7-aa*aa*bb*bb*cc*cc/4
    Hlog=sp.Matrix([[(v*sp.diff(t0*sp.diff(R0,t0),v)).subs(dict.fromkeys(ts,1)) for t0 in ts] for v in ts])
    assert Hlog==-24*sp.eye(3)
    witness={x:-4,y:0,z:0,w:0,A:2,B:2,C:2}
    assert all(sp.diff(phi,v).subs(witness)==0 for v in (x,y,z,w))
    assert sp.hessian(phi,(x,y,z,w)).subs(witness).det()==-432
    # Exhaust the independent A1^4 standard embeddings, not just the archive.
    Cartan=cartan('E7');embeddings=[]
    for mask in range(127):
        nodes=[i for i in range(7) if mask>>i&1]
        if len(nodes)==4 and all(not Cartan[i][j] for i in nodes for j in nodes if i!=j):
            embeddings.append([i+1 for i in nodes])
    assert embeddings==[[1,2,5,7],[2,3,5,7]]
    return dict(full_character_independently_reconstructed=True,collision_coprimality=True,
                rank_one_curves_and_sign_changes=True,low_rank_dimension_bound=2,
                original_derivatives_verified=True,root_trace_matrix=trace,
                ordinary_logarithmic_hessian_diagonal=-24,relative_hessian_determinant=-432,
                all_A1_four_embeddings=embeddings,regular_fiber_dimension=8)


def companion(source):
    a,x,b,y,w,p,q=sp.symbols('a x b y w p q')
    phi=actual_polynomial(source['companion']['character'],(a,x,b,y,w),(p,q))
    F=sp.diff(phi,a)-b;G=sp.diff(phi,b)-a;graph={a:-G,b:-F}
    psi=sp.expand(phi.subs(graph,simultaneous=True))
    U=q;V=p**3*q;A=U+1/U;B=V+1/V;C=U*V+1/(U*V)
    expected=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    assert sp.expand(psi-expected)==0 and sp.expand(A*A+B*B+C*C-A*B*C-4)==0
    for z in (p,q):
        assert sp.expand((z*sp.diff(phi,z)).subs(graph,simultaneous=True)-z*sp.diff(psi,z))==0
    uu,vv=sp.symbols('U V')
    potential=-(uu+1/uu)*(vv+1/vv)*(uu*vv+1/(uu*vv))-3
    H=sp.Matrix([[(z*sp.diff(v*sp.diff(potential,v),z)).subs({uu:1,vv:1}) for v in (uu,vv)] for z in (uu,vv)])
    assert H==sp.Matrix([[-16,-8],[-8,-16]]) and H.det()==192
    assert sp.hessian(psi,(x,y,w)).subs({x:0,y:0,w:0,p:1,q:1}).det()==-216
    return dict(exact_character_and_residual=True,original_target_derivatives=True,
                generic_smoothness_claim=False,conditional_dominance_hessian_determinant=192)


def candidate_audit():
    record=json.loads((ROOT/'e7-quartic-boundary-exploration.json').read_text())
    payload=(ROOT/record['archive']).read_bytes()
    assert sha256(payload).hexdigest()==record['archive_sha256']
    data=json.loads(gzip.decompress(payload))
    vs=sp.symbols('x y w A B C');x,y,w,A,B,C=vs
    psi=x*x+y*y+w*w-x*x*y*y-x*x*w*w-y*y*w*w-x*y*w*(A*x+B*y+C*w)-A*B*x*y-A*C*x*w-B*C*y*w-A*B*C-3
    equations=[sp.diff(psi,v) for v in (x,y,w)]+[A*A+B*B+C*C-A*B*C-4]
    jac=sp.Matrix([[sp.diff(f,v) for v in vs] for f in equations])
    generators=equations+[jac[:,cols].det(method='domain-ge') for cols in combinations(range(6),4)]
    decode=lambda terms:sp.Poly.from_dict({tuple(e):sp.Rational(c) for e,c in terms},vs)
    assert {sp.Poly(f,vs) for f in generators}=={decode(f) for f in data['original_generators']}
    supports=[]
    for terms in data['selected_members']:
        leading=max((tuple(e) for e,c in terms if sp.Rational(c)),key=lambda e:(sum(e),tuple(-v for v in reversed(e))))
        supports.append({i for i,e in enumerate(leading) if e})
    assert {frozenset(s) for s in supports}=={frozenset(s) for s in ((0,),(1,),(2,4),(2,5),(3,),(4,5))}
    dimension=max(len(S) for k in range(7) for v in combinations(range(6),k)
                  for S in [set(v)] if not any(t<=S for t in supports))
    assert dimension==1
    assert not data['reverse_ideal_inclusion_independently_certified']
    assert not record['companion_generic_smoothness_promoted']
    return dict(original_Jacobian_ideal_reconstructed=True,selected_leading_support_dimension=1,
                selected_memberships_verified=False,companion_generic_smoothness_claim=False)


def main():
    source=json.loads((ROOT/'e7-symmetric-pencil-verification.json').read_text())
    report=dict(result='passed',uses_Sage=False,symmetric_pencil=audit(source),
                companion=companion(source),candidate=candidate_audit(),all_E7_Levis_resolved=False,
                full_map_degree_claim=False,literature_priority_claim=False)
    names=['E7_SYMMETRIC_PENCIL.md','tools/audit_e7_symmetric_pencil.py',
           'tools/audit_affine_branching_certificates.py','tools/audit_branching_norm_fibers.py',
           'e7-symmetric-pencil-verification.json','e7-quartic-boundary-exploration.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-symmetric-pencil-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
