"""Independent native-root and SymPy audit of every E8 A2+A1^2 boundary."""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip,json
import sympy as s
from audit_affine_branching_certificates import cartan,root_data
from audit_branching_norm_fibers import reconstruct

ROOT=Path(__file__).resolve().parents[1]

def main():
    rows=json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    row=next(r['character'] for r in rows if r['character']['type']=='E8' and r['character']['nodes']==[1,2,3,5])
    reconstruct(row) # Native roots and subset fundamental weights, no Sage.
    A,B,x,y,u,S,T,V=s.symbols('A B x y u S T V')
    FA=x*y+S*x/u+V*y/u+T+u**-2
    FB=x*y+u*V*x+u*S*y+T+u*u
    H=x*x+y*y+T*x*y+(u*S+V/u)*x+(S/u+u*V)*y+S*V-3
    phi=A*B+A*FA+B*FB+H
    z,q,r,v=s.symbols('z q r v');p=z*z
    ws=[z**3*q*r*v,z**3*q*r,z**3*q,z**-9*q**-3*r**-2*v**-1]
    substitution={u:z**3,S:sum(ws),T:sum(ws[i]*ws[j] for i,j in combinations(range(4),2)),V:sum(1/w for w in ws)}
    actual=s.expand(sum(s.Rational(t['coefficient'])*s.prod(w**e for w,e in zip((A/p,x,p*B,y),t['exponent']))*
                        s.prod(w**e for w,e in zip((p,q,r,v),t['charge'])) for t in row['terms']))
    assert s.expand(actual-phi.subs(substitution,simultaneous=True))==0
    residual=s.expand(H-FA*FB); graph={A:-FB,B:-FA}
    assert s.expand(phi.subs(graph,simultaneous=True)-residual)==0
    for w in (x,y,u,S,T,V):
        assert s.expand(s.diff(phi,w).subs(graph,simultaneous=True)-s.diff(residual,w))==0
    k=u*V+S/u;m=u*S+V/u;c0=1-S*V;v0=u*u+u**-2
    alpha=c0-y*y-k*y
    beta=-m*y*y-(T+v0+S*S+V*V)*y-T*k
    gamma=c0*y*y-T*m*y+S*V-4-T*T-v0*T
    assert s.expand(residual-alpha*x*x-beta*x-gamma)==0
    assert s.gcd(s.expand(u*alpha),s.expand(u*u*beta))==1
    P=s.expand(4*alpha**2*s.diff(gamma,y)-2*alpha*beta*s.diff(beta,y)+beta**2*s.diff(alpha,y))
    assert s.Poly(P,y).degree()==5
    assert s.expand(s.Poly(P,y).LC()-2*((u*S-V/u)**2+4))==0
    coeffT2=s.expand(P).coeff(T,2)
    assert s.expand(coeffT2.coeff(y,2)+k)==0
    assert s.expand(coeffT2.coeff(y,1)+2*(k*k+c0))==0
    assert s.gcd(S+u*u*V,c0)==1
    spec={u:1,S:0,T:2,V:0}
    ps=s.Poly(P.subs(spec),y)
    assert ps.as_expr()==s.expand(8*y*(y*y-3)*(y*y+1)) and ps.gcd(ps.diff()).degree()==0
    # Exact elimination of the derivative equation; both complex components
    # and their k=0 subloci remain in these identities.
    y0=-k/2;delta=u*V-S/u
    assert s.expand(alpha.subs(y,y0)-(delta*delta+4)/4)==0
    assert s.expand(beta.subs(y,y0)-k*(S*S+V*V+v0*(2-S*V)-2*T)/4)==0
    for eps in (-1,1):
        sigma=u*u*V-2*eps*s.I*u
        b_boundary=s.expand(beta.subs(y,y0).subs(S,sigma))
        assert s.expand(s.diff(b_boundary,T)+u*V-eps*s.I)==0
        # k=0 is expressly included: it is a divisor in this component,
        # not an excluded source open.
        assert s.expand(b_boundary.subs(V,eps*s.I/u))==0
    # Identify the six branch walls as ACTUAL E8 roots centralizing the Levi.
    C=s.Matrix(cartan('E8'));positive=root_data(cartan('E8'))[1]
    commuting=[r for r in positive if all(r[j]==0 for j in (0,1,2,4))]
    assert len(commuting)==6
    HC=s.Matrix(row['central_cocharacters']).T*s.diag(2,1,1,1)
    M=s.Matrix([[3,0,0,0],[3,1,1,1],[3,1,1,0],[3,1,0,0]])
    assert M.det()==-3
    weights=[s.Matrix(t) for t in [(3,1,1,1),(3,1,1,0),(3,1,0,0),(-9,-3,-2,-1)]]
    charges={tuple(s.Matrix(r).T*HC) for r in commuting}
    all_charges=charges|{tuple(-v for v in r) for r in charges}
    assert all_charges=={tuple(weights[i]-weights[j]) for i in range(4) for j in range(4) if i!=j}
    directions=[]
    for root in commuting:
        d=M*(HC.T*C*HC).inv()*HC.T*s.Matrix(root)
        assert d[0]==0
        vector4=list(d[1:])+[-sum(d[1:])]
        assert sorted(vector4)==[-1,0,0,1]
        i,j=vector4.index(1),vector4.index(-1)
        assert tuple(s.Matrix(root).T*HC)==tuple(weights[i]-weights[j])
        directions.append(list(map(int,d)))
    w=s.symbols('w0:4')
    invariants=[sum(w),sum(w[i]*w[j] for i,j in combinations(range(4),2)),sum(1/a for a in w)]
    for i,j in combinations(range(4),2):
        divided=[1,sum(w[k] for k in range(4) if k not in (i,j)),1/(w[i]*w[j])]
        for f,g in zip(invariants,divided):
            assert s.cancel(w[i]*s.diff(f,w[i])-w[j]*s.diff(f,w[j])-(w[i]-w[j])*g)==0
    w1,w2,w3=s.symbols('w1 w2 w3');w4=1/(w1*w2*w3)
    ww=[w1,w2,w3,w4]
    inv=[sum(ww),sum(ww[i]*ww[j] for i,j in combinations(range(4),2)),sum(1/a for a in ww)]
    J=s.Matrix([[a*s.diff(f,a) for a in (w1,w2,w3)] for f in inv])
    vand=s.prod(ww[i]-ww[j] for i,j in combinations(range(4),2))
    assert s.cancel(J.det()/vand) in (-1,1)
    pt={x:-1,y:-1,u:1,S:4,T:6,V:4}
    assert all(s.diff(residual,a).subs(pt)==0 for a in (x,y))
    Hr=s.hessian(residual,(x,y)).subs(pt)
    assert Hr.det()==112
    D=lambda f:u*s.diff(f,u)
    assert all(s.diff(D(residual),a).subs(pt)==0 for a in (x,y))
    assert D(D(residual)).subs(pt)==56
    assert [s.diff(residual,a).subs(pt) for a in (S,T,V)]==[2,1,2]
    Hinv=[]
    for f in inv:
        Hinv.append(s.Matrix([[a*s.diff(b*s.diff(f,b),a).subs({w1:1,w2:1,w3:1}) for b in (w1,w2,w3)] for a in (w1,w2,w3)]).subs({w1:1,w2:1,w3:1}))
    block=2*Hinv[0]+Hinv[1]+2*Hinv[2]
    assert block==6*(s.eye(3)+s.ones(3))
    assert 56*block.det()==48384
    result=dict(result='passed',independent_full_root_character=True,all_elimination_chain_rules=True,
                gcd_alpha_beta_is_one=True,primitive_quintic=True,squarefree_specialization=True,
                two_complex_boundary_components_checked=True,k_zero_intersections_retained=True,
                exceptional_base_dimension_bound=2,full_boundary_dimension_bound=3,
                actual_commuting_E8_roots=6,central_root_directions=directions,
                all_six_ramification_wall_factorizations=True,Vandermonde_jacobian=True,
                residual_hessian_determinant=112,central_schur_determinant=48384,
                nonempty_smooth_generic_regular_fiber_dimension=10,
                scope='Exact independent root, coefficient, boundary and ramification checks; geometric implications are in the written proof. No full-map degree or literature novelty assertion.')
    names=['E8_A2_A1_PAIR.md','tools/audit_e8_a2_a1_pair.py','e8-a2-a1-pair-verification.json']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                            for n in names}
    (ROOT/'e8-a2-a1-pair-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
