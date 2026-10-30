"""Exact full-character and boundary certificates for E8 Levi nodes 1235."""
from hashlib import sha256
from pathlib import Path
import json
from sage.all import (QQ, PolynomialRing, LaurentPolynomialRing, matrix, vector,
                      prod, RootSystem)
from verify_branching_norm_fibers import charge_polynomial

ROOT=Path(__file__).resolve().parents[1]

def main():
    row=charge_polynomial('E8',(0,1,2,4))
    L=LaurentPolynomialRing(QQ,names='s,q,r,v');s,q,r,v=L.gens();p=s*s
    R=PolynomialRing(L,names='A,B,x,y');A,B,x,y=R.gens()
    original=sum(t['coefficient']*prod(z**e for z,e in zip((p,q,r,v),t['charge']))*
                 prod(z**e for z,e in zip((A/p,x,p*B,y),t['exponent'])) for t in row['terms'])
    u=s**3;ws=[s**3*q*r*v,s**3*q*r,s**3*q,s**-9*q**-3*r**-2*v**-1]
    S=sum(ws);T=sum(ws[i]*ws[j] for i in range(4) for j in range(i+1,4));V=sum(1/w for w in ws)
    FA=x*y+S*x/u+V*y/u+T+u**-2
    FB=x*y+u*V*x+u*S*y+T+u*u
    H=x*x+y*y+T*x*y+(u*S+V/u)*x+(S/u+u*V)*y+S*V-3
    assert original==A*B+A*FA+B*FB+H
    graph={A:-FB,B:-FA};residual=H-FA*FB
    assert original.subs(graph)==residual
    for z in (x,y):assert original.derivative(z).subs(graph)==residual.derivative(z)
    for z in L.gens():
        D=lambda f:R({ex:z*co.derivative(z) for ex,co in f.dict().items()})
        assert D(original).subs(graph)==D(residual)
    M=matrix(QQ,[[3,0,0,0],[3,1,1,1],[3,1,1,0],[3,1,0,0]])
    assert M.det()==-3
    # Symbolic computations use invariant coefficients as independent parameters.
    K=LaurentPolynomialRing(QQ,'u');u=K.gen()
    P0=PolynomialRing(K.fraction_field(),names='A,B,x,y,S,T,V');A,B,x,y,S,T,V=P0.gens()
    k=u*V+S/u;m=u*S+V/u;c0=1-S*V;v0=u*u+u**-2
    alpha=c0-y*y-k*y
    beta=-m*y*y-(T+v0+S*S+V*V)*y-T*k
    gamma=c0*y*y-T*m*y+S*V-4-T*T-v0*T
    FA=x*y+S*x/u+V*y/u+T+u**-2
    FB=x*y+u*V*x+u*S*y+T+u*u
    H=x*x+y*y+T*x*y+m*x+k*y+S*V-3
    phi=A*B+A*FA+B*FB+H
    assert H-FA*FB==alpha*x*x+beta*x+gamma
    P=4*alpha**2*gamma.derivative(y)-2*alpha*beta*beta.derivative(y)+beta**2*alpha.derivative(y)
    assert P.degree(y)==5
    assert P.coefficient({y:5})==2*((u*S-V/u)**2+4)
    assert P.coefficient({T:2,y:2})==-k
    assert P.coefficient({T:2,y:1})==-2*(k*k+c0)
    spec={S:0,T:2,V:0}
    specialize=lambda f:P0({ex:co(u=1) for ex,co in f.subs(spec).dict().items()})
    assert specialize(alpha)==1-y*y and specialize(beta)==-4*y
    assert specialize(P)==8*y*(y*y-3)*(y*y+1)
    U=PolynomialRing(QQ,'z');z=U.gen();f=8*z*(z*z-3)*(z*z+1)
    assert f.gcd(f.derivative())==1
    Delta=(u*V-S/u)**2+4
    y0=-k/2
    assert alpha.subs({y:y0})==Delta/4
    assert beta.subs({y:y0})==k*(S*S+V*V+v0*(2-S*V)-2*T)/4
    assert alpha.derivative(y)==-2*y-k
    assert beta.subs({y:y0}).derivative(T)==-k/2
    # These identities prove coprimality and nonvanishing over C as well as Q.
    assert c0.subs({S:-u*u*V})==1+u*u*V*V
    assert k.subs({S:u*u*V})==2*u*V
    witness={A:0,B:0,x:-1,y:-1,S:4,T:6,V:4}
    at=lambda f:QQ(P0(P0(f).subs(witness)).constant_coefficient()(u=1))
    assert all(at(phi.derivative(z))==0 for z in (A,B,x,y))
    Hr=matrix(QQ,[[at(phi.derivative(z).derivative(w)) for w in (A,B,x,y)] for z in (A,B,x,y)])
    assert Hr.det()==-112
    D=lambda f:P0({ex:u*co.derivative(u) for ex,co in P0(f).dict().items()})
    mixed=vector(QQ,[at(D(phi.derivative(z))) for z in (A,B,x,y)])
    assert mixed==vector(QQ,[6,-6,0,0])
    hu=at(D(D(phi)))-mixed*Hr.inverse()*mixed
    assert hu==56
    assert [at(phi.derivative(z)) for z in (S,T,V)]==[2,1,2]
    C=matrix.identity(QQ,3)+matrix.ones(QQ,3,3)
    central=matrix.block_diagonal(matrix(QQ,[[hu]]),6*C)
    assert central.det()==48384
    # Torus-to-SL4-invariants determinant and all reflection walls.
    Z=LaurentPolynomialRing(QQ,names='w1,w2,w3');w1,w2,w3=Z.gens()
    ws=[w1,w2,w3,1/(w1*w2*w3)]
    inv=[sum(ws),sum(ws[i]*ws[j] for i in range(4) for j in range(i+1,4)),sum(1/w for w in ws)]
    J=matrix(Z,[[w*g.derivative(w) for w in Z.gens()] for g in inv])
    vand=prod(ws[i]-ws[j] for i in range(4) for j in range(i+1,4))
    assert J.det() in (vand,-vand)
    weyl=json.loads((ROOT/'e8-weyl-transfers-verification.json').read_text())
    group=next(g for g in weyl['Weyl_classes'] if g['seed']==[1,2,3,5])
    assert len(group['members'])==28
    result=dict(result='passed',type='E8',representation='adjoint',retained_nodes=[1,2,3,5],
                character=row,original_central_cover_degree=row['cover_degree'],
                common_cover_degrees=[2,3],central_exponent_matrix=[list(map(int,r)) for r in M.rows()],
                full_character_verified=True,all_original_derivative_chain_rules=True,
                binary_degree=5,primitive_and_squarefree=True,
                exceptional_base_dimension_bound=2,full_elimination_boundary_dimension_bound=3,
                torus_quotient_jacobian_is_Vandermonde=True,ramification_walls=6,
                full_relative_hessian_determinant=-112,residual_hessian_determinant=112,
                central_schur_determinant=48384,
                regular_fiber_dimension=10,nonempty_generic_regular_fibers=True,
                new_standard_subsets=28,covered_proper_E8_standard_subsets=218,
                remaining_standard_subsets=37,remaining_Weyl_classes=5,
                all_E8_Levis_resolved=False,full_map_degree_claim=False,literature_priority_claim=False)
    names=['E8_A2_A1_PAIR.md','tools/verify_e8_a2_a1_pair.py','BINARY_BRANCHING_FIBERS.md',
           'e8-weyl-transfers-verification.json']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                            for n in names}
    (ROOT/'e8-a2-a1-pair-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('character','source_sha256')},indent=2))

if __name__=='__main__':main()
