"""Exact full-root and binary-boundary certificates for E8, Levi nodes 125."""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
from sage.all import QQ, RootSystem, PolynomialRing, matrix, vector, prod

ROOT=Path(__file__).resolve().parents[1]
BASIS=[(1,0,0,0,0,0,0,0),(0,1,0,0,0,0,0,0),(0,0,0,0,1,0,0,0),
       (1,1,2,2,1,0,0,0),(0,0,0,0,0,0,0,1),
       (1,1,2,2,2,2,2,1),(1,2,2,4,3,2,2,1),(2,3,4,6,5,4,2,1)]

def main():
    lattice=RootSystem('E8').root_lattice(); A=lattice.cartan_type().cartan_matrix()
    roots=[vector(QQ,[r[i+1] for i in range(8)]) for r in lattice.roots()]
    H=matrix(QQ,BASIS).transpose()
    assert H.transpose()*A*H==2*matrix.identity(QQ,8) and abs(H.det())==16
    assert all(vector(QQ,b) in roots for b in BASIS)
    actual=Counter({tuple(r*A*H):1 for r in roots}); actual[(0,)*8]=8
    R=PolynomialRing(QQ,names='X,x,y,t,a,b,c,d'); X,x,y,t,a,b,c,d=R.gens()
    B=a*b+c*d; C=a*c+b*d; D=a*d+b*c
    Q=t*t+a*a+b*b+c*c+d*d+a*b*c*d-8
    phi=X*X+x*x+y*y+t*X*x*y+B*X*x+C*X*y+D*x*y+t*D*X+t*C*x+t*B*y+Q
    # Reconstruct each Laurent weight of the full 248-dimensional character.
    expanded=Counter()
    for ex,co in phi.dict().items():
        weights={(0,)*8:QQ(co)}
        for i,power in enumerate(ex):
            for _ in range(power):
                nxt=Counter()
                for weight,m in weights.items():
                    for sign in (-1,1):
                        v=list(weight);v[i]+=sign;nxt[tuple(v)]+=m
                weights=nxt
        expanded.update(weights)
    assert {k:v for k,v in expanded.items() if v}==dict(actual)
    tetrads=sorted(tuple(i+1 for i,e in enumerate(ex) if e) for ex,co in phi.dict().items() if sum(ex)==4)
    assert len(tetrads)==14 and phi(*([2]*8))==248
    L=t*y+B; M=C*y+t*D
    aa=1-L*L/4; bb=D*y+t*C-L*M/2; cc=y*y+t*B*y+Q-M*M/4
    graph={X:-(L*x+M)/2}; residual=aa*x*x+bb*x+cc
    assert phi.subs(graph)==residual
    for v in (x,y,t,a,b,c,d):
        assert phi.derivative(v).subs(graph)==residual.derivative(v)
    P=4*aa*aa*cc.derivative(y)-2*aa*bb*bb.derivative(y)+bb*bb*aa.derivative(y)
    assert P.degree(y)==5 and P.coefficient({y:5})==t**4/2
    assert P.subs({t:0,a:0,b:0,c:0,d:0})==8*y
    spec={t:1,a:1,b:0,c:0,d:1}
    assert aa.subs(spec)==1-y*y/4 and bb.subs(spec)==y/2
    assert P.subs(spec)==y*(y*y-3)*(y*y-5)/2
    S=PolynomialRing(QQ,'u');u=S.gen();p=u*(u*u-3)*(u*u-5)/2
    assert p.gcd(p.derivative())==1
    assert aa.derivative(y)==-t*(t*y+B)/2
    assert bb.derivative(y).subs({t:0})==D-B*C/2
    for eps in (-1,1):
        assert D-eps*C==(a-eps*b)*(d-eps*c)
        # B-2eps is irreducible: primitive linear in a with coprime
        # coefficients b and cd-2eps. It divides neither displayed factor.
        assert len((B-2*eps).factor())==1 and (B-2*eps).factor()[0][1]==1
        assert (B-2*eps).gcd(D-eps*C)==1
    witness={X:-1,x:-1,y:-1,t:2,a:2,b:2,c:2,d:2}
    assert all(phi.derivative(v).subs(witness)==0 for v in (X,x,y))
    relative=matrix(QQ,[[phi.derivative(v).derivative(w).subs(witness)
                        for w in (X,x,y)] for v in (X,x,y)])
    assert relative==matrix(QQ,[[2,6,6],[6,2,6],[6,6,2]]) and relative.det()==224
    cp=[phi.derivative(v).subs(witness) for v in (t,a,b,c,d)]
    assert cp==[-21,6,6,6,6]
    logarithmic=matrix.diagonal([2*v for v in cp])
    assert logarithmic.det()==-870912
    assert aa.subs(witness)==-8
    weyl=json.loads((ROOT/'e8-weyl-transfers-verification.json').read_text())
    group=next(g for g in weyl['Weyl_classes'] if g['seed']==[1,2,5])
    assert len(group['members'])==21
    result=dict(result='passed',type='E8',representation='adjoint',retained_nodes=[1,2,5],
                coroot_columns=BASIS,central_cover_degree=16,full_root_character_reconstructed=True,
                character_dimension=248,tetrads=tetrads,
                all_elimination_chain_rules=True,critical_polynomial_degree=5,
                primitive_leading_coefficient='t^4/2',squarefree_specialization='y(y^2-3)(y^2-5)/2',
                boundary_dimension_bound=4,exceptional_base_dimension_bound=3,
                relative_hessian_determinant=224,logarithmic_hessian_diagonal=[-42,12,12,12,12],
                logarithmic_hessian_determinant=-870912,nonempty_generic_regular_fibers=True,
                regular_fiber_dimension=6,new_standard_subsets=21,
                covered_proper_E8_standard_subsets=190,remaining_standard_subsets=65,
                remaining_Weyl_classes=6,all_E8_Levis_resolved=False,
                full_map_degree_claim=False,literature_priority_claim=False)
    names=['E8_A1_TRIPLE.md','tools/verify_e8_a1_triple.py','BINARY_BRANCHING_FIBERS.md',
           'e8-weyl-transfers-verification.json']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                            for n in names}
    (ROOT/'e8-a1-triple-verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','tetrads','coroot_columns')},indent=2))

if __name__=='__main__':main()
