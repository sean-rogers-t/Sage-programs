"""Sage-free Euclidean-root and exact SymPy audit of the E8 A1^3 proof."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def main():
    # Construct E8 directly, independently of Sage and the character archive.
    roots=[]
    for i,j in combinations(range(8),2):
        for a,b in product((-2,2),repeat=2):
            v=[0]*8;v[i]=a;v[j]=b;roots.append(tuple(v))
    roots += [v for v in product((-1,1),repeat=8) if sum(x<0 for x in v)%2==0]
    simple=[(1,-1,-1,-1,-1,-1,-1,1),(2,2,0,0,0,0,0,0),(-2,2,0,0,0,0,0,0)]
    for i in range(1,6):
        v=[0]*8;v[i]=-2;v[i+1]=2;simple.append(tuple(v))
    E=s.Matrix.hstack(*(s.Matrix(v) for v in simple))
    # Recover the completion by searching roots perpendicular to nodes 125.
    coef=sorted(tuple(E.inv()*s.Matrix(r)) for r in roots if min(E.inv()*s.Matrix(r))>=0)
    A=E.T*E/4
    basis=[tuple(int(i==j) for i in range(8)) for j in (0,1,4)]
    inner=lambda u,v:(s.Matrix(u).T*A*s.Matrix(v))[0]
    for r in coef:
        if all(inner(r,b)==0 for b in basis):basis.append(r)
    assert len(basis)==8
    # Fourth direction in character order is the unique root completing
    # the tetrad through all three retained roots.
    labels=[tuple(inner(tuple(E.inv()*s.Matrix(r)),b) for b in basis) for r in roots]
    supports={tuple(i for i,v in enumerate(row) if v) for row in labels if sum(bool(v) for v in row)==4}
    fourth=next(iter(set(next(v for v in supports if {0,1,2}<=set(v)))-{0,1,2}))
    order=[0,1,2,fourth]+[i for i in range(3,8) if i!=fourth]
    basis=[basis[i] for i in order]
    H=s.Matrix.hstack(*(s.Matrix(v) for v in basis))
    assert H.T*A*H==2*s.eye(8) and abs(H.det())==16
    labels=Counter(tuple((s.Matrix(r).T*E*H/4)) for r in roots)
    tetrads={tuple(i for i,v in enumerate(row) if v) for row in labels if sum(bool(v) for v in row)==4}
    rebuilt=Counter()
    for i in range(8):
        for sign in (-2,2):
            v=[0]*8;v[i]=sign;rebuilt[tuple(v)]+=1
    for indices in tetrads:
        for signs in product((-1,1),repeat=4):
            v=[0]*8
            for i,sign in zip(indices,signs):v[i]=sign
            rebuilt[tuple(v)]+=1
    assert rebuilt==labels and len(labels)==240 and len(tetrads)==14
    X,x,y,t,a,b,c,d=s.symbols('X x y t a b c d');vs=(X,x,y,t,a,b,c,d)
    derived=sum(v*v for v in vs)-8+sum(s.prod(vs[i] for i in indices) for indices in tetrads)
    B=a*b+c*d;C=a*c+b*d;D=a*d+b*c
    Q=t*t+a*a+b*b+c*c+d*d+a*b*c*d-8
    phi=X*X+x*x+y*y+t*X*x*y+B*X*x+C*X*y+D*x*y+t*D*X+t*C*x+t*B*y+Q
    assert s.expand(derived-phi)==0 and phi.subs(dict.fromkeys(vs,2))==248
    L=t*y+B;M=C*y+t*D;graph=-(L*x+M)/2
    residual=s.expand(phi.subs(X,graph))
    aa=1-L*L/4;bb=D*y+t*C-L*M/2;cc=y*y+t*B*y+Q-M*M/4
    assert s.expand(residual-aa*x*x-bb*x-cc)==0
    for v in vs[1:]:
        assert s.expand(s.diff(phi,v).subs(X,graph)-s.diff(residual,v))==0
    assert s.gcd(aa,bb)==1
    P=s.expand(4*aa*aa*s.diff(cc,y)-2*aa*bb*s.diff(bb,y)+bb*bb*s.diff(aa,y))
    assert s.degree(P,y)==5 and s.Poly(P,y).LC()==t**4/2
    assert P.subs({t:0,a:0,b:0,c:0,d:0})==8*y
    spec={t:1,a:1,b:0,c:0,d:1}
    p=s.Poly(P.subs(spec),y)
    assert p.degree()==5 and p.gcd(p.diff()).degree()==0
    assert s.factor(p.as_expr())==y*(y*y-5)*(y*y-3)/2
    # Exact ideal identity, stronger than a set-theoretic boundary inference.
    assert s.expand(t*aa-L*s.diff(aa,y)/2-t)==0
    assert s.expand(aa.subs(t,0)+(B*B-4)/4)==0
    for eps in (-1,1):
        assert s.expand(D-eps*C-(a-eps*b)*(d-eps*c))==0
        assert s.gcd(b,c*d-2*eps)==1  # primitive linear B-2eps in a
        assert s.gcd(B-2*eps,D-eps*C)==1
    pt={X:-1,x:-1,y:-1,t:2,a:2,b:2,c:2,d:2}
    assert all(s.diff(phi,v).subs(pt)==0 for v in (X,x,y))
    Hr=s.hessian(phi,(X,x,y)).subs(pt)
    assert Hr.det()==224
    # At z_i=1, d(z_i+z_i^-1)=0 and d_log^2(z_i+z_i^-1)=2.
    Hc=s.diag(*(2*s.diff(phi,v).subs(pt) for v in (t,a,b,c,d)))
    assert Hc==s.diag(-42,12,12,12,12) and Hc.det()==-870912
    source=json.loads((ROOT/'e8-a1-triple-verification.json').read_text())
    assert [list(map(int,b)) for b in basis]==source['coroot_columns']
    assert source['relative_hessian_determinant']==int(Hr.det())
    assert source['logarithmic_hessian_determinant']==int(Hc.det())
    result=dict(result='passed',full_Euclidean_root_character=True,
                orthogonal_completion_independently_recovered=True,central_cover_degree=16,
                exact_binary_elimination_and_chain_rules=True,coprime_binary_coefficients=True,
                primitive_and_separable_degree_five_polynomial=True,
                entire_boundary_control=True,exact_boundary_ideal_identity=True,
                exceptional_base_dimension_bound=3,full_boundary_dimension_bound=4,
                relative_hessian_determinant=224,central_logarithmic_hessian_determinant=-870912,
                smooth_nonempty_generic_regular_fiber_dimension=6,
                scope='Independent exact root and polynomial inputs. Scheme geometry and genericity are proved in the note; no full Springer degree or literature novelty assertion.')
    names=['E8_A1_TRIPLE.md','tools/audit_e8_a1_triple.py','e8-a1-triple-verification.json']
    result['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                            for n in names}
    (ROOT/'e8-a1-triple-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
