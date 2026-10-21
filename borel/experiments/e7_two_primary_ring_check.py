"""Plain Python/SymPy check of the structural identities in the E7 proof."""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def main():
    C=s.Matrix([[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],
      [0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]])
    w=s.eye(7)[:,6];V=[C[:,i] for i in range(6)]
    XB=s.Matrix.hstack(*V,w);QB=C.inv()*s.Matrix.hstack(*V,2*w)
    assert XB.det()==QB.det()==3 and all(s.denom(c)==1 for c in QB)
    assert all(w[i]==0 for i in range(6))
    eye=s.eye(7);chain=(0,2,3,4,5)
    vs=[eye[:,6]+sum((eye[:,j] for j in chain[i:]),s.zeros(7,1)) for i in range(6)]
    u=eye[:,1]+sum(vs[:3],s.zeros(7,1));E=sum(vs,s.zeros(7,1))
    root_model=s.Matrix.hstack(*vs,u)
    zz=C*(u-E)/2;weight_model=s.Matrix.hstack(*(C*v for v in vs),zz)
    assert abs(root_model.det())==abs(weight_model.det())==1
    assert C*(3*u-E)/2==w
    for node in chain:
        reflection=lambda v:v-(C.row(node)*v)[0]*eye[:,node]
        assert sorted(vs.index(reflection(v)) for v in vs)==list(range(6))
        assert reflection(u)==u and reflection((u-E)/2)==(u-E)/2
    alpha=s.symbols('a1:8')
    X=[[5,0,4,3,2,1,0],[-1,0,4,3,2,1,0],[-1,0,-2,3,2,1,0],
       [-1,0,-2,-3,2,1,0],[-1,0,-2,-3,-4,1,0],[-1,0,-2,-3,-4,-5,0]]
    z=[-3,-6,-6,-9,-6,-3,0]
    assert [list(6*v-E) for v in vs]==X and list(3*E-6*u)==z
    xs=[sum(a*c for a,c in zip(alpha,row)) for row in X];Z=sum(a*c for a,c in zip(alpha,z))
    linear=[x+sign*Z for x in xs for sign in (-1,1)]+[-x-y for x,y in combinations(xs,2)]
    B=s.Poly(-sum(f*f for f in linear)/144,*alpha)
    assert all(s.denom(c)==1 for c in B.coeffs()) and B.degree(alpha[6])==0
    E6=C[:6,:6];three_inverse=3*E6.inv()
    assert all(s.denom(c)==1 for c in three_inverse)
    assert all(three_inverse[i,i]%2==0 for i in range(6))
    a=s.Matrix(alpha[:6]);F=(a.T*E6.inv()*a)[0]
    w_root=sum(C.inv()[6,j]*alpha[j] for j in range(7))
    full=(s.Matrix(alpha).T*C.inv()*s.Matrix(alpha))[0]
    assert s.expand(B.as_expr()+s.Rational(3,2)*F)==0
    assert s.expand(F+s.Rational(2,3)*w_root*w_root-full)==0
    h,ss,t=s.symbols('h s t');S=2*ss;T=2*t
    relations=[ss*ss-2*h*t,2*ss*t-9*h**4*ss*ss+6*h**9*ss-h**14,
               t*t+10*h**3*ss**3-9*h**8*ss**2+2*h**13*ss]
    G=s.groebner(relations,h,ss,t,domain=s.QQ)
    equations=[S*S-4*h*T,S*T-2*h**14+6*h**9*S-18*h**5*T,
      T*T+40*h**18-116*h**13*S+324*h**9*T,
      h**14*S-s.Rational(98,285)*h**19-s.Rational(14,5)*h**10*T,
      h**14*T-s.Rational(238,1311)*h**23,h**28]
    assert all(G.reduce(f)[1]==0 for f in equations)
    q=s.symbols('q');bracket=lambda n:sum(q**i for i in range(n))
    poincare=s.cancel((bracket(10)*bracket(14)*bracket(18))/(bracket(5)*bracket(9)))
    towers=bracket(28)+(q**5+q**9)*bracket(14)
    assert s.expand(poincare-towers)==0 and towers.subs(q,1)==56
    multiplicities=q*sum(q**(2*i) for i in range(14))+(q**6+q**10)*sum(q**(2*i) for i in range(7))
    assert multiplicities.subs(q,1)==28
    output=dict(result='passed',engine='Plain Python/SymPy; no Sage or project imports',
                root_complement_index=3,weight_complement_index=3,
                root_permutation_basis_determinant=int(root_model.det()),
                weight_permutation_basis_determinant=int(weight_model.det()),permutation_reflections_checked=5,
                root_quadratic_lift=str(B.as_expr()),ring_identities=6,
                root_quadratic_from_invariant_form=True,
                graded_basis_count=56,total_relative_factors_2=28,
                multiplicity_coefficients=[int(s.expand(multiplicities).coeff(q,i)) for i in range(28)],
                refinement_scope='All 64 retained subsets of {1,2,3,4,5,6}; exponent 2 only, not their full multiplicities.',
                source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/e7-two-primary-ring-check.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
