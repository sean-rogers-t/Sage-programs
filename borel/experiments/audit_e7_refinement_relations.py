"""Plain Python/SymPy audit of the uniform E7 refinement truncation.

No Sage or project imports. The independent sextic input uses 27 E6 forms;
the primary certificate instead uses the orbit of 56 E7 weights.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from math import gcd
from functools import reduce
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
def coeffs(degrees,bound):
 out=[1]+[0]*bound
 for degree in degrees:
  for j in range(degree,bound+1):out[j]+=out[j-degree]
 return out

def levi_degrees(I):
 edges=((1,3),(3,4),(2,4),(4,5),(5,6));left=set(I);ds=[];types=[]
 while left:
  component={min(left)};pending=list(component)
  while pending:
   v=pending.pop()
   for e in edges:
    if v in e:
     u=e[1] if e[0]==v else e[0]
     if u in left and u not in component:component.add(u);pending.append(u)
  left-=component;n=len(component)
  branching=any(sum(v in e and set(e)<=component for e in edges)==3 for v in component)
  if not branching:ds+=list(range(2,n+2));types.append('A'+str(n))
  elif n==6:ds +=[2,5,6,8,9,12];types.append('E6')
  else:ds+=list(range(2,2*n,2))+[n];types.append('D'+str(n))
 return sorted(ds),'+'.join(sorted(types)) or 'trivial'

def main():
 A=s.symbols('a1:7');w=s.symbols('w');variables=A+(w,)
 X=[[5,0,4,3,2,1],[-1,0,4,3,2,1],[-1,0,-2,3,2,1],[-1,0,-2,-3,2,1],[-1,0,-2,-3,-4,1],[-1,0,-2,-3,-4,-5]]
 xs=[sum(a*c for a,c in zip(A,row)) for row in X];z=sum(a*c for a,c in zip(A,[-3,-6,-6,-9,-6,-3]))
 forms=[x+sign*z for x in xs for sign in (-1,1)]+[-x-y for x,y in combinations(xs,2)]
 psi={d:s.Poly(sum(f**d for f in forms),*A).as_expr() for d in (2,5,6,8)}
 B=s.Poly(-psi[2]/144,*A);F=s.Poly((psi[6]+117504*B.as_expr()**3-12*w*psi[5])/3456,*variables)
 assert all(s.denom(c)==1 for c in F.coeffs()) and reduce(gcd,map(int,F.coeffs()))==1
 even=s.Poly(F.as_expr().subs(w,0),*A);odd=s.Poly(s.diff(F.as_expr(),w),*A)
 val=lambda n:(abs(int(n))&-abs(int(n))).bit_length()-1
 assert min(map(val,even.coeffs()))==0 and min(map(val,odd.coeffs()))==2
 primary=json.loads((ROOT/'experiments/e7-refinement-relations.json').read_text());assert primary['result']=='passed'
 def decode(rows):return s.Poly.from_dict({tuple(e):s.Rational(c) for e,c in rows},variables)
 G=decode(primary['sextic']);assert F==G or F==-G
 assert decode(primary['quadratic'])==s.Poly(w*w-B.as_expr(),*variables)
 # Reconstruct an octic independently from the printed E6 psi images.
 O=s.Poly(psi[8]+1568*w*B.as_expr()*psi[5]-9796608*B.as_expr()**4,*variables)
 g=reduce(gcd,map(int,O.coeffs()));O=s.Poly(O.as_expr()/g,*variables)
 L=s.Poly(B.as_expr()*F.as_expr(),*variables)
 for step in range(8):
  if all((O.as_dict().get(e,0)+L.as_dict().get(e,0))%2==0 for e in set(O.as_dict())|set(L.as_dict())):O=s.Poly((O.as_expr()+L.as_expr())/2,*variables)
  else:break
 H=decode(primary['octic']);ld=L.as_dict();od=O.as_dict();hd=H.as_dict();keys=set(ld)|set(od)|set(hd)
 pivot=next(e for e in keys if od.get(e,0))
 second=next(e for e in keys if od.get(pivot,0)*ld.get(e,0)-od.get(e,0)*ld.get(pivot,0))
 mat=s.Matrix([[od.get(e,0),ld.get(e,0)] for e in (pivot,second)])
 aa,bb=mat.inv()*s.Matrix([hd.get(e,0) for e in (pivot,second)])
 assert all(hd.get(e,0)==aa*od.get(e,0)+bb*ld.get(e,0) for e in keys)
 assert s.denom(aa)%2 and s.numer(aa)%2 and s.denom(bb)%2
 assert all(c%2==0 for e,c in hd.items() if e[6])
 # Independently certify that the degree-6 and degree-8 relation multiples
 # have primitive spans through degree 9, after monic quadratic reduction.
 def comps(n,k=6):
  if k==1:yield (n,);return
  for j in range(n+1):
   for rest in comps(n-j,k-1):yield (j,)+rest
 bd=B.as_dict()
 def mul(f,g):
  out={}
  for e,c in f.items():
   for u,v in g.items():
    k=tuple(a+b for a,b in zip(e,u));out[k]=out.get(k,0)+c*v
  return {e:c for e,c in out.items() if c}
 def quad_reduce(f):
  out={}
  for e,c in f.items():
   factor={tuple(e[:6]):c}
   for k in range(e[6]//2):factor=mul(factor,bd)
   for u,v in factor.items():
    key=u+(e[6]%2,);out[key]=out.get(key,0)+v
  return {e:c for e,c in out.items() if c}
 kernel_checks=[]
 for d in range(6,10):
  columns=[]
  for degree,polynomial in ((6,F),(8,H)):
   for parity in (0,1):
    if d-degree-parity<0:continue
    for e in comps(d-degree-parity):
     product={tuple(a+b for a,b in zip(u,e+(parity,))):c for u,c in polynomial.as_dict().items()}
     columns.append(quad_reduce(product))
  pivots={}
  for e in set().union(*(set(f) for f in columns)):
   row=sum((int(f.get(e,0))%2)<<j for j,f in enumerate(columns))
   while row:
    j=row.bit_length()-1
    if j in pivots:row^=pivots[j]
    else:pivots[j]=row;break
  assert len(pivots)==len(columns)
  kernel_checks.append(dict(degree=d,columns=len(columns),rank=len(pivots)))
 C=s.Matrix([[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],
             [0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
 assert s.Poly(-s.Rational(3,2)*(s.Matrix(A).T*C.inv()*s.Matrix(A))[0],*A)==B
 for i in range(6):
  images={A[j]:A[j]-C[j,i]*A[i] for j in range(6)}
  assert s.Poly(F.as_expr().subs(images,simultaneous=True)-F.as_expr(),*variables).is_zero
  assert s.Poly(B.as_expr().subs(images,simultaneous=True)-B.as_expr(),*A).is_zero
 # The E6 longest-element symmetry exchanges the two D5 retained sets.
 C7=C.row_join(s.Matrix([0,0,0,0,0,-1])).col_join(s.Matrix([[0,0,0,0,0,-1,2]]))
 XB=C7[:,:6].row_join(s.eye(7)[:,6]);tau=[5,1,4,3,2,0]
 J=s.zeros(7)
 for i,j in enumerate(tau):J[j,i]=-1
 J[6,6]=1
 symmetry=XB*J*XB.inv();root_symmetry=C7.inv()*symmetry*C7
 assert abs(symmetry.det())==abs(root_symmetry.det())==1
 assert all(s.denom(c)==1 for c in list(symmetry)+list(root_symmetry))
 assert symmetry*C7*symmetry.T==C7
 # Independent graph classification and rank-series enumeration for every I.
 rows=[]
 for mask in range(64):
  I=[i+1 for i in range(6) if mask>>i&1];ds,kind=levi_degrees(I)
  hs=coeffs([1]*(6-len(I))+ds,8)
  m=[0]+[hs[d-1]-(hs[d-7] if d>=7 else 0)-(hs[d-9] if d>=9 else 0) for d in range(1,10)]
  rows.append(dict(retained_nodes=I,levi_type=kind,invariant_degrees=ds,multiplicity_2=m))
 assert next(r for r in rows if not r['retained_nodes'])['multiplicity_2']==[0,1,6,21,56,126,252,461,786,1265]
 assert rows[-1]['multiplicity_2']==[0,1,0,1,0,1,1,1,1,1]
 assert next(r for r in rows if r['retained_nodes']==[1,3,4,5,6])['multiplicity_2']==[0,1,1,2,3,5,7,10,13,17]
 bad=[r['retained_nodes'] for r in rows if any(x.startswith(('D','E')) for x in r['levi_type'].split('+'))]
 assert len(bad)==4 and all(set((2,3,4,5))<=set(I) for I in bad)
 out=dict(result='passed',engine='Plain Python/SymPy, independent 27-form sextic and graph enumeration',refinements=rows,degree_bound=9,kernel_primitivity=kernel_checks,independent_octic_transition=[str(aa),str(bb)],exact_sextic_agrees_up_to_sign=True,odd_coefficient_valuation=2,good_levi_count=60,d5_symmetry_integral=True,bad_levi_sets=bad,scope='Uniform proved low-degree formula only; all-degree multiplicities remain open beyond P7.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7-refinement-relations.json')})
 (ROOT/'experiments/e7-refinement-independent.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED: 64 refinements, degrees 0..9; sextic matches independent 56-weight route')
if __name__=='__main__':main()
