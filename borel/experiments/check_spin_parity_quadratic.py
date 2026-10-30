"""Independent quadratic-straightening check, without Schur-Q diagonal lattices or Pieri matrices."""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import *
from functools import lru_cache
from itertools import combinations
from time import time

def case(n):
 t=time(); top=n*(n+1)//2
 bydeg=[[] for _ in range(top+1)]
 for k in range(n+1):
  for x in combinations(range(1,n+1),k):bydeg[sum(x)].append(x)
 index=[{x:i for i,x in enumerate(b)} for b in bydeg]
 @lru_cache(None)
 def nf(x):
  for j in range(len(x)-1):
   if x[j]==x[j+1]:
    i=x[j]; rest=x[:j]+x[j+2:]; out={}
    terms=[]
    if 2*i<=n:terms.append(((-1)**(i+1),(2*i,)))
    for k in range(1,i):
     if i+k<=n:terms.append((-2*(-1)**k,(i-k,i+k)))
    for c,ys in terms:
     for key,v in nf(tuple(sorted(rest+ys))):out[key]=out.get(key,0)+c*v
    return tuple((k,v) for k,v in out.items() if v)
  return ((x,1),)
 @lru_cache(None)
 def col(j,d,i):
  if d+j>top:return vector(ZZ,0)
  out=vector(ZZ,len(bydeg[d+j])); factor=1 if j==1 else 2
  for key,v in nf(tuple(sorted(bydeg[d][i]+(j,)))):out[index[d+j][key]]+=factor*v
  return out
 def mul(j,d,v):
  if d+j>top:return vector(ZZ,0)
  out=vector(ZZ,len(bydeg[d+j]))
  for i,c in enumerate(v):
   if c:out+=c*col(j,d,i)
  return out
 H=[matrix(ZZ,[[1]])]; delta=[matrix(GF(2),0,1)]
 for d in range(1,top+2):
  aa=[]; bb=[]
  for j in range(1,min(n,d)+1):
   for k,b in enumerate(H[d-j].rows()):
    aa.append(mul(j,d-j,b))
    db=vector(ZZ,len(bydeg[d-1]) if d-1<=top else 0)
    if j==1:db+=b
    if d-j:
     prev=vector(ZZ,[ZZ(c) for c in delta[d-j].column(k)])*H[d-j-1]
     db+=mul(j,d-j-1,prev)
    bb.append(db)
  width=len(bydeg[d]) if d<=top else 0
  A=matrix(ZZ,aa) if aa else matrix(ZZ,0,width)
  h=A.row_module().basis_matrix() if width else matrix(ZZ,0,0)
  T=h.change_ring(QQ).solve_left(A.change_ring(QQ)).transpose() if h.nrows() else matrix(ZZ,0,A.nrows())
  T=matrix(GF(2),T)
  B=matrix(ZZ,bb) if bb else matrix(ZZ,0,len(bydeg[d-1]))
  C=matrix(GF(2),H[d-1].change_ring(QQ).solve_left(B.change_ring(QQ)).transpose()) if H[d-1].nrows() else matrix(GF(2),0,A.nrows())
  if (C*T.right_kernel().basis_matrix().transpose()).is_zero():
   dd=T.transpose().solve_right(C.transpose()).transpose()
   H.append(h);delta.append(dd)
  else:
   print('FAIL',n,d,'time',round(time()-t,2),flush=True);return d
 print('STABLE',n,'time',round(time()-t,2),flush=True);return None

if __name__ == '__main__':
 root=Path(__file__).resolve().parents[1]
 cases=[dict(rank=n,first_failure_degree=case(n)) for n in (3,4,7,8,11,12)]
 primary=json.loads((root/'experiments/spin-parity-recursion.json').read_text())
 expected={c['rank']:c['first_failure_degree'] for c in primary['cases']}
 assert all(c['first_failure_degree']==expected[c['rank']] for c in cases)
 source=Path(__file__)
 report=dict(result='passed',date='2026-10-22',cases=cases,
  method='Integral quadratic Schur-P straightening in squarefree special-class monomials; generate image by z and 2*tau_j; Hermite and derivative relation matrices',
  independence='Does not use Schur-Q diagonal image lattices, one-box Pieri matrices, polynomial saturation or saved image bases. Shares Sage exact arithmetic and the classical integral Grassmannian presentation.',
  source_sha256={'experiments/check_spin_parity_quadratic.py':sha256(source.read_bytes().replace(b"\r\n",b"\n")).hexdigest()})
 (root/'experiments/spin-parity-quadratic-check.json').write_text(json.dumps(report,indent=2)+'\n')
