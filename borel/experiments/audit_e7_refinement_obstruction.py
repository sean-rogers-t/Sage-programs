"""Independent sparse-integer audit of the E6 invariant-ring obstruction.

No Sage, SymPy, or project imports. Moments use multinomial coefficients;
primitivity uses bit elimination rather than a computer algebra matrix rank.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from math import factorial,prod,comb
import json
ROOT=Path(__file__).resolve().parents[1]
ZERO=(0,)*6

def compositions(n,k=6):
 if k==1:yield (n,);return
 for j in range(n+1):
  for rest in compositions(n-j,k-1):yield (j,)+rest

def weighted(ds,n):
 if not ds:
  if n==0:yield ()
  return
 for j in range(n//ds[0]+1):
  for rest in weighted(ds[1:],n-j*ds[0]):yield (j,)+rest

def multiply(f,g):
 out={}
 for e,c in f.items():
  for u,v in g.items():
   k=tuple(a+b for a,b in zip(e,u));out[k]=out.get(k,0)+c*v
 return {e:c for e,c in out.items() if c}

def power(f,n):
 out={ZERO:1}
 while n:
  if n&1:out=multiply(out,f)
  n//=2
  if n:f=multiply(f,f)
 return out

def half_sum(f,g):
 keys=set(f)|set(g);out={e:f.get(e,0)+g.get(e,0) for e in keys}
 assert all(c%2==0 for c in out.values())
 return {e:c//2 for e,c in out.items() if c}

def rank_mod2(polys):
 pivots={}
 for e in set().union(*(set(f) for f in polys)):
  row=sum((f.get(e,0)%2)<<j for j,f in enumerate(polys))
  while row:
   p=row.bit_length()-1
   if p in pivots:row^=pivots[p]
   else:pivots[p]=row;break
 return len(pivots)

def main():
 X=[[5,0,4,3,2,1],[-1,0,4,3,2,1],[-1,0,-2,3,2,1],[-1,0,-2,-3,2,1],[-1,0,-2,-3,-4,1],[-1,0,-2,-3,-4,-5]];Z=[-3,-6,-6,-9,-6,-3]
 ls=[tuple(x+sign*z for x,z in zip(row,Z)) for row in X for sign in (-1,1)]
 ls +=[tuple(-x-y for x,y in zip(u,v)) for u,v in combinations(X,2)]
 assert len(ls)==27 and all(c%2==0 for row in ls for c in row)
 weights=[tuple(c//2 for c in row) for row in ls]
 degrees=(2,5,6,8,9,12);contents=(144,5760,3456,9216,387072,663552)
 lookup=[[[c**j for j in range(13)] for c in row] for row in ls];fs={}
 for d,g in zip(degrees,contents):
  f={}
  for e in compositions(d):
   mult=factorial(d)//prod(factorial(j) for j in e)
   c=mult*sum(prod(pows[j][e[j]] for j in range(6)) for pows in lookup)
   assert c%g==0
   if c:f[e]=c//g
  fs[d]=f
 fs[8]=half_sum(fs[8],power(fs[2],4));fs[12]=half_sum(fs[12],power(fs[6],2))
 checks=[]
 for d in degrees:
  exps=list(weighted(degrees,d));polys=[]
  for e in exps:
   p={ZERO:1}
   for k,j in zip(degrees,e):
    if j:p=multiply(p,power(fs[k],j))
   polys.append(p)
  rank=rank_mod2(polys);assert rank==len(polys);checks.append(dict(degree=d,rank=rank))
 vals=[fs[d].get((0,0,0,0,0,d),0) for d in degrees];assert all(v%2==0 for v in vals)
 C=[[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]]
 for i in range(6):
  changed=[]
  for row in weights:
   v=list(row);v[i]-=sum(row[j]*C[j][i] for j in range(6));changed.append(tuple(v))
  assert set(changed)==set(weights)
 wvals=[row[-1] for row in weights];assert sum(v%2 for v in wvals)==16
 es=[1]+[0]*16
 for v in wvals:
  for j in range(16,0,-1):es[j]+=v*es[j-1]
 assert es[16]%2==1
 primary_path=ROOT/'experiments/e7-refinement-invariant-obstruction.json';primary=json.loads(primary_path.read_text())
 assert primary['result']=='passed' and primary['basic_values']==vals and primary['chern_invariant_value']==es[16]
 assert [r['rank'] for r in checks]==[r['primitive_mod2_rank'] for r in primary['primitive_basic_degree_certificates']]
 result=dict(result='passed',engine='Plain Python integer multinomial expansions, sparse convolution and bit elimination; no Sage, SymPy, or project imports',basic_values=vals,basic_degree_ranks=checks,chern_invariant_degree=16,chern_invariant_value=es[16],permutation_reflections=6,scope='Independent certificate for the source-ring obstruction. No assertion that the proposed all-degree image formula is false.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path)})
 (ROOT/'experiments/e7-refinement-invariant-obstruction-audit.json').write_text(json.dumps(result,indent=2)+'\n');print('PASSED: independent integer source primitivity and degree-16 witness')
if __name__=='__main__':main()
