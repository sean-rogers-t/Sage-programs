"""Independent Schur/Pieri image checks for the odd-block type A theorem.

Native Python only. Root generators are centered elementary functions in
each of two odd blocks. Compute their full mod-two image by multiplication
closure in the integral Schur basis, not by invariant-kernel dimensions.
The proved quadratic lift makes this reduction sufficient for relative
2-primary factors. This does not check odd-primary factors of SL(n).
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from math import comb
from functools import lru_cache
import json
ROOT=Path(__file__).resolve().parents[1]

def partitions(k,b,cap=None):
 if k==0:yield ();return
 if cap is None:cap=b
 for j in range(cap,-1,-1):
  for tail in partitions(k-1,b,j):yield (j,)+tail

def transpose(lam,b):return tuple(sum(x>=j for x in lam) for j in range(1,b+1))

def check(n,k):
 b=n-k;assert n%4==2 and k%2 and b%2
 basis={d:[] for d in range(k*b+1)}
 for lam in partitions(k,b):basis[sum(lam)].append(lam)
 index={d:{lam:i for i,lam in enumerate(v)} for d,v in basis.items()}
 @lru_cache(None)
 def vertical(lam,l,width):
  out=[]
  for rows in combinations(range(len(lam)),l):
   mu=tuple(x+int(i in rows) for i,x in enumerate(lam))
   if mu[0]<=width and all(mu[i]>=mu[i+1] for i in range(len(mu)-1)):out.append(mu)
  return tuple(out)
 def strip(lam,l,horizontal):
  if not horizontal:return vertical(lam,l,b)
  return tuple(transpose(mu,k) for mu in vertical(transpose(lam,b),l,k))
 def pieri(v,d,l,horizontal=False):
  if d+l>k*b:return 0
  out=0
  while v:
   i=(v&-v).bit_length()-1;v^=1<<i
   for mu in strip(basis[d][i],l,horizontal):out^=1<<index[d+l][mu]
  return out
 def centered(v,d,j,size,horizontal):
  out=0
  for l in range(j+1):
   if comb(size-l,j-l)%2:
    value=pieri(v,d,l,horizontal)
    for s in range(j-l):value=pieri(value,d+l+s,1)
    out^=value
  return out
 generators=[(j,size,horizontal) for size,horizontal in ((k,False),(b,True)) for j in range(2,size+1)]
 images={0:[1]};rows=[]
 for d in range(k*b+1):
  if d:
   pivots={}
   for j,size,horizontal in generators:
    for value in images.get(d-j,[]):
     value=centered(value,d-j,j,size,horizontal)
     while value:
      p=value.bit_length()-1
      if p in pivots:value^=pivots[p]
      else:pivots[p]=value;break
   images[d]=list(pivots.values())
  rows.append(dict(degree=d,cohomology_rank=len(basis[d]),root_image_rank_mod2=len(images[d]),multiplicity_2=len(basis[d])-len(images[d])))
 # Division of q*P by 1+q via coefficient recurrence, after the image calculation.
 predicted=[]
 for d,row in enumerate(rows):
  predicted.append((len(basis[d-1]) if d else 0)-(predicted[-1] if d else 0))
 assert [r['multiplicity_2'] for r in rows]==predicted
 assert sum(predicted)*2==comb(n,k)
 return dict(n=n,k=k,top_degree=k*b,rank=comb(n,k),degrees=rows)

def compositions(n):
 if not n:yield ();return
 for j in range(1,n+1):
  for rest in compositions(n-j):yield (j,)+rest

def main():
 cases=[check(n,k) for n,k in ((6,1),(6,3),(10,1),(10,3),(10,5))]
 inventories={}
 for n in (6,10):
  cs=list(compositions(n));odd=[c for c in cs if any(j%2 for j in c)]
  inventories[str(n)]=dict(total_standard_parabolics=len(cs),odd_block_family=len(odd),all_even_proper=len(cs)-len(odd)-1)
  assert len(odd)==2**(n-1)-2**(n//2-1)
 # Generic quadratic lifting alone does not imply the desired splitting:
 # R[a,h]/(h-a,a^2) has h^2=0 and M=N, so d(h)=1,d(a)=0 cannot descend.
 out=dict(result='passed',engine='Native Python partitions, vertical/horizontal Pieri and bit elimination; no Sage/runtime/table imports',cases=cases,inventories=inventories,scope='Complete 2-primary root-image calculations on five odd-block Grassmannians. The arbitrary-rank theorem is proved in the note; these are independent finite checks, not its proof. No odd-primary or all-even-block formula asserted.',source_sha256={'experiments/'+Path(__file__).name:sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
 (ROOT/'experiments/type-a-quadratic-derivation-audit.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASSED',[(c['n'],c['k'],c['rank']) for c in cases],inventories)
if __name__=='__main__':main()
