"""Independent exact 2-adic audit of bounded E7 full-flag kernel lattices.

No Sage or project imports. Integer multinomial expansions use the 27 E6
forms to construct the 56 E7 weights. Saturation witnesses are replayed
modulo 2^64, tracking the loss of precision after every exact halving.
Ranks use Python bit elimination, independently of the Sage matrix engine.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from math import comb,factorial,gcd,prod
from functools import reduce
from time import perf_counter
import gzip,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
ZERO=(0,)*6

def compositions(n,k=6):
 if n<0:return
 if k==1:yield (n,);return
 for j in range(n+1):
  for rest in compositions(n-j,k-1):yield (j,)+rest

def multiply(f,g):
 out={}
 for e,c in f.items():
  for u,v in g.items():
   k=tuple(a+b for a,b in zip(e,u));out[k]=out.get(k,0)+c*v
 return {e:c for e,c in out.items() if c}

def rank2(A):
 packed=np.packbits(np.bitwise_and(A,1).astype(np.uint8),axis=1,bitorder='little');pivots={}
 for packed_row in packed:
  row=int.from_bytes(packed_row.tobytes(),'little')
  while row:
   k=row.bit_length()-1
   if k in pivots:row^=pivots[k]
   else:pivots[k]=row;break
 return len(pivots)

def main():
 start=perf_counter();report_path=ROOT/'experiments/e7-full-flag-kernel.json';primary=json.loads(report_path.read_text());assert primary['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 maximum=max(r['degree'] for r in primary['degrees']);needed=[d for d in (6,8,10,12) if d<=maximum]
 X=[[5,0,4,3,2,1],[-1,0,4,3,2,1],[-1,0,-2,3,2,1],[-1,0,-2,-3,2,1],[-1,0,-2,-3,-4,1],[-1,0,-2,-3,-4,-5]];Z=[-3,-6,-6,-9,-6,-3]
 ls=[tuple(x+sign*z for x,z in zip(row,Z)) for row in X for sign in (-1,1)]+[tuple(-x-y for x,y in zip(u,v)) for u,v in combinations(X,2)]
 assert all(c%2==0 for row in ls for c in row);ell=[tuple(c//2 for c in row) for row in ls]
 # The 27-form construction proposes weights +/-w and +/-(ell-w)/3.
 C=[[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]]
 highest=(0,0,0,0,0,0,1);weights={highest,tuple(-c for c in highest)}
 for row in ell:
  numerator=tuple(sum(C[j][i]*row[i] for i in range(6))-int(j==6) for j in range(7));assert all(c%3==0 for c in numerator)
  v=tuple(c//3 for c in numerator);weights.add(v);weights.add(tuple(-c for c in v))
 assert len(weights)==56;orbit={highest};todo=[highest]
 while todo:
  v=todo.pop()
  for i in range(7):
   t=tuple(v[j]-v[i]*C[i][j] for j in range(7));assert t in weights
   if t not in orbit:orbit.add(t);todo.append(t)
 assert orbit==weights
 lookup=[[[c**j for j in range(maximum+1)] for c in row] for row in ell];moments={0:{ZERO:27}}
 for k in range(1,max(needed)+1):
  f={}
  for e in compositions(k):
   c=(factorial(k)//prod(factorial(j) for j in e))*sum(prod(p[j][e[j]] for j in range(6)) for p in lookup)
   if c:f[e]=c
  moments[k]=f
 assert all(c%36==0 for c in moments[2].values());B={e:-c//36 for e,c in moments[2].items()}
 quad={e+(0,):-c for e,c in B.items()};quad[(0,)*6+(2,)]=1
 assert quad=={tuple(e):int(c) for e,c in cert['quadratic']}
 bpowers=[{ZERO:1}]
 for j in range(max(needed)//2):bpowers.append(multiply(bpowers[-1],B))
 relations={}
 for d in needed:
  out={e+(0,):2*3**d*c for e,c in bpowers[d//2].items()}
  for k in range(d+1):
   j=d-k;f=multiply(moments[k],bpowers[j//2]);scalar=2*comb(d,k)*(-1)**j
   for e,c in f.items():
    key=e+(j%2,);out[key]=out.get(key,0)+scalar*c
  out={e:c for e,c in out.items() if c};g=reduce(gcd,out.values());f={e:c//g for e,c in out.items()}
  assert f=={tuple(e):int(c) for e,c in cert['relations'][str(d)]},d
  relations[d]=f
  print('independent 27-form moment',d,len(f),round(perf_counter()-start,2),flush=True)
 def polynomial_column(rel,monomial):
  out={}
  for e,c in rel.items():
   exponent=tuple(a+b for a,b in zip(e,monomial));k=exponent[6]
   if k<2:out[exponent]=out.get(exponent,0)+c
   else:
    for u,v in B.items():
     key=tuple(a+b for a,b in zip(exponent[:6],u))+(0,);out[key]=out.get(key,0)+c*v
  return {e:c for e,c in out.items() if c}
 pc=[1]
 for d in (2,6,8,10,12,14,18):
  new=[0]*(len(pc)+d-1)
  for i,c in enumerate(pc):
   for j in range(d):new[i+j]+=c
  pc=new
 results=[]
 for case,saved in zip(cert['degrees'],primary['degrees']):
  d=case['degree'];support=[e+(p,) for p in (0,1) for e in compositions(d-p)];index={e:i for i,e in enumerate(support)}
  labels=[(degree,e+(p,)) for degree in relations for p in (0,1) for e in compositions(d-degree-p)]
  assert labels==[(k,tuple(e)) for k,e in case['initial_columns']]
  A=np.zeros((len(support),len(labels)),dtype=np.uint64);mask=(1<<64)-1
  for j,(degree,exponent) in enumerate(labels):
   for e,c in polynomial_column(relations[degree],exponent).items():A[index[e],j]=c&mask
  precision=64
  for stage in case['saturation']:
   positions=stage['replace_columns'];combinations0=stage['combinations'];assert len(positions)==len(combinations0)
   for i,indices in enumerate(combinations0):
    assert [int(p in indices) for p in positions]==[int(j==i) for j in range(len(positions))]
   changes=[]
   for indices in combinations0:
    value=A[:,indices].sum(axis=1,dtype=np.uint64);assert not np.any(value&np.uint64(1))
    changes.append(value>>np.uint64(1))
   for j,value in zip(positions,changes):A[:,j]=value
   precision-=1;assert precision>0
  rank=rank2(A);assert rank==len(labels)==len(support)-pc[d]
  oddrows=[i for i,e in enumerate(support) if e[6]];odd_rank=rank2(A[oddrows,:]);m=comb(d+4,5)-odd_rank
  assert rank==saved['kernel_rank'] and odd_rank==saved['odd_relation_rank'] and m==saved['multiplicity_2']
  results.append(dict(degree=d,kernel_rank=rank,odd_relation_rank=odd_rank,multiplicity_2=m,minimum_exact_binary_precision=precision))
  print('independent exact degree',d,'kernel rank',rank,'multiplicity',m,round(perf_counter()-start,2),flush=True)
 output=dict(result='passed',engine='Plain Python integer moments and NumPy 2^64 coefficient congruences; own Python bit ranks; no Sage or project imports',degrees=results,seconds=round(perf_counter()-start,3),scope='Rebuilds all initial relations and coefficients, replays every integral saturation witness, and checks full and odd ranks in the recorded degrees only.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),report_path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-full-flag-kernel-audit.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED',output['seconds'],flush=True)
if __name__=='__main__':main()
