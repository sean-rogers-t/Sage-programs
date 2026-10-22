"""Independent exact audit of E7 quadratic-derivation kernel stability.

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
 start=perf_counter();report_path=ROOT/'experiments/e7-quadratic-derivation.json';primary=json.loads(report_path.read_text());assert primary['result']=='passed'
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
 oddvals={str(d):min((abs(c)&-abs(c)).bit_length()-1 for e,c in f.items() if e[-1]) for d,f in relations.items()}
 assert oddvals=={'6':2,'8':1,'10':2,'12':2}
 # Audit the lattice hypothesis and canonical mod-two derivative in all 112 cases.
 # SymPy is used only for seven-dimensional rational matrix arithmetic here.
 from sympy import Matrix,eye
 CM=Matrix(C);lattice_checks=[]
 for i,n in ((1,7),(4,15),(6,3)):
  K=[j for j in range(7) if j!=i];basis=CM[:,K].row_join(eye(7)[:,i])
  assert abs(int(basis.det()))==n
  quad=-(n*CM.extract(K,K).inv()) # twice the quadratic B_i coefficients
  assert all(c.q==1 for c in quad) and all(int(quad[j,j])%2==0 for j in range(6))
  # In fundamental-weight coordinates, every omitted odd node gives the same chi.
  deriv=(Matrix([[0]*6+[1]])*basis.inv())
  assert [int(c.p)%2 for c in deriv]==[0,1,0,0,1,0,1]
  lattice_checks.append(dict(omitted_node=i+1,odd_index=n,quadratic_integral=True,canonical_derivative=True))
 family=[mask for mask in range(127) if any(not(mask>>j&1) for j in (1,4,6))]
 assert len(family)==112
 results=[];previous_support=[];previous_pivots={}
 def pivots_of(columns):
  pivots={}
  for row in columns:
   while row:
    k=row.bit_length()-1
    if k in pivots:row^=pivots[k]
    else:pivots[k]=row;break
  return pivots
 def packed_columns(A):
  return [int.from_bytes(row.tobytes(),'little') for row in np.packbits((A.T&np.uint64(1)).astype(np.uint8),axis=1,bitorder='little')]
 for case,saved in zip(cert['degrees'],primary['degrees']):
  d=case['degree'];support=[e+(p,) for p in (0,1) for e in compositions(d-p)];index={e:i for i,e in enumerate(support)}
  assert d>10 or not case['saturation']
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
  oddrows=[i for i,e in enumerate(support) if e[6]];evenrows=[i for i,e in enumerate(support) if not e[6]]
  odd_rank=rank2(A[oddrows,:]);even_rank=rank2(A[evenrows,:]);assert odd_rank+even_rank==rank
  if d:
   ix={e:i for i,e in enumerate(previous_support)}
   derivative=np.zeros((len(previous_support),len(labels)),dtype=np.uint64)
   for i in oddrows:derivative[ix[support[i][:-1]+(0,)],:]=A[i,:]
   for value in packed_columns(derivative):
    while value:
     k=value.bit_length()-1;assert k in previous_pivots,(d,k)
     value^=previous_pivots[k]
  assert (rank,even_rank,odd_rank)==(saved['kernel_rank'],saved['even_projection_rank'],saved['odd_projection_rank'])
  previous_support=support;previous_pivots=pivots_of(packed_columns(A))
  results.append(dict(degree=d,kernel_rank=rank,even_projection_rank=even_rank,odd_projection_rank=odd_rank,derivative_preserves_full_kernel=True,minimum_exact_binary_precision=precision))
  print('independent derivation',d,'kernel rank',rank,round(perf_counter()-start,2),flush=True)
 output=dict(result='passed',maximum_degree=maximum,uniform_family_maximum_degree=10,uniform_family_size=112,primitive_relation_odd_valuations=oddvals,lattice_checks=lattice_checks,engine='Plain Python 27-form integer moments, NumPy 2^64 congruences and own bit elimination; SymPy only for lattice matrices; no Sage or project imports',degrees=results,seconds=round(perf_counter()-start,3),scope='Fresh 27-form weight construction, replayed exact saturation witnesses and derivative membership in the primitive kernel, degrees 0..13; no relative multiplicity input. Lattice splitting, primitive relation parity and unique unsaturated relation multiples through 10 certify the uniform invariant-source argument in the note.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),report_path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-quadratic-derivation-audit.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED',output['seconds'],flush=True)
if __name__=='__main__':main()
