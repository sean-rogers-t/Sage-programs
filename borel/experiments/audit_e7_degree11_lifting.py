"""Independent degree-eleven half-sum and invariant-lifting audit.

No Sage or runtime imports. Integer multinomial expansions use 27 E6
forms to construct the 56 E7 weights. The cubic, quintic and exact half-sum
are rebuilt, with Python bit ranks and a linear-ideal restriction algorithm
independent of the primary Sage matrix and Groebner engines.
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
 start=perf_counter();report_path=ROOT/'experiments/e7-degree11-lifting.json';primary=json.loads(report_path.read_text());assert primary['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 maximum=11;needed=[6,8,10]
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
 # END OF SHARED 27-FORM CONSTRUCTION; all new checks below use integers.
 def parity(f):return {e:1 for e,c in f.items() if c%2}
 def add(f,g):
  out=dict(f)
  for e,c in g.items():out[e]=out.get(e,0)+c
  return {e:c for e,c in out.items() if c}
 def frobenius_derivative(f,power):
  out={}
  for e,c in f.items():
   for i,k in enumerate(e):
    if k%2:
     u=list(e);u[i]+=power-1;u=tuple(u);out[u]=out.get(u,0)+c
  return parity(out)
 g3={e+(0,):c for e,c in frobenius_derivative(B,2).items()}
 g5={e+(0,):c for e,c in frobenius_derivative(B,4).items()}
 assert g3=={tuple(e):int(c) for e,c in cert['cubic']}
 assert g5=={tuple(e):int(c) for e,c in cert['quintic']}
 b7={e+(0,):c for e,c in B.items()}
 assert parity(relations[6])==parity(multiply(g3,g3))
 assert parity(relations[8])==parity(multiply(g3,g5))
 odd=lambda f:{e[:-1]+(0,):c for e,c in f.items() if e[-1]}
 c5=odd(relations[6]);c7=odd(relations[8]);c9=odd(relations[10])
 assert all(c%4==0 for c in c5.values()) and all(c%2==0 for c in c7.values()) and all(c%4==0 for c in c9.values())
 assert parity({e:c//4 for e,c in c5.items()})==parity(add(g5,multiply(b7,g3)))
 assert parity({e:c//2 for e,c in c7.items()})==parity(add(multiply(b7,g5),multiply(multiply(b7,b7),g3)))
 numerator=add(multiply(relations[6],g5),multiply(relations[8],g3));assert all(c%2==0 for c in numerator.values())
 H={e:c//2 for e,c in numerator.items() if c}
 assert H=={tuple(e):int(c) for e,c in cert['half_sum']}
 assert parity(odd(H))==parity(add(multiply(b7,relations[8]),multiply(multiply(b7,b7),relations[6])))
 support=[e+(p,) for p in (0,1) for e in compositions(11-p)];ix={e:i for i,e in enumerate(support)}
 labels=[(d,e+(p,)) for d in relations for p in (0,1) for e in compositions(11-d-p)]
 A=np.zeros((len(support),len(labels)+1),dtype=np.uint8)
 for j,(d,e) in enumerate(labels):
  for u,c in polynomial_column(relations[d],e).items():A[ix[u],j]=c%2
 for e,c in H.items():A[ix[e],-1]=c%2
 r0=rank2(A[:,:-1]);r1=rank2(A);assert (r0,r1)==(461,462)
 syzygy=np.array([g5.get(e,0) if d==6 else g3.get(e,0) if d==8 else 0 for d,e in labels],dtype=np.uint64)
 assert syzygy.any() and not np.any((A[:,:-1]@syzygy)%2)
 pc=[1]
 for d in (2,6,8,10,12,14,18):
  new=[0]*(len(pc)+d-1)
  for i,c in enumerate(pc):
   for j in range(d):new[i+j]+=c
  pc=new
 assert len(labels)==len(support)-pc[11]==462
 # Independent linear-ideal quotient, not the primary Sage Groebner reducer.
 from sympy import Matrix
 C6=Matrix(C)[:6,:6];ci=3*C6.inv();assert C6.det()==3
 def rank_bits(columns):
  piv={}
  for x in columns:
   while x:
    j=x.bit_length()-1
    if j in piv:x^=piv[j]
    else:piv[j]=x;break
  return len(piv)
 def restrict(f,ls):
  piv={}
  for x in ls:
   while x:
    j=x.bit_length()-1
    if j in piv:x^=piv[j]
    else:piv[j]=x;break
  for j in sorted(piv):
   for k in piv:
    if k>j and piv[k]>>j&1:piv[k]^=piv[j]
  images=[]
  for j in range(6):
   bits=(piv[j]^(1<<j)) if j in piv else 1<<j
   images.append({tuple(int(k==i) for k in range(6)):1 for i in range(6) if bits>>i&1})
  out={}
  for e,c in f.items():
   term={ZERO:c}
   for i,k in enumerate(e):
    for z in range(k):term=parity(multiply(term,images[i]))
   out=parity(add(out,term))
  return out
 bad=[]
 for I in ((1,2,3,4,5,6),(1,2,3,4,5),(2,3,4,5,6),(2,3,4,5)):
  omitted=[j for j in range(6) if j+1 not in I]
  ls=[sum((int(ci[i,j])%2)<<i for i in range(6)) for j in omitted]
  assert rank_bits(ls)==len(omitted)
  restriction=restrict({e[:-1]:c for e,c in g3.items()},ls);assert restriction
  bad.append(dict(retained_nodes=list(I),cubic_restriction=[[list(e),c] for e,c in sorted(restriction.items())],cubic_does_not_lift=True))
 family=[];good=[]
 for mask in range(127):
  I=[j for j in range(7) if mask>>j&1]
  if any(j not in I for j in (1,4,6)):family.append(tuple(j+1 for j in I))
  is_A=all(sum(C[i][j]==-1 for j in I)<=2 for i in I)
  cols=[sum((C[i][j]%2)<<i for i in range(7)) for j in I]
  if is_A and rank_bits(cols)==len(I):good.append(tuple(j+1 for j in I))
 assert len(family)==112 and len(good)==108
 assert set(family)-set(good)=={tuple(c['retained_nodes']) for c in bad}
 output=dict(result='passed',degree=11,uniform_family_size=112,initial_kernel_rank_mod2=r0,primitive_kernel_rank_mod2=r1,mod2_syzygy_dimension=1,quadratic_Frobenius_identities_verified=True,half_sum_derivative_identity_verified=True,good_Levi_count=108,exceptional_noncubic_lifting_cases=bad,seconds=round(perf_counter()-start,3),engine='Native 27-form integer moments, exact integer half-sum, own bit ranks and linear-ideal restriction; SymPy only for 6x6 inverse; no Sage or runtime imports',scope='Independently verifies the coefficient identities, primitive full kernel, one-dimensional saturation defect and four cubic nonlifting obstructions used in the uniform degree-eleven proof. Reads no relative Smith table and does not prove an all-degree E7 structural result.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),report_path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-degree11-lifting-audit.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED',output['seconds'],flush=True)
if __name__=='__main__':main()
