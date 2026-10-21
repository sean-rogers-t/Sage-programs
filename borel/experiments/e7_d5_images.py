"""Full D5-refinement root/weight images from the coefficient-complete source."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from math import factorial
from itertools import combinations
import json,gzip
from sage.all import QQ,ZZ,PolynomialRing,TermOrder,SymmetricFunctions,Partitions,prod,binomial,matrix,vector,CartanMatrix,pari,GF,inverse_mod
from e7_d5_source import context,powers,WEIGHTS
ROOT=Path(__file__).resolve().parents[1]
DEGREES=(1,1,4,5,6,8)
def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def main():
 start=perf_counter();P=PolynomialRing(QQ,6,names=('h','w','v4','v5','v6','v8'),order=TermOrder('wdegrevlex',DEGREES));h,w,v4,v5,v6,v8=P.gens();e1=-3*h*h-QQ(8)/3*w*w;es=[P.one(),e1,v4,v6,v8,v5*v5]
 SF=SymmetricFunctions(QQ);MM=SF.monomial();EE=SF.elementary()
 def symmetric(part):
  return P(sum(c*prod(es[j] if j<=5 else P.zero() for j in mu) for mu,c in EE(MM[part])))
 spin={}
 for k in range(19):
  if k%2==0:n=k//2;offset=0
  elif k>=5:n=(k-5)//2;offset=1
  else:spin[k]=P.zero();continue
  spin[k]=P(sum((QQ(factorial(k))/prod(factorial(2*j+offset) for j in list(part)+[0]*(5-len(part))))*symmetric(part) for part in Partitions(n,max_length=5))*(v5 if offset else 1))
 ps=[P(5)]
 for d in range(1,10):ps.append(P(sum((-1)**(j-1)*es[j]*ps[d-j] for j in range(1,min(d,6)))+((-1)**(d-1)*d*es[d] if d<=5 else 0)))
 # The 27 E6 weights: h; -h/2 +/- x_i/2; h/4 + sum(eps_i x_i)/4, product eps=1.
 moments={}
 for d in (2,6,8,10,12,14,18):
  f=2*w**d+2*(h-w/3)**d
  f+=4*sum(binomial(d,k)*(-h/2-w/3)**(d-k)*ps[k//2]/2**k for k in range(0,d+1,2))
  f+=32*sum(binomial(d,k)*(h/4-w/3)**(d-k)*spin[k]/4**k for k in range(d+1))
  moments[d]=P(f)
 assert moments[2]==0
 relations=[moments[d] for d in (6,8,10,12,14,18)]
 G=P.ideal(relations).groebner_basis();assert P.ideal(G).vector_space_dimension()==1512
 print('D5 quotient',len(G),'relations; dimension1512',round(perf_counter()-start,2),flush=True)
 leading=[tuple(f.lm().exponents()[0]) for f in G];zero=(0,)*6;todo=[zero];standard={zero}
 while todo:
  e=todo.pop()
  for j in range(6):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   if u in standard or any(all(a<=b for a,b in zip(l,u)) for l in leading):continue
   standard.add(u);todo.append(u)
 normal={d:sorted(e for e in standard if sum(a*b for a,b in zip(e,DEGREES))==d) for d in range(44)}
 assert len(standard)==1512 and sum(map(len,normal.values()))==1512
 assert max(sum(a*b for a,b in zip(e,DEGREES)) for e in standard)==43
 # Recover each new invariant as a rational expression in the classical six generators.
 T,orig=context();source=json.loads((ROOT/'experiments/e7-d5-source.json').read_text());new=[(T('h'),1)]+[(T(f),r['degree']) for r in source['degrees'] for f in r['new_generators']]
 target=[h,e1,v4/2,v5,v6/2,v8];mapped=[];preimages=[]
 for f,d in new:
  exps=list(powers((1,2,4,5,6,8),d));fs=[prod(g**k for (g,_),k in zip(orig,e)) for e in exps]
  support=list(powers(WEIGHTS,d));columns=[{tuple(e):c for e,c in g.dict().items()} for g in fs];D={tuple(e):c for e,c in f.dict().items()}
  A=matrix(QQ,[[g.get(e,0) for g in columns] for e in support]);piv=list(A.transpose().pivots());v=A.matrix_from_rows(piv).inverse()*vector(QQ,[D.get(support[i],0) for i in piv]);assert A*v==vector(QQ,[D.get(e,0) for e in support])
  image=P(sum(c*prod(g**k for g,k in zip(target,e)) for c,e in zip(v,exps)));mapped.append((image,d));preimages.append(dict(degree=d,polynomial=str(f),rational_exponents=[list(e) for e in exps],rational_coefficients=list(map(str,v))))
 assert len(mapped)==10
 gen={'N':[(w,1)]+mapped,'M':[(2*w,1)]+mapped}
 mono={d:[P({e:1}) for e in normal[d]] for d in normal};index={d:{e:i for i,e in enumerate(normal[d])} for d in normal}
 from functools import lru_cache
 reducers=[]
 for f in G:
  lead=tuple(f.lm().exponents()[0]);c=f.monomial_coefficient(P({lead:1}));reducers.append((lead,{tuple(e):-v/c for e,v in f.dict().items() if tuple(e)!=lead}))
 @lru_cache(None)
 def normal_monomial(e):
  d=sum(a*b for a,b in zip(e,DEGREES))
  if d>43:return {}
  if e in index[d]:return {e:QQ.one()}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,k in normal_monomial(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,QQ.zero())+c*k
    return {u:c for u,c in out.items() if c}
  raise AssertionError('nonstandard monomial without divisor')
 def coords(f,d):
  out=[QQ.zero()]*len(normal[d])
  for e,c in P(f).dict().items():
   for u,k in normal_monomial(tuple(e)).items():out[index[d][u]]+=c*k
  return vector(QQ,out)
 basis={k:{0:matrix(QQ,1,1,[1])} for k in gen};records=[]
 from functools import lru_cache
 @lru_cache(None)
 def multiplication(i,d):
  g,k=gen['N'][i];return matrix(QQ,[coords(g*m,d) for m in mono[d-k]]).transpose()
 def local_basis(A):
  k=max([0]+[int(c.denominator().valuation(2)) for c in A.list()]);A0=2**k*A
  for precision in (32,64,128,256,512):
   modulus=ZZ(2)**precision;inverses={}
   values=[]
   for c in A0.list():
    den=c.denominator()%modulus
    if den not in inverses:inverses[den]=inverse_mod(den,modulus)
    values.append(ZZ(c.numerator()*inverses[den])%modulus)
   Z=matrix(ZZ,A.nrows(),A.ncols(),values);H=matrix(ZZ,pari(Z).mathnfmodid(modulus)).change_ring(QQ)/2**k
   C=H.inverse()*A
   if C.denominator()%2 and matrix(GF(2),C).rank()==H.nrows():return H
  raise RuntimeError('modular image precision bound')
 for d in range(1,44):
  for kind,gs in gen.items():
   blocks=[(2 if kind=='M' and i==0 else 1)*multiplication(i,d)*basis[kind][d-k] for i,(g,k) in enumerate(gs) if k<=d]
   A=blocks[0]
   for block in blocks[1:]:A=A.augment(block)
   H=local_basis(A);assert H.is_square() and H.nrows()==len(normal[d]);basis[kind][d]=H
  B=basis['N'][d].inverse()*basis['M'][d];den=B.denominator();assert den%2
  vals=[int(c.valuation(2)) for c in matrix(ZZ,den*B).elementary_divisors()];assert max(vals,default=0)<=1
  m=vals.count(1);records.append(dict(degree=d,rank=len(normal[d]),multiplicity_2=m,weight_basis=[[str(c) for c in row] for row in basis['N'][d]],root_basis=[[str(c) for c in row] for row in basis['M'][d]]))
  print('D5 image',d,'rank',len(normal[d]),'factors2',m,round(perf_counter()-start,2),flush=True)
 cert=dict(relations=[str(f) for f in relations],groebner_basis=[str(f) for f in G],normal_exponents={str(d):[list(e) for e in es] for d,es in normal.items()},image_degrees=records,source_generators={kind:[[str(f),d] for f,d in gs] for kind,gs in gen.items()},preimages=preimages)
 archive=ROOT/'experiments/e7-d5-images-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0))
 q=PolynomialRing(ZZ,'q').gen();WG=prod(sum(q**i for i in range(d)) for d in (2,6,8,10,12,14,18));WI=prod(sum(q**i for i in range(d)) for d in (2,4,5,6,8));cohom=WG//WI;prediction=q*cohom//(1+q);coeff=[0]+[r['multiplicity_2'] for r in records]
 agrees=all(coeff[d]==prediction[d] for d in range(44));print('FORMULA AGREES',agrees,flush=True)
 report=dict(result='passed',coefficient_ring='Z_(2)',retained_nodes=[2,3,4,5,6],dimension=43,total_rank=1512,multiplicity_coefficients=coeff,half_poincare_prediction=[int(prediction[d]) for d in range(44)],prediction_agrees=agrees,certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in [Path(__file__).resolve(),ROOT/'experiments/e7_d5_source.py',ROOT/'experiments/e7-d5-source.json']})
 (ROOT/'experiments/e7-d5-images.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
