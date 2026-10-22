"""D4+A1 images in a rank-ten quadric tower over the audited P6 quotient.

Avoid a new 7560-dimensional Groebner elimination. The complete D4 weight
source is a known 34-generator module over the D5 source, with common t.
The root algebra below is only a lower image until source completeness
or a separate upper-image comparison proves equality.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
from math import comb
import gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,matrix,prod,pari,inverse_mod,block_diagonal_matrix
ROOT=Path(__file__).resolve().parents[1]
WS=(1,2,4,5,6,8);FD=(0,1,2,3,4,5,6,7,8,4)

def main():
 start=perf_counter();bp=ROOT/'experiments/e7-d5a1-images.json';br=json.loads(bp.read_text());ba=ROOT/br['certificate_archive'];bc=json.loads(gzip.decompress(ba.read_bytes()));tp=ROOT/'experiments/e7-d4-tower-certificates.json.gz';tc=json.loads(gzip.decompress(tp.read_bytes()));sp=ROOT/'experiments/e7-d4a1-sources.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];sc=json.loads(gzip.decompress(sa.read_bytes()));P=PolynomialRing(QQ,6,names=('h','t','v4','p5','v6','v8'),order=TermOrder('wdegrevlex',WS));h,t,E2,E5,E3,E4=P.gens();E1=8*t-3*h*h;G=[P(f) for f in bc['groebner_basis']];norm={int(d):[tuple(e) for e in ee] for d,ee in bc['normal_exponents'].items()};ix={d:{e:i for i,e in enumerate(ee)} for d,ee in norm.items()};NH={0:matrix(QQ,1,1,[1])}
 for row in bc['degrees']:NH[row['degree']]=matrix(QQ,[[QQ(c) for c in line] for line in row['image_bases']['N']])
 NI={d:H.inverse() for d,H in NH.items()};reducers=[]
 for f in G:
  e=tuple(f.lm().exponents()[0]);reducers.append((e,{tuple(u):-c/f.lc() for u,c in f.dict().items() if tuple(u)!=e}))
 @lru_cache(None)
 def mon(e):
  if sum(a*b for a,b in zip(e,WS))>42:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,b in mon(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,QQ.zero())+c*b
    return {e:c for e,c in out.items() if c}
  return {e:QQ.one()}
 def red(f):
  out={}
  for e,c in P(f).dict().items():
   for u,b in mon(tuple(e)).items():out[u]=out.get(u,QQ.zero())+c*b
  return P(out)
 @lru_cache(None)
 def fp(i,j):
  if j>=2:
   out={}
   for shift,c in ((0,E4),(2,-E3),(4,E2),(6,-E1),(8,P.one())):
    for k,v in fp(i+shift,j-2).items():out[k]=out.get(k,P.zero())+c*v
   return out
  if i and j:return {k:E5*v for k,v in fp(i-1,0).items()}
  if i>=9:
   out={}
   for ii,jj,c in ((0,1,E5),(1,0,-E4),(3,0,E3),(5,0,-E2),(7,0,E1)):
    for k,v in fp(i-9+ii,jj).items():out[k]=out.get(k,P.zero())+c*v
   return out
  return {9 if j else i:P.one()}
 def tmul(a,b):
  out={}
  for i,f in a.items():
   for j,g in b.items():
    for k,c in fp((i if i<9 else 0)+(j if j<9 else 0),int(i==9)+int(j==9)).items():out[k]=out.get(k,P.zero())+f*g*c
  return {i:red(f) for i,f in out.items() if f}
 def tpower(f,k):
  out={0:P.one()}
  for _ in range(k):out=tmul(out,f)
  return out
 def addterms(parts):
  out={}
  for c,f in parts:
   for i,g in f.items():out[i]=out.get(i,P.zero())+c*g
  return {i:red(f) for i,f in out.items() if f}
 # The old universal D4 formulas contain the old complement only via
 # w^2. Substitution w^2=-3t gives the new D5+A1 quadratic relation.
 O=PolynomialRing(QQ,8,names=('h','w','v4','v5','v6','v8','t','r'))
 def oldpoly(text):
  out={}
  for e,c in O(text).dict().items():
   assert e[1]%2==0;base=prod(g**k for g,k in zip(P.gens(),(e[0],e[1]//2,e[2],e[3],e[4],e[5])))*(-3)**(e[1]//2)
   for i,f in fp(e[6],e[7]).items():out[i]=out.get(i,P.zero())+c*base*f
  return {i:red(f) for i,f in out.items() if f}
 ninputs=[oldpoly(row['tower_polynomial']) for row in tc['source_preimages']];reference=tc['reference_quartic_generator'];gref=ninputs[reference];c=gref[9];assert c.is_constant() and c;refs=[tpower({0:h/2,1:P(QQ(1)/2)},i) for i in range(9)]+[gref]
 def toref(f):
  q=f.get(9,P.zero())/c;out=[f.get(i,P.zero())-q*gref.get(i,P.zero()) for i in range(9)]+[q];res=[P.zero() for _ in range(10)];res[9]=red(q)
  for i,g in enumerate(out[:9]):
   for j in range(i+1):res[j]+=comb(i,j)*2**j*(-h)**(i-j)*g
  return list(map(red,res))
 # Reconstruct from the universal source polynomials, not the old
 # reduced module coefficients: the old Borel GB mixes odd w powers.
 nm=[]
 for e in tc['module_exponents']:
  ff={0:P.one()}
  for i,k in enumerate(e):
   if k:ff=tmul(ff,tpower(ninputs[i],k))
  degree=sum(k*row['degree'] for k,row in zip(e,tc['source_preimages']));nm.append((degree,toref(ff)))
 raw=[{0:h},{0:h/2,1:P(QQ(1)/2)},{0:t},{0:E1,2:-P.one()},{0:E2,2:-E1,4:P.one()},{9:P.one()},{0:E3,2:-E2,4:E1,6:-P.one()}];raw=[{i:f/QQ(c0) for i,f in ff.items()} for ff,c0 in zip(raw,sc['normalization_factors'])]
 @lru_cache(None)
 def classical(e):
  if not any(e):return {0:P.one()}
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return tmul(raw[j],classical(tuple(u)))
 mg=[]
 for row in sc['degrees']:
  for j in row['new_generator_indices']:
   f=addterms((QQ(line[j]),classical(tuple(e))) for e,line in zip(row['classical_exponents'],row['rational_coordinates']) if QQ(line[j]));mg.append((row['degree'],f))
 products=[[toref(tmul(f,r)) for r in refs] for d,f in mg];Nproducts=[(d,cc) for d,cc in nm]
 @lru_cache(1024)
 def small(text,bd):
  f=P(text);k=sum(a*b for a,b in zip(f.exponents()[0],WS));dd=bd+k;A=matrix(QQ,len(norm[dd]),len(norm[bd]))
  for j,e in enumerate(norm[bd]):
   for u,c0 in f.dict().items():
    for v,b0 in mon(tuple(a+b for a,b in zip(u,e))).items():A[ix[dd][v],j]+=c0*b0
  return NI[dd]*A*NH[bd]
 def positions(d):
  out={};n=0
  for i,k in enumerate(FD):
   if 0<=d-k<=42:out[i]=n;n+=len(norm[d-k])
  return out,n
 @lru_cache(None)
 def operator(j,d):
  k=mg[j][0];out,n=positions(d);ins,m=positions(d-k);A=matrix(QQ,n,m)
  for ii,off in ins.items():
   bd=d-k-FD[ii]
   for oi,f in enumerate(products[j][ii]):
    if f and oi in out:A.set_block(out[oi],off,small(str(f),bd))
  return A
 def bitrank(vectors):
  piv={}
  for v in vectors:
   while v:
    k=v.bit_length()-1
    if k in piv:v^=piv[k]
    else:piv[k]=v;break
  return len(piv)
 def rational_mod2_rank(A):return bitrank([sum((int(c0.numerator())&1)<<j for j,c0 in enumerate(line)) for line in A])
 def lb(A):
  if not A.nrows():return matrix(QQ,0,0)
  k=max([0]+[int(c.denominator().valuation(2)) for c in A.list()]);scaled=2**k*A
  for precision in (64,128,256):
   mod=ZZ(2)**precision;Z=matrix(ZZ,A.nrows(),A.ncols(),[ZZ(c.numerator()*inverse_mod(c.denominator(),mod))%mod for c in scaled.list()]);H=matrix(QQ,pari(Z).mathnfmodid(mod))/2**k;C=H.inverse()*A
   if C.denominator()%2 and rational_mod2_rank(C)==H.nrows():return H
  raise AssertionError('local precision')
 classical.cache_clear()
 def matrix_digest(A):
  digest=sha256()
  for line in A:digest.update(('|'.join(map(str,line))+'\n').encode())
  return digest.hexdigest()
 print('D4+A1 quadric tower ready',len(mg),round(perf_counter()-start,2),flush=True);NB={0:matrix(QQ,1,1,[1])};MB={0:matrix(QQ,1,1,[1])};records=[];certrows=[]
 for d in range(1,9):
  offsets,n=positions(d);blocks=[]
  for k,cc in Nproducts:
   bd=d-k
   if not 0<=bd<=42:continue
   A=matrix(QQ,n,len(norm[bd]))
   for i,f in enumerate(cc):
    if f and i in offsets:A.set_block(offsets[i],0,small(str(f),bd))
   blocks.append(A)
  NA=blocks[0]
  for A in blocks[1:]:NA=NA.augment(A)
  NB[d]=lb(NA);blocks=[operator(j,d)*MB[d-k] for j,(k,f) in enumerate(mg) if k<=d];MA=blocks[0]
  for A in blocks[1:]:MA=MA.augment(A)
  MB[d]=lb(MA);T=NB[d].inverse()*MB[d];assert T.denominator()%2;assert (2*T.inverse()).denominator()%2;mm=n-rational_mod2_rank(T);counts={'1':mm} if mm else {};full=d<=sr['max_degree'] or not counts;r=dict(degree=d,rank=n,smith_multiplicities=counts,full_root_image_proved=full);records.append(r);certrows.append(dict(**r,weight_basis=[[str(c0) for c0 in line] for line in NB[d]],root_lower_basis=[[str(c0) for c0 in line] for line in MB[d]],weight_columns_sha256=matrix_digest(NA),root_lower_columns_sha256=matrix_digest(MA)));print('D4+A1 tower image',r,round(perf_counter()-start,2),flush=True);operator.cache_clear();small.cache_clear()
 archive=ROOT/'experiments/e7-d4a1-tower-pilot-certificates.json.gz';cert=dict(root_generator_degrees=[d for d,f in mg],root_generators=[{str(i):str(f) for i,f in ff.items()} for d,ff in mg],reference_generators=[{str(i):str(f) for i,f in ff.items()} for ff in refs],root_multiplication_components=[[[str(f) for f in cc] for cc in pp] for pp in products],weight_source_modules=[dict(degree=d,components=list(map(str,cc))) for d,cc in nm],fiber_degrees=list(FD),degrees=certrows);archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));report=dict(result='passed',complete=False,degree_bound=8,top_polynomial_degree=50,total_rank=7560,degrees=records,seconds=round(perf_counter()-start,3),scope='Full weight source using the complete D4 finite module and common A1 norm variable. Root algebra generated by full source components through degree 16; higher nonzero counts are bounds until a full-source or upper-image proof.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),bp,sp)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (ba,tp,sa)});(ROOT/'experiments/e7-d4a1-tower-pilot.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
