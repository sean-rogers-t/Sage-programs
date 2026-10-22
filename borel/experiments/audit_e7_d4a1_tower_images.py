"""Independent native verification of complete D4+A1 Borel images.

Reconstruct the quadric tower and every source generator. Verify full
weight bases by local primitive ranks and full root images modulo four.
Only addition and multiplication are used over composite residue rings.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb
from time import perf_counter
import argparse,gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mat,fmpq_mpoly_ctx,nmod_mat,fmpz_mod_mat,fmpz_mod_ctx
ROOT=Path(__file__).resolve().parents[1]
WS=(1,2,4,5,6,8);FD=(0,1,2,3,4,5,6,7,8,4)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');args=ap.parse_args();start=perf_counter();rp=ROOT/('experiments/e7-d4a1-tower-'+('pilot' if args.pilot else 'final')+'.json');r=json.loads(rp.read_text());archive=ROOT/r['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==r['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()));bp=ROOT/'experiments/e7-d5a1-images.json';br=json.loads(bp.read_text());ba=ROOT/br['certificate_archive'];bc=json.loads(gzip.decompress(ba.read_bytes()));tp=ROOT/'experiments/e7-d4-tower-certificates.json.gz';tc=json.loads(gzip.decompress(tp.read_bytes()));sp=ROOT/'experiments/e7-d4a1-sources.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];sc=json.loads(gzip.decompress(sa.read_bytes()));sap=ROOT/'experiments/e7-d4a1-initial-source-audit.json';gap=ROOT/'experiments/e7-d4a1-generators-audit.json';assert json.loads(sap.read_text())['full_source_complete_all_degrees'] and json.loads(gap.read_text())['result']=='passed';P=fmpq_mpoly_ctx.get(('h','t','v4','p5','v6','v8'));h,t,E2,E5,E3,E4=P.gens();E1=8*t-3*h*h;zero=P.constant(0);one=P.constant(1)
 def parse(text,ctx=P):
  names=dict(zip(ctx.names(),range(len(ctx.names()))));out={}
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(names)
   for factor in term.lstrip('-').split('*'):
    name,_,power=factor.partition('^')
    if name in names:e[names[name]]+=int(power or 1)
    else:c*=fmpq(factor)
   e=tuple(e);out[e]=out.get(e,fmpq(0))+c
  return ctx.from_dict(out)
 def key(e):return (sum(a*b for a,b in zip(e,WS)),tuple(-a for a in reversed(e)))
 reducers=[]
 for text in bc['groebner_basis']:
  f=parse(text).to_dict();e=max(f,key=key);reducers.append((e,{u:-c/f[e] for u,c in f.items() if u!=e}))
 @lru_cache(None)
 def normal_mon(e):
  if key(e)[0]>42:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,b in normal_mon(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,fmpq(0))+c*b
    return {e:c for e,c in out.items() if c}
  return {e:fmpq(1)}
 def red(f):
  out={}
  for e,c in f.to_dict().items():
   for u,b in normal_mon(e).items():out[u]=out.get(u,fmpq(0))+c*b
  return P.from_dict(out)
 @lru_cache(None)
 def fiber(i,j):
  if j>=2:
   out={}
   for shift,c in ((0,E4),(2,-E3),(4,E2),(6,-E1),(8,one)):
    for k,v in fiber(i+shift,j-2).items():out[k]=out.get(k,zero)+c*v
   return out
  if i and j:return {k:E5*v for k,v in fiber(i-1,0).items()}
  if i>=9:
   out={}
   for ii,jj,c in ((0,1,E5),(1,0,-E4),(3,0,E3),(5,0,-E2),(7,0,E1)):
    for k,v in fiber(i-9+ii,jj).items():out[k]=out.get(k,zero)+c*v
   return out
  return {9 if j else i:one}
 def mul(a,b):
  out={}
  for i,f in a.items():
   for j,g in b.items():
    for k,c in fiber((i if i<9 else 0)+(j if j<9 else 0),int(i==9)+int(j==9)).items():out[k]=out.get(k,zero)+f*g*c
  return {i:red(f) for i,f in out.items() if f}
 def power(f,k):
  ans={0:one}
  for _ in range(k):ans=mul(ans,f)
  return ans
 def add(parts):
  out={}
  for c,ff in parts:
   for i,f in ff.items():out[i]=out.get(i,zero)+c*f
  return {i:red(f) for i,f in out.items() if f}
 O=fmpq_mpoly_ctx.get(('h','w','v4','v5','v6','v8','t','r'))
 def universal(text):
  out={}
  for e,c in parse(text,O).to_dict().items():
   assert e[1]%2==0;ff=P.from_dict({(e[0],e[1]//2,e[2],e[3],e[4],e[5]):c*(-3)**(e[1]//2)})
   for i,g in fiber(e[6],e[7]).items():out[i]=out.get(i,zero)+ff*g
  return {i:red(f) for i,f in out.items() if f}
 ni=[universal(row['tower_polynomial']) for row in tc['source_preimages']];gref=ni[tc['reference_quartic_generator']];assert set(gref[9].to_dict())=={(0,)*6};quartic=gref[9].to_dict()[(0,)*6];refs=[power({0:h/2,1:one/2},i) for i in range(9)]+[gref]
 def coordinates(ff):
  q=ff.get(9,zero)/quartic;out=[zero]*10;out[9]=red(q)
  for i in range(9):
   f=ff.get(i,zero)-q*gref.get(i,zero)
   for j in range(i+1):out[j]+=comb(i,j)*2**j*(-h)**(i-j)*f
  return list(map(red,out))
 modules=[]
 for exp,row in zip(tc['module_exponents'],cert['weight_source_modules']):
  ff={0:one}
  for i,k in enumerate(exp):
   if k:ff=mul(ff,power(ni[i],k))
  cc=coordinates(ff);assert cc==[parse(f) for f in row['components']];modules.append((row['degree'],cc))
 raw=[{0:h},{0:h/2,1:one/2},{0:t},{0:E1,2:-one},{0:E2,2:-E1,4:one},{9:one},{0:E3,2:-E2,4:E1,6:-one}];raw=[{i:f/fmpq(c) for i,f in ff.items()} for ff,c in zip(raw,sc['normalization_factors'])]
 @lru_cache(None)
 def classical(e):
  if not any(e):return {0:one}
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return mul(raw[j],classical(tuple(u)))
 mg=[]
 for row in sc['degrees']:
  for j in row['new_generator_indices']:mg.append((row['degree'],add((fmpq(line[j]),classical(tuple(e))) for e,line in zip(row['classical_exponents'],row['rational_coordinates']) if fmpq(line[j]))))
 assert len(mg)==18
 products=[]
 for j,(k,ff) in enumerate(mg):
  assert ff=={int(i):parse(f) for i,f in cert['root_generators'][j].items()};pp=[coordinates(mul(ff,rr)) for rr in refs];assert pp==[[parse(f) for f in row] for row in cert['root_multiplication_components'][j]];products.append(pp)
 assert refs==[{int(i):parse(f) for i,f in row.items()} for row in cert['reference_generators']];classical.cache_clear();print('independent D4+A1 complete source/tower reconstruction',round(perf_counter()-start,2),flush=True)
 norm={int(d):[tuple(e) for e in ee] for d,ee in bc['normal_exponents'].items()};ix={d:{e:i for i,e in enumerate(ee)} for d,ee in norm.items()};BN={0:fmpq_mat([[1]])}
 for row in bc['degrees']:BN[row['degree']]=fmpq_mat([[fmpq(c) for c in line] for line in row['image_bases']['N']])
 BI={d:H.inv() for d,H in BN.items()}
 @lru_cache(1024)
 def small(text,bd):
  f=parse(text);k=key(next(iter(f.to_dict())))[0];dd=bd+k;A=fmpq_mat(len(norm[dd]),len(norm[bd]))
  for j,e in enumerate(norm[bd]):
   out={}
   for u,c in f.to_dict().items():
    for v,b in normal_mon(tuple(a+b for a,b in zip(e,u))).items():out[v]=out.get(v,fmpq(0))+c*b
   for v,c in out.items():A[ix[dd][v],j]=c
  return BI[dd]*A*BN[bd]
 def pos(d):
  offsets={};n=0
  for i,k in enumerate(FD):
   if 0<=d-k<=42:offsets[i]=n;n+=len(norm[d-k])
  return offsets,n
 def setblock(A,i,j,B):
  for rr in range(B.nrows()):
   for cc in range(B.ncols()):A[i+rr,j+cc]=B[rr,cc]
 def local(A):return all(int(c.denominator)%2 for c in A.entries())
 def valuation(n):n=int(n);return (n&-n).bit_length()-1
 def denominator_exponent(A):return max([0]+[valuation(c.denominator) for c in A.entries()])
 @lru_cache(None)
 def modular_context(precision):return fmpz_mod_ctx(1<<precision)
 def residue(A,scale,precision):
  assert precision>0;mod=1<<precision;inverses={};values=[]
  for c in A.entries():
   den=int(c.denominator);v=valuation(den);assert scale>=v;odd=den>>v
   if odd not in inverses:inverses[odd]=pow(odd,-1,mod)
   values.append((int(c.numerator)*inverses[odd]<<(scale-v))%mod)
  return nmod_mat(A.nrows(),A.ncols(),values,mod) if precision<63 else fmpz_mod_mat(A.nrows(),A.ncols(),values,modular_context(precision))
 def rank2(vectors):
  piv={}
  for v in vectors:
   while v:
    k=v.bit_length()-1
    if k in piv:v^=piv[k]
    else:piv[k]=v;break
  return piv
 def ranks4(columns,n):
  # Independent row operations on the reconstructed mod-four image.
  low=[0]*n;high=[0]*n
  for j,(a,b) in enumerate(columns):
   for rows,bits in ((low,a),(high,b)):
    while bits:
     bit=bits&-bits;rows[bit.bit_length()-1]|=1<<j;bits-=bit
  piv={};even=[]
  for a,b in zip(low,high):
   while a:
    k=(a&-a).bit_length()-1
    if k not in piv:
     if b>>k&1:b^=a
     piv[k]=(a,b);break
    x,y=piv[k]
    b^=y^((a&x) if b>>k&1 else ((~a)&x));a^=x
   else:even.append(b)
  ep=[]
  for b in even:
   for k,(a,_) in sorted(piv.items()):
    if b>>k&1:b^=a
   ep.append(b)
  return len(piv),len(rank2(ep))
 NB={0:fmpq_mat([[1]])};MB={0:fmpq_mat([[1]])};records=[]
 for row in cert['degrees']:
  d=row['degree'];offset,n=pos(d);HN=fmpq_mat([[fmpq(c) for c in line] for line in row['weight_basis']]);HM=fmpq_mat([[fmpq(c) for c in line] for line in row['root_lower_basis']]);HI=HN.inv();assert local(HI);blocks=[]
  for k,cc in modules:
   bd=d-k
   if not 0<=bd<=42:continue
   A=fmpq_mat(n,len(norm[bd]))
   for i,f in enumerate(cc):
    if f and i in offset:setblock(A,offset[i],0,small(str(f),bd))
   blocks.append(A)
  total=sum(A.ncols() for A in blocks);NA=fmpq_mat(n,total);column=0
  for A in blocks:setblock(NA,0,column,A);column+=A.ncols()
  scale=denominator_exponent(NA);C=residue(HI,0,scale+1)*residue(NA,scale,scale+1);mask=(1<<scale)-1;assert all(int(c)&mask==0 for c in C.entries());assert len(rank2([sum(((int(C[i,j])>>scale)&1)<<j for j in range(C.ncols())) for i in range(n)]))==n
  R=HI*HM;assert local(R) and local(2*R.inv());rcols=[sum((int(R[i,j].numerator)&1)<<i for i in range(n)) for j in range(n)];rpiv=rank2(rcols);actual=[]
  for gi,(k,ff) in enumerate(mg):
   if k>d:continue
   oldoff,oldn=pos(d-k);OP=fmpq_mat(n,oldn)
   for ii,off in oldoff.items():
    bd=d-k-FD[ii]
    for oi,f in enumerate(products[gi][ii]):
     if f and oi in offset:setblock(OP,offset[oi],off,small(str(f),bd))
   a=denominator_exponent(OP);b=denominator_exponent(MB[d-k]);p=a+b+2;C=residue(HI,0,p)*residue(OP,a,p)*residue(MB[d-k],b,p);mask=(1<<(a+b))-1;assert all(int(c)&mask==0 for c in C.entries())
   for j in range(C.ncols()):
    lo=hi=0
    for i in range(n):
     v=(int(C[i,j])>>(a+b))&3;lo|=(v&1)<<i;hi|=(v>>1)<<i
    v=lo
    while v:
     pivot=v.bit_length()-1;assert pivot in rpiv;v^=rpiv[pivot]
    actual.append((lo,hi))
  r0,r1=ranks4(actual,n);assert r0+r1==n and r0==len(rpiv);m=n-r0;assert row['smith_multiplicities']==({'1':m} if m else {});NB[d]=HN;MB[d]=HM;records.append(dict(degree=d,rank=n,multiplicity_2=m,full_weight_image_verified=True,full_root_image_verified=True));small.cache_clear();print('independent complete D4+A1 images',d,m,round(perf_counter()-start,2),flush=True)
 assert len(records)==(8 if args.pilot else 50)
 if not args.pilot:assert sum(row['rank'] for row in records)+1==7560
 out=dict(result='passed',complete=not args.pilot,degree_bound=8 if args.pilot else 50,retained_nodes=[2,3,4,5,7],levi_type='D4+A1',dimension=50,total_rank=7560,degrees=records,multiplicity_coefficients=[0]+[row['multiplicity_2'] for row in records],total_nonunit_factors=sum(row['multiplicity_2'] for row in records),all_nonunit_factors=2,seconds=round(perf_counter()-start,3),engine='Native FLINT reconstruction of the exact quadric tower, all source images and 34 full weight-source modules. Exact full weight primitive ranks and independent mod-four root-lattice verification with row operations. Uses only ring addition/multiplication over composite residues. No Sage or primary experiment imports.',scope='Full-source images in every recorded degree; pilot mode stops at degree eight. The separately audited all-degree root source theorem supplies generator completeness. The integral 34-generator weight source module is a shared audited input; integral fiber freeness is not assumed.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),rp,bp,sp,sap,gap)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (archive,ba,tp,sa)});(ROOT/('experiments/e7-d4a1-tower-'+('pilot-' if args.pilot else '')+'images-audit.json')).write_text(json.dumps(out,indent=2)+'\n');print('PASSED D4+A1 image audit',out['seconds'],flush=True)

if __name__=='__main__':main()
