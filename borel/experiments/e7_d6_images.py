"""E7/P1 rational quotient and exact full D6 image comparisons."""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from math import factorial
import argparse,gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,SymmetricFunctions,Partitions,prod,binomial,matrix,pari,inverse_mod
ROOT=Path(__file__).resolve().parents[1]
WS=(1,4,6,6,8,10)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=12);args=ap.parse_args();start=perf_counter()
 P=PolynomialRing(QQ,6,names=('h','v4','v6','p6','v8','v10'),order=TermOrder('wdegrevlex',WS));h,v4,v6,p6,v8,v10=P.gens()
 S=PolynomialRing(QQ,7,names=('h','v2','v4','v6','p6','v8','v10'));hh,v2,ee2,ee3,pp,ee4,ee5=S.gens();es=[S.one(),v2,ee2,ee3,ee4,ee5,pp*pp]
 SF=SymmetricFunctions(QQ);MM=SF.monomial();EE=SF.elementary()
 def symmetric(part):return S(sum(c*prod(es[j] if j<=6 else S.zero() for j in mu) for mu,c in EE(MM[part])))
 ps=[S(6)]
 for d in range(1,10):ps.append(S(sum((-1)**(j-1)*es[j]*ps[d-j] for j in range(1,min(d,7)))+((-1)**(d-1)*d*es[d] if d<=6 else 0)))
 moments={}
 for d in (2,6,8,10,12,14,18):
  spin=S.zero()
  for offset in (0,1):
   if d<6*offset:continue
   spin+=sum((QQ(factorial(d))/prod(factorial(2*j+offset) for j in list(part)+[0]*(6-len(part))))*symmetric(part) for part in Partitions((d-6*offset)//2,max_length=6))*(pp if offset else 1)
  f=4*sum(binomial(d,k)*hh**(d-k)*ps[k//2] for k in range(0,d+1,2))/2**d+32*spin/4**d;moments[d]=S(f)
 substitute=S.hom([h,-2*h*h,v4,v6,p6,v8,v10],P);assert substitute(moments[2])==0
 rel=[substitute(moments[d]) for d in (6,8,10,12,14,18)];G=P.ideal(rel).groebner_basis();assert P.ideal(G).vector_space_dimension()==126
 leads=[tuple(g.lm().exponents()[0]) for g in G];found={(0,)*6};todo=list(found)
 while todo:
  e=todo.pop()
  for j in range(6):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   if u not in found and not any(all(a<=b for a,b in zip(l,u)) for l in leads):found.add(u);todo.append(u)
 normal={d:sorted(e for e in found if sum(a*b for a,b in zip(e,WS))==d) for d in range(34)};assert len(found)==126 and max(sum(a*b for a,b in zip(e,WS)) for e in found)==33
 print('D6 quotient rank126, GB',len(G),'seconds',round(perf_counter()-start,2),flush=True)
 srcpath=ROOT/'experiments/e7-d6-sources.json';src=json.loads(srcpath.read_text());assert args.max_degree<=src['max_degree'];archive=ROOT/src['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==src['certificate_archive_sha256'];sources=json.loads(gzip.decompress(archive.read_bytes()))
 target=[h,-2*h*h,v4,v6,p6,v8,v10];records=[];certrows=[];history={'M':{0:matrix(QQ,1,1,[1])}};root_generators=[];witness=None
 ds=sources['M']['original_divisors']
 for row in sources['M']['degrees']:
  C=matrix(QQ,row['rational_coordinates']);mons=[prod((g/k)**a for g,k,a in zip(target,ds,e)) for e in row['classical_exponents']]
  for j in row['new_generator_indices']:root_generators.append((P(sum(mons[i]*C[i,j] for i in range(len(mons)))).reduce(G),row['degree']))
 def local_basis(A):
  k=max([0]+[int(c.denominator().valuation(2)) for c in A.list()]);scaled=2**k*A
  for prec in (64,128,256):
   mod=ZZ(2)**prec;invs={};vals=[]
   for c in scaled.list():
    den=c.denominator()%mod
    if den not in invs:invs[den]=inverse_mod(den,mod)
    vals.append(ZZ(c.numerator()*invs[den])%mod)
   Z=matrix(ZZ,A.nrows(),A.ncols(),vals);H=matrix(QQ,pari(Z).mathnfmodid(mod))/2**k;C=H.inverse()*A
   if C.denominator()%2 and matrix(GF(2),C).rank()==A.nrows():return H
  raise AssertionError('image precision bound')
 for d in range(1,args.max_degree+1):
  images={};bases={};input_columns={}
  for kind in ('N','M'):
   if d<=len(sources[kind]['degrees']):
    entry=sources[kind]['degrees'][d-1];ds=sources[kind]['original_divisors'];C=matrix(QQ,entry['rational_coordinates']);mon=[prod((g/k)**a for g,k,a in zip(target,ds,e)) for e in entry['classical_exponents']]
    fs=[P(sum(mon[i]*C[i,j] for i in range(len(mon)))).reduce(G) for j in range(C.ncols())]
   else:
    assert kind=='M';fs=[]
    for g,k in root_generators:
     if k<=d:
      for col in history['M'][d-k].columns():fs.append((g*P({e:c for e,c in zip(normal[d-k],col) if c})).reduce(G))
   idx={e:i for i,e in enumerate(normal[d])};A=matrix(QQ,len(normal[d]),len(fs),{(idx[tuple(e)],j):c for j,f in enumerate(fs) for e,c in f.dict().items()});H=local_basis(A);images[kind]=A;bases[kind]=H
   if kind=='M':history['M'][d]=H
  T=bases['N'].inverse()*bases['M'];den=T.denominator();assert den%2;vals=[int(c.valuation(2)) for c in matrix(ZZ,den*T).elementary_divisors()];assert all(0<=v<=d for v in vals)
  if d==8:
   V=bases['M'].inverse()*images['N']
   for j,col in enumerate(V.columns()):
    if col.denominator()%2==0:
     assert (2*col).denominator()%2
     witness=dict(degree=d,source_column=j,integral_weight_invariant=sources['N']['degrees'][d-1]['basis_polynomials'][j],borel_normal_form=str(P({e:c for e,c in zip(normal[d],images['N'].column(j)) if c})),root_basis_coordinates=list(map(str,col)),nonzero_relative_class_of_order=2);break
   assert witness is not None
  counts={str(j):vals.count(j) for j in sorted(set(vals)) if j}
  if d>len(sources['M']['degrees']):assert not counts # lower root image = full weight image, hence full root image by inclusion.
  row=dict(degree=d,rank=len(normal[d]),smith_multiplicities=counts,root_image_completeness='complete source component' if d<=len(sources['M']['degrees']) else 'root-generator lower image equals full weight image');records.append(row);certrows.append(dict(**row,source_image_columns={k:[[str(c) for c in row] for row in A] for k,A in images.items()},image_bases={k:[[str(c) for c in row] for row in H] for k,H in bases.items()}));print('D6 images',row,'seconds',round(perf_counter()-start,2),flush=True)
 cert=dict(full_moments={str(d):str(f) for d,f in moments.items()},relations=list(map(str,rel)),groebner_basis=list(map(str,G)),normal_exponents={str(d):[list(e) for e in es] for d,es in normal.items()},root_generators=[[str(f),k] for f,k in root_generators],first_nonroot_weight_image=witness,degrees=certrows)
 path=ROOT/'experiments/e7-d6-images-certificates.json.gz';path.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0))
 report=dict(result='passed',scope='Full image comparison in every recorded degree. N uses complete components; M uses complete components through 27, then a lower-image/full-weight-image equality. No unbounded source-generation assumption.',max_degree=args.max_degree,complete=args.max_degree==33,source_degree_bounds=src.get('source_degree_bounds'),retained_nodes=[2,3,4,5,6,7],dimension=33,total_rank=126,first_nonroot_weight_image=witness,degrees=records,seconds=round(perf_counter()-start,3),certificate_archive=path.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(path.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),srcpath)})
 (ROOT/'experiments/e7-d6-images.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
