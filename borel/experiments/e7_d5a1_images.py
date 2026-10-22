"""Full weight image and bounded full-root comparison for E7/P6.

The full weight source is the previously certified D5 algebra tensored
with R[t]. Higher root-generator images are lower bounds only, unless
they fill the complete weight image. No full-root claim is inferred from
a truncated algebra generating set.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import factorial
from time import perf_counter
import gzip,json
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,SymmetricFunctions,Partitions,prod,binomial,matrix,pari,inverse_mod
ROOT=Path(__file__).resolve().parents[1]
WS=(1,2,4,5,6,8)

def parse(text,P):
 names={name:i for i,name in enumerate(P.variable_names())};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=QQ(-1 if term.startswith('-') else 1);e=[0]*P.ngens()
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in names:e[names[name]]+=int(power) if power else 1
   else:c*=QQ(f)
  out[tuple(e)]=out.get(tuple(e),0)+c
 return P(out)

def main():
 start=perf_counter();P=PolynomialRing(QQ,6,names=('h','t','v4','p5','v6','v8'),order=TermOrder('wdegrevlex',WS));h,t,v4,p5,v6,v8=P.gens()
 S=PolynomialRing(QQ,7,names=('h','t','v2','v4','p5','v6','v8'));hh,tt,E1,E2,pf,E3,E4=S.gens();es=[S.one(),E1,E2,E3,E4,pf*pf];SF=SymmetricFunctions(QQ);MM=SF.monomial();EE=SF.elementary()
 def symmetric(part):return S(sum(c*prod(es[j] if j<=5 else S.zero() for j in mu) for mu,c in EE(MM[part])))
 spin={}
 for k in range(19):
  offset=k%2
  if k<5*offset:spin[k]=S.zero();continue
  spin[k]=S(sum((QQ(factorial(k))/prod(factorial(2*j+offset) for j in list(part)+[0]*(5-len(part))))*symmetric(part) for part in Partitions((k-5*offset)//2,max_length=5)))*(pf if offset else 1)
 ps=[S(5)]
 for k in range(1,10):ps.append(sum((-1)**(j-1)*es[j]*ps[k-j] for j in range(1,min(k,6)))+((-1)**(k-1)*k*es[k] if k<=5 else 0))
 pair=[S(2),hh]
 for k in range(2,19):pair.append(hh*pair[-1]-tt*pair[-2])
 moments={d:S(2*pair[d]+4*sum(binomial(d,k)*(hh*hh/4-tt)**((d-k)//2)*ps[k//2]/2**k for k in range(0,d+1,2))+32*sum(binomial(d,k)*hh**(d-k)*spin[k] for k in range(d+1))/4**d) for d in (2,6,8,10,12,14,18)}
 assert moments[2]==3*(E1+3*hh*hh-8*tt);phi=S.hom([h,t,8*t-3*h*h,v4,p5,v6,v8],P);relations=[phi(moments[d]) for d in (6,8,10,12,14,18)];G=P.ideal(relations).groebner_basis();assert P.ideal(G).vector_space_dimension()==756
 leads=[tuple(g.lm().exponents()[0]) for g in G];todo=[(0,)*6];found=set(todo)
 while todo:
  e=todo.pop()
  for j in range(6):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   if u not in found and not any(all(a<=b for a,b in zip(l,u)) for l in leads):found.add(u);todo.append(u)
 normal={d:sorted(e for e in found if sum(a*b for a,b in zip(e,WS))==d) for d in range(43)};assert len(found)==756 and max(sum(a*b for a,b in zip(e,WS)) for e in found)==42
 print('quotient',len(G),'GB elements; 756 standard monomials; degree42',round(perf_counter()-start,2),flush=True)
 oldpath=ROOT/'experiments/e7-d5-images.json';old=json.loads(oldpath.read_text());oldarchive=ROOT/old['certificate_archive'];assert sha256(oldarchive.read_bytes()).hexdigest()==old['certificate_archive_sha256'];oldcert=json.loads(gzip.decompress(oldarchive.read_bytes()));oldtarget=[h,8*t-3*h*h,v4/2,p5,v6/2,v8]
 generators={'N':[(t,2)],'M':[]}
 for row in oldcert['preimages']:generators['N'].append((P(sum(QQ(c)*prod(g**k for g,k in zip(oldtarget,e)) for c,e in zip(row['rational_coefficients'],row['rational_exponents']))).reduce(G),row['degree']))
 srcpath=ROOT/'experiments/e7-d5a1-sources-22.json';src=json.loads(srcpath.read_text());sourcearchive=ROOT/src['certificate_archive'];assert sha256(sourcearchive.read_bytes()).hexdigest()==src['certificate_archive_sha256'];sources=json.loads(gzip.decompress(sourcearchive.read_bytes()));target=[h,t,8*t-3*h*h,v4,p5,v6,v8];norms=list(map(QQ,sources['normalization_factors']))
 for row in sources['degrees']:
  C=matrix(QQ,row['rational_coordinates']);mons=[prod((g/k)**a for g,k,a in zip(target,norms,e)) for e in row['classical_exponents']]
  for j in row['new_generator_indices']:generators['M'].append((P(sum(mons[i]*C[i,j] for i in range(len(mons)))).reduce(G),row['degree']))
 normpath=ROOT/'experiments/e7-d5a1-norm.json';norm=json.loads(normpath.read_text());normarchive=ROOT/norm['certificate_archive'];assert sha256(normarchive.read_bytes()).hexdigest()==norm['certificate_archive_sha256'];normdata=json.loads(gzip.decompress(normarchive.read_bytes()));generators['M'].append((phi(parse(normdata['norm_classical'],S)).reduce(G),32))
 bases={kind:{0:matrix(QQ,1,1,[1])} for kind in generators};records=[];certrows=[]
 @lru_cache(None)
 def multiplication(kind,i,d):
  g,k=generators[kind][i];ix={e:j for j,e in enumerate(normal[d])};columns=[]
  for e in normal[d-k]:
   f=(g*P({e:1})).reduce(G);columns.append([f.monomial_coefficient(P({u:1})) for u in normal[d]])
  return matrix(QQ,columns).transpose()
 def local_basis(A):
  k=max([0]+[int(c.denominator().valuation(2)) for c in A.list()]);scaled=2**k*A
  for precision in (64,128,256):
   mod=ZZ(2)**precision;Z=matrix(ZZ,A.nrows(),A.ncols(),[ZZ(c.numerator()*inverse_mod(c.denominator(),mod))%mod for c in scaled.list()]);H=matrix(QQ,pari(Z).mathnfmodid(mod))/2**k;C=H.inverse()*A
   if C.denominator()%2 and matrix(GF(2),C).rank()==H.nrows():return H
  raise AssertionError('local image precision')
 for d in range(1,43):
  columns={}
  for kind,gs in generators.items():
   if kind=='M' and d<=src['max_degree']:
    row=sources['degrees'][d-1];C=matrix(QQ,row['rational_coordinates']);mons=[prod((g/k)**a for g,k,a in zip(target,norms,e)).reduce(G) for e in row['classical_exponents']];F=matrix(QQ,len(normal[d]),len(mons),{(normal[d].index(tuple(e)),j):c for j,f in enumerate(mons) for e,c in f.dict().items()});A=F*C
   else:
    blocks=[multiplication(kind,i,d)*bases[kind][d-k] for i,(g,k) in enumerate(gs) if k<=d];A=blocks[0]
    for block in blocks[1:]:A=A.augment(block)
   H=local_basis(A);bases[kind][d]=H;columns[kind]=A
  T=bases['N'][d].inverse()*bases['M'][d];assert T.denominator()%2;vals=[int(c.valuation(2)) for c in matrix(ZZ,T.denominator()*T).elementary_divisors()];counts={str(j):vals.count(j) for j in sorted(set(vals)) if j};status='full source comparison' if d<=src['max_degree'] else 'full images equal by sandwich' if not counts else 'upper bound using lower root image only'
  row=dict(degree=d,rank=len(normal[d]),smith_multiplicities=counts,status=status);records.append(row);certrows.append(dict(**row,image_columns={kind:[[str(c) for c in line] for line in A] for kind,A in columns.items()},image_bases={kind:[[str(c) for c in line] for line in bases[kind][d]] for kind in bases}));print('images',row,round(perf_counter()-start,2),flush=True)
 cert=dict(full_moments={str(d):str(f) for d,f in moments.items()},relations=list(map(str,relations)),groebner_basis=list(map(str,G)),normal_exponents={str(d):[list(e) for e in es] for d,es in normal.items()},generators={kind:[[str(g),k] for g,k in gs] for kind,gs in generators.items()},degrees=certrows);archive=ROOT/'experiments/e7-d5a1-images-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));report=dict(result='passed',dimension=42,total_rank=756,full_root_source_degree_bound=src['max_degree'],complete=all(r['status']!='upper bound using lower root image only' for r in records),degrees=records,seconds=round(perf_counter()-start,3),scope='Full weight image in all degrees using the established D5 source and the A1 norm variable. Full root comparisons only in the stated source range or where the lower root image fills the full weight image; other rows are upper bounds, not full-source Smith factors.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),srcpath,oldpath,normpath,ROOT/'experiments/e7-d5-source-algebra.json',ROOT/'experiments/e7-d5-source-audit.json')},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (sourcearchive,oldarchive,normarchive)});(ROOT/'experiments/e7-d5a1-images.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
