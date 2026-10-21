"""Complete D4 images by an integral finite-module calculation over D5.

No primitive integral fiber basis or freeness of invariant sources is assumed.
The rational quadric tower supplies reference lattices. All 34 integral
module generators are retained before computing their extension indices.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
from math import comb
import argparse,gzip,json,zlib
from sage.all import QQ,ZZ,GF,PolynomialRing,TermOrder,matrix,vector,prod,pari,inverse_mod
from e7_d4_source import context,powers,WEIGHTS
ROOT=Path(__file__).resolve().parents[1]
WS=(1,1,4,5,6,8)
FIBER=(0,1,2,3,4,5,6,7,8,4)
def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def encode(f):return [[list(e),str(c)] for e,c in f.dict().items()]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=51);ap.add_argument('--resume',action='store_true');args=ap.parse_args();start=perf_counter()
 basepath=ROOT/'experiments/e7-d5-images.json';prior=json.loads(basepath.read_text());archive=ROOT/prior['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==prior['certificate_archive_sha256'];base=json.loads(gzip.decompress(archive.read_bytes()))
 P=PolynomialRing(QQ,6,names=('h','w','v4','v5','v6','v8'),order=TermOrder('wdegrevlex',WS));h,w,E2,E5,E3,E4=P.gens();E1=-3*h*h-QQ(8)/3*w*w
 G=[P(f) for f in base['groebner_basis']];normal={int(d):[tuple(e) for e in es] for d,es in base['normal_exponents'].items()};idx={d:{e:i for i,e in enumerate(es)} for d,es in normal.items()}
 reducers=[]
 for f in G:
  e=tuple(f.lm().exponents()[0]);c=f.monomial_coefficient(P({e:1}));reducers.append((e,{tuple(u):-v/c for u,v in f.dict().items() if tuple(u)!=e}))
 @lru_cache(None)
 def norm(e):
  d=sum(a*b for a,b in zip(e,WS))
  if d>43:return {}
  if e in idx[d]:return {e:QQ.one()}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,k in norm(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,QQ.zero())+c*k
    return {e:c for e,c in out.items() if c}
  raise AssertionError('missing reducer')
 def reduce_base(f):
  out={}
  for e,c in P(f).dict().items():
   for u,v in norm(tuple(e)).items():out[u]=out.get(u,QQ.zero())+c*v
  return P(out)
 @lru_cache(None)
 def fiber_power(i,j):
  if j>=2:
   out={}
   for shift,c in ((0,E4),(2,-E3),(4,E2),(6,-E1),(8,P.one())):
    for k,v in fiber_power(i+shift,j-2).items():out[k]=out.get(k,P.zero())+c*v
   return out
  if i and j:return {k:E5*v for k,v in fiber_power(i-1,0).items()}
  if i>=9:
   out={}
   for ii,jj,c in ((0,1,E5),(1,0,-E4),(3,0,E3),(5,0,-E2),(7,0,E1)):
    for k,v in fiber_power(i-9+ii,jj).items():out[k]=out.get(k,P.zero())+c*v
   return out
  return {9 if j else i:P.one()}
 def tmul(a,b):
  out={}
  for i,f in a.items():
   for j,g in b.items():
    ii=(i if i<9 else 0)+(j if j<9 else 0);jj=int(i==9)+int(j==9)
    for k,c in fiber_power(ii,jj).items():out[k]=out.get(k,P.zero())+f*g*c
  return {i:reduce_base(f) for i,f in out.items() if f}
 S=PolynomialRing(QQ,8,names=('h','w','v4','v5','v6','v8','t','r'));hs,ww,e2,e5,e3,e4,t,r=S.gens();embed=P.hom(S.gens()[:6],S);e1=embed(E1)
 raw=[hs,(t+hs)/2,e1-t*t,e2-t*t*e1+t**4,r,e3-t*t*e2+t**4*e1-t**6]
 T,original=context();source=json.loads((ROOT/'experiments/e7-d4-source-algebra.json').read_text());inputs=[];preimages=[]
 lead_values=(1,None,4,6,1,4)
 for i,((f,d),v) in enumerate(zip(original,lead_values)):
  divisor=QQ.one() if i in (0,1) else QQ(v)/f.monomial_coefficient(T.gen(0)**d);raw[i]/=divisor
 for text,d in zip(source['integral_generators'],source['generator_degrees']):
  f=T(text);exps=list(powers((1,1,2,4,4,6),d));fs=[prod(g**k for (g,_),k in zip(original,e)) for e in exps];support=list(powers(WEIGHTS,d));dicts=[{tuple(e):c for e,c in g.dict().items()} for g in fs];fd={tuple(e):c for e,c in f.dict().items()}
  A=matrix(QQ,[[g.get(e,0) for g in dicts] for e in support]);piv=list(A.transpose().pivots());cs=A.matrix_from_rows(piv).inverse()*vector(QQ,[fd.get(support[i],0) for i in piv]);assert A*cs==vector(QQ,[fd.get(e,0) for e in support]);q=S(sum(c*prod(g**k for g,k in zip(raw,e)) for c,e in zip(cs,exps)));out={}
  for exp,c in q.dict().items():
   f=P({tuple(exp[:6]):c})
   for k,g in fiber_power(exp[6],exp[7]).items():out[k]=out.get(k,P.zero())+f*g
  inputs.append({i:reduce_base(f) for i,f in out.items() if f});preimages.append(dict(degree=d,source=text,exponents=[list(e) for e in exps],coefficients=list(map(str,cs)),tower_polynomial=str(q)))
 assert len(inputs)==10
 reference=next(i for i in (3,4) if inputs[i].get(9,0));c=P(inputs[reference][9]);assert c.is_constant() and c
 # Relative quotient monomials give a complete 34-generator B_D5 module.
 U=PolynomialRing(GF(2),10,names=source['elimination_variables'][3:],order=TermOrder('wdegrevlex',source['generator_degrees']));QG=[U(f) for f in source['relative_quotient_groebner']];leads=[tuple(f.lm().exponents()[0]) for f in QG];zero=(0,)*10;found={zero};todo=[zero]
 while todo:
  e=todo.pop()
  for j in range(10):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   if u not in found and not any(all(a<=b for a,b in zip(l,u)) for l in leads):found.add(u);todo.append(u)
 assert len(found)==34
 module_exponents=sorted(found,key=lambda e:(sum(a*b for a,b in zip(e,source['generator_degrees'])),e))
 @lru_cache(None)
 def module_poly(e):
  if not any(e):return {0:P.one()}
  i=next(i for i,k in enumerate(e) if k);u=list(e);u[i]-=1;return tmul(inputs[i],module_poly(tuple(u)))
 modules=[]
 for e in module_exponents:
  old=module_poly(e);rs=old.get(9,P.zero())/c;coeff=[old.get(i,P.zero())-rs*inputs[reference].get(i,P.zero()) for i in range(9)];new=[P.zero() for i in range(10)];new[9]=reduce_base(rs)
  for i,f in enumerate(coeff):
   for j in range(i+1):new[j]+=comb(i,j)*2**j*(-h)**(i-j)*f
  new=list(map(reduce_base,new));degree=sum(a*b for a,b in zip(e,source['generator_degrees']));modules.append((new,degree))
 tower=dict(base_archive=prior['certificate_archive'],base_archive_sha256=prior['certificate_archive_sha256'],source_preimages=preimages,module_exponents=[list(e) for e in module_exponents],reference_quartic_generator=reference,fiber_degrees=list(FIBER),module_coefficients=[dict(degree=d,coefficients=[encode(f) for f in fs]) for fs,d in modules])
 towerpath=ROOT/'experiments/e7-d4-tower-certificates.json.gz';towerpath.write_bytes(gzip.compress(json.dumps(tower,separators=(',',':')).encode(),mtime=0))
 print('tower ready; module generators34; reference',reference,'seconds',round(perf_counter()-start,2),flush=True)
 H={kind:{0:matrix(QQ,1,1,[1])} for kind in ('N','M')}
 for row in base['image_degrees']:
  for kind,key in (('N','weight_basis'),('M','root_basis')):H[kind][row['degree']]=matrix(QQ,row[key])
 HI={kind:{d:A.inverse() for d,A in rows.items()} for kind,rows in H.items()}
 @lru_cache(None)
 def variable_operator(j,degree):
  target=degree+WS[j]
  entries={}
  for column,e in enumerate(normal[degree]):
   u=tuple(a+int(i==j) for i,a in enumerate(e))
   for v,a in norm(u).items():
    if a:entries[idx[target][v],column]=a
  return matrix(QQ,len(normal[target]),len(normal[degree]),entries)
 @lru_cache(None)
 def monomial_operator(e,degree):
  if not any(e):return matrix.identity(QQ,len(normal[degree]))
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return monomial_operator(tuple(u),degree+WS[j])*variable_operator(j,degree)
 @lru_cache(None)
 def polynomial_operator(terms,degree):
  target=degree+sum(a*b for a,b in zip(terms[0][0],WS));out=matrix(QQ,len(normal[target]),len(normal[degree]))
  for e,c in terms:out+=QQ(c)*monomial_operator(e,degree)
  return out
 def multiplication(i,j,degree):return polynomial_operator(tuple((tuple(e),str(c)) for e,c in sorted(modules[i][0][j].dict().items())),degree)
 output=ROOT/'experiments/e7-d4-images.json';certpath=ROOT/'experiments/e7-d4-index-certificates.jsonl.gz';report=dict(result='running',retained_nodes=[2,3,4,5],dimension=51,total_rank=15120,degrees=[])
 saved=json.loads(output.read_text()) if args.resume and output.exists() else None
 prefix=[]
 if saved:
  report['degrees']=saved['degrees'];assert [r['degree'] for r in report['degrees']]==list(range(len(report['degrees'])))
  prefix=[line for line in zlib.decompressobj(31).decompress(certpath.read_bytes()).decode().splitlines() if json.loads(line)['degree']<len(report['degrees'])];assert len(prefix)==2*len(report['degrees'])
  print('resuming after degree',len(report['degrees'])-1,flush=True)
 stream=gzip.open(certpath,'wt',encoding='utf-8')
 for line in prefix:stream.write(line+'\n')
 for degree in range(len(report['degrees']),args.max_degree+1):
  rows=[(j,degree-k) for j,k in enumerate(FIBER) if 0<=degree-k<=43];offset={};n=0
  for j,d in rows:offset[j]=n;n+=len(normal[d])
  record=dict(degree=degree,rank=n,indices={})
  for kind in ('N','M'):
   blocks=[]
   for i,(fs,k) in enumerate(modules):
    prev=degree-k
    if not 0<=prev<=43:continue
    B=matrix(QQ,n,len(normal[prev]))
    for j,tar in rows:
     if not fs[j]:continue
     V=HI[kind][tar]*multiplication(i,j,prev)*H[kind][prev];B.set_block(offset[j],0,V)
    blocks.append(B)
   A=blocks[0]
   for B in blocks[1:]:A=A.augment(B)
   k=int(A.denominator().valuation(2));entry=dict(degree=degree,kind=kind,denominator_2_exponent=k,rows=n,columns=A.ncols())
   if k:
    modulus=ZZ(2)**k;scaled=2**k*A;invs={};values=[]
    for x in scaled.list():
     den=x.denominator()%modulus
     if den not in invs:invs[den]=inverse_mod(den,modulus)
     values.append(ZZ(x.numerator()*invs[den])%modulus)
    C=matrix(ZZ,n,A.ncols(),values);L=matrix(ZZ,pari(C).mathnfmodid(modulus));val=sum(int(L[i,i].valuation(2)) for i in range(n));extension=k*n-val;entry['residues']=[[int(x) for x in row] for row in C];entry['hnf_diagonal']=[int(L[i,i]) for i in range(n)]
   else:extension=0
   entry['extension_index_2_valuation']=extension;stream.write(json.dumps(entry,separators=(',',':'))+'\n');stream.flush();record['indices'][kind]=extension
  predicted=sum(prior['multiplicity_coefficients'][degree-k] for k in FIBER if 0<=degree-k<=43);actual=predicted+record['indices']['N']-record['indices']['M'];record.update(reference_multiplicity=predicted,multiplicity_2=actual,prediction_agrees=actual==predicted);report['degrees'].append(record)
  print('D4',record,'seconds',round(perf_counter()-start,2),flush=True)
  output.write_text(json.dumps(report,indent=2)+'\n')
 stream.close();report.update(result='passed',complete=args.max_degree==51,seconds=round(perf_counter()-start,3),tower_archive=towerpath.relative_to(ROOT).as_posix(),tower_archive_sha256=sha256(towerpath.read_bytes()).hexdigest(),certificate_archive=certpath.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(certpath.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in (Path(__file__).resolve(),ROOT/'experiments/e7_d4_source.py',basepath,ROOT/'experiments/e7-d4-source-algebra.json')});output.write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
