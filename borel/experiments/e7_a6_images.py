"""Full 2-local root/weight image comparison on E7/P2 (Levi A6).

A common A6 root complement reduces the full-flag comparison to a
576-dimensional quotient. Both invariant sources are proved polynomial
at 2, so multiplication closure computes their complete images.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import gzip,json
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,prod,binomial,TermOrder
ROOT=Path(__file__).resolve().parents[1]
DEGREES=(1,3,4,5,6,7)

def powers(weights,total):
 if total<0:return
 if not weights:
  if total==0:yield ()
  return
 for e in range(total//weights[0]+1):
  for rest in powers(weights[1:],total-e*weights[0]):yield (e,)+rest

def main():
 start=perf_counter();P=PolynomialRing(QQ,6,names=('c1','c3','c4','c5','c6','c7'),order=TermOrder('wdegrevlex',DEGREES))
 c1,c3,c4,c5,c6,c7=P.gens();es=[P.one(),c1,QQ(4)/9*c1*c1,c3,c4,c5,c6,c7];t=c1/3
 ps=[P(7)]
 for k in range(1,19):
  ps.append(P(sum((-1)**(i-1)*es[i]*ps[k-i] for i in range(1,min(k,8)))+((-1)**(k-1)*k*es[k] if k<=7 else 0)))
 relations=[]
 for d in (6,8,10,12,14,18):
  f=2*ps[d]+sum(binomial(d,k)*(-t)**(d-k)*(sum(binomial(k,j)*ps[j]*ps[k-j] for j in range(k+1))-2**k*ps[k]) for k in range(d+1))
  relations.append(P(f))
 print('A6 moments built',round(perf_counter()-start,2),flush=True)
 ideal=P.ideal(relations);G=ideal.groebner_basis();assert ideal.vector_space_dimension()==576
 print('A6 Groebner basis',len(G),round(perf_counter()-start,2),flush=True)
 normal={d:[e for e in powers(DEGREES,d) if prod(x**k for x,k in zip(P.gens(),e)).reduce(G)==prod(x**k for x,k in zip(P.gens(),e))] for d in range(43)}
 assert sum(map(len,normal.values()))==576
 assert all(not prod(x**k for x,k in zip(P.gens(),e)).reduce(G) for e in powers(DEGREES,43))
 mu=c1/7
 common=[P(sum(binomial(7-j,k-j)*(-mu)**(k-j)*es[j] for j in range(k+1))) for k in range(2,8)]
 assert common[0]==c1*c1/63
 generators={'N':[(es[i],i) for i in range(1,8)],'M':[(2*mu,1)]+[(f,d) for f,d in zip(common,range(2,8))]}
 monomials={d:[P({e:1}) for e in normal[d]] for d in normal}
 indices={d:{e:i for i,e in enumerate(normal[d])} for d in normal}
 def coordinates(f,d):
  f=P(f).reduce(G);assert all(tuple(e) in indices[d] for e in f.dict())
  return vector(QQ,[f.monomial_coefficient(m) for m in monomials[d]])
 bases={kind:{0:matrix(QQ,1,1,[1])} for kind in generators};polys={kind:{0:[P.one()]} for kind in generators};records=[]
 for degree in range(1,43):
  for kind,gs in generators.items():
   columns=[]
   for g,k in gs:
    if k<=degree:columns.extend(coordinates(g*f,degree) for f in polys[kind][degree-k])
   A=matrix(QQ,columns).transpose();den=A.denominator();H=matrix(ZZ,den*A).transpose().hermite_form(include_zero_rows=False).transpose()/den
   assert H.is_square() and H.nrows()==len(normal[degree])
   bases[kind][degree]=H;polys[kind][degree]=[sum(c*m for c,m in zip(col,monomials[degree])) for col in H.columns()]
  transition=bases['N'][degree].inverse()*bases['M'][degree];den=transition.denominator();assert den%2
  smith=matrix(ZZ,den*transition).elementary_divisors();valuations=[int(v.valuation(2)) for v in smith];assert max(valuations,default=0)<=1
  m=valuations.count(1);records.append(dict(degree=degree,rank=len(normal[degree]),multiplicity_2=m,weight_basis=[[str(c) for c in r] for r in bases['N'][degree]],root_basis=[[str(c) for c in r] for r in bases['M'][degree]]))
  print('A6 exact',degree,'rank',len(normal[degree]),'factors 2',m,round(perf_counter()-start,2),flush=True)
 cert=dict(relations=[str(f) for f in relations],groebner_basis=[str(f) for f in G],normal_exponents={str(d):[list(e) for e in exps] for d,exps in normal.items()},image_degrees=records,source_generators={kind:[[str(f),d] for f,d in gs] for kind,gs in generators.items()})
 archive=ROOT/'experiments/e7-a6-images-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0))
 coeff=[0]+[r['multiplicity_2'] for r in records];q=PolynomialRing(ZZ,'q').gen();WG=prod(sum(q**i for i in range(d)) for d in (2,6,8,10,12,14,18));WA=prod(sum(q**i for i in range(d)) for d in range(2,8));target=WG//WA;prediction=(q*target)//(1+q)
 assert all(coeff[d]==prediction[d] for d in range(43))
 out=dict(result='passed',parabolic='E7/P2',retained_nodes=[1,3,4,5,6,7],dimension=42,total_rank=576,coefficient_ring='Z_(2)',multiplicity_coefficients=coeff,prediction_agrees=True,full_flag_polynomial=str(prediction*WA),seconds=round(perf_counter()-start,3),certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),scope='Complete A6 parabolic images from full invariant sources; classical common-fiber transfer gives all-degree full-flag factors. Independent audit is required separately.',source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
 (ROOT/'experiments/e7-a6-images.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],flush=True)
if __name__=='__main__':main()
