"""Plain Python audit of complete E7/D5 local image bases.

No Sage or project imports. Uses a sparse rational weighted reducer,
Buchberger checks, Gaussian inverses, and mod-2 bit ranks. The separate
adjoint-moment audit checks the rational ideal via a different input.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from functools import lru_cache
from time import perf_counter
import gzip,json
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
WS=(1,1,4,5,6,8)
ZERO=(0,)*6

def rank2(rows):
 pivots={}
 for row0 in rows:
  row=sum((int(c.numerator)%2)<<j for j,c in enumerate(row0))
  while row:
   k=row.bit_length()-1
   if k in pivots:row^=pivots[k]
   else:pivots[k]=row;break
 return len(pivots)

def inverse(A):
 n=len(A);B=[list(row)+[mpq(int(i==j)) for j in range(n)] for i,row in enumerate(A)]
 for k in range(n):
  pivot=next(j for j in range(k,n) if B[j][k]);B[k],B[pivot]=B[pivot],B[k]
  c=B[k][k];B[k]=[v/c for v in B[k]]
  for j in range(n):
   if j!=k and B[j][k]:
    c=B[j][k];B[j]=[v-c*u for v,u in zip(B[j],B[k])]
 return [row[n:] for row in B]

def matmul(A,B):
 return [[sum((a*b for a,b in zip(row,col)),mpq(0)) for col in zip(*B)] for row in A]

def local(A):return all(int(c.denominator)%2 for row in A for c in row)

def main():
 start=perf_counter();primary_path=ROOT/'experiments/e7-d5-images.json'
 primary=json.loads(primary_path.read_text());assert primary['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 variables=('h','w','v4','v5','v6','v8');ix={v:i for i,v in enumerate(variables)}
 def parse(text):
  out={}
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   coefficient=mpq(-1 if term.startswith('-') else 1);term=term.lstrip('-');exponents=[0]*6
   for factor in term.split('*'):
    name,_,power=factor.partition('^')
    if name in ix:exponents[ix[name]]+=int(power) if power else 1
    else:coefficient*=mpq(factor)
   e=tuple(exponents);out[e]=out.get(e,mpq(0))+coefficient
  return {e:c for e,c in out.items() if c}
 def key(e):return (sum(w*j for w,j in zip(WS,e)),tuple(-j for j in reversed(e)))
 G=[parse(f) for f in cert['groebner_basis']];leads=[max(f,key=key) for f in G]
 assert all(len({key(e)[0] for e in f})==1 for f in G)
 reducers=[(e,{u:-c/f[e] for u,c in f.items() if u!=e}) for e,f in zip(leads,G)]
 @lru_cache(None)
 def normal_monomial(e):
  # Every monomial above degree 43 lies in the leading ideal, as verified
  # below by exhaustive enumeration of its 1512 standard monomials.
  if key(e)[0]>43:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     exponent=tuple(a+b for a,b in zip(u,shift))
     for v,d in normal_monomial(exponent).items():out[v]=out.get(v,mpq(0))+c*d
    return {u:c for u,c in out.items() if c}
  return {e:mpq(1)}
 def normal(f):
  out={}
  for e,c in f.items():
   for u,d in normal_monomial(e).items():out[u]=out.get(u,mpq(0))+c*d
  return {u:c for u,c in out.items() if c}
 # Leading-ideal enumeration is independent of the recorded monomial list.
 standard={ZERO};pending=[ZERO]
 while pending:
  e=pending.pop()
  for j in range(6):
   v=list(e);v[j]+=1;v=tuple(v)
   if v not in standard and not any(all(a>=b for a,b in zip(v,u)) for u in leads):standard.add(v);pending.append(v)
 assert len(standard)==1512 and max(key(e)[0] for e in standard)==43
 bydegree={d:[tuple(e) for e in cert['normal_exponents'][str(d)]] for d in range(44)}
 assert set().union(*(set(es) for es in bydegree.values()))==standard
 assert all(key(e)[0]==d for d,es in bydegree.items() for e in es)
 # Buchberger's criterion, skipping relatively prime leading monomials.
 pairs=0;chain_pairs=0;product_pairs=0;high_pairs=0;resolved=set();queue=[]
 for i in range(len(G)):
  for j in range(i):
   a,b=leads[i],leads[j];lcm=tuple(max(x,y) for x,y in zip(a,b));queue.append((key(lcm)[0],i,j,lcm))
 for degree,i,j,lcm in sorted(queue):
  a,b=leads[i],leads[j];pair=(i,j)
  if not any(x and y for x,y in zip(a,b)):product_pairs+=1;resolved.add(pair);continue
  if degree>43:high_pairs+=1;resolved.add(pair);continue
  # Buchberger chain criterion: both smaller-lcm pairs already have
  # certified standard representations, hence so does this S-pair.
  witness=next((k for k,c in enumerate(leads) if k not in (i,j)
     and all(x<=y for x,y in zip(c,lcm))
     and tuple(max(x,y) for x,y in zip(a,c))!=lcm
     and tuple(max(x,y) for x,y in zip(b,c))!=lcm
     and (max(i,k),min(i,k)) in resolved and (max(j,k),min(j,k)) in resolved),None)
  if witness is not None:chain_pairs+=1;resolved.add(pair);continue
  shifts=[tuple(x-y for x,y in zip(lcm,lead)) for lead in (a,b)];poly={}
  for f,lead,shift,sign in ((G[i],a,shifts[0],1),(G[j],b,shifts[1],-1)):
   for e,c in f.items():
    u=tuple(x+y for x,y in zip(e,shift));poly[u]=poly.get(u,mpq(0))+sign*c/f[lead]
  assert not normal(poly),(i,j)
  pairs+=1;resolved.add(pair)
 print('independent Buchberger/standard basis',pairs,'pairs',round(perf_counter()-start,2),flush=True)
 assert all(not normal(parse(f)) for f in cert['relations'])
 # Source completeness and the fresh-weight rational-ideal identity are
 # audited separately. This verifier independently checks all image arithmetic.
 sources={kind:[(parse(f),d) for f,d in gs] for kind,gs in cert['source_generators'].items()}
 assert len(sources['N'])==len(sources['M'])==11
 assert sources['N'][1:]==sources['M'][1:]
 assert sources['N'][0]==(parse('w'),1) and sources['M'][0]==(parse('2*w'),1)
 proposed={kind:{0:[[mpq(1)]]} for kind in sources};inverses={kind:{0:[[mpq(1)]]} for kind in sources};records=[]
 @lru_cache(None)
 def multiplication(kind,index,degree):
  f,k=sources[kind][index];target=bydegree[degree];ix={e:i for i,e in enumerate(target)};columns=[]
  for e in bydegree[degree-k]:
   out=normal({tuple(a+b for a,b in zip(u,e)):c for u,c in f.items()});col=[mpq(0)]*len(target)
   for u,c in out.items():col[ix[u]]=c
   columns.append(col)
  return list(map(list,zip(*columns)))
 for entry in cert['image_degrees']:
  degree=entry['degree'];n=len(bydegree[degree])
  for kind in sources:
   H=[[mpq(c) for c in row] for row in entry['weight_basis' if kind=='N' else 'root_basis']];HI=inverse(H);proposed[kind][degree]=H;inverses[kind][degree]=HI
   products=[]
   for i,(f,k) in enumerate(sources[kind]):
    if k>degree:continue
    values=matmul(multiplication(kind,i,degree),proposed[kind][degree-k]);coords=matmul(HI,values);assert local(coords)
    products.extend(list(col) for col in zip(*coords))
   assert rank2(list(map(list,zip(*products))))==n
  transition=matmul(inverses['N'][degree],proposed['M'][degree]);assert local(transition)
  twice_inverse=[[2*c for c in row] for row in inverse(transition)];assert local(twice_inverse)
  m=n-rank2(transition);assert m==entry['multiplicity_2']==primary['multiplicity_coefficients'][degree]
  records.append(dict(degree=degree,rank=n,multiplicity_2=m,both_image_closures_and_primitive_generation=True,twice_weight_image_in_root_image=True))
  if degree%5==0 or degree==43:print('independent complete images degree',degree,round(perf_counter()-start,2),flush=True)
 output=dict(result='passed',engine='Plain Python/gmpy2 sparse weighted reduction and Buchberger checks; own rational inverses and bit ranks. A separate native parser reads canonical polynomial text. No Sage or project imports.',degree_rows=records,total_rank=1512,dimension=43,buchberger_pairs_checked=pairs,buchberger_chain_pairs=chain_pairs,buchberger_product_pairs=product_pairs,buchberger_above_top_degree_pairs=high_pairs,seconds=round(perf_counter()-start,3),scope='Rational standard basis, every proposed image basis by closure and primitive generation, and every local comparison. Full-source completeness and the fresh-weight ideal identity are checked in separate audits.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-d5-images-audit.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED',output['seconds'],flush=True)
if __name__=='__main__':main()
