"""Plain Python audit of E7/P6 full-component and bounded image bases.

No Sage or project imports. Uses a sparse rational weighted reducer,
Buchberger checks, Gaussian inverses, and mod-2 bit ranks. The separate
fresh-weight moment audit checks the rational ideal via a different input.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from functools import lru_cache
from time import perf_counter
import gzip,json
from gmpy2 import mpq
ROOT=Path(__file__).resolve().parents[1]
WS=(1,2,4,5,6,8)
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
 start=perf_counter();primary_path=ROOT/'experiments/e7-d5a1-images.json'
 primary=json.loads(primary_path.read_text());assert primary['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 variables=('h','t','v4','p5','v6','v8');ix={v:i for i,v in enumerate(variables)}
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
  # The leading ideal has no standard monomials above degree 42, as
  # verified by exhaustive enumeration of its 756 standard monomials.
  if key(e)[0]>42:return {}
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
 assert len(standard)==756 and max(key(e)[0] for e in standard)==42
 bydegree={d:[tuple(e) for e in cert['normal_exponents'][str(d)]] for d in range(43)}
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
  if degree>42:high_pairs+=1;resolved.add(pair);continue
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
 # Rebuild every source generator image independently of the primary.
 sourcepath=ROOT/'experiments/e7-d5a1-sources-22.json';source=json.loads(sourcepath.read_text());sourcearchive=ROOT/source['certificate_archive'];assert sha256(sourcearchive.read_bytes()).hexdigest()==source['certificate_archive_sha256'];sources=json.loads(gzip.decompress(sourcearchive.read_bytes()))
 oldpath=ROOT/'experiments/e7-d5-images.json';old=json.loads(oldpath.read_text());oldarchive=ROOT/old['certificate_archive'];assert sha256(oldarchive.read_bytes()).hexdigest()==old['certificate_archive_sha256'];oldcert=json.loads(gzip.decompress(oldarchive.read_bytes()))
 normpath=ROOT/'experiments/e7-d5a1-norm.json';norm=json.loads(normpath.read_text());normarchive=ROOT/norm['certificate_archive'];assert sha256(normarchive.read_bytes()).hexdigest()==norm['certificate_archive_sha256'];normcert=json.loads(gzip.decompress(normarchive.read_bytes()))
 def multiply(f,g):
  out={}
  for e,a in f.items():
   for u,b in g.items():
    v=tuple(x+y for x,y in zip(e,u));out[v]=out.get(v,mpq(0))+a*b
  return {e:c for e,c in out.items() if c}
 targets=dict(M=[parse(f) for f in ('h','t','8*t-3*h^2','v4','p5','v6','v8')],N=[parse(f) for f in ('h','8*t-3*h^2','1/2*v4','p5','1/2*v6','v8')])
 factors=list(map(mpq,sources['normalization_factors']));targets['raw']=targets['M'];targets['M']=[{e:c/k for e,c in f.items()} for f,k in zip(targets['M'],factors)]
 @lru_cache(None)
 def rational_monomial(kind,e):
  if not any(e):return {ZERO:mpq(1)}
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return normal(multiply(targets[kind][j],rational_monomial(kind,tuple(u))))
 def combine(mons,coefficients):
  out={}
  for g,c in zip(mons,coefficients):
   c=mpq(c)
   if c:
    for e,a in g.items():out[e]=out.get(e,mpq(0))+c*a
  return {e:c for e,c in out.items() if c}
 generators=dict(N=[(parse('t'),2)],M=[])
 for row in oldcert['preimages']:generators['N'].append((combine([rational_monomial('N',tuple(e)) for e in row['rational_exponents']],row['rational_coefficients']),row['degree']))
 for row in sources['degrees']:
  mons=[rational_monomial('M',tuple(e)) for e in row['classical_exponents']]
  for j in row['new_generator_indices']:generators['M'].append((combine(mons,[line[j] for line in row['rational_coordinates']]),row['degree']))
 # Parse the seven-variable classical norm separately; substitute the
 # already checked quadratic relation by exact rational monomial reduction.
 normterms={};names7={v:i for i,v in enumerate(('h','t','v2','v4','p5','v6','v8'))}
 for term in normcert['norm_classical'].replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=mpq(-1 if term.startswith('-') else 1);e=[0]*7
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in names7:e[names7[name]]+=int(power) if power else 1
   else:c*=mpq(f)
  normterms[tuple(e)]=c
 generators['M'].append((combine([rational_monomial('raw',e) for e in normterms],normterms.values()),32))
 assert all(generators[kind]==[(parse(f),d) for f,d in cert['generators'][kind]] for kind in generators)
 histories={kind:{0:[[mpq(1)]]} for kind in generators};records=[]
 for entry in cert['degrees']:
  d=entry['degree'];n=len(bydegree[d]);bases={};inverses={}
  for kind,gs in generators.items():
   if kind=='M' and d<=source['max_degree']:
    row=sources['degrees'][d-1];mons=[rational_monomial(kind,tuple(e)) for e in row['classical_exponents']];F=[[f.get(e,mpq(0)) for f in mons] for e in bydegree[d]];C=[[mpq(c) for c in line] for line in row['rational_coordinates']];A=matmul(F,C)
   else:
    columns=[]
    for g,k in gs:
     if k<=d:
      for col in zip(*histories[kind][d-k]):
       f=normal(multiply(g,{e:c for e,c in zip(bydegree[d-k],col) if c}));columns.append([f.get(e,mpq(0)) for e in bydegree[d]])
    A=list(map(list,zip(*columns)))
   assert A==[[mpq(c) for c in line] for line in entry['image_columns'][kind]],(kind,d,'source image')
   H=[[mpq(c) for c in line] for line in entry['image_bases'][kind]];HI=inverse(H);coords=matmul(HI,A);assert local(coords) and rank2(coords)==n,(kind,d,'primitive generation');bases[kind]=H;inverses[kind]=HI;histories[kind][d]=H
  if d in (21,22):
   lower=[]
   for g,k in generators['M']:
    if k<=20:
     for col in zip(*histories['M'][d-k]):
      f=normal(multiply(g,{e:c for e,c in zip(bydegree[d-k],col) if c}));lower.append([f.get(e,mpq(0)) for e in bydegree[d]])
   coords=matmul(inverses['M'],list(map(list,zip(*lower))));assert local(coords) and rank2(coords)==n,(d,'new generator does not enlarge Borel image')
  T=matmul(inverses['N'],bases['M']);assert local(T);assert local([[2*c for c in row] for row in inverse(T)]),(d,'twice full weight image')
  m=n-rank2(T);assert ({'1':m} if m else {})==entry['smith_multiplicities'];status='full source comparison' if d<=source['max_degree'] else 'full images equal by sandwich' if not m else 'upper bound using lower root image only';assert status==entry['status']
  records.append(dict(degree=d,rank=n,multiplicity_2=m,full_relative_multiplicity_proved=d<=source['max_degree'] or m==0,twice_full_weight_image_in_root_image=True,new_source_generator_image_already_in_lower_image=d in (21,22),status=status))
  if d%5==0 or d==42:print('independent images',d,round(perf_counter()-start,2),flush=True)
 out=dict(result='passed',engine='Native Python/gmpy2: separate sparse rational reducer, Buchberger criteria, standard-monomial enumeration, all source images reconstructed, rational inverses and mod-two ranks. No Sage or project imports.',dimension=42,total_rank=756,full_root_source_degree_bound=source['max_degree'],all_degree_exponent_divides_2=True,full_multiplicities_proved=primary['complete'],degree_rows=records,buchberger_pairs_checked=pairs,buchberger_chain_pairs=chain_pairs,buchberger_product_pairs=product_pairs,buchberger_above_top_degree_pairs=high_pairs,seconds=round(perf_counter()-start,3),scope='All-degree exponent bound from 2N contained in an explicitly generated lower root image. Multiplicities above the source bound are upper bounds unless the lower root and full weight images coincide. Shared full D5 source and fresh-weight rational ideal are checked separately.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path,sourcepath,oldpath,normpath)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (archive,sourcearchive,oldarchive,normarchive)});(ROOT/'experiments/e7-d5a1-images-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],flush=True)

if __name__=='__main__':main()
