"""Independent D4 tower and image-index audit, without Sage/project imports.

Requires gmpy2 and python-flint. Set PYTHON_FLINT_PATH when the latter is
installed in an isolated directory. Previously audited D5 image bases are
shared inputs; all new D4 substitutions and extension matrices are rebuilt.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb
from time import perf_counter,sleep
import argparse,gzip,json,os,sys,zlib
from gmpy2 import mpq
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq_mat,fmpq
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
WS=(1,1,4,5,6,8);ZERO=(0,)*6;FIBER=(0,1,2,3,4,5,6,7,8,4)
sys.setrecursionlimit(20000)

def add(f,g,c=1):
 out=dict(f)
 for e,v in g.items():out[e]=out.get(e,mpq(0))+c*v
 return {e:v for e,v in out.items() if v}
def mul(f,g):
 out={}
 for e,a in f.items():
  for u,b in g.items():v=tuple(x+y for x,y in zip(e,u));out[v]=out.get(v,mpq(0))+a*b
 return {e:v for e,v in out.items() if v}
def scale(f,c):return {e:v*c for e,v in f.items() if v*c}
def power(f,k,n=6):
 out={(0,)*n:mpq(1)}
 while k:
  if k%2:out=mul(out,f)
  f=mul(f,f);k//=2
 return out
def parse(text,variables):
 ix={v:i for i,v in enumerate(variables)};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=mpq(-1 if term.startswith('-') else 1);term=term.lstrip('-');e=[0]*len(variables)
  for factor in term.split('*'):
   name,_,k=factor.partition('^')
   if name in ix:e[ix[name]]+=int(k) if k else 1
   else:c*=mpq(factor)
  e=tuple(e);out[e]=out.get(e,mpq(0))+c
 return {e:c for e,c in out.items() if c}
def valuation(n):
 n=abs(int(n));assert n;return (n&-n).bit_length()-1
def decode(f):return {tuple(e):mpq(c) for e,c in f}

def extension_index(residues,k):
 """Smith elimination over Z/(2^k), with exact bounded integer entries.

 A unit pivot contributes k to the image length, then is removed by
 invertible row/column operations. If no unit remains, divide every
 entry by two and decrement k. This computes the image cardinality,
 not a rank over a field incorrectly applied to the composite ring.
 """
 assert 0<k<=30 # Products of reduced residues fit signed int64 exactly.
 A=np.array(residues,dtype=np.int64);extension=0;pivots={}
 while k and A.size:
  positions=np.argwhere(A%2)
  if not len(positions):A=A//2;k-=1;continue
  r,c=positions[0];A[[0,r]]=A[[r,0]];A[:,[0,c]]=A[:,[c,0]];mod=2**k
  row=(A[0,1:]*pow(int(A[0,0]),-1,mod))%mod
  A=(A[1:,1:]-A[1:,0,None]*row[None,:])%mod
  extension+=k;pivots[k]=pivots.get(k,0)+1
 return extension,pivots

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--max-degree',type=int,default=51);ap.add_argument('--follow',action='store_true');args=ap.parse_args();assert 0<=args.max_degree<=51;full=args.max_degree==51
 start=perf_counter();primary_path=ROOT/'experiments/e7-d4-images.json';primary=json.loads(primary_path.read_text());assert not full or args.follow or (primary['result']=='passed' and primary['complete'])
 towerpath=ROOT/'experiments/e7-d4-tower-certificates.json.gz';towerhash=sha256(towerpath.read_bytes()).hexdigest();assert not full or args.follow or towerhash==primary['tower_archive_sha256'];tower=json.loads(gzip.decompress(towerpath.read_bytes()))
 basepath=ROOT/tower['base_archive'];assert sha256(basepath.read_bytes()).hexdigest()==tower['base_archive_sha256'];base=json.loads(gzip.decompress(basepath.read_bytes()))
 auditpath=ROOT/'experiments/e7-d5-images-audit.json';audit=json.loads(auditpath.read_text());assert audit['result']=='passed' and audit['raw_sha256'][tower['base_archive']]==tower['base_archive_sha256']
 sourcepath=ROOT/'experiments/e7-d4-source-audit.json';sourceaudit=json.loads(sourcepath.read_text());assert sourceaudit['full_source_complete_all_degrees'] and sourceaudit['raw_sha256'][towerpath.relative_to(ROOT).as_posix()]==towerhash
 source=json.loads((ROOT/'experiments/e7-d4-source-algebra.json').read_text())
 names=('h','w','v4','v5','v6','v8');G=[parse(f,names) for f in base['groebner_basis']]
 key=lambda e:(sum(a*b for a,b in zip(e,WS)),tuple(-x for x in reversed(e)))
 reducers=[]
 for f in G:
  l=max(f,key=key);reducers.append((l,{e:-c/f[l] for e,c in f.items() if e!=l}))
 normal={int(d):list(map(tuple,es)) for d,es in base['normal_exponents'].items()};ix={d:{e:i for i,e in enumerate(es)} for d,es in normal.items()}
 @lru_cache(None)
 def nm(e):
  d=key(e)[0]
  if d>43:return {}
  if e in ix[d]:return {e:mpq(1)}
  for l,tail in reducers:
   if all(a>=b for a,b in zip(e,l)):
    shift=tuple(a-b for a,b in zip(e,l));out={}
    for u,c in tail.items():
     for v,k in nm(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,mpq(0))+c*k
    return {e:c for e,c in out.items() if c}
  raise AssertionError('missing reducer')
 def reduce_base(f):
  out={}
  for e,c in f.items():
   for u,k in nm(e).items():out[u]=out.get(u,mpq(0))+c*k
  return {e:c for e,c in out.items() if c}
 one={ZERO:mpq(1)};h=parse('h',names);E2=parse('v4',names);E5=parse('v5',names);E3=parse('v6',names);E4=parse('v8',names);E1=parse('-3*h^2-8/3*w^2',names)
 @lru_cache(None)
 def fp(i,j):
  if j>=2:
   out={}
   for shift,c in ((0,E4),(2,scale(E3,-1)),(4,E2),(6,scale(E1,-1)),(8,one)):
    for k,v in fp(i+shift,j-2).items():out[k]=add(out.get(k,{}),mul(c,v))
   return out
  if i and j:return {k:mul(E5,v) for k,v in fp(i-1,0).items()}
  if i>=9:
   out={}
   for ii,jj,c in ((0,1,E5),(1,0,scale(E4,-1)),(3,0,E3),(5,0,scale(E2,-1)),(7,0,E1)):
    for k,v in fp(i-9+ii,jj).items():out[k]=add(out.get(k,{}),mul(c,v))
   return out
  return {9 if j else i:one}
 def tmul(a,b):
  out={}
  for i,f in a.items():
   for j,g in b.items():
    for k,c in fp((i if i<9 else 0)+(j if j<9 else 0),int(i==9)+int(j==9)).items():out[k]=add(out.get(k,{}),mul(mul(f,g),c))
  return {i:q for i,f in out.items() if (q:=reduce_base(f))}
 # Verify the rational quadric tower by three S-polynomials over a fresh
 # polynomial coefficient ring (no base quotient relations are needed).
 import sympy as s
 tt,rr,e1,e2,e3,e4,e5=s.symbols('t r E1 E2 E3 E4 E5')
 rel=[tt*rr-e5,rr**2-e4+e3*tt**2-e2*tt**4+e1*tt**6-tt**8,tt**9-e5*rr+e4*tt-e3*tt**3+e2*tt**5-e1*tt**7]
 # Weighted grevlex, weights (4,1) on (r,t), gives the monic leads
 # tr,r^2,t^9. The third pair has coprime leads; these are the other
 # two S-pair reductions over the polynomial coefficient ring itself.
 f1,f2,f3=rel;assert s.expand(rr*f1-tt*f2-f3)==0
 sp=s.Poly(s.expand(tt**8*f1-rr*f3-e5*f2),rr,tt);res=0
 for (a,b),coefficient in sp.terms():
  assert a<=1 and (a==0 or b>=1)
  res+=coefficient*(e5*tt**(b-1) if a else tt**b)
 assert s.expand(res)==0
 # Independently reconstruct all ten source substitutions before reducing.
 n8=names+('t','r');rawtext=['h','1/2*t+1/2*h','-3*h^2-8/3*w^2-t^2','v4+3*h^2*t^2+8/3*w^2*t^2+t^4','r','v6-t^2*v4-3*h^2*t^4-8/3*w^2*t^4-t^6']
 raw=[scale(parse(f,n8),mpq(1,k)) for f,k in zip(rawtext,sourceaudit['rational_generator_divisors'])]
 @lru_cache(None)
 def rawmon(e):
  if not any(e):return {(0,)*8:mpq(1)}
  i=next(i for i,k in enumerate(e) if k);u=list(e);u[i]-=1;return mul(raw[i],rawmon(tuple(u)))
 inputs=[]
 for record in tower['source_preimages']:
  f={}
  for e,c in zip(record['exponents'],record['coefficients']):f=add(f,rawmon(tuple(e)),mpq(c))
  assert f==parse(record['tower_polynomial'],n8)
  out={}
  for e,c in f.items():
   for k,g in fp(e[6],e[7]).items():out[k]=add(out.get(k,{}),mul({e[:6]:c},g))
  inputs.append({i:q for i,f in out.items() if (q:=reduce_base(f))})
 ref=tower['reference_quartic_generator'];assert set(inputs[ref][9])=={ZERO};c=inputs[ref][9][ZERO]
 @lru_cache(None)
 def module(e):
  if not any(e):return {0:one}
  i=next(i for i,k in enumerate(e) if k);u=list(e);u[i]-=1;return tmul(inputs[i],module(tuple(u)))
 modules=[]
 for e,entry in zip(tower['module_exponents'],tower['module_coefficients']):
  e=tuple(e);old=module(e);rs=scale(old.get(9,{}),1/c);coeff=[add(old.get(i,{}),mul(rs,inputs[ref].get(i,{})),-1) for i in range(9)];new=[{} for i in range(10)];new[9]=reduce_base(rs)
  for i,f in enumerate(coeff):
   for j in range(i+1):new[j]=add(new[j],mul(power(h,i-j),f),comb(i,j)*2**j*(-1)**(i-j))
  new=list(map(reduce_base,new));assert new==[decode(f) for f in entry['coefficients']]
  assert entry['degree']==sum(a*b for a,b in zip(e,source['generator_degrees']));modules.append((new,entry['degree']))
 print('independent ten substitutions and 34 module coefficients',round(perf_counter()-start,2),flush=True)
 H={k:{0:fmpq_mat([[1]])} for k in ('N','M')}
 for row in base['image_degrees']:
  for k,name in (('N','weight_basis'),('M','root_basis')):H[k][row['degree']]=fmpq_mat(row[name])
 HI={k:{d:A.inv() for d,A in hs.items()} for k,hs in H.items()}
 @lru_cache(None)
 def variable_operator(j,d):
  tar=d+WS[j];A=[['0']*len(normal[d]) for _ in normal[tar]]
  for col,e in enumerate(normal[d]):
   u=tuple(a+int(k==j) for k,a in enumerate(e))
   for v,c in nm(u).items():A[ix[tar][v]][col]=str(c)
  return fmpq_mat(A)
 @lru_cache(4096)
 def monomial_operator(e,d):
  if not any(e):return fmpq_mat([[int(i==j) for j in normal[d]] for i in normal[d]])
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
  return monomial_operator(tuple(u),d+WS[j])*variable_operator(j,d)
 @lru_cache(512)
 def polynomial_operator(terms,d):
  degree=key(terms[0][0])[0];A=fmpq_mat(len(normal[d+degree]),len(normal[d]))
  for e,c in terms:A+=fmpq(c)*monomial_operator(e,d)
  return A
 def multiplication(i,j,d):
  # Compose the six variable operators, rather than copying the primary
  # code's direct normal reduction of each large polynomial product.
  return polynomial_operator(tuple((e,str(c)) for e,c in sorted(modules[i][0][j].items())),d)
 certpath=ROOT/'experiments/e7-d4-index-certificates.jsonl.gz';assert not full or args.follow or sha256(certpath.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];stream=iter(zlib.decompressobj(31).decompress(certpath.read_bytes()).decode().splitlines());records=[];checked_lines=[]
 for degree in range(args.max_degree+1):
  if args.follow:
   while True:
    try:
     latest=json.loads(primary_path.read_text());lines=zlib.decompressobj(31).decompress(certpath.read_bytes()).decode().splitlines()
     if len(latest['degrees'])>degree and len(lines)>=2*degree+2:break
    except (json.JSONDecodeError,zlib.error):pass
    sleep(2)
   primary=latest;stream=iter(lines[2*degree:2*degree+2])
  rows=[(j,degree-k) for j,k in enumerate(FIBER) if 0<=degree-k<=43];n=sum(len(normal[d]) for j,d in rows);out=dict(degree=degree,rank=n,indices={})
  for kind in ('N','M'):
   line=next(stream);checked_lines.append(line);cert=json.loads(line);assert (cert['degree'],cert['kind'],cert['rows'])==(degree,kind,n);k=cert['denominator_2_exponent'];mod=2**k;actual=[[] for _ in range(n)];found_k=0
   for i,(fs,birth) in enumerate(modules):
    prev=degree-birth
    if not 0<=prev<=43:continue
    width=len(normal[prev]);offset=0
    for j,tar in rows:
     if fs[j]:
      A=HI[kind][tar]*multiplication(i,j,prev)*H[kind][prev];num,den=A.numer_denom();a=valuation(den);found_k=max(found_k,a);assert a<=k;inv=pow(int(den)//(2**a),-1,mod) if k else 0;factor=2**(k-a)*inv
      for r0 in range(len(normal[tar])):actual[offset+r0].extend((int(num[r0,col])*factor)%mod for col in range(width))
     else:
      for r0 in range(len(normal[tar])):actual[offset+r0].extend([0]*width)
     offset+=len(normal[tar])
   assert found_k==k and len(actual)==n and len(actual[0])==cert['columns']
   if k:
    assert actual==cert['residues'],(degree,kind,'independent matrix reconstruction')
    extension,pivots=extension_index(actual,k)
   else:extension=0
   assert extension==cert['extension_index_2_valuation']==primary['degrees'][degree]['indices'][kind];out['indices'][kind]=extension
  assert out['indices']['N']==out['indices']['M'];records.append(out)
  if degree%4==0 or degree==51:print('independent full D4 indices',degree,round(perf_counter()-start,2),flush=True)
 if not full:print('PASSED independent prefix through',args.max_degree,'seconds',round(perf_counter()-start,3),flush=True);return
 if args.follow:
  while True:
   try:
    primary=json.loads(primary_path.read_text())
    if primary['result']=='passed' and primary['complete']:break
   except json.JSONDecodeError:pass
   sleep(2)
  assert sha256(certpath.read_bytes()).hexdigest()==primary['certificate_archive_sha256'] and towerhash==primary['tower_archive_sha256']
  assert gzip.decompress(certpath.read_bytes()).decode().splitlines()==checked_lines
 assert next(stream,None) is None;assert sum(r['rank'] for r in records)==15120
 report=dict(result='passed',engine='No Sage/project imports. Own sparse rational torus substitutions and normal reductions; FLINT rational block matrices. Own unit-pivot Smith elimination over Z/(2^k), using bounded exact NumPy int64 arithmetic, independently replays PARI modular-column HNF. SymPy checks both nontrivial rational quadric S-pairs.',shared_inputs='Previously independently audited D5 rational quotient and complete image bases; ten new D4 preimages audited separately by the native source audit.',degree_rows=records,total_rank=15120,dimension=51,all_extension_indices_agree=True,seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path,sourcepath,auditpath)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (basepath,towerpath,certpath)})
 (ROOT/'experiments/e7-d4-images-audit.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED independent complete D4 image indices',report['seconds'],flush=True)

if __name__=='__main__':main()
