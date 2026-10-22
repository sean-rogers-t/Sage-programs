"""Exact full-component source audit using FLINT, without Sage/project imports.

Set PYTHON_FLINT_PATH if python-flint is installed outside the usual path.
The proof is componentwise: integral invariant columns, full rational rank,
and primitive reduction modulo two imply the entire local invariant lattice.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb,gcd,lcm
from time import perf_counter
import gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mpoly_ctx,fmpz_mat
ROOT=Path(__file__).resolve().parents[1]
WS=(1,1,2,3,4,5,6)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def parse(text,ctx):
 names=ctx.names();indices={str(v):j for j,v in enumerate(names)};terms={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(names)
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in indices:e[indices[name]]+=int(power) if power else 1
   else:c*=fmpq(f)
  e=tuple(e);terms[e]=terms.get(e,fmpq(0))+c
 return ctx.from_dict(terms)

def rank2(columns):
 piv={}
 for col in columns:
  a=sum((int(v)&1)<<j for j,v in enumerate(col))
  while a:
   k=a.bit_length()-1
   if k in piv:a^=piv[k]
   else:piv[k]=a;break
 return len(piv)

def main():
 start=perf_counter();path=ROOT/'experiments/e7-d6-sources.json';source=json.loads(path.read_text());archive=ROOT/source['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==source['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()));records=[]
 ctx=fmpq_mpoly_ctx.get(('z','c1','c2','c3','c4','c5','c6'));z,*cs=ctx.gens();one=ctx.constant(1)
 for kind in ('N','M'):
  bound=source.get('source_degree_bounds',{}).get(kind,source['max_degree']);assert [r['degree'] for r in cert[kind]['degrees']]==list(range(1,bound+1))
  h=cs[0]-2*z if kind=='N' else cs[0]+2*z;shift=-h if kind=='N' else -cs[0]-z;es=[one]+cs;ps=[6*one]
  for k in range(1,11):ps.append(sum(((-1)**(j-1)*es[j]*ps[k-j] for j in range(1,min(k,7))),ctx.constant(0))+((-1)**(k-1)*k*es[k] if k<=6 else ctx.constant(0)))
  xp={k:sum((comb(k,j)*2**j*shift**(k-j)*ps[j] for j in range(k+1)),ctx.constant(0)) for k in (2,4,6,8,10)};sq=[one]
  for k in range(1,6):sq.append(sum(((-1)**(j-1)*sq[k-j]*xp[2*j] for j in range(1,k+1)),ctx.constant(0))/k)
  eu=sum((2**j*shift**(6-j)*es[j] for j in range(7)),ctx.constant(0));raw=[h,sq[1],sq[2],sq[3],eu,sq[4],sq[5]]
  ds=[gcd(*(int(c) for c in f.to_dict().values())) for f in raw];orig=[f/k for f,k in zip(raw,ds)]
  assert ds==cert[kind]['original_divisors'];assert orig==[parse(f,ctx) for f in cert[kind]['original_generators']]
  @lru_cache(None)
  def mon(e):
   if not any(e):return one
   j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1
   return orig[j]*mon(tuple(u))
  for row in cert[kind]['degrees']:
   d=row['degree'];exps=list(powers((1,2,4,6,6,8,10),d));support=list(powers(WS,d));n=len(exps);assert [list(e) for e in exps]==row['classical_exponents']
   F=[mon(e).to_dict() for e in exps];A=[parse(f,ctx).to_dict() for f in row['basis_polynomials']] if 'basis_polynomials' in row else None;assert all(c.denominator==1 for f in F+(A or []) for c in f.values());assert all(set(f)<=set(support) for f in (A or []))
   C=[[fmpq(c) for c in line] for line in row['rational_coordinates']];D=lcm(*(int(c.denominator) for line in C for c in line));assert D>0 and D&(D-1)==0
   B=fmpz_mat([[int(c*D) for c in line] for line in C]);assert B.nrows()==B.ncols()==n
   # Chunk rows to keep peak memory bounded on the higher components.
   bits=[0]*n
   for k in range(0,len(support),128):
    ss=support[k:k+128];left=fmpz_mat([[int(f.get(e,0)) for f in F] for e in ss])*B
    if A is not None:
     right=fmpz_mat([[D*int(f.get(e,0)) for f in A] for e in ss]);assert left==right,(kind,d,k)
    for i in range(len(ss)):
     for j in range(n):
      a=int(left[i,j]);assert a%D==0,(kind,d,k+i,j);bits[j]|=((a//D)&1)<<(k+i)
   piv={}
   for a in bits:
    while a:
     k=a.bit_length()-1
     if k in piv:a^=piv[k]
     else:piv[k]=a;break
   assert len(piv)==n
   records.append(dict(source=kind,degree=d,rank=n,integral_invariant_identity=True,rational_basis=True,primitive_mod2=True))
   print('independent full source',kind,d,n,round(perf_counter()-start,2),flush=True)
 out=dict(result='passed',engine='Native Python/python-flint: independently constructed Newton identities, exact integer matrix products and own GF(2) bit ranks. The full-rank identity F C=A implies the rational basis change is nonsingular. No Sage or project imports.',max_degree=source['max_degree'],source_degree_bounds=source.get('source_degree_bounds'),degrees=records,seconds=round(perf_counter()-start,3),scope='Full localized invariant components in every recorded degree; no unbounded algebra-generation claim.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-d6-sources-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED',out['seconds'],flush=True)

if __name__=='__main__':main()
