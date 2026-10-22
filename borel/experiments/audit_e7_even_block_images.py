"""Native independent reconstruction of the finite-module E7 images.

Use separate rational normal forms and FLINT matrices, then row rather
than column mod-four elimination. Also verify the actual root-coordinate
formulas and every finite-module image before computing any indices.
"""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from time import perf_counter
import argparse,gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mat
from gmpy2 import mpq
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
BW=(2,3,3,4,5,6)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--report',required=True);args=ap.parse_args();start=perf_counter();rp=ROOT/'experiments'/args.report;r=json.loads(rp.read_text());archive=ROOT/r['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==r['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()));kind=cert['blocks'];aux=('s','t') if kind=='42' else ('s2','t2','s','t');aw=(1,2) if kind=='42' else (1,2,1,2);names=aux+('q2','q3','e3','e4','e5','e6');naux=len(aux);ws=aw+BW;top=r['top_polynomial_degree'];nvars=len(names)
 bp=ROOT/'experiments/e7-a5a2-images.json';br=json.loads(bp.read_text());ba=ROOT/br['certificate_archive'];bc=json.loads(gzip.decompress(ba.read_bytes()));sp=ROOT/'experiments/a5-even-block-sources.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];source=next(c for c in json.loads(gzip.decompress(sa.read_bytes()))['cases'] if c['blocks']==kind);assert json.loads((ROOT/'experiments/a5-even-block-sources-audit.json').read_text())['full_sources_complete_all_degrees']
 def parse(text,nn):
  ix=dict(zip(nn,range(len(nn))));out={}
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   c=mpq(-1 if term.startswith('-') else 1);e=[0]*len(nn)
   for f in term.lstrip('-').split('*'):
    name,_,power=f.partition('^')
    if name in ix:e[ix[name]]+=int(power) if power else 1
    else:c*=mpq(f)
   e=tuple(e);out[e]=out.get(e,mpq(0))+c
  return {e:c for e,c in out.items() if c}
 def key(e):
  out=[]
  for k in range(naux//2):
   part=e[2*k:2*k+2];out.extend((part[0]+2*part[1],tuple(-a for a in reversed(part))))
  out.extend((sum(a*b for a,b in zip(e[naux:],BW)),tuple(-a for a in reversed(e[naux:]))));return tuple(out)
 GG=[parse(f,names) for f in cert['fiber_groebner_basis']]+[{(0,)*naux+e:c for e,c in parse(f,names[naux:]).items()} for f in bc['groebner_basis']];reducers=[]
 for f in GG:
  e=max(f,key=key);reducers.append((e,{u:-c/f[e] for u,c in f.items() if u!=e}))
 @lru_cache(None)
 def normal_mon(e):
  if sum(a*b for a,b in zip(e[naux:],BW))>45:return {}
  for lead,tail in reducers:
   if all(a>=b for a,b in zip(e,lead)):
    shift=tuple(a-b for a,b in zip(e,lead));out={}
    for u,c in tail.items():
     for v,b in normal_mon(tuple(a+b for a,b in zip(u,shift))).items():out[v]=out.get(v,mpq(0))+c*b
    return {e:c for e,c in out.items() if c}
  return {e:mpq(1)}
 def normal(f):
  out={}
  for e,c in f.items():
   for u,b in normal_mon(e).items():out[u]=out.get(u,mpq(0))+c*b
  return {e:c for e,c in out.items() if c}
 def multiply(f,g):
  out={}
  for e,c in f.items():
   for u,b in g.items():
    v=tuple(a+b for a,b in zip(e,u));out[v]=out.get(v,mpq(0))+c*b
  return normal(out)
 # Independently check the nine root expressions BEFORE Borel, in
 # universal block elementary coefficients and a quadratic w algebra.
 syms=S.symbols(' '.join(names));env=dict(zip(names,syms));z,w=S.symbols('z w')
 if kind=='42':
  s,t,a2,a3,a4=S.symbols('ss tt a2 a3 a4');poly=S.Poly((1+s*z+t*z*z)*(1-s*z+a2*z*z+a3*z**3+a4*z**4),z);ev=[poly.nth(j) for j in range(7)];a=s-2*w;b=a2+3*s*w+6*w*w;c3=a3-2*w*a2-3*s*w*w-4*w**3;c4=a4-w*a3+w*w*a2+s*w**3+w**4;u=-9*s;v=c3-a*b+a**3;c=c4-a*a*b/3;subs={env['s']:s,env['t']:t};hh=s;tt=t
 else:
  s1,s2,t1,t2,t3=S.symbols('ss1 ss2 tt1 tt2 tt3');s3=-s1-s2;poly=S.Poly((1+s1*z+t1*z*z)*(1+s2*z+t2*z*z)*(1+s3*z+t3*z*z),z);ev=[poly.nth(j) for j in range(7)];u=2*s1+s2;v=s1+2*s2;a=s3-2*w;b=t1-t3-u*w;c=t2-t3-v*w;subs={env['s2']:s2,env['t2']:t2,env['s']:s3,env['t']:t3};hh=s3;tt=t3
 subs.update({env['q2']:-ev[2],env['q3']:0,**{env['e'+str(j)]:ev[j] for j in range(3,7)}});actual=[u,v,a*a,b*b-u*a*b,c*c-v*a*c,2*b-u*a,2*c-v*a,v*b-u*c,2*b*c-a*(v*b+u*c)]
 for f,text in zip(actual,cert['partial_root_generators']):
  ff=S.sympify(text.replace('^','**'),locals=env);assert not ff.has(env['q3']);assert S.expand(S.rem(S.expand(f),w*w-hh*w+tt,w)-ff.subs(subs,simultaneous=True))==0
 # Product/chain checks for the small monic fiber lists, in the full
 # coefficient polynomial ring. Their standard monomials are precisely
 # the classical Grassmannian Schur-module staircase.
 nf=len(cert['fiber_groebner_basis']);leads=[p[0] for p in reducers]
 for i in range(nf):
  assert not any(leads[i][naux:]) and GG[i][leads[i]]==1
  for j in range(i):
   if not any(a and b for a,b in zip(leads[i],leads[j])):continue
   lcm=tuple(max(a,b) for a,b in zip(leads[i],leads[j]));f={}
   for k,sign in ((i,1),(j,-1)):
    shift=tuple(a-b for a,b in zip(lcm,leads[k]))
    for e,c in GG[k].items():
     u=tuple(a+b for a,b in zip(e,shift));f[u]=f.get(u,mpq(0))+sign*c
   assert not normal(f)
 for f in cert['fiber_relations']:assert not normal(parse(f,names))
 fiber=[tuple(e) for e in cert['fiber_exponents']];assert len(fiber)==source['rational_module_rank'] and all(not any(all(a<=b for a,b in zip(l[:naux],e)) for l in leads[:nf]) for e in fiber);fd=cert['fiber_degrees'];inv=[parse(f,names) for f in cert['partial_root_generators']]
 @lru_cache(None)
 def module(e):
  if not any(e):return {(0,)*nvars:mpq(1)}
  j=next(i for i,k in enumerate(e) if k);u=list(e);u[j]-=1;return multiply(inv[j],module(tuple(u)))
 for e,row in zip(source['module_generator_exponents'],cert['module_components']):
  f={}
  for ii,text in row['components'].items():
   for u,c in parse(text,names[naux:]).items():f[fiber[int(ii)]+u]=c
  assert module(tuple(e))==f
 print('independent source/fiber/module reconstruction',kind,round(perf_counter()-start,2),flush=True)
 norm={int(d):[tuple(e) for e in ee] for d,ee in bc['normal_exponents'].items()};ix={d:{e:i for i,e in enumerate(ee)} for d,ee in norm.items()};NB={0:fmpq_mat([[1]])};MB={0:fmpq_mat([[1]])}
 for row in bc['degrees']:
  d=row['degree'];n=len(norm[d]);NB[d]=fmpq_mat(n,n,[fmpq(c) for line in row['image_bases']['N'] for c in line]);MB[d]=fmpq_mat(n,n,[fmpq(c) for line in row['image_bases']['M'] for c in line])
 NI={d:A.inv() for d,A in NB.items()};components=[(row['degree'],{int(i):parse(text,names[naux:]) for i,text in row['components'].items()}) for row in cert['module_components']]
 @lru_cache(None)
 def block(items,bd):
  f=dict(items);k=sum(a*b for a,b in zip(next(iter(f)),BW));outd=bd+k;A=fmpq_mat(len(norm[outd]),len(norm[bd]))
  for j,e in enumerate(norm[bd]):
   out={}
   for u,c in f.items():
    for v,b in normal_mon((0,)*naux+tuple(a+b for a,b in zip(u,e))).items():
     assert not any(v[:naux]);vv=v[naux:];out[vv]=out.get(vv,mpq(0))+c*b
   for v,c in out.items():A[ix[outd][v],j]=fmpq(str(c))
  T=NI[outd]*A*MB[bd];assert all(int(c.denominator)%2 for c in T.entries());return [[int(T[i,j].numerator)*pow(int(T[i,j].denominator),-1,4)%4 for i in range(T.nrows())] for j in range(T.ncols())]
 def row_ranks(columns,n):
  # Transpose the reconstructed Z/4 columns into bit-plane rows.
  low=[0]*n;high=[0]*n
  for j,(lo,hi) in enumerate(columns):
   for values,bits in ((low,lo),(high,hi)):
    while bits:
     bit=bits&-bits;i=bit.bit_length()-1;values[i]|=1<<j;bits-=bit
  rows=list(zip(low,high));piv={};even=[]
  for lo,hi in rows:
   while lo:
    bit=lo&-lo;k=bit.bit_length()-1
    if k not in piv:
     if hi&bit:hi^=lo
     piv[k]=(lo,hi);break
    a,b=piv[k];borrow=(~lo)&a
    if hi&bit:hi=hi^b^(lo&a)
    else:hi=hi^b^borrow
    lo^=a
   else:even.append(hi)
  ep={}
  for hi in even:
   for k,(lo,_) in sorted(piv.items()):
    if hi>>k&1:hi^=lo
   while hi:
    k=(hi&-hi).bit_length()-1
    if k in ep:hi^=ep[k]
    else:ep[k]=hi;break
  return len(piv),len(ep)
 records=[];hashes=[]
 for d in range(top+1):
  row=cert['degrees'][d] if cert['degrees'] else None
  offsets={};n=0
  for fi,k in enumerate(fd):
   if 0<=d-k<=45:offsets[fi]=n;n+=len(norm[d-k])
  columns=[]
  for md,comp in components:
   bd=d-md
   if not 0<=bd<=45:continue
   cc=[[0,0] for _ in norm[bd]]
   for fi,f in comp.items():
    if fi not in offsets or not norm[d-fd[fi]]:continue
    for j,values in enumerate(block(tuple(sorted(f.items())),bd)):
     for i,c in enumerate(values):cc[j][0]|=(c&1)<<(i+offsets[fi]);cc[j][1]|=(c>>1)<<(i+offsets[fi])
   columns.extend(cc)
  digest=sha256(json.dumps(columns,separators=(',',':')).encode()).hexdigest();hashes.append(digest);r0,r1=row_ranks(columns,n);assert r0+r1==n
  if row is not None:assert digest==cert['mod_four_column_sha256'][d] and (r0,r1)==(row['mod_four_unit_rank'],row['mod_four_even_rank']) and row['multiplicity_2']==n-r0
  records.append(dict(degree=d,rank=n,multiplicity_2=n-r0,mod_four_unit_rank=r0,mod_four_even_rank=r1,full_source_image_verified=True));print('independent even-block images',kind,d,round(perf_counter()-start,2),flush=True)
 out=dict(result='passed',complete=True,blocks=kind,total_rank=sum(row['rank'] for row in records),top_polynomial_degree=top,degrees=records,mod_four_column_sha256=hashes,primary_comparison_checked=bool(cert['degrees']),seconds=round(perf_counter()-start,3),engine='Native SymPy pre-Borel root-coordinate identities; separate sparse gmpy2 normal forms, complete module-image reconstruction, FLINT rational multiplication maps and row-based mod-four elimination. No Sage or primary experiment imports.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),rp,bp,sp)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (archive,ba,sa)});(ROOT/('experiments/e7-even-block-'+kind+'-images-audit.json')).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
