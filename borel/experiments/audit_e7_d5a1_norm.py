"""Fresh root-character orbit and independent FLINT orbit-norm audit."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations,product
from math import comb,prod
from random import Random
from time import perf_counter
import gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mpoly_ctx
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def parse(text,ctx):
 ix={str(v):j for j,v in enumerate(ctx.names())};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(ix)
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in ix:e[ix[name]]+=int(power) if power else 1
   else:c*=fmpq(f)
  e=tuple(e);out[e]=out.get(e,fmpq(0))+c
 return ctx.from_dict(out)

def powers(ws,d):
 if not ws:
  if not d:yield ()
  return
 for k in range(d//ws[0]+1):
  for e in powers(ws[1:],d-k*ws[0]):yield (k,)+e

def main():
 start=perf_counter();path=ROOT/'experiments/e7-d5a1-norm.json';primary=json.loads(path.read_text());archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];data=json.loads(gzip.decompress(archive.read_bytes()))
 C=2*s.eye(7)
 for i,j in ((1,3),(3,4),(2,4),(4,5),(5,6),(6,7)):C[i-1,j-1]=C[j-1,i-1]=-1
 om=[s.eye(7).col(i) for i in range(7)];h,w,z=om[5],om[6],om[4];us=[om[0],om[2]-om[0],om[3]-om[2],om[1]+om[4]-om[3],om[4]-om[1]];ys=[2*u-h for u in us];TI=s.Matrix.hstack(h,*ys,w).inv()
 def orbit(v,nodes):
  found={tuple(v)};todo=list(found)
  while todo:
   v=todo.pop()
   for i in nodes:
    u=tuple(v[j]-v[i]*C[j,i] for j in range(7))
    if u not in found:found.add(u);todo.append(u)
  return found
 Z=orbit(z,range(5));RP=orbit(w-z,(0,1,2,3,4,6));assert len(Z)==16 and len(RP)==32;assert all(c.q==1 for v in RP for c in C.inv()*s.Matrix(v));assert RP=={tuple(a) for v in Z for a in (w-s.Matrix(v),h-w-s.Matrix(v))}
 forms=[tuple(fmpq(int(c.p),int(c.q)) for c in TI*s.Matrix(v)) for v in Z]
 expected={(fmpq(5,4),)+tuple(fmpq(v,4) for v in signs)+(fmpq(0),) for signs in product((-1,1),repeat=5) if prod(signs)==1};assert set(forms)==expected
 ctx=fmpq_mpoly_ctx.get(('h','t','v2','v4','p5','v6','v8'));hh,t,E1,E2,pf,E3,E4=ctx.gens();one=ctx.constant(1);zero=ctx.constant(0);ps=[parse(f,ctx) for f in data['spin_power_sums']];assert ps[0]==16*one
 # Independent exact unisolvence for every spin power-sum identity.
 ws=(1,2,4,5,6,8);exps={d:list(powers(ws,d)) for d in range(1,17)};ranks={d:{} for d in exps};rng=Random(3141592653);samples=0;terms=[f.to_dict() for f in ps]
 while any(len(ranks[d])<len(exps[d]) for d in exps):
  samples+=1;assert samples<1000;hv=rng.randrange(-31,32);yy=[rng.randrange(-31,32) for _ in range(5)];sq=[v*v for v in yy];ee=[sum(prod(sq[j] for j in ix) for ix in combinations(range(5),k)) for k in range(1,5)];vals=(hv,0,ee[0],ee[1],prod(yy),ee[2],ee[3]);v6=(vals[0],*vals[2:]);spinvals=[sum(c*v for c,v in zip(f,(hv,*yy,0))) for f in forms]
  for d in exps:
   assert sum(c*prod(v**k for v,k in zip(vals,e)) for e,c in terms[d].items())==sum(v**d for v in spinvals),(samples,d)
   assert all(e[1]==0 and sum(a*b for a,b in zip(e,(1,2,2,4,5,6,8)))==d for e in terms[d]);row=[prod(pow(v%101,k,101) for v,k in zip(v6,e))%101 for e in exps[d]]
   for j,old in sorted(ranks[d].items()):
    if row[j]:
     a=row[j];row=[(x-a*y)%101 for x,y in zip(row,old)]
   j=next((j for j,v in enumerate(row) if v),None)
   if j is not None:
    a=pow(row[j],-1,101);ranks[d][j]=[a*v%101 for v in row]
 elementary=[one]
 for k in range(1,17):elementary.append(sum(((-1)**(j-1)*elementary[k-j]*ps[j] for j in range(1,k+1)),zero)/k)
 assert elementary==[parse(f,ctx) for f in data['spin_elementary']];A=B=zero;rA=zero;rB=one
 for k in range(17):
  c=(-1)**(16-k)*elementary[16-k];A+=c*rA;B+=c*rB;rA,rB=hh*rA+rB,-t*rA
 assert [A,B]==[parse(f,ctx) for f in data['quadratic_remainder']];F=A*A*t+A*B*hh+B*B;assert F==parse(data['norm_classical'],ctx)
 P=fmpq_mpoly_ctx.get(('h','z','c2','c3','c4','c5','r'));h,z,c2,c3,c4,c5,r=P.gens();one=P.constant(1);zero=P.constant(0);es=[one,z,c2,c3,c4,c5];newton=[5*one]
 for k in range(1,9):newton.append(sum(((-1)**(j-1)*es[j]*newton[k-j] for j in range(1,min(k,6))),zero)+((-1)**(k-1)*k*es[k] if k<=5 else zero))
 xp={k:sum((comb(k,j)*2**j*(-h)**(k-j)*newton[j] for j in range(k+1)),zero) for k in (2,4,6,8)};sq=[one]
 for k in range(1,5):sq.append(sum(((-1)**(j-1)*sq[k-j]*xp[2*j] for j in range(1,k+1)),zero)/k)
 pf=sum((2**j*(-h)**(5-j)*es[j] for j in range(6)),zero);expanded=F.compose(h,r+h*z/2-z*z/4,sq[1],sq[2],pf,sq[3],sq[4]);assert expanded==parse(data['norm_root_coefficients'],P);assert all(c.denominator==1 for c in expanded.to_dict().values());mod2=P.from_dict({e:int(c)&1 for e,c in expanded.to_dict().items() if int(c)&1});assert mod2==parse(data['norm_mod2'],P);assert {e:c for e,c in mod2.to_dict().items() if not any(e[:6])}=={(0,0,0,0,0,0,16):one.to_dict()[(0,)*7]}
 out=dict(result='passed',engine='Native SymPy/python-flint: fresh root-character orbit, exact spin moment identities by unisolvence, Newton recurrence, quadratic resultant, and full coefficient integrality. No Sage or project imports.',root_character='omega7-omega5',root_orbit=32,spin_orbit=16,degree=32,monic_last_axis='r^16',unisolvent_samples=samples,seconds=round(perf_counter()-start,3),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()});(ROOT/'experiments/e7-d5a1-norm-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED norm audit',out['seconds'],flush=True)

if __name__=='__main__':main()
