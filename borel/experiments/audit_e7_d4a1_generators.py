"""Fresh D4+A1 coefficient-lattice and actual integral-generator audit."""
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from math import comb,gcd,lcm
from time import perf_counter
import gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq,fmpq_mpoly_ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1]

def main():
 start=perf_counter();sp=ROOT/'experiments/e7-d4a1-sources.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];assert sha256(sa.read_bytes()).hexdigest()==sr['certificate_archive_sha256'];sc=json.loads(gzip.decompress(sa.read_bytes()));P=fmpq_mpoly_ctx.get(('h','y','b','c2','c3','c4','r'));h,y,b,c2,c3,c4,r=P.gens();one=P.constant(1);zero=P.constant(0)
 def parse(text):
  names=dict(zip(P.names(),range(7)));out={}
  for term in text.replace(' ','').replace('-','+-').split('+'):
   if not term:continue
   c=fmpq(-1 if term.startswith('-') else 1);e=[0]*7
   for f in term.lstrip('-').split('*'):
    n,_,k=f.partition('^')
    if n in names:e[names[n]]+=int(k) if k else 1
    else:c*=fmpq(f)
   e=tuple(e);out[e]=out.get(e,fmpq(0))+c
  return P.from_dict(out)
 es=[one,b-y,c2,c3,c4];ps=[4*one]
 for k in range(1,7):ps.append(sum(((-1)**(j-1)*es[j]*ps[k-j] for j in range(1,min(k,5))),zero)+((-1)**(k-1)*k*es[k] if k<=4 else zero))
 xp={k:sum((comb(k,j)*2**j*(-h)**(k-j)*ps[j] for j in range(k+1)),zero) for k in (2,4,6)};sq=[one]
 for k in range(1,4):sq.append(sum(((-1)**(j-1)*sq[k-j]*xp[2*j] for j in range(1,k+1)),zero)/k)
 pf=sum((2**j*(-h)**(4-j)*es[j] for j in range(5)),zero);raw=[h,y,r+h*b/2-b*b/4,sq[1],sq[2],pf,sq[3]];factors=[fmpq(gcd(*(int(c.numerator) for c in f.to_dict().values())),lcm(*(int(c.denominator) for c in f.to_dict().values()))) for f in raw];orig=[f/c for f,c in zip(raw,factors)];assert orig==list(map(parse,sc['original_generators'])) and list(map(str,factors))==sc['normalization_factors']
 @lru_cache(None)
 def mon(e):
  if not any(e):return one
  i=next(i for i,k in enumerate(e) if k);u=list(e);u[i]-=1;return orig[i]*mon(tuple(u))
 verified=[]
 for row in sc['degrees']:
  for j in row['new_generator_indices']:
   f=sum((mon(tuple(e))*fmpq(line[j]) for e,line in zip(row['classical_exponents'],row['rational_coordinates']) if fmpq(line[j])),zero);assert all(c.denominator==1 for c in f.to_dict().values());text,d=sr['generators'][len(verified)];assert f==parse(text) and d==row['degree'];verified.append(dict(degree=d,integral_actual_invariant=True,terms=len(f.to_dict())))
 assert len(verified)==18
 C=2*S.eye(7)
 for i,j in ((1,3),(3,4),(2,4),(4,5),(5,6),(6,7)):C[i-1,j-1]=C[j-1,i-1]=-1
 om=[S.eye(7).col(i) for i in range(7)];hh,yy,zz,ww=om[5],om[0],om[4],om[6];us=[om[2]-om[0],om[3]-om[2],om[1]+om[4]-om[3],om[4]-om[1]];pp=ww-zz;QB=S.Matrix.hstack(hh,yy,2*zz,*us[:3],pp);XB=S.Matrix.hstack(hh,yy,zz,*us[:3],ww);assert abs(XB.det())==1 and abs((C.inv()*QB).det())==1 and all(S.denom(c)==1 for c in C.inv()*QB);assert sum(us,S.zeros(7,1))==2*zz-yy
 def reflection(v,i):return v-v[i-1]*C.col(i-1)
 for i in (2,3,4):assert set(tuple(reflection(u,i)) for u in us)==set(map(tuple,us)) and all(reflection(v,i)==v for v in (hh,yy,zz,pp,ww))
 ys=[2*u-hh for u in us]
 for i in (2,3,4,5):
  actual={tuple(reflection(v,i)) for v in ys};assert actual<={tuple(e*v) for v in ys for e in (-1,1)}
  signs=[]
  for v in ys:
   vv=reflection(v,i);signs.append(1 if tuple(vv) in set(map(tuple,ys)) else -1)
  assert signs.count(-1)%2==0
 assert all(reflection(v,7)==v for v in (hh,yy,zz,*us)) and reflection(pp,7)==hh-2*zz-pp and reflection(ww,7)==hh-ww
 out=dict(result='passed',actual_root_basis_determinant=1,actual_weight_basis_determinant=1,full_coefficient_rings_verified=True,coefficient_rings=dict(M='R[h,y,b,c2,c3,c4,r], sum(u)=b-y, r=(w-z)(h-2z-(w-z))',N='R[h,y,z,c2,c3,c4,t], sum(u)=2z-y, t=w(h-w)'),rational_invariant_degrees=[1,1,2,2,4,4,6],generators=verified,seconds=round(perf_counter()-start,3),scope='Each recorded root generator is an actual integral invariant in the exact root coefficient ring. This audit alone does not prove that they generate the full invariant algebra.',engine='Native FLINT Newton identities and exact rational source preimages; fresh SymPy Cartan lattice determinants and S4, even-sign and A1 actions. No Sage or primary experiment imports.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),sp)},raw_sha256={sa.relative_to(ROOT).as_posix():sha256(sa.read_bytes()).hexdigest()});(ROOT/'experiments/e7-d4a1-generators-audit.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED D4+A1 actual integral generators',out['seconds'],flush=True)

if __name__=='__main__':main()
