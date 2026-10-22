"""Primitive integral Chern coefficients of the 32-character root orbit.

Compute the complete characteristic polynomial, not just its constant
term. Every coefficient, and its coefficientwise primitive normalization,
is an actual integral D5+A1 root invariant. No generation claim is made.
"""
from pathlib import Path
from hashlib import sha256
from math import comb,gcd
from time import perf_counter
import gzip,json,os,sys
if os.environ.get('PYTHON_FLINT_PATH'):sys.path.insert(0,os.environ['PYTHON_FLINT_PATH'])
from flint import fmpq_mpoly_ctx
ROOT=Path(__file__).resolve().parents[1]

def parse(text,P):
 from flint import fmpq
 ix={str(v):j for j,v in enumerate(P.names())};out={}
 for term in text.replace(' ','').replace('-','+-').split('+'):
  if not term:continue
  c=fmpq(-1 if term.startswith('-') else 1);e=[0]*len(ix)
  for f in term.lstrip('-').split('*'):
   name,_,power=f.partition('^')
   if name in ix:e[ix[name]]+=int(power) if power else 1
   else:c*=fmpq(f)
  e=tuple(e);out[e]=out.get(e,fmpq(0))+c
 return P.from_dict(out)

def main():
 start=perf_counter();np=ROOT/'experiments/e7-d5a1-norm.json';nr=json.loads(np.read_text());na=ROOT/nr['certificate_archive'];assert sha256(na.read_bytes()).hexdigest()==nr['certificate_archive_sha256'];nc=json.loads(gzip.decompress(na.read_bytes()));sp=ROOT/'experiments/e7-d5a1-sources.json';sr=json.loads(sp.read_text());sa=ROOT/sr['certificate_archive'];assert sha256(sa.read_bytes()).hexdigest()==sr['certificate_archive_sha256'];sc=json.loads(gzip.decompress(sa.read_bytes()));S=fmpq_mpoly_ctx.get(('h','t','v2','v4','p5','v6','v8'));h,t,*_=S.gens();one=S.constant(1);zero=S.constant(0);es=[parse(f,S) for f in nc['spin_elementary']];WA=[zero,one];WB=[one,zero]
 for k in range(2,17):WA.append(h*WA[-1]-t*WA[-2]);WB.append(h*WB[-1]-t*WB[-2])
 A=[sum(((-1)**(16-k+i)*comb(k,i)*es[16-k]*WA[k-i] for k in range(i,17)),zero) for i in range(17)];B=[sum(((-1)**(16-k+i)*comb(k,i)*es[16-k]*WB[k-i] for k in range(i,17)),zero) for i in range(17)];F=[sum((t*A[j]*A[k]+h*A[j]*B[k]+B[j]*B[k] for j in range(17) for k in range(17) if j+k==i),zero) for i in range(33)];assert F[0]==parse(nc['norm_classical'],S) and F[32]==one
 P=fmpq_mpoly_ctx.get(('h','z','c2','c3','c4','c5','r'));orig=[parse(f,P) for f in sc['original_generators']];raw=[f*__import__('flint').fmpq(k) for f,k in zip(orig,sc['normalization_factors'])];rows=[]
 for shift in (0,1):
  # Changing p to p+h shifts the characteristic polynomial U -> U-h.
  coeffs=F if not shift else [sum((comb(j,i)*(-h)**(j-i)*F[j] for j in range(i,33)),zero) for i in range(33)]
  for d in range(1,33):
   f=(-1)**d*coeffs[32-d];expanded=f.compose(*raw);assert all(c.denominator==1 for c in expanded.to_dict().values());div=gcd(*(int(c) for c in expanded.to_dict().values()));div=abs(div);assert div;f=f/div;expanded=expanded/div;mod=P.from_dict({e:int(c)&1 for e,c in expanded.to_dict().items() if int(c)&1});rows.append(dict(orbit_shift_by_h=shift,degree=d,primitive_divisor=str(div),classical_polynomial=str(f),root_polynomial_mod2=str(mod),root_terms=len(expanded.to_dict()),mod2_terms=len(mod.to_dict())));print('orbit Chern',shift,d,'terms',len(expanded.to_dict()),'divisor',div,round(perf_counter()-start,2),flush=True)
 archive=ROOT/'experiments/e7-d5a1-orbit-chern-certificates.json.gz';archive.write_bytes(gzip.compress(json.dumps(dict(quadratic_remainder_A=list(map(str,A)),quadratic_remainder_B=list(map(str,B)),characteristic_coefficients=list(map(str,F)),generators=rows),separators=(',',':')).encode(),mtime=0));report=dict(result='passed',orbits=2,orbit_size=32,generators=64,max_degree=32,seconds=round(perf_counter()-start,3),scope='Explicit coefficient-integral root invariants only; no source or Borel-image completeness claim.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),np,sp)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (na,sa)});(ROOT/'experiments/e7-d5a1-orbit-chern.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
