"""The degree-32 root-orbit norm for E7/(D5 A1).

For z_lambda=(5h+sum eps_i y_i)/4, product eps=1, take
prod_lambda (t-h*z_lambda+z_lambda^2). It is the orbit norm of
w-z in the actual root lattice. Newton identities and a quadratic
resultant compute it without expanding all 32 linear factors.
"""
from pathlib import Path
from hashlib import sha256
from math import factorial
from time import perf_counter
import gzip,json,sys
from sage.all import QQ,GF,PolynomialRing,TermOrder,SymmetricFunctions,Partitions,prod,binomial
from e7_d5a1_sources import context,CLASSICAL
ROOT=Path(__file__).resolve().parents[1]

def main():
 start=perf_counter();S=PolynomialRing(QQ,7,names=('h','t','v2','v4','p5','v6','v8'),order=TermOrder('wdegrevlex',CLASSICAL));h,t,E1,E2,pf,E3,E4=S.gens();es=[S.one(),E1,E2,E3,E4,pf*pf];SF=SymmetricFunctions(QQ);MM=SF.monomial();EE=SF.elementary()
 def symmetric(part):return S(sum(c*prod(es[j] if j<=5 else S.zero() for j in mu) for mu,c in EE(MM[part])))
 spin={}
 for k in range(17):
  offset=k%2
  spin[k]=S.zero() if k<5*offset else S(sum((QQ(factorial(k))/prod(factorial(2*j+offset) for j in list(part)+[0]*(5-len(part))))*symmetric(part) for part in Partitions((k-5*offset)//2,max_length=5)))*(pf if offset else 1)
 ps=[S(16)]
 for k in range(1,17):ps.append(16*sum(binomial(k,j)*(5*h)**(k-j)*spin[j] for j in range(k+1))/4**k)
 elementary=[S.one()]
 for k in range(1,17):elementary.append(sum((-1)**(j-1)*elementary[k-j]*ps[j] for j in range(1,k+1))/k)
 A=S.zero();B=S.zero();rA=S.zero();rB=S.one()
 for k in range(17):
  c=(-1)**(16-k)*elementary[16-k];A+=c*rA;B+=c*rB;rA,rB=h*rA+rB,-t*rA
 F=S(A*A*t+A*B*h+B*B);P,orig,factors=context();sub=S.hom([f*k for f,k in zip(orig,factors)],P);expanded=sub(F);assert expanded.denominator()==1;P2=P.change_ring(GF(2));mod2=P2(expanded);r=P2.gen(6);assert mod2.subs({P2.gen(j):0 for j in range(6)})==r**16
 archive=ROOT/'experiments/e7-d5a1-norm-certificates.json.gz';cert=dict(spin_power_sums=list(map(str,ps)),spin_elementary=list(map(str,elementary)),quadratic_remainder=[str(A),str(B)],norm_classical=str(F),norm_root_coefficients=str(expanded),norm_mod2=str(mod2));archive.write_bytes(gzip.compress(json.dumps(cert,separators=(',',':')).encode(),mtime=0));report=dict(result='passed',degree=32,orbit_size=32,rational_terms=len(F.dict()),root_coefficient_terms=len(expanded.dict()),mod2_terms=len(mod2.dict()),monic_on_last_axis='r^16',seconds=round(perf_counter()-start,3),scope='Explicit integral root-orbit norm; this alone does not establish a complete invariant generating set.',certificate_archive=archive.relative_to(ROOT).as_posix(),certificate_archive_sha256=sha256(archive.read_bytes()).hexdigest(),source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7_d5a1_sources.py')});(ROOT/'experiments/e7-d5a1-norm.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)

if __name__=='__main__':main()
