"""Independent adjoint-root construction of the rational E7/A6 Borel ideal.

Girard's partition formula replaces the primary Newton recurrence; the
126 roots replace its 56 minuscule weights. Singular arithmetic is shared.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
from math import factorial
import gzip,json
from sage.all import QQ,ZZ,PolynomialRing,TermOrder,prod,binomial,multinomial
ROOT=Path(__file__).resolve().parents[1]
DEGREES=(1,3,4,5,6,7)
def powers(ds,n):
 if not ds:
  if n==0:yield ()
  return
 for j in range(n//ds[0]+1):
  for e in powers(ds[1:],n-j*ds[0]):yield (j,)+e

def main():
 start=perf_counter();primary_path=ROOT/'experiments/e7-a6-images.json';primary=json.loads(primary_path.read_text());assert primary['result']=='passed'
 archive=ROOT/primary['certificate_archive'];assert sha256(archive.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];cert=json.loads(gzip.decompress(archive.read_bytes()))
 P=PolynomialRing(QQ,6,names=('c1','c3','c4','c5','c6','c7'),order=TermOrder('wdegrevlex',DEGREES));c1,c3,c4,c5,c6,c7=P.gens();es=[P.one(),c1,QQ(4)/9*c1*c1,c3,c4,c5,c6,c7];t=c1/3
 ps=[P(7)]
 for k in range(1,19):
  f=P.zero()
  for e in powers(tuple(range(1,8)),k):
   length=sum(e);coefficient=QQ((-1)**(k-length)*k*factorial(length-1))/prod(factorial(a) for a in e)
   f+=coefficient*prod(g**j for g,j in zip(es[1:],e))
  ps.append(P(f))
 triples=[]
 for k in range(19):
  all0=sum(multinomial(i,j,k-i-j)*ps[i]*ps[j]*ps[k-i-j] for i in range(k+1) for j in range(k-i+1))
  doubles=sum(binomial(k,j)*2**j*ps[j]*ps[k-j] for j in range(k+1))
  triples.append(P((all0-3*doubles+2*3**k*ps[k])/6))
 roots=[]
 for d in (2,6,8,10,12,14,18):
  differences=sum((-1)**j*binomial(d,j)*ps[j]*ps[d-j] for j in range(d+1))
  shifts=2*sum(binomial(d,k)*(-t)**(d-k)*(triples[k]+ps[k]) for k in range(d+1))
  roots.append(P(differences+shifts))
 assert roots[0]==0
 ideal=P.ideal(roots[1:]);G=ideal.groebner_basis();assert ideal.vector_space_dimension()==576
 old=[P(f) for f in cert['groebner_basis']]
 assert all(not f.reduce(G) for f in old) and all(not f.reduce(old) for f in G)
 # Also check every primary moment belongs to the independently rebuilt ideal.
 assert all(not P(f).reduce(G) for f in cert['relations'])
 output=dict(result='passed',engine='Fresh 126-root power sums using Girard partitions; same Sage/Singular polynomial arithmetic backend as primary',dimension=576,groebner_basis_elements=len(G),primary_ideal_equals_adjoint_ideal=True,seconds=round(perf_counter()-start,3),relations=[str(f) for f in roots[1:]],groebner_basis=[str(f) for f in G],scope='Independent construction of the rational ideal only; the separate plain-Python audit checks full integral images and multiplicities.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),primary_path)},raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
 (ROOT/'experiments/e7-a6-ideal-audit.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED',output['seconds'],flush=True)
if __name__=='__main__':main()
