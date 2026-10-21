"""Direct full-source BGG test on the A5 refinement of E7/P7."""
from sage.all import QQ,ZZ,PolynomialRing,matrix,prod
from itertools import combinations
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import argparse,json
from sage_borel import BorelMap
from e7_two_primary_basis import powers
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--start-degree',type=int,default=0);ap.add_argument('--max-degree',type=int,default=9);args=ap.parse_args()
 start=perf_counter();b=BorelMap('E7',(2,7));a=b.simple_roots;chain=(1,3,4,5,6)
 vs=[a[7]+sum(a[j] for j in chain[i:]) for i in range(6)];u=a[2]+sum(vs[:3]);z=(u-sum(vs))/2
 es=[sum(prod(sub) for sub in combinations(vs,k)) for k in range(1,7)]
 gs={'M':es+[u],'N':es+[z]};weights=(1,2,3,4,5,6,1);out=[]
 for d in range(args.start_degree,args.max_degree+1):
  labels=b.basis(d);exps=list(powers(weights,d));images={}
  for kind,inputs in gs.items():
   cs=[]
   for e in exps:
    f=prod(g**k for g,k in zip(inputs,e));v=b.expand(f);cs.append([v.get(t,0) for t in labels])
   images[kind]=matrix(QQ,cs).transpose();assert images[kind].denominator()==1
  H=matrix(ZZ,images['N']).transpose().hermite_form(include_zero_rows=False).transpose();assert H.is_square()
  A=H.change_ring(QQ).inverse()*images['M'];assert A.denominator()==1
  snf=matrix(ZZ,A).elementary_divisors();m=sum(int(n.valuation(2)) for n in snf);assert all(n.valuation(2)<=1 for n in snf)
  # A5 in the six-dimensional common complement: degrees 1,2,3,4,5,6.
  def count(k):return len(list(powers((1,2,3,4,5,6),k))) if k>=0 else 0
  expected=count(d-1)-count(d-7)-count(d-9)
  if d<=9:assert m==expected,(d,m,expected)
  out.append(dict(degree=d,target_rank=len(labels),source_monomials=len(exps),multiplicity_2=m,proved_truncation_expected=expected))
  print('A5 direct',out[-1],round(perf_counter()-start,2),flush=True)
 result=dict(result='passed',retained_nodes=[1,3,4,5,6],degree_interval=[args.start_degree,args.max_degree],degrees=out,seconds=round(perf_counter()-start,3),scope='Direct full-source BGG/Hermite/Smith calculations in degrees in the recorded interval only; no transfer formula supplies the image columns.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'src/sage_borel/core.py',ROOT/'experiments/e7_two_primary_basis.py')})
 (ROOT/('experiments/e7-refinement-bgg-check'+('-higher' if args.start_degree else '')+'.json')).write_text(json.dumps(result,indent=2)+'\n')
 print('PASSED',result['seconds'],flush=True)
if __name__=='__main__':main()
