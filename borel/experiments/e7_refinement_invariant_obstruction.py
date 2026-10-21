"""Certificate for an integral-invariant obstruction to a naive common basis.

The six rational basic E6 degrees do not supply six integral generators at 2.
A degree-16 Chern invariant separates their common mod-2 zero from the origin.
This is a source-ring obstruction, not a counterexample to image multiplicities.
"""
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,prod,gcd
from itertools import combinations
from hashlib import sha256
from pathlib import Path
from time import perf_counter
import json
from e7_two_primary_basis import powers
ROOT=Path(__file__).resolve().parents[1]
def main():
 start=perf_counter();R=PolynomialRing(QQ,6,'a');a=R.gens();X=[[5,0,4,3,2,1],[-1,0,4,3,2,1],[-1,0,-2,3,2,1],[-1,0,-2,-3,2,1],[-1,0,-2,-3,-4,1],[-1,0,-2,-3,-4,-5]]
 xs=[sum(c*t for c,t in zip(row,a)) for row in X];z=sum(c*t for c,t in zip([-3,-6,-6,-9,-6,-3],a));ls=[v+sg*z for v in xs for sg in (-1,1)]+[-v-u for v,u in combinations(xs,2)]
 C=matrix(ZZ,[[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],[0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]])
 weights=[f/2 for f in ls];assert all(f.denominator()==1 for f in weights)
 for i in range(6):
  new=[a[j]-C[j,i]*a[i] for j in range(6)]
  assert set(f(*new) for f in weights)==set(weights)
 ds=(2,5,6,8,9,12);contents=(144,5760,3456,9216,387072,663552)
 ps={d:sum(f**d for f in ls) for d in ds};fs={d:ps[d]/g for d,g in zip(ds,contents)}
 assert all(gcd(ps[d].coefficients())==g for d,g in zip(ds,contents))
 fs[8]=(fs[8]+fs[2]**4)/2;fs[12]=(fs[12]+fs[6]**2)/2
 assert all(f.denominator()==1 for f in fs.values())
 records=[]
 for d in ds:
  exps=list(powers(ds,d));polys=[prod(fs[k]**j for k,j in zip(ds,e)) for e in exps];maps=[f.dict() for f in polys];support=sorted(set().union(*(set(t) for t in maps)))
  M=matrix(GF(2),[[t.get(e,0) for t in maps] for e in support]);assert M.rank()==len(exps)
  records.append(dict(degree=d,rational_invariant_monomials=[list(e) for e in exps],coefficient_rows=len(support),primitive_mod2_rank=int(M.rank())))
 point=(0,0,0,0,0,1);vals=[int(fs[d](*point)) for d in ds];assert all(v%2==0 for v in vals)
 wvals=[ZZ(f(*point)) for f in weights];assert sum(v%2 for v in wvals)==16
 coeff=[ZZ(1)]+[ZZ(0)]*16
 for v in wvals:
  for j in range(16,0,-1):coeff[j]+=v*coeff[j-1]
 assert coeff[16]%2==1
 result=dict(result='passed',basic_degrees=list(ds),basic_values=vals,point=list(point),primitive_basic_degree_certificates=records,weight_forms=[[int(f.monomial_coefficient(x)) for x in a] for f in weights],permutation_reflections=6,chern_invariant_degree=16,chern_invariant_value=int(coeff[16]),odd_weight_values=16,seconds=round(perf_counter()-start,3),scope='The six normalized rational-degree generators do not generate the full 2-local E6 invariant ring. This obstructs that polynomial common-basis proof, not the conjectured E7 image transfer formula.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),ROOT/'experiments/e7_two_primary_basis.py')})
 (ROOT/'experiments/e7-refinement-invariant-obstruction.json').write_text(json.dumps(result,indent=2)+'\n');print('PASSED',result['seconds'],'degree-16 witness',coeff[16],flush=True)
if __name__=='__main__':main()
