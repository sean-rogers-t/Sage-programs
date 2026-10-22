"""Bind complete relative image formulas for all 127 proper E7 parabolics."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from collections import Counter
import json
from e7_even_block_formula import trim,add,mul,bracket
ROOT=Path(__file__).resolve().parents[1]
EDGES=((1,3),(3,4),(2,4),(4,5),(5,6),(6,7))

def divide(a,b):
 a=trim(a);b=trim(b);assert b[-1]==1;out=[0]*max(1,len(a)-len(b)+1)
 while len(a)>=len(b) and any(a):
  k=len(a)-len(b);c=a[-1];out[k]=c
  for j,v in enumerate(b):a[k+j]-=c*v
  a=trim(a)
 assert not any(a);return trim(out)
def poincare(degrees):
 out=[1]
 for d in degrees:out=mul(out,bracket(d))
 return out
def levi_degrees(nodes):
 todo=set(nodes);out=[];types=[]
 while todo:
  component={todo.pop()};pending=list(component)
  while pending:
   i=pending.pop()
   for a,b in EDGES:
    j=b if i==a else a if i==b else None
    if j in todo:todo.remove(j);component.add(j);pending.append(j)
  adjacency={i:{b if a==i else a for a,b in EDGES if (a==i and b in component) or (b==i and a in component)} for i in component};n=len(component);branch=[i for i,neighbors in adjacency.items() if len(neighbors)==3]
  if not branch:degrees=list(range(2,n+2));types.append('A'+str(n))
  else:
   assert len(branch)==1;center=branch[0];arms=[]
   for first in adjacency[center]:
    length=1;previous=center;current=first
    while adjacency[current]-{previous}:
     nxt=(adjacency[current]-{previous}).pop();previous,current=current,nxt;length+=1
    arms.append(length)
   arms.sort()
   if arms[:2]==[1,1]:degrees=list(range(2,2*n,2))+[n];types.append('D'+str(n))
   else:assert arms==[1,2,2] and n==6;degrees=[2,5,6,8,9,12];types.append('E6')
  out+=degrees
 return sorted(out),'+'.join(sorted(types)) if types else 'trivial'

def main():
 paths=[]
 def read(name):
  p=ROOT/'experiments'/name;c=json.loads(p.read_text());assert c['result']=='passed',name;paths.append(p)
  for filename,expected in c.get('source_sha256',{}).items():assert sha256((ROOT/filename).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected,(name,filename)
  for filename,expected in c.get('raw_sha256',{}).items():assert sha256((ROOT/filename).read_bytes()).hexdigest()==expected,(name,filename)
  return c
 old=read('e7-d4-summary.json');assert old['proper_parabolics_settled']==112;d6=read('e7-d6-summary.json');d5=read('e7-d5a1-complete-summary.json');d4=read('e7-d4a1-tower-images-audit.json');assert d4['complete'];assert read('e7-d4a1-weight-transfer-audit.json')['full_weight_source_transfer_verified'];formula=read('e7-even-block-formula.json');read('e7-a5-fiber-transfers-audit.json');read('e7-remaining-conjugacy-audit.json');T0=formula['base_multiplicity_polynomial'];T42=formula['cases'][0]['multiplicity_coefficients'];T222=formula['cases'][1]['multiplicity_coefficients'];partial=bracket(3);full=mul(bracket(2),partial)
 tables={
  (2,5,7):('A5-even-blocks-222/full-A2',mul(T222,full)),
  (1,2,5,7):('A5-even-blocks-222/partial-A2',mul(T222,partial)),
  (2,3,5,7):('A5-even-blocks-222/partial-A2/conjugate',mul(T222,partial)),
  (2,4,5,7):('A5-even-blocks-42/full-A2',mul(T42,full)),
  (2,5,6,7):('A5-even-blocks-42/full-A2/conjugate',mul(T42,full)),
  (1,2,3,5,7):('A5-even-blocks-222',T222),
  (1,2,4,5,7):('A5-even-blocks-42/partial-A2',mul(T42,partial)),
  (1,2,5,6,7):('A5-even-blocks-42/partial-A2/conjugate',mul(T42,partial)),
  (2,3,5,6,7):('A5-even-blocks-42/partial-A2/conjugate',mul(T42,partial)),
  (2,3,4,5,7):('D4+A1',d4['multiplicity_coefficients']),
  (2,4,5,6,7):('A5/full-A2',mul(T0,full)),
  (1,2,3,4,5,7):('D5+A1',d5['multiplicity_coefficients']),
  (1,2,3,5,6,7):('A5-even-blocks-42',T42),
  (1,2,4,5,6,7):('A5/partial-A2',mul(T0,partial)),
  (2,3,4,5,6,7):('D6',d6['multiplicity_coefficients'])}
 assert len(tables)==15;PW=poincare([2,6,8,10,12,14,18]);records=[];gaps=Counter();families=Counter()
 for count in range(7):
  for I in combinations(range(1,8),count):
   ds,typ=levi_degrees(I);PI=divide(PW,poincare(ds));dim=63-sum(d-1 for d in ds);assert len(PI)==dim+1 and PI==PI[::-1] and all(c>=0 for c in PI)
   if {2,5,7}<=set(I):family,T=tables[I]
   else:family='uniform-112';T=divide([0]+PW,mul([1,1],poincare(ds)))
   T=trim(T);assert len(T)<=len(PI);T=T+[0]*(len(PI)-len(T));assert all(0<=m<=n for m,n in zip(T,PI));first=next(d for d,m in enumerate(T) if m);gaps[first]+=1;families['uniform-112' if family=='uniform-112' else 'remaining-15']+=1;records.append(dict(retained_nodes=list(I),omitted_nodes=[i for i in range(1,8) if i not in I],levi_type=typ,dimension=dim,schubert_rank=sum(PI),cohomology_ranks=PI,multiplicity_coefficients=T,total_factors_2=sum(T),first_relative_degree=first,proof_family=family))
 assert len(records)==127 and families==dict([('uniform-112',112),('remaining-15',15)]);assert gaps=={1:112,2:12,4:2,8:1};out=dict(result='passed',complete=True,ambient_type='simply connected E7',coefficient_ring='Z_(2), determining the integral relative group since Q=X at odd primes',proper_standard_parabolics=127,improper_G_over_G='trivial comparison',all_nonunit_relative_factors=2,relative_group='N_I^d/M_I^d = (Z/2)^(m_I,d)',root_source='Full (Sym^d Q)^(W_I)',weight_source='Full (Sym^d X)^(W_I)',polynomial_degree_convention='d means cohomology degree 2d',first_gap_counts=dict(sorted(gaps.items())),uniform_formula_parabolics=112,remaining_formula_parabolics=15,all_remaining_Weyl_types=11,parabolics=records,Springer_normalizations=dict(minimal_56='J56^d=12^d M_I^d',adjoint_133='J133^d=36^d N_I^d',comparison='(J56^d+J133^d)/J56^d is isomorphic to N_I^d/M_I^d'),scope='Complete relative full integral root/weight Borel-image comparison for every proper standard E7 parabolic. Ambient cohomology/image cokernels, other groups, geometric folding explanations and literature priority are separate. Classical Weyl/Poincare inputs are distinguished from the repository source and image deductions.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve(),ROOT/'experiments/e7_even_block_formula.py']+paths});(ROOT/'experiments/e7-all-parabolics-summary.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED all 127 proper E7 relative image formulas',dict(gaps),flush=True)

if __name__=='__main__':main()
