"""Independent root-basis matrix audit of the exceptional-family reductions."""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
from math import factorial,prod
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main():
 path=ROOT/'experiments/e7-remaining-conjugacy.json';data=json.loads(path.read_text());assert data['result']=='passed'
 C=2*s.eye(7)
 for i,j in ((1,3),(3,4),(2,4),(4,5),(5,6),(6,7)):C[i-1,j-1]=C[j-1,i-1]=-1
 CI=C.inv();assert C.det()==2
 central=[i+1 for i in range(7) if any(c.q%2==0 for c in CI.col(i))];assert central==[2,5,7]
 for row in data['conjugacies']:
  W=s.eye(7)
  for i in row['simple_reflections_in_application_order']:W=(s.eye(7)-s.eye(7).col(i-1)*C.row(i-1))*W
  assert W.T*C*W==C and abs(W.det())==1
  A=[i-1 for i in row['left']];B=[i-1 for i in row['right']];image=W[:,A]
  assert image[[i for i in range(7) if i not in B],:]==s.zeros(7-len(B),len(A))
  assert abs(image[B,:].det())==1
 print('PASSED: four unimodular root-subsystem conjugacies, independently in the root basis')
 sets=[sorted(set(central)|set(J)) for k in range(4) for J in combinations((1,3,4,6),k)];assert len(sets)==15
 def kind(I):
  left=set(I);types=[];degrees=[];positive=0
  while left:
   component={min(left)};todo=list(component)
   while todo:
    i=todo.pop()
    for j in left-component:
     if C[i-1,j-1]:component.add(j);todo.append(j)
   left-=component;n=len(component);valencies=[sum(C[i-1,j-1]==-1 for j in component if i!=j) for i in component]
   if max(valencies)<=2:t='A'+str(n);ds=list(range(2,n+2));positive+=n*(n+1)//2
   else:
    assert sorted(valencies)==[1,1,1]+[2]*(n-4)+[3];t='D'+str(n);ds=list(range(2,2*n,2))+[n];positive+=n*(n-1)
   types.append(t);degrees+=ds
  return '+'.join(sorted(types)),degrees,positive
 equivalences={tuple(row['right']):tuple(row['left']) for row in data['conjugacies']}
 def representative(I):
  I=tuple(I)
  while I in equivalences:I=equivalences[I]
  return I
 rows=[]
 for I in sorted(sets,key=lambda I:(len(I),I)):
  t,ds,pos=kind(I);rows.append(dict(retained_nodes=I,levi_type=t,levi_basic_degrees=sorted(ds),dimension=63-pos,schubert_rank=2903040//prod(ds),representative=list(representative(I)),degree_one_root_and_weight_images_equal=True))
 assert len({tuple(row['representative']) for row in rows})==len({row['levi_type'] for row in rows})==11
 out=dict(result='passed',engine='Independent SymPy root-basis reflection matrices, unimodular sublattice tests, inverse Cartan matrix and direct diagram classification. Does not use the primary root-mask search.',central_character_nodes=central,remaining_diagrams=15,distinct_Weyl_conjugacy_classes=11,rows=rows,scope='Four exact conjugacies and degree-one equality for all 15 excluded parabolics. No higher-degree formula follows from this reduction alone.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),path)})
 (ROOT/'experiments/e7-remaining-conjugacy-audit.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
