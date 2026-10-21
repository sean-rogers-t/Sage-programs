"""Plain Python check of E7 common-Levi transfers and their exact scope."""
from pathlib import Path
from hashlib import sha256
import json
ROOT=Path(__file__).resolve().parents[1]
C=[[2,0,-1,0,0,0,0],[0,2,0,-1,0,0,0],[-1,0,2,-1,0,0,0],[0,-1,-1,2,-1,0,0],[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,-1],[0,0,0,0,0,-1,2]]
EDGES=((1,3),(3,4),(2,4),(4,5),(5,6),(6,7))
def multiply(f,g):
 h=[0]*(len(f)+len(g)-1)
 for i,a in enumerate(f):
  for j,b in enumerate(g):h[i+j]+=a*b
 return h

def divide(f,g):
 f=list(f);out=[0]*(len(f)-len(g)+1)
 for k in range(len(out)-1,-1,-1):
  assert f[k+len(g)-1]%g[-1]==0;c=f[k+len(g)-1]//g[-1];out[k]=c
  for j,b in enumerate(g):f[k+j]-=c*b
 assert not any(f)
 while len(out)>1 and out[-1]==0:out.pop()
 return out

def brackets(ds):
 p=[1]
 for d in ds:p=multiply(p,[1]*d)
 return p

def rank2(rows):
 pivots={}
 for row in rows:
  while row:
   k=row.bit_length()-1
   if k in pivots:row^=pivots[k]
   else:pivots[k]=row;break
 return len(pivots)

def main():
 names=('e7-a6-images.json','e7-a6-images-audit.json','e7-a6-ideal-audit.json','e7-full-flag-kernel.json','e7-full-flag-kernel-audit.json','e7-two-primary-ring-check.json');data={n:json.loads((ROOT/'experiments'/n).read_text()) for n in names};assert all(d['result']=='passed' for d in data.values())
 WA=brackets(range(2,8));WG=brackets((2,6,8,10,12,14,18));full=multiply(data[names[0]]['multiplicity_coefficients'],WA);candidate=[0]+brackets((6,8,10,12,14,18));assert full==candidate
 assert sum(full)==1451520 and full[-1]==1
 for case in data['e7-full-flag-kernel.json']['degrees']:assert full[case['degree']]==case['multiplicity_2']
 rows=[]
 for mask in range(128):
  I=[i+1 for i in range(7) if mask>>i&1];S=set(I);left=set(I);degrees=[];components=[];type_a=True
  while left:
   component={min(left)};todo=list(component)
   while todo:
    node=todo.pop()
    for edge in EDGES:
     if node in edge:
      other=edge[1] if edge[0]==node else edge[0]
      if other in left and other not in component:component.add(other);todo.append(other)
   left-=component;components.append(sorted(component))
   if any(sum(node in e and set(e)<=component for e in EDGES)>2 for node in component):type_a=False
   degrees+=list(range(2,len(component)+2))
  rank=rank2([sum((C[r][i-1]%2)<<j for j,i in enumerate(I)) for r in range(7)])
  if not type_a or rank!=len(I):continue
  PI=brackets(degrees);m=divide(full,PI);target=divide(WG,PI)
  assert multiply(m,[1,1])==[0]+target
  assert all(c>=0 for c in m) and all(c<=target[d] for d,c in enumerate(m))
  rows.append(dict(retained_nodes=I,components=components,multiplicity_coefficients=m,cohomology_ranks=target))
 assert len(rows)==108
 within=[r for r in rows if 7 not in r['retained_nodes']];assert len(within)==60
 E6=brackets((2,5,6,8,9,12));p7=divide(full,E6)
 assert p7==data['e7-two-primary-ring-check.json']['multiplicity_coefficients']
 output=dict(result='passed',full_flag_multiplicity_coefficients=full,good_levi_parabolics=rows,proper_parabolics_by_transfer=108,additional_P7_separate_basis_proof=True,proper_parabolics_settled=109,P7_refinements_settled=61,P7_refinements_remaining=[[2,3,4,5],[1,2,3,4,5],[2,3,4,5,6]],scope='All-degree formulas for type A Levi diagrams with full mod-2 Cartan column rank; P7 is added by its separate basis theorem. No formula is asserted for other excluded parabolics.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+[ROOT/'experiments'/n for n in names]})
 (ROOT/'experiments/e7-full-flag-transfer.json').write_text(json.dumps(output,indent=2)+'\n');print('PASSED: full flag, 108 good-Levi parabolics, P7 separately; 61/64 refinements')
if __name__=='__main__':main()
