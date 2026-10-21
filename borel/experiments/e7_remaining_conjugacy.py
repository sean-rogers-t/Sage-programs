"""Exact Weyl-conjugacy reductions among the 15 excluded E7 Levi systems."""
from pathlib import Path
from hashlib import sha256
from collections import deque
from time import perf_counter
import json
ROOT=Path(__file__).resolve().parents[1]
C=((2,0,-1,0,0,0,0),(0,2,0,-1,0,0,0),(-1,0,2,-1,0,0,0),(0,-1,-1,2,-1,0,0),(0,0,0,-1,2,-1,0),(0,0,0,0,-1,2,-1),(0,0,0,0,0,-1,2))
PAIRS=(([1,2,5,7],[2,3,5,7]),([2,4,5,7],[2,5,6,7]),([1,2,4,5,7],[1,2,5,6,7]),([1,2,4,5,7],[2,3,5,6,7]))

def reflect(v,i):return tuple(v[j]-v[i]*C[i][j] for j in range(7))
def positive(v):return v if next(x for x in v if x)>0 else tuple(-x for x in v)

def main():
 start=perf_counter();roots=set(C);todo=list(roots)
 while todo:
  v=todo.pop()
  for i in range(7):
   u=reflect(v,i)
   if u not in roots:roots.add(u);todo.append(u)
 assert len(roots)==126;lines=sorted({positive(v) for v in roots});indices={v:j for j,v in enumerate(lines)};perms=[[indices[positive(reflect(v,i))] for v in lines] for i in range(7)]
 def subsystem(I):
  r={positive(C[i-1]) for i in I};todo=list(r)
  while todo:
   v=todo.pop()
   for i in I:
    u=positive(reflect(v,i-1))
    if u not in r:r.add(u);todo.append(u)
  return sum(1<<indices[v] for v in r)
 def act(mask,i):
  result=0
  while mask:
   b=mask&-mask;mask-=b;result|=1<<perms[i][b.bit_length()-1]
  return result
 records=[]
 for left,right in PAIRS:
  initial=subsystem(left);target=subsystem(right);previous={initial:None};queue=deque([initial])
  while target not in previous and queue:
   a=queue.popleft()
   for i in range(7):
    b=act(a,i)
    if b not in previous:previous[b]=(a,i);queue.append(b)
  assert target in previous,(left,right,len(previous));word=[];b=target
  while previous[b] is not None:b,i=previous[b];word.append(i+1)
  word.reverse();b=initial
  for i in word:b=act(b,i-1)
  assert b==target
  records.append(dict(left=left,right=right,simple_reflections_in_application_order=word,root_lines=initial.bit_count(),search_states=len(previous)))
  print('conjugacy',left,right,'word',word,'states',len(previous),flush=True)
 out=dict(result='passed',engine='Native integer E7 root reflections and exact root-subsystem bit masks; no Sage or project imports.',positive_roots=len(lines),conjugacies=records,remaining_diagrams=15,distinct_representatives_after_these_reductions=11,scope='Exact Weyl conjugacies give graded isomorphisms of root/weight Borel-image comparisons; they do not determine the common unknown multiplicity polynomials.',seconds=round(perf_counter()-start,3),source_sha256={Path(__file__).relative_to(ROOT).as_posix():sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
 (ROOT/'experiments/e7-remaining-conjugacy.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
