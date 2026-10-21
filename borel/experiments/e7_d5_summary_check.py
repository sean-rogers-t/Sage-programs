"""Integer-only D5 formula, overlap and completion-scope audit."""
from pathlib import Path
from hashlib import sha256
import json
ROOT=Path(__file__).resolve().parents[1]
def multiply(f,g):
 h=[0]*(len(f)+len(g)-1)
 for i,a in enumerate(f):
  for j,b in enumerate(g):h[i+j]+=a*b
 return h

def main():
 names=('e7-d5-images.json','e7-d5-images-audit.json','e7-d5-source-audit.json','e7-d5-moments-audit.json','e7-d5-kernel.json','e7-refinement-independent.json','e7-full-flag-transfer.json');data={n:json.loads((ROOT/'experiments'/n).read_text()) for n in names};assert all(d['result']=='passed' for d in data.values())
 f=[0,1]
 for g in ([1]+[0]*4+[1],[1,0,0,0,1,0,0,0,1],[int(i%2==0) for i in range(13)],[1]*18):f=multiply(f,g)
 p=data[names[0]];a=data[names[1]];assert f==p['multiplicity_coefficients']==[0]+[r['multiplicity_2'] for r in a['degree_rows']]
 assert len(f)==44 and sum(f)==756 and max(f)==38 and f[-1]==1
 ranks=[int(i>0)*f[i]+(f[i+1] if i<43 else 0) for i in range(44)];assert sum(ranks)==1512
 assert ranks[1:]==[r['rank'] for r in a['degree_rows']] and ranks==list(reversed(ranks))
 assert all(r['multiplicity_2']==f[r['degree']] for r in data[names[4]]['degrees'])
 for I in ([1,2,3,4,5],[2,3,4,5,6]):
  row=next(r for r in data[names[5]]['refinements'] if r['retained_nodes']==I);assert row['multiplicity_2']==f[:10]
 assert data[names[5]]['d5_symmetry_integral']
 baseline=data[names[6]];remaining=baseline['P7_refinements_remaining'];assert remaining==[[2,3,4,5],[1,2,3,4,5],[2,3,4,5,6]]
 report=dict(result='passed',formula='q(1+q^5)(1+q^4+q^8)[7]_(q^2)[18]_q',multiplicity_coefficients=f,cohomology_ranks=ranks,total_relative_factors_2=756,proper_parabolics_settled=baseline['proper_parabolics_settled']+2,P7_refinements_settled=baseline['P7_refinements_settled']+2,P7_refinements_remaining=[[2,3,4,5]],scope='Both D5 cases follow from the full-source and image proofs plus the previously certified E6 opposition. D4 remains; other excluded E7 parabolics and ambient H/N are separate.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+[ROOT/'experiments'/n for n in names]})
 assert report['proper_parabolics_settled']==111 and report['P7_refinements_settled']==63
 (ROOT/'experiments/e7-d5-summary.json').write_text(json.dumps(report,indent=2)+'\n');print('PASSED: exact D5 formula, both overlap checks, 111 proper parabolics, 63/64 P7 refinements')
if __name__=='__main__':main()
