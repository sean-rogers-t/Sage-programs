"""Bind the completed D6 comparison and the remaining-family reductions."""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s
ROOT=Path(__file__).resolve().parents[1]

def main():
 names=('e7-d6-sources.json','e7-d6-sources-audit.json','e7-d6-images.json','e7-d6-images-audit.json','e7-d6-moments-audit.json','e7-remaining-conjugacy.json','e7-remaining-conjugacy-audit.json','e7-d5a1-degree-two.json','e7-d5a1-degree-two-audit.json');data=[json.loads((ROOT/'experiments'/n).read_text()) for n in names];assert all(d['result']=='passed' for d in data)
 src,sa,im,ia,mom,conj,ca,d5a1,d5audit=data;assert src['source_degree_bounds']==sa['source_degree_bounds']==dict(N=33,M=27);assert im['complete'] and im['max_degree']==ia['max_degree']==33
 assert [r['degree'] for r in im['degrees']]==[r['degree'] for r in ia['degree_rows']]==list(range(1,34));assert len(sa['degrees'])==60;assert mom['basic_invariant_algebra_complete'] and mom['weight_orbit']==56
 coefficients=[0]
 for row,checked in zip(im['degrees'],ia['degree_rows']):
  assert set(row['smith_multiplicities'])<=set(('1',));m=row['smith_multiplicities'].get('1',0);assert m==checked['multiplicity_2'];coefficients.append(m)
  assert checked['twice_weight_image_in_root_image']
  if row['degree']>27:assert m==0 and row['root_image_completeness']=='root-generator lower image equals full weight image'
 old=json.loads((ROOT/'experiments/e7-d4-summary.json').read_text());remaining=[I for I in old['other_proper_parabolics_remaining'] if I!=[2,3,4,5,6,7]];assert len(remaining)==14
 assert conj['distinct_representatives_after_these_reductions']==ca['distinct_Weyl_conjugacy_classes']==11
 q=s.symbols('q');T=sum(c*q**d for d,c in enumerate(coefficients));cohom=(1+q**4+q**8)*(1+q**6+q**12)*sum(q**i for i in range(14));assert s.Poly(cohom,q).degree()==33 and cohom.subs(q,1)==126
 assert [sum(1 for r in sa['degrees'] if r['source']==kind) for kind in ('N','M')]==[33,27]
 assert d5a1['relative_smith_2_valuations']==[0,1] and d5audit['first_nonzero_relative_degree']==2
 expected=q**8*(1+q)*(sum(q**(2*i) for i in range(10))+q**4*sum(q**(2*i) for i in range(6))+q**12);assert s.expand(T-expected)==0 and sum(coefficients)==34
 out=dict(result='passed',retained_nodes=[2,3,4,5,6,7],parabolic='E7/P1',levi_type='D6',dimension=33,schubert_rank=126,multiplicity_coefficients=coefficients,multiplicity_polynomial=str(s.factor(T)),total_relative_factors_2=sum(coefficients),first_nonzero_degree=next(d for d,c in enumerate(coefficients) if c),last_nonzero_degree=max(d for d,c in enumerate(coefficients) if c),relative_exponent=2,source_degree_bounds=src['source_degree_bounds'],bounded_N_generator_degrees=src['sources']['N']['generator_degrees'],proper_parabolics_settled=113,other_proper_parabolics_remaining=remaining,remaining_distinct_Weyl_conjugacy_classes=10,full_unbounded_source_generation_claim=False,scope='Computer-assisted all-degree relative image theorem on E7/P1 only, with full components in the required bounds and equality sandwich above degree 27. Four independently checked Weyl conjugacies reduce the other 14 diagrams to ten distinct cases; their formulas remain unresolved here. Classical Borel and Weyl-invariant inputs are separated. No ambient-image classification, structural explanation or novelty claim.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+[ROOT/'experiments'/n for n in names]+[ROOT/'experiments/e7-d4-summary.json']})
 out['partial_D5_A1']=dict(retained_nodes=[1,2,3,4,5,7],first_nonzero_relative_degree=2,relative_group_in_degree_two='Z/2',higher_multiplicities_unresolved=True);out['first_gap_not_determined_by_degree_one_equality']=True
 (ROOT/'experiments/e7-d6-summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
