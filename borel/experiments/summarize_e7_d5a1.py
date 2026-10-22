"""Bind the proved E7/P6 exponent, exact prefix and unresolved counts."""
from pathlib import Path
from hashlib import sha256
import json
ROOT=Path(__file__).resolve().parents[1]

def digest(p):return sha256(p.read_bytes() if p.suffix=='.gz' else p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
 names=['e7-d5a1-sources.json','e7-d5a1-sources-audit.json','e7-d5a1-sources-22.json','e7-d5a1-sources-22-audit.json','e7-d5a1-norm.json','e7-d5a1-norm-audit.json','e7-d5a1-images.json','e7-d5a1-images-audit.json','e7-d5a1-moments-audit.json','e7-d5a1-source-algebra.json','e7-d5a1-free-source.json'];paths=[ROOT/'experiments'/name for name in names];checks=[json.loads(p.read_text()) for p in paths]
 for p,c in zip(paths,checks):
  assert c['result']==('incomplete' if 'source-algebra' in p.name or 'free-source.json' in p.name else 'passed'),p
  for name,h in c.get('source_sha256',{}).items():assert digest(ROOT/name)==h,(p,name)
  for name,h in c.get('raw_sha256',{}).items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(p,name)
  if 'certificate_archive' in c:assert sha256((ROOT/c['certificate_archive']).read_bytes()).hexdigest()==c['certificate_archive_sha256']
 images=checks[6];audit=checks[7];assert images['dimension']==audit['dimension']==42;assert images['total_rank']==audit['total_rank']==756;assert not images['complete'] and not audit['full_multiplicities_proved'];assert audit['all_degree_exponent_divides_2'];assert len(audit['degree_rows'])==42
 exact=[0]+[r['multiplicity_2'] if r['full_relative_multiplicity_proved'] else None for r in audit['degree_rows']];upper=[0]+[r['multiplicity_2'] for r in audit['degree_rows']];assert exact[:23]==[0,0,1,1,1,2,3,4,4,5,7,9,9,11,12,15,14,15,17,19,17,17,19];assert exact[41]==0;assert all(r['twice_full_weight_image_in_root_image'] for r in audit['degree_rows']);assert all(r['new_source_generator_image_already_in_lower_image'] for r in audit['degree_rows'] if r['degree'] in (21,22));assert checks[3]['indecomposable_checks']==[dict(degree=21,full_component_rank=548,products_of_lower_generators_rank=547,necessary_new_generators=1),dict(degree=22,full_component_rank=698,products_of_lower_generators_rank=697,necessary_new_generators=1)];assert checks[9]['first_missing_degree']==21
 out=dict(result='passed',parabolic='E7/P6',retained_nodes=[1,2,3,4,5,7],levi_type='D5+A1',coefficient_ring='Z_(2)',root_source='Full invariants of Sym Q',weight_source='Full invariants of Sym X = previously audited D5 algebra A[t]',dimension=42,schubert_rank=756,all_degree_exponent_divides_2=True,first_relative_degree=2,complete_multiplicity_formula=False,exact_multiplicities=exact,multiplicity_upper_bounds=upper,unknown_multiplicity_degrees=[d for d,m in enumerate(exact) if m is None],full_root_source_degree_bound=22,necessary_new_source_generator_degrees=[21,22],new_source_generators_enlarge_borel_image_in_their_degrees=False,proper_E7_parabolics_with_complete_formulas=113,proper_E7_diagrams_without_complete_formulas=14,remaining_distinct_Weyl_conjugacy_classes=10,scope='All-degree exponent theorem and exact prefix, not a complete all-degree multiplicity formula, ambient cohomology cokernel, or originality claim. Earlier full D5 source is a shared audited input.',source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in [Path(__file__).resolve(),*paths]});(ROOT/'experiments/e7-d5a1-summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
