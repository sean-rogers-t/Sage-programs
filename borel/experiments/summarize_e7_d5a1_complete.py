"""Promote the audited P6 lower image using the complete root source."""
from pathlib import Path
from hashlib import sha256
from collections import Counter
import json,gzip
ROOT=Path(__file__).resolve().parents[1]

def main():
 paths=[ROOT/'experiments'/name for name in ('e7-d5a1-images.json','e7-d5a1-images-audit.json','e7-d5a1-initial-source.json','e7-d5a1-initial-source-audit.json','e7-d5a1-extra23-audit.json','e7-d5a1-norm-audit.json','e7-d5a1-sources-22-audit.json')];checks=[]
 for p in paths:
  c=json.loads(p.read_text());assert c['result']=='passed',p
  for name,expected in c.get('source_sha256',{}).items():assert sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==expected,(p,name)
  for name,expected in c.get('raw_sha256',{}).items():assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,(p,name)
  if 'certificate_archive' in c:assert sha256((ROOT/c['certificate_archive']).read_bytes()).hexdigest()==c['certificate_archive_sha256']
  checks.append(c)
 images,audit,source,source_audit,extra,norm,seed=checks;assert images['total_rank']==audit['total_rank']==756;assert audit['all_degree_exponent_divides_2'];assert source['full_source_complete_all_degrees'] and source_audit['full_source_complete_all_degrees'];assert len(source_audit['generator_degrees'])==25 and source_audit['finite_polynomial_witness_degree']==32;assert extra['borel_image_in_constructed_root_image'] and extra['degree']==23;assert norm['degree']==32
 # The old lower algebra already includes every <=22 generator and the
 # norm. The sole additional degree-23 source generator has its Borel
 # image in that lower algebra. Hence the lower algebra IS the full image.
 coefficients=[0]+[row['multiplicity_2'] for row in audit['degree_rows']];assert len(coefficients)==43;out=dict(result='passed',complete_multiplicity_formula=True,parabolic='E7/P6',retained_nodes=[1,2,3,4,5,7],levi_type='D5+A1',coefficient_ring='Z_(2)',dimension=42,schubert_rank=756,full_root_source_complete_all_degrees=True,full_root_source_generator_degrees=source_audit['generator_degrees'],full_root_source_generators=25,source_witness_top_degree=32,all_degree_exponent_divides_2=True,first_relative_degree=2,multiplicity_coefficients=coefficients,total_nonunit_factors=sum(coefficients),all_nonunit_factors=2,unknown_multiplicity_degrees=[],old_constructed_root_image_equals_full_root_image=True,degree23_source_generator_borel_redundant=True,scope='Complete relative full root/weight Borel images in every degree for E7/P6. Source completeness follows from an independently checked all-degree initial-monoid Hilbert identity. The degree-23 correction is necessary in the source and redundant in the Borel image. Ambient cohomology/image and novelty claims are separate.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+paths});(ROOT/'experiments/e7-d5a1-complete-summary.json').write_text(json.dumps(out,indent=2)+'\n');print('PASSED complete E7/P6 formula',sum(coefficients),coefficients,flush=True)

if __name__=='__main__':main()
