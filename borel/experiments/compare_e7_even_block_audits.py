"""Bind the two independent complete finite-module image calculations."""
from pathlib import Path
from hashlib import sha256
import argparse,gzip,json
ROOT=Path(__file__).resolve().parents[1]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--blocks',required=True,choices=['42','222']);args=ap.parse_args();kind=args.blocks;pp=ROOT/('experiments/e7-even-block-fast-'+kind+'.json');primary=json.loads(pp.read_text());pa=ROOT/primary['certificate_archive'];assert sha256(pa.read_bytes()).hexdigest()==primary['certificate_archive_sha256'];pc=json.loads(gzip.decompress(pa.read_bytes()));ip=ROOT/('experiments/e7-even-block-'+kind+'-images-audit.json');independent=json.loads(ip.read_text());assert independent['result']=='passed' and independent['complete'];prep=next(ROOT/p for p in independent['source_sha256'] if 'prepared-'+kind+'.json' in p);pr=json.loads(prep.read_text());prep_a=ROOT/pr['certificate_archive'];prep_c=json.loads(gzip.decompress(prep_a.read_bytes()))
 for key in ('blocks','fiber_exponents','fiber_degrees','fiber_relations','fiber_groebner_basis','partial_root_generators','module_components'):assert pc[key]==prep_c[key]
 assert pc['mod_four_column_sha256']==independent['mod_four_column_sha256'];assert len(primary['degrees'])==len(independent['degrees'])
 for a,b in zip(primary['degrees'],independent['degrees']):assert all(a[key]==b[key] for key in ('degree','rank','multiplicity_2','mod_four_unit_rank','mod_four_even_rank'))
 counts=[row['multiplicity_2'] for row in primary['degrees']];out=dict(result='passed',complete=True,blocks=kind,multiplicity_coefficients=counts,total_nonunit_factors=sum(counts),all_nonunit_factors=2,total_rank=primary['total_rank'],top_polynomial_degree=primary['top_polynomial_degree'],primary_and_independent_all_columns_agree=True,scope='Both complete calculations independently reconstruct the finite source module and compute every image coordinate modulo four, using column and row elimination respectively. Their source modules and common fiber bases are separately proved.',source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),pp,ip,prep)},raw_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in (pa,prep_a)});(ROOT/('experiments/e7-even-block-'+kind+'-comparison-audit.json')).write_text(json.dumps(out,indent=2)+'\n');print('PASSED two independent complete image calculations',kind,sum(counts),flush=True)

if __name__=='__main__':main()
