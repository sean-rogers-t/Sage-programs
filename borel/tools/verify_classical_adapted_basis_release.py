"""Bind the finite classical-basis checks to their exact research sources.

This is a consistency audit, not an arbitrary-rank proof or an external review.
"""
from hashlib import sha256
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
PAIRS=(
    ('type-b-adapted-basis-verification.json','experiments/verify_type_b_adapted_basis.py'),
    ('type-b-adapted-basis-audit.json','experiments/audit_type_b_adapted_basis.py'),
    ('type-d-adapted-basis-verification.json','experiments/verify_type_d_adapted_basis.py'),
    ('type-d-adapted-basis-audit.json','experiments/audit_type_d_adapted_basis.py'),
    ('classical-kernel-generator-test.json','experiments/test_classical_kernel_generators.py'),
    ('classical-divided-generator-test.json','experiments/test_classical_divided_generators.py'))

def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
    records=[];bound=[]
    for name,script in PAIRS:
        data=json.loads((ROOT/name).read_text(encoding='utf-8'))
        assert data['status'] in ('passed','completed')
        assert data['source_sha256']==sha256((ROOT/script).read_bytes()).hexdigest(),script
        records.append(data);bound.extend((name,script))
    b,ba,d,da,k,g=records
    assert len(b['actual_embeddings'])==len(ba['embeddings'])==96
    assert len(d['constructed'])==da['constructed_matches']==124
    assert len(d['actual_D'])==len(da['actual_D'])==28
    assert d['multiplication']==da['multiplication']
    assert da['input_sha256']==sha256((ROOT/PAIRS[2][0]).read_bytes()).hexdigest()
    assert len(k['cases'])==len(g['cases'])==9
    assert all(c['all_equal'] for c in k['cases'])
    assert all(r['spans'] for c in g['cases'] for r in c['degrees'])
    for typ,count,source in [('b',26,'experiments/type_b_multiplicity_search_results.json'),
                             ('d',28,'type-d-adapted-basis-verification.json')]:
        name=f'type-{typ}-basis-chevalley-audit.json'
        script=f'experiments/audit_type_{typ}_basis_chevalley.py'
        record=json.loads((ROOT/name).read_text(encoding='utf-8'))
        assert record['result']=='passed' and len(record['actual_'+typ.upper()])==count
        assert all(r['exact_match'] for r in record['actual_'+typ.upper()])
        assert record['source_sha256']==sha256((ROOT/script).read_bytes()).hexdigest()
        assert record['input_sha256']==sha256((ROOT/source).read_bytes()).hexdigest()
        bound.extend((name,script))
    two_name='classical-two-head-generator-audit.json'
    two_script='experiments/audit_classical_two_head_generators.py'
    two=json.loads((ROOT/two_name).read_text(encoding='utf-8'))
    assert two['result']=='passed' and len(two['cases'])==15
    assert two['source_sha256']==sha256((ROOT/two_script).read_bytes()).hexdigest()
    assert sum(c['monic_S_pairs'] for c in two['cases'])==283
    bound.extend((two_name,two_script,'CLASSICAL_TWO_HEAD_GENERATORS.md'))
    bound.extend(('CLASSICAL_ADAPTED_SMITH_BASIS.md','TYPE_D_EXPLICIT_INTEGRAL_BASIS.md',
        'CLASSICAL_SATURATION_GENERATORS.md','TYPE_B_GENERAL_SMITH.md',
        'TYPE_B_FILTRATION.md','TYPE_B_DEFECT_COMPARISON.md','TYPE_D_SQUARE_THEOREM.md',
        'TYPE_D_EVEN_GLUING_THEOREM.md','TYPE_D_HEAD_EULER_B_COMPARISON.md',
        'TYPE_D_FROBENIUS.md','experiments/verify_type_b_general_smith.py',
        'experiments/verify_type_b_filtration.py','experiments/type_b_multiplicity_search_results.json',
        'experiments/verify_type_d_candidate.py','experiments/verify_type_d_square_comparison.py',
        'experiments/type_d_cartesian_frontier.py',Path(__file__).relative_to(ROOT).as_posix()))
    result=dict(result='passed',date='2026-10-29',B_actual_embeddings=96,
        D_constructed_components=124,D_actual_embeddings=28,D_independent_product_pairs=7,
        independent_Chevalley_B_components=26,independent_Chevalley_D_components=28,
        two_head_symbolic_parameters=15,two_head_integer_S_pairs=283,
        two_head_all_rank_presentation_proved=True,
        finite_kernel_ideal_pairs=9,arbitrary_rank_kernel_presentation_proved=False,
        scope='Exact finite arithmetic and source consistency. Uniform basis/denominator proofs are in the notes; intermediate kernel-generator equality remains conjectural.',
        source_sha256={name:digest(ROOT/name) for name in sorted(set(bound))})
    (ROOT/'classical-adapted-basis-release.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'},indent=2))

if __name__=='__main__':main()
