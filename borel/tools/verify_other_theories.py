"""Validate the exploration reports, their source bindings and symbolic SL3 blocks.

This checker does not rerun every experiment or certify mathematical
originality. Each report names its own independent method and finite scope.
"""
from hashlib import sha256
import json
from pathlib import Path
import re
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
REPORTS=(
    'other-theories-k-pilot.json','other-theories-k-audit.json',
    'other-theories-sl3-classification.json','other-theories-cotangent-pilot.json',
    'other-theories-height-pilot.json','other-theories-symplectic-k-pilot.json',
    'other-theories-equivariant-k-pilot.json','other-theories-symplectic-equivariant-pilot.json',
    'other-theories-affine-pilot.json','other-theories-affine-c-audit.json',
    'other-theories-quantum-pilot.json','other-theories-folding-pilot.json',
    'other-theories-scalar-check.json')
NOTES=tuple('OTHER_THEORIES_'+name+'.md' for name in (
    'PROGRAM','FINITE','SL3','TYPE_A','HEIGHT','COTANGENT','EQUIVARIANT_K','QUANTUM','AFFINE','FOLDING'))


def digest(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def main():
    reports={};bound=[]
    for name in REPORTS:
        path=ROOT/'experiments'/name;record=json.loads(path.read_text())
        assert record['result']=='passed',name
        for key in ('source_sha256','raw_sha256'):
            for relative,expected in record.get(key,{}).items():
                assert digest(ROOT/relative)==expected,(name,relative)
                bound.append(relative)
        reports[name]=record
    k=reports['other-theories-k-audit.json'];sl3=reports['other-theories-sl3-classification.json']
    assert len(k['cases'])==6 and k['G2_Davis_primitive_comparison']['verified']
    assert k['G2_adjoint_not_contained_in_minimal']['nonintegral_coordinate']=='53/576'
    assert sl3['representation_cases']==41 and sl3['proper_flag_images_checked']==123
    h=reports['other-theories-height-pilot.json']
    assert len(h['cases'])==5 and sum(len(c['checks']) for c in h['cases'])==320
    assert len(h['independent_height_two_checks'])==64
    q=reports['other-theories-quantum-pilot.json']
    assert q['independent_degree_gcd_cases']==25991 and q['multiplication_integrality_checks']==173879
    a=reports['other-theories-affine-pilot.json']
    assert a['multinomial_gcd_vs_prime_power_cases']==6231
    assert a['type_C_weighted_gcd_vs_prime_power_cases']==4824
    ca=reports['other-theories-affine-c-audit.json']
    assert [(c['rank'],c['degree_through']) for c in ca['cases']]==[(2,20),(3,12)]
    eq=reports['other-theories-equivariant-k-pilot.json']
    assert len(eq['rank_two_determinants'])==3
    assert sum(len(c['points']) for c in eq['rank_two_determinants'])==9
    assert [r['cosets'] for r in eq['independent_coset_orbit_counts'][-3:]]==[27,56,240]
    assert len(reports['other-theories-symplectic-equivariant-pilot.json']['cases'])==3
    sk=reports['other-theories-symplectic-k-pilot.json']
    assert len(sk['cases'])==63 and sk['inverse_coefficient_denominator_checks']==1000
    folding=reports['other-theories-folding-pilot.json']
    assert [r['quotient_of_full_fixed_source_by_restricted_source'] for r in folding['cases']]==[[1,2],[1,1,2,2],[1,3]]
    scalar=reports['other-theories-scalar-check.json']
    assert scalar['global_Laurent_column_scaling_verified'] and scalar['reference_equivariant_matrix_is_unit']
    assert scalar['diagonal_character_exponents']==[0,1,1,2,2,3]
    # Symbolic all-parameter integral equivalences behind the SL3 theorem.
    s,c=sp.symbols('s c');C=sp.Matrix([[1,1],[-1,2]])
    E=sp.Matrix([[-2,1],[-1,2]]);D=sp.Matrix([[0,1],[-1,1]])
    M=(s*C).row_join(sp.zeros(2)).col_join((c*E).row_join(3*s*s*D))
    left=sp.diag(sp.Matrix([[1,0],[1,1]]),D.inv())
    right=sp.diag(sp.Matrix([[1,-1],[0,1]]),sp.eye(2))
    transformed=left*M*right
    transformed[3,:]=transformed[3,:]-2*transformed[2,:]
    transformed[:,2]=transformed[:,2]+2*transformed[:,3]
    permutation=[0,2,1,3];transformed=transformed[permutation,permutation]
    expected=sp.diag(sp.Matrix([[s,0],[-c,3*s*s]]),3*sp.Matrix([[s,0],[c,s*s]]))
    assert transformed==expected
    B=sp.Matrix([[-1,-1],[0,1]])
    P=(s*sp.eye(2)).row_join(sp.zeros(2)).col_join((c*B).row_join(s*s*sp.eye(2)))
    P=sp.diag(sp.eye(2),B.inv())*P*sp.diag(sp.eye(2),B)
    assert P[permutation,permutation]==sp.diag(sp.Matrix([[s,0],[c,s*s]]),sp.Matrix([[s,0],[c,s*s]]))
    u=sp.Symbol('u',nonzero=True)
    lower=sp.Matrix([[1,0],[1,1]]);upper=sp.diag(u,1/u)
    theta=lambda M:M-sp.trace(M)*sp.eye(2)/2
    assert sp.simplify((theta(lower*upper)-theta(lower))[1,0])==u-1
    ledger=(ROOT/'RESEARCH_LEDGER.md').read_text(encoding='utf-8')
    ids=re.findall(r'^\| ((?:BASE|SRC|IMG|GEO|CASE|LIMIT|PART|API|OPEN)-[0-9]+) \|',ledger,re.M)
    assert len(ids)==len(set(ids))==178
    assert sum(i.startswith('OPEN') for i in ids)==22
    paths=[ROOT/name for name in NOTES]+[ROOT/'experiments'/name for name in REPORTS]
    paths += [ROOT/name for name in set(bound)]+[Path(__file__)]
    report=dict(result='passed',date='2026-10-23',reports_checked=len(REPORTS),
                proof_notes_and_program=len(NOTES),source_bindings_checked=len(bound),
                symbolic_all_parameter_SL3_integral_block_reductions=True,
                constant_loop_flag_descent_counterexample_verified=True,
                ledger_entries=178,open_questions=22,
                scope='Report consistency, normalized source/raw hashes and exact symbolic SL3 integral reductions; not an external proof review, priority certificate or new runtime release.',
                source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(set(paths))})
    (ROOT/'other-theories-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
