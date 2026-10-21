"""Bind the degree derivation, finite corroboration, tests and installed wheel."""
from hashlib import sha256
from pathlib import Path
import json
import zipfile
from sage.all import GF, PolynomialRing
import sage_borel
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    raw=path.read_bytes()
    return sha256(raw if path.suffix=='.gz' else raw.replace(b'\r\n',b'\n')).hexdigest()


def read(name):
    data=json.loads((ROOT/name).read_text())
    assert data['result']=='passed',name
    return data


def main():
    assert staircase_length([(2,0),(0,2)])==4
    assert staircase_length([(3,0),(1,1),(0,4)])==6
    assert staircase_length([(0,0)])==0
    version=sage_borel.__version__
    regression=read('character-verification.json')
    docs=read('documentation-verification.json')
    wheel=read('geometry-wheel-verification.json')
    geometry=read('experiments/exceptional-degree-geometry.json')
    e6_checks=[read('experiments/'+name) for name in (
        'e6-remaining-a5.json','e6-remaining-a2a2a1.json',
        'e6-remaining-a5-audit.json','e6-remaining-a2a2a1-audit.json',
        'e6-remaining-bgg-check.json','e6-remaining-structure.json',
        'e6-independent-rational-base.json')]
    assert len(e6_checks[2]['degree_rows'])==len(e6_checks[2]['source_rows'])==22
    assert len(e6_checks[3]['degree_rows'])==len(e6_checks[3]['source_rows'])==30
    assert len(e6_checks[4]['cases'])==20
    e7_checks=[read('experiments/'+name) for name in (
        'e7-two-primary-basis.json','e7-two-primary-basis-audit.json',
        'e7-two-primary-schubert-check.json','e7-two-primary-ring-check.json')]
    assert len(e7_checks[0]['degrees'])==len(e7_checks[1]['degrees'])==28
    assert e7_checks[2]['fresh_bgg_generators']==7 and e7_checks[2]['total_basis_classes']==56
    assert e7_checks[3]['ring_identities']==6 and e7_checks[3]['total_relative_factors_2']==28
    coefficients=e7_checks[3]['multiplicity_coefficients']
    assert [r['multiplicity_2'] for r in e7_checks[0]['degrees']]==coefficients
    assert [r['multiplicity_2'] for r in e7_checks[1]['degrees']]==coefficients
    assert sum(coefficients)==28
    assert all(r['multiplicity_4']==r['multiplicity_8']==0 for r in e7_checks[0]['degrees'])
    assert sum(r['source_certificates']['M']['rank'] for r in e7_checks[0]['degrees'])==1967
    assert sum(r['source_certificates']['N']['rank'] for r in e7_checks[0]['degrees'])==1967
    e7_refinements=[read('experiments/'+name) for name in (
        'e7-refinement-relations.json','e7-refinement-independent.json',
        'e7-refinement-bgg-check.json','e7-refinement-invariant-obstruction.json',
        'e7-refinement-invariant-obstruction-audit.json')]
    assert e7_refinements[1]['degree_bound']==9 and len(e7_refinements[1]['refinements'])==64
    assert [r['columns'] for r in e7_refinements[0]['kernel_primitivity']]==[1,7,28,84]
    assert all(r['columns']==r['primitive_mod2_rank'] for r in e7_refinements[0]['kernel_primitivity'])
    assert e7_refinements[2]['degree_interval']==[0,9] and len(e7_refinements[2]['degrees'])==10
    row=next(r for r in e7_refinements[1]['refinements'] if r['retained_nodes']==[1,3,4,5,6])
    assert row['multiplicity_2']==[r['multiplicity_2'] for r in e7_refinements[2]['degrees']]
    assert e7_refinements[3]['chern_invariant_value']==e7_refinements[4]['chern_invariant_value']==-696735
    assert e7_refinements[3]['basic_values']==e7_refinements[4]['basic_values']==[2,-4,88,954,-340,55780]
    e7_full=[read('experiments/'+name) for name in (
        'e7-a6-images.json','e7-a6-ideal-audit.json','e7-a6-images-audit.json',
        'e7-full-flag-kernel.json','e7-full-flag-kernel-audit.json','e7-full-flag-transfer.json')]
    assert e7_full[0]['dimension']==e7_full[2]['dimension']==42
    assert e7_full[0]['total_rank']==e7_full[2]['total_rank']==576
    assert len(e7_full[2]['degree_rows'])==42 and e7_full[2]['buchberger_pairs_checked']==9731
    assert e7_full[1]['primary_ideal_equals_adjoint_ideal']
    assert [r['degree'] for r in e7_full[3]['degrees']]==[10,11,12,13]
    assert [r['multiplicity_2'] for r in e7_full[3]['degrees']]==[1940,2855,4054,5578]
    assert [r['multiplicity_2'] for r in e7_full[4]['degrees']]==[1940,2855,4054,5578]
    assert e7_full[5]['proper_parabolics_settled']==109 and e7_full[5]['P7_refinements_settled']==61
    assert len(e7_full[5]['good_levi_parabolics'])==108
    assert sum(e7_full[5]['full_flag_multiplicity_coefficients'])==1451520
    e7_d5=[read('experiments/'+name) for name in (
        'e7-d5-source-algebra.json','e7-d5-images.json','e7-d5-source-audit.json',
        'e7-d5-moments-audit.json','e7-d5-images-audit.json','e7-d5-kernel.json','e7-d5-summary.json')]
    bounded_source=json.loads((ROOT/'experiments/e7-d5-source.json').read_text())
    assert bounded_source['source_sha256']==digest(ROOT/'experiments/e7_d5_source.py')
    assert e7_d5[0]['generator_degrees']==[1,2,4,5,6,8,8,9,12,16]
    assert e7_d5[1]['prediction_agrees'] and e7_d5[1]['dimension']==43 and e7_d5[1]['total_rank']==1512
    assert sum(e7_d5[1]['multiplicity_coefficients'])==756
    assert e7_d5[2]['full_source_complete_all_degrees'] and e7_d5[2]['integral_invariant_generators']==10
    assert e7_d5[3]['basic_invariant_algebra_complete'] and e7_d5[3]['weight_orbit']==56
    assert len(e7_d5[4]['degree_rows'])==43 and e7_d5[4]['total_rank']==1512
    assert [0]+[r['multiplicity_2'] for r in e7_d5[4]['degree_rows']]==e7_d5[1]['multiplicity_coefficients']
    assert all(r['multiplicity_2']==e7_d5[1]['multiplicity_coefficients'][r['degree']] for r in e7_d5[5]['degrees'])
    assert e7_d5[6]['proper_parabolics_settled']==111 and e7_d5[6]['P7_refinements_settled']==63
    for check in e6_checks+e7_checks+e7_refinements+e7_full+e7_d5:
        for name,expected in check.get('source_sha256',{}).items():
            path=ROOT/name
            if not path.exists():path=ROOT/'experiments'/name
            assert digest(path)==expected,name
        for name,expected in check.get('raw_sha256',{}).items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
        if 'certificate_archive' in check:
            assert sha256((ROOT/check['certificate_archive']).read_bytes()).hexdigest()==check['certificate_archive_sha256']
    assert regression['version']==wheel['version']==version=='0.20.0'
    for name,expected in regression['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    wheel_path=ROOT/'dist'/wheel['wheel']
    assert sha256(wheel_path.read_bytes()).hexdigest()==wheel['wheel_sha256']
    runtime=list((ROOT/'src/sage_borel').glob('*.py'))
    runtime+=list((ROOT/'src/sage_borel/data').glob('*.json'))
    with zipfile.ZipFile(wheel_path) as archive:
        for path in runtime:
            name=path.relative_to(ROOT/'src').as_posix()
            assert archive.read(name).replace(b'\r\n',b'\n')==path.read_bytes().replace(b'\r\n',b'\n'),name
    degrees={}
    for entry in geometry['degree_deductions']:
        key=(entry['cartan_type'],entry['representation'])
        c=sage_borel.exceptional_springer_degree(key[0],representation=key[1])
        assert (c.normalized_volume,c.boundary_contribution,c.central_cover,c.degree)==(
            entry['normalized_volume'],entry['boundary_correction'],entry['central_cover'],entry['degree'])
        degrees[key]=c.degree
    fibers=[]
    fiber_paths=sorted((ROOT/'experiments').glob('e[67]-*-fiber-*.json'))
    for path in fiber_paths:
        f=json.loads(path.read_text())
        key=(f['cartan_type'],f['representation'])
        assert f['status']=='specialized-fiber-only'
        assert f['fiber_dimension']==0 and f['fiber_length']==degrees[key],path
        # Recompute the length from the recorded leading ideal. This checks
        # the report's staircase count, not the Groebner-basis derivation.
        n=(f['base_rank'] if 'module_leading_terms' in f else
           len(f['leading_exponents'][0]) if 'leading_exponents' in f else len(f['targets']))
        R=PolynomialRing(GF(f['characteristic']),n,names=['x%d'%(i+1) for i in range(n)])
        if 'module_leading_terms' in f:
            ideals=[[] for _ in range(f['module_rank'])]
            for row in f['module_leading_terms']:
                ideals[row[0]-1].append(R({tuple(row[1:]):1}))
            components=[int(R.ideal(gens).vector_space_dimension()) for gens in ideals]
            assert components==f['module_component_lengths']
            for gens,component_length in zip(ideals,components):
                assert staircase_length([tuple(g.exponents()[0]) for g in gens])==component_length
            assert f['preflight_lengths']==[1,4]
            if 'tower_reduction' in f:
                assert f['rank_four_preflight_length']==4
                assert f['tower_reduction']['first_target_matrix_identity']
                assert f['tower_reduction']['second_target_matrix_identity']
                assert f['tower_reduction']['commuting_regular_representation']
            length=sum(components)
        else:
            monomials=([R({tuple(e):1}) for e in f['leading_exponents']]
                       if 'leading_exponents' in f else [R(m) for m in f['leading_monomials']])
            length=int(R.ideal(monomials).vector_space_dimension())
            assert staircase_length([tuple(g.exponents()[0]) for g in monomials])==length
        cover=(f.get('center_chart') or {}).get('cover_degree',1)
        if f.get('center_chart'):
            assert f['center_chart']['excluded_original_divisor_empty']
            assert f['center_chart']['excluded_quotient_divisor_empty']
        assert cover*length==f['fiber_length'],path
        fibers.append(dict(file=path.relative_to(ROOT).as_posix(),cartan_type=key[0],
            representation=key[1],prime=f['characteristic'],seed=f['seed'],
            length=f['fiber_length'],engine=f['engine']))
    paths=fiber_paths+runtime+[ROOT/name for name in (
        'EXCEPTIONAL_DEGREES.md','experiments/exceptional_degree_geometry.py',
        'experiments/exceptional-degree-geometry.json','experiments/exceptional_degree_fibers.py',
        'experiments/exceptional-degree-bounded-attempts.json',
        'experiments/e7_adjoint_module_fiber.py',
        'experiments/e7_adjoint_shifted_fiber.py',
        'experiments/e7_adjoint_tower_fiber.py',
        'experiments/e6_matrix_oracle.py','experiments/e7_matrix_oracle.py',
        'README.md','../README.md','E6_GEOMETRY.md','E7_GEOMETRY.md',
        'RESEARCH_LEDGER.md','AGENTS.md',
        'LOCAL_IMAGE_BOUND.md','EQUIVARIANT_IMAGE_PILOT.md','LOCAL_IMAGE_AUDIT.md',
        'E6_ROOT_WEIGHT_COMPARISON.md','E6_ROOT_WEIGHT_AUDIT.md',
        'E6_REMAINING_CASES.md','E6_REMAINING_AUDIT.md','E6_REMAINING_TABLES.md',
        'E7_ROOT_WEIGHT_COMPARISON.md','E7_ROOT_WEIGHT_AUDIT.md',
        'E7_REFINEMENT_MULTIPLICITIES.md',
        'E7_FULL_FLAG_COMPARISON.md','E7_D5_REFINEMENT.md',
        'experiments/e7_d5_source.py','experiments/e7_d5_source_algebra.py','experiments/e7_d5_images.py',
        'experiments/e7_d5_kernel.py','experiments/e7_d5_summary_check.py','experiments/e7-d5-summary.json','experiments/audit_e7_d5_source.py','experiments/audit_e7_d5_images.py',
        'experiments/audit_e7_d5_moments.py','experiments/e7-d5-source.json',
        'experiments/e7-d5-source-algebra.json','experiments/e7-d5-images.json','experiments/e7-d5-kernel.json',
        'experiments/e7-d5-source-audit.json','experiments/e7-d5-images-audit.json',
        'experiments/e7-d5-moments-audit.json','experiments/e7-d5-images-certificates.json.gz',
        'experiments/e7-d5-kernel-certificates.json.gz',
        'experiments/e7_full_flag_kernel.py','experiments/audit_e7_full_flag_kernel.py',
        'experiments/e7_a6_images.py','experiments/audit_e7_a6_ideal.py','experiments/audit_e7_a6_images.py',
        'experiments/e7_full_flag_transfer_check.py',
        'experiments/e7-full-flag-kernel.json','experiments/e7-full-flag-kernel-audit.json',
        'experiments/e7-a6-images.json','experiments/e7-a6-ideal-audit.json','experiments/e7-a6-images-audit.json',
        'experiments/e7-full-flag-transfer.json','experiments/e7-a6-images-certificates.json.gz',
        'experiments/e7-full-flag-kernel-certificates.json.gz',
        'experiments/e7_refinement_relations.py','experiments/audit_e7_refinement_relations.py',
        'experiments/e7_refinement_bgg_check.py','experiments/e7_refinement_invariant_obstruction.py',
        'experiments/audit_e7_refinement_obstruction.py',
        'experiments/e7-refinement-relations.json','experiments/e7-refinement-independent.json',
        'experiments/e7-refinement-bgg-check.json','experiments/e7-refinement-invariant-obstruction.json',
        'experiments/e7-refinement-invariant-obstruction-audit.json',
        'E7_REPRESENTATION_COMPARISON.md',
        'experiments/e7_two_primary_basis.py','experiments/audit_e7_two_primary_basis.py',
        'experiments/e7_two_primary_schubert_check.py','experiments/e7_two_primary_ring_check.py',
        'experiments/e7-two-primary-basis.json','experiments/e7-two-primary-basis-audit.json',
        'experiments/e7-two-primary-schubert-check.json','experiments/e7-two-primary-ring-check.json',
        'experiments/e7-character-form-results.json','experiments/e7-adjoint-comparison-results.json',
        'tests/test_exceptional.py',
        'experiments/e6_remaining_images.py','experiments/audit_e6_remaining_images.py',
        'experiments/build_e6_independent_base.py','experiments/e6_remaining_structure_check.py',
        'experiments/e6_remaining_bgg_check.py','experiments/package_e6_remaining.py',
        'experiments/summarize_e6_remaining.py',
        'experiments/e6-independent-rational-base.json','experiments/e6-remaining-structure.json',
        'experiments/e6-remaining-a5.json','experiments/e6-remaining-a2a2a1.json',
        'experiments/e6-remaining-a2a2a1-certificates.json.gz',
        'experiments/e6-remaining-a5-audit.json','experiments/e6-remaining-a2a2a1-audit.json',
        'experiments/e6-remaining-bgg-check.json',
        'E6_REPRESENTATION_COMPARISON.md',
        'experiments/e6_root_weight_comparison.py',
        'experiments/verify_e6_root_weight_independent.py',
        'experiments/e6_root_weight_refinement_check.py',
        'experiments/e6-root-weight-comparison.json',
        'experiments/e6-root-weight-independent.json',
        'experiments/e6-root-weight-refinements.json',
        'experiments/local_image_equivariant_pilot.py',
        'experiments/verify_local_image_independent.py',
        'experiments/local-image-equivariant-pilot.json',
        'experiments/local-image-independent.json',
        '../articles/research-program-proposal/main.tex',
        '../articles/research-program-proposal/README.md',
        'PROOF_STATUS.md','EXTENSIONS.md','SPRINGER_GEOMETRY.md',
        'tools/verify_documentation.py','tools/verify_geometry_wheel.py',
        'tools/monomial_staircase.py',
        'tests/test_exceptional_degrees.py','character-verification.json',
        'documentation-verification.json','geometry-wheel-verification.json')]+[Path(__file__).resolve()]
    report=dict(result='passed',version=version,
        characteristic_zero_proof='EXCEPTIONAL_DEGREES.md; toric localized top Chern and explicit boundary Segre classes',
        degree_deductions=geometry['degree_deductions'],specialized_fibers=fibers,
        full_fiber_cases_not_yet_corroborated=[list(key) for key in degrees
            if key not in {(f['cartan_type'],f['representation']) for f in fibers}],
        finite_check_scope='Specialized lengths corroborate the separate geometric proof; no inference of characteristic-zero generic degree from specialization.',
        standard_monomial_counts_rechecked='Singular vdim and independent integer exponent-slice recursion',
        regression_tests=regression['test_methods_run'],documentation_examples=docs['python_blocks_executed'],
        local_links_and_anchors=docs['local_links_and_anchors_checked'],isolated_wheel_tests=wheel['tests'],
        wheel_sha256=wheel['wheel_sha256'],wheel_sources_and_data_match=True,
        e6_root_weight_completion=dict(proper_parabolics=63,previous_transfer_family=60,
            new_base_degrees_independently_audited=52,direct_bgg_degree_cases=20,
            third_case='Proved common (1,z) basis; multiplicity polynomials multiplied by 1+q.',
            scope='Relative root/weight image factors; ambient H/N and literature priority are separate.'),
        e7_root_weight_basis=dict(parabolic='E7/P7',basis_elements=56,full_sources=2,all_degrees=28,
            relative_factors='Only 2; polynomial q[14]_(q2)+(q6+q10)[7]_(q2).',
            source_monomials=1967,refinement_exponent=2,refinement_count=64,
            scope='Full P7 multiplicities; exponent bound only on the other refinements. Shared Sage arithmetic in the independent source replay; separate BGG/Chevalley and plain Python structural checks.'),
        e7_refinement_progress=dict(refinements=64,multiplicity_degree_bound=9,
            independent_relation_constructions=2,direct_bgg_degrees=10,all_degree_unknown_base_polynomials=1,all_degree_refinements_settled=63,
            source_obstruction_degree=16,scope='Uniform degree-nine formula; A6/P7 and D5 source proofs settle 63/64 refinements in all degrees. D4 remains; no novelty claim.'),
        e7_full_flag_comparison=dict(formula='q P_W/(1+q)',intermediary='E7/P2 with A6 Levi',
            intermediary_rank=576,intermediary_dimension=42,intermediary_all_degrees_audited=True,
            full_flag_total_factors_2=1451520,proper_good_levi_parabolics=108,additional_P7_basis_case=1,
            proper_parabolics_settled=109,P7_refinements_settled=61,
            scope='Complete computer-assisted A6 image proof and classical common-fiber transfer. A6-transfer scope only; the D5 pair is completed separately below, and D4 remains. No novelty or ambient H/N classification claim.'),
        e7_d5_refinement=dict(source_generators=10,source_complete_all_degrees=True,quotient_rank=1512,dimension=43,
            independent_image_degrees=43,total_relative_factors_2=756,proper_parabolics_settled=111,P7_refinements_settled=63,
            scope='Complete computer-assisted D5 full-source comparison, with separate native-Python source, weight and image audits. E6 opposition gives the second D5 case; D4 remains.'),
        hash_semantics='SHA256 after CRLF to LF normalization for text; raw bytes for .gz archives',
        audit_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(set(paths))})
    (ROOT/'release-0.20-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='audit_sha256'},indent=2))


if __name__=='__main__':
    main()
