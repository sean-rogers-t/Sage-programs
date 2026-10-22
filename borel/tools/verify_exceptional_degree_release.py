"""Bind the degree derivation, finite corroboration, tests and installed wheel."""
from hashlib import sha256
from pathlib import Path
import json
import re
import zipfile
from sage.all import GF, PolynomialRing
import sage_borel
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    raw=path.read_bytes()
    return sha256(raw if path.suffix in ('.gz','.pdf') else raw.replace(b'\r\n',b'\n')).hexdigest()


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
    other=read('other-theories-verification.json')
    classical_equivariant=read('equivariant-verification.json')
    assert classical_equivariant['full_regression_methods']==regression['test_methods_run']
    assert classical_equivariant['symbolic_rank_three_unimodular_reduction']
    assert classical_equivariant['symbolic_doubling_unimodular_certificates']==list(range(1,9))
    assert classical_equivariant['symbolic_rank_four_minor_ideal_certificates']
    assert classical_equivariant['initial_column_minor_theorem']['symbolic_identities']==94
    assert classical_equivariant['initial_column_minor_theorem']['original_inclusion_checks']==72
    assert classical_equivariant['hankel_balancing']['universal_quartic_integer_identities']==70
    assert classical_equivariant['hankel_balancing']['original_inclusion_checks']==68
    assert classical_equivariant['type_A']['fixed_points']==10
    assert classical_equivariant['type_C']['index_exponent']==153
    assert classical_equivariant['type_C']['symbolic_matrix_avoided']
    for name,expected in classical_equivariant['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    assert other['reports_checked']==13 and other['proof_notes_and_program']==10
    assert other['symbolic_all_parameter_SL3_integral_block_reductions']
    for name,expected in other['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    other_notes={name for name in other['source_sha256'] if name.startswith('OTHER_THEORIES_') and name.endswith('.md')}
    assert other_notes <= {r['file'] for r in docs['documents']}
    exposition_prefix='../articles/classical-equivariant-exposition/'
    exposition=read(exposition_prefix+'verification.json')
    exposition_examples=read(exposition_prefix+'examples-verification.json')
    assert exposition['pages']==12 and exposition['numbered_equations']==43
    assert exposition['page_layout']=='SinglePage' and exposition['manuscript_pages_per_pdf_page']==1
    assert exposition['bibliography_items']==13 and exposition['characters_outside_page']==0
    assert exposition['overfull_boxes_or_undefined_references']==0
    assert len(exposition_examples['cases'])==6
    assert digest(ROOT/(exposition_prefix+'main.tex'))==exposition['source_sha256']==exposition_examples['tex_sha256']
    assert digest(ROOT/(exposition_prefix+'verify_examples.py'))==exposition_examples['source_sha256']
    assert digest(ROOT/(exposition_prefix+'output/pdf/Classical_Springer_Borel_equivariant_exposition.pdf'))==exposition['pdf_sha256']
    assert exposition_prefix+'README.md' in {r['file'] for r in docs['documents']}
    assert {'E6_RELATIVE_IMAGE_API.md','G2_FULL_FLAG.md','E7_QUADRATIC_DERIVATION.md','E7_DEGREE11_LIFTING.md','CENTRAL_PARITY_KERNEL.md','SPIN_PARITY_EXTENSION.md','SPIN_PARITY_RECURSION.md','SPIN_PARITY_BOUNDARY.md'} <= {r['file'] for r in docs['documents']}
    wheel=read('geometry-wheel-verification.json')
    relative_api=read('relative-image-verification.json')
    assert relative_api['version']==version
    assert relative_api['proper_standard_parabolics']==190
    assert relative_api['E7']['cohomology_polynomials_reconstructed']==127
    assert relative_api['E7']['multiplicity_polynomials_reconstructed']==125
    assert relative_api['E7']['multiplicity_polynomials_from_independent_image_audits']==2
    assert relative_api['E7']['fresh_character_degrees']==[0,1,2,3,6,10]
    assert relative_api['E6']['proper_standard_parabolics']==63
    assert relative_api['E6']['nonunit_factors']==[3,9,27]
    assert relative_api['E6']['fresh_a5_full_integer_invariant_kernel_degrees']==[3]
    for name,expected in relative_api['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    e6_api_summary=read('experiments/e6-all-parabolics-summary.json')
    assert e6_api_summary['complete'] and len(e6_api_summary['parabolics'])==63
    g2_full=read('experiments/g2-full-flag.json')
    g2_full_audit=read('experiments/g2-full-flag-audit.json')
    assert g2_full['complete'] and g2_full_audit['complete']
    assert g2_full_audit['checked_monomial_columns']==72
    assert g2_full['cases'][0]['degrees'][6]['smith_factors']==[3456]
    for record in (e6_api_summary,g2_full,g2_full_audit):
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    structural=[read('experiments/'+name) for name in (
        'e7-quadratic-derivation.json','e7-quadratic-derivation-audit.json',
        'type-a-quadratic-derivation-audit.json')]
    assert structural[0]['maximum_degree']==structural[1]['maximum_degree']==13
    assert structural[1]['uniform_family_maximum_degree']==10
    assert structural[1]['uniform_family_size']==112
    assert structural[1]['primitive_relation_odd_valuations']=={'6':2,'8':1,'10':2,'12':2}
    for report0 in structural[:2]:
        assert [r['degree'] for r in report0['degrees']]==list(range(14))
        assert all(r['derivative_preserves_full_kernel'] for r in report0['degrees'])
        assert report0['degrees'][13]['kernel_rank']==1716
        assert report0['degrees'][13]['even_projection_rank']==1106
    assert [c['rank'] for c in structural[2]['cases']]==[6,20,10,120,252]
    for record in structural:
        for key in ('source_sha256','raw_sha256'):
            for name,expected in record.get(key,{}).items():
                assert digest(ROOT/name)==expected,name
    assert digest(ROOT/structural[0]['certificate_archive'])==structural[0]['certificate_archive_sha256']
    lifting11=[read('experiments/'+name) for name in (
        'e7-degree11-lifting.json','e7-degree11-lifting-audit.json')]
    for record in lifting11:
        assert record['degree']==11 and record['uniform_family_size']==112
        assert record['initial_kernel_rank_mod2']==461
        assert record['primitive_kernel_rank_mod2']==462
        assert record['mod2_syzygy_dimension']==1 and record['good_Levi_count']==108
        assert record['quadratic_Frobenius_identities_verified']
        assert record['half_sum_derivative_identity_verified']
        assert len(record['exceptional_noncubic_lifting_cases'])==4
        assert all(c['cubic_does_not_lift'] for c in record['exceptional_noncubic_lifting_cases'])
        for key in ('source_sha256','raw_sha256'):
            for name,expected in record.get(key,{}).items():
                assert digest(ROOT/name)==expected,name
    assert lifting11[0]['saturation_quotient_order']==2
    assert digest(ROOT/lifting11[0]['certificate_archive'])==lifting11[0]['certificate_archive_sha256']
    central=[read('experiments/'+name) for name in (
        'central-parity-kernel.json','central-parity-kernel-audit.json')]
    assert central[0]['Spin7']['degrees'][4]['actual_integral_kernel_preserved'] is False
    assert central[0]['Spin7']['degrees'][3]['actual_integral_kernel_preserved'] is True
    assert central[1]['Spin7']['integral_kernel_preservation'] is False
    assert central[1]['Spin7']['root_lattice_index']==2
    assert central[1]['Spin7']['monic_division_remainder']=='v^3'
    assert len(central[0]['type_C_checks'])==3
    assert all(row.get('actual_integral_kernel_preserved',True)
               and row.get('actual_source_preserved',True)
               for case in central[0]['type_C_checks'] for row in case['degrees'])
    for record in central:
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    spin=[read('experiments/'+name) for name in (
        'spin-parity-extension.json','spin-parity-extension-audit.json')]
    for record in spin:
        assert [c['rank'] for c in record['cases']]==list(range(3,9))
        assert [c['entire_kernel_preserved'] for c in record['cases']]==[False,True,True,True,False,False]
        assert [c['first_failure_degree'] for c in record['cases']]==[4,None,None,None,11,16]
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    assert spin[1]['Spin9_closed_presentation_checked']
    assert all(c['finite_flatness_certificate'] and c['reduced_quotient_length']==2**c['rank']
               for c in spin[1]['cases'])
    spin_recursive=[read('experiments/'+name) for name in (
        'spin-parity-recursion.json','spin-parity-recursion-audit.json',
        'spin-parity-quadratic-check.json')]
    for record in spin_recursive[:2]:
        assert [c['rank'] for c in record['cases']]==list(range(1,13))
        assert [c['first_failure_degree'] for c in record['cases']]==[None,None,4,None,None,None,11,16,None,None,None,34]
        assert record['arbitrary_rank_closed_classification_complete'] is False
    assert spin_recursive[1]['lowering_proof']['symbolic_commutator_coefficients']
    assert spin_recursive[1]['lowering_proof']['strict_partitions_checked']==1022
    assert spin_recursive[1]['lowering_proof']['infinite_quadratic_derivation_identities_checked']==9
    assert [c['rank'] for c in spin_recursive[2]['cases']]==[3,4,7,8,11,12]
    assert [c['first_failure_degree'] for c in spin_recursive[2]['cases']]==[4,None,11,16,None,34]
    for record in spin_recursive:
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    assert digest(ROOT/spin_recursive[0]['certificate_archive'])==spin_recursive[0]['certificate_archive_sha256']
    spin_boundary=[read('experiments/'+name) for name in (
        'spin-parity-boundary.json','spin-parity-boundary-audit.json',
        'spin-parity-boundary-rank15.json','spin-parity-boundary-rank15-audit.json',
        'spin-parity-even-source-check.json')]
    for record in spin_boundary[:2]:
        assert [c['rank'] for c in record['cases']]==[3,4,7,8,11,12]
        assert [c['first_failure_degree'] for c in record['cases']]==[4,None,11,16,None,34]
    for record in spin_boundary[2:4]:
        assert len(record['cases'])==1
        assert record['cases'][0]['rank']==15
        assert record['cases'][0]['first_failure_degree']==37
        assert record['cases'][0]['tested_through_degree']==37
    for record in spin_boundary:
        assert record['closed_numerical_rank_classification_proved'] is False
        for key in ('source_sha256','raw_sha256'):
            for name,expected in record.get(key,{}).items():
                assert digest(ROOT/name)==expected,name
    for record in (spin_boundary[0],spin_boundary[2]):
        assert digest(ROOT/record['certificate_archive'])==record['certificate_archive_sha256']
    for record in (spin_boundary[1],spin_boundary[3]):
        assert record['independent_derivative_factorization']
        assert record['schur_P_pfaffian_lowering_checks']==94
    assert [c['rank'] for c in spin_boundary[4]['cases']]==[3,4,7,8,11,12,15]
    for c in spin_boundary[4]['cases']:
        if c['rank']==15:
            assert not c['complete_index_polynomial'] and c['first_failure_degree']==37
        else:
            assert c['complete_index_polynomial']
            assert sum(r['obstruction_rank'] for r in c['degree_rows'])==2**c['rank']-2*c['total_index_exponent']
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
    e7_d4=[read('experiments/'+name) for name in (
        'e7-d4-source-algebra.json','e7-d4-images.json','e7-d4-source-audit.json',
        'e7-d4-images-audit.json','e7-d4-summary.json')]
    d4_bounded=json.loads((ROOT/'experiments/e7-d4-source.json').read_text())
    assert d4_bounded['source_sha256']==digest(ROOT/'experiments/e7_d4_source.py')
    assert e7_d4[0]['generator_degrees']==[1,1,2,4,4,5,6,8,8,8]
    assert e7_d4[0]['relative_quotient_dimension']==e7_d4[2]['relative_module_generators']==34
    assert e7_d4[2]['full_source_complete_all_degrees'] and e7_d4[1]['complete']
    assert e7_d4[1]['dimension']==e7_d4[3]['dimension']==51
    assert e7_d4[1]['total_rank']==e7_d4[3]['total_rank']==15120
    assert len(e7_d4[1]['degrees'])==len(e7_d4[3]['degree_rows'])==52
    assert e7_d4[4]['total_relative_factors_2']==7560
    assert e7_d4[4]['proper_parabolics_settled']==112 and e7_d4[4]['P7_refinements_settled']==64
    assert len(e7_d4[4]['other_proper_parabolics_remaining'])==15
    assert [r['weight_and_root_split_odd_index'] for r in e7_d4[4]['quadratic_lifts']]==[7,15,3]
    assert all(r['all_coefficients_integral'] and r['all_simple_reflections_checked'] for r in e7_d4[4]['quadratic_lifts'])
    assert [r['indices'] for r in e7_d4[1]['degrees']]==[r['indices'] for r in e7_d4[3]['degree_rows']]
    assert all(r['indices']['N']==r['indices']['M'] and r['prediction_agrees'] for r in e7_d4[1]['degrees'])
    e7_d6=[read('experiments/'+name) for name in (
        'e7-d6-sources.json','e7-d6-sources-audit.json','e7-d6-images.json',
        'e7-d6-images-audit.json','e7-d6-moments-audit.json','e7-d6-summary.json',
        'e7-remaining-conjugacy.json','e7-remaining-conjugacy-audit.json',
        'e7-d5a1-degree-two.json','e7-d5a1-degree-two-audit.json')]
    assert e7_d6[0]['source_degree_bounds']==e7_d6[1]['source_degree_bounds']==dict(N=33,M=27)
    assert len(e7_d6[1]['degrees'])==60 and e7_d6[2]['complete']
    assert e7_d6[2]['max_degree']==e7_d6[3]['max_degree']==33
    assert len(e7_d6[2]['degrees'])==len(e7_d6[3]['degree_rows'])==33
    assert e7_d6[4]['basic_invariant_algebra_complete'] and e7_d6[4]['weight_orbit']==56
    assert e7_d6[5]['proper_parabolics_settled']==113
    assert len(e7_d6[5]['other_proper_parabolics_remaining'])==14
    assert e7_d6[5]['remaining_distinct_Weyl_conjugacy_classes']==10
    assert e7_d6[5]['total_relative_factors_2']==34
    assert e7_d6[8]['relative_smith_2_valuations']==[0,1] and e7_d6[9]['first_nonzero_relative_degree']==2
    assert e7_d6[6]['distinct_representatives_after_these_reductions']==e7_d6[7]['distinct_Weyl_conjugacy_classes']==11
    assert e7_d6[5]['multiplicity_coefficients']==[0]+[r['multiplicity_2'] for r in e7_d6[3]['degree_rows']]
    assert all(r['smith_multiplicities']==({'1':m} if m else {}) for r,m in zip(e7_d6[2]['degrees'],e7_d6[5]['multiplicity_coefficients'][1:]))
    e7_d5a1=[read('experiments/'+name) for name in (
        'e7-d5a1-sources.json','e7-d5a1-sources-audit.json',
        'e7-d5a1-sources-22.json','e7-d5a1-sources-22-audit.json',
        'e7-d5a1-norm.json','e7-d5a1-norm-audit.json',
        'e7-d5a1-images.json','e7-d5a1-images-audit.json',
        'e7-d5a1-moments-audit.json','e7-d5a1-summary.json')]
    e7_d5a1_attempts=[json.loads((ROOT/'experiments'/name).read_text()) for name in (
        'e7-d5a1-source-algebra.json','e7-d5a1-free-source.json')]
    assert all(c['result']=='incomplete' and not c['full_source_complete_all_degrees'] for c in e7_d5a1_attempts)
    assert e7_d5a1_attempts[0]['first_missing_degree']==21
    assert e7_d5a1[0]['max_degree']==e7_d5a1[1]['max_degree']==20
    assert e7_d5a1[2]['max_degree']==e7_d5a1[3]['max_degree']==22
    assert [r['necessary_new_generators'] for r in e7_d5a1[3]['indecomposable_checks']]==[1,1]
    assert e7_d5a1[4]['degree']==e7_d5a1[5]['degree']==32
    assert e7_d5a1[6]['dimension']==e7_d5a1[7]['dimension']==42
    assert e7_d5a1[6]['total_rank']==e7_d5a1[7]['total_rank']==756
    assert e7_d5a1[6]['full_root_source_degree_bound']==e7_d5a1[7]['full_root_source_degree_bound']==22
    assert not e7_d5a1[6]['complete'] and not e7_d5a1[7]['full_multiplicities_proved']
    assert e7_d5a1[7]['all_degree_exponent_divides_2'] and len(e7_d5a1[7]['degree_rows'])==42
    assert all(r['twice_full_weight_image_in_root_image'] for r in e7_d5a1[7]['degree_rows'])
    assert e7_d5a1[8]['samples']==333 and e7_d5a1[8]['basic_invariant_algebra_complete']
    assert e7_d5a1[9]['exact_multiplicities'][21:23]==[17,19]
    assert e7_d5a1[9]['unknown_multiplicity_degrees']==list(range(23,41))+[42]
    assert e7_d5a1[9]['proper_E7_parabolics_with_complete_formulas']==113
    # Historical bounded reports above remain reproducible. The source
    # Hilbert certificates and complete image audits below supersede their
    # conservative incomplete flags; none are promoted by a cutoff alone.
    completion_patterns=(
        'a5*','audit_a5*','audit_e7_a5*','e7-a5*','e7_a5*',
        '*even_block*','e7-even-block*','*d4a1*',
        '*d5a1_initial_source*','e7-d5a1-initial-source*',
        '*d5a1_extra23*','e7-d5a1-extra23*',
        '*d5a1_orbit_chern*','e7-d5a1-orbit-chern*',
        '*d5a1_refinement_bound*','e7-d5a1-refinement-bound*',
        'audit_e7_initial_sources.py','summarize_e7_d5a1_complete.py',
        'e7-d5a1-complete-summary.json','summarize_e7_all_parabolics.py',
        'e7-all-parabolics-summary.json')
    completion_paths=sorted({p for pattern in completion_patterns
        for p in (ROOT/'experiments').glob(pattern) if p.is_file()
        and (p.suffix in ('.py','.json') or p.name.endswith('.json.gz'))})
    completion_checks=[]
    for path in completion_paths:
        if path.suffix!='.json':continue
        check=json.loads(path.read_text())
        assert check['result'] in ('passed','prepared'),path
        if check['result']=='prepared':assert not check['complete']
        completion_checks.append(check)
    complete_e7=read('experiments/e7-all-parabolics-summary.json')
    complete_d5a1=read('experiments/e7-d5a1-complete-summary.json')
    complete_d4a1=read('experiments/e7-d4a1-tower-images-audit.json')
    assert complete_e7['complete'] and complete_e7['proper_standard_parabolics']==127
    assert complete_e7['all_nonunit_relative_factors']==2
    assert complete_e7['first_gap_counts']=={'1':112,'2':12,'4':2,'8':1}
    assert len({tuple(row['retained_nodes']) for row in complete_e7['parabolics']})==127
    assert complete_d5a1['complete_multiplicity_formula']
    assert complete_d5a1['full_root_source_generators']==25
    assert complete_d5a1['unknown_multiplicity_degrees']==[]
    assert complete_d5a1['total_nonunit_factors']==366
    assert complete_d4a1['complete'] and complete_d4a1['degree_bound']==50
    assert complete_d4a1['total_rank']==7560
    assert complete_d4a1['total_nonunit_factors']==3618
    assert read('experiments/e7-d4a1-weight-transfer-audit.json')['full_weight_source_transfer_verified']
    for levi,count in (('d5a1',25),('d4a1',18)):
        check=read('experiments/e7-'+levi+'-initial-source-audit.json')
        assert check['full_source_complete_all_degrees']
        assert check['entire_hilbert_series_verified']
        assert len(check['generator_degrees'])==count
    for blocks,total in (('42',4911),('222',29982)):
        check=read('experiments/e7-even-block-'+blocks+'-comparison-audit.json')
        assert check['complete'] and check['total_nonunit_factors']==total
    assert read('experiments/e7-even-block-formula.json')['uniform_even_block_formula_verified']
    # Bind the printed coefficient lists as well as their proof data.
    def numeric_blocks(note):
        return [list(map(int,re.findall(r'\d+',block))) for block in
            re.findall(r'^```\n([0-9,\s.]+)\n```',(ROOT/note).read_text(),re.MULTILINE)]
    assert numeric_blocks('E7_D4_A1_COMPARISON.md')[0]==complete_d4a1['multiplicity_coefficients']
    p6_rows=re.findall(r'^\| \d+ through \d+ \| ([0-9, ]+) \|',
        (ROOT/'E7_D5_A1_COMPARISON.md').read_text(),re.MULTILINE)
    assert [int(value) for row in p6_rows for value in row.split(',')]==complete_d5a1['multiplicity_coefficients']
    a5_lists=numeric_blocks('E7_A5_FAMILY.md')
    a5_formula=read('experiments/e7-even-block-formula.json')
    assert a5_lists[0]==a5_formula['base_multiplicity_polynomial']
    assert a5_lists[1]==a5_formula['cases'][0]['multiplicity_coefficients']
    for check in e6_checks+e7_checks+e7_refinements+e7_full+e7_d5+e7_d4+e7_d6+e7_d5a1+e7_d5a1_attempts:
        for name,expected in check.get('source_sha256',{}).items():
            path=ROOT/name
            if not path.exists():path=ROOT/'experiments'/name
            assert digest(path)==expected,name
        for name,expected in check.get('raw_sha256',{}).items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
        if 'certificate_archive' in check:
            assert sha256((ROOT/check['certificate_archive']).read_bytes()).hexdigest()==check['certificate_archive_sha256']
        if 'tower_archive' in check:
            assert sha256((ROOT/check['tower_archive']).read_bytes()).hexdigest()==check['tower_archive_sha256']
        if 'seed_certificate_archive' in check:
            assert sha256((ROOT/check['seed_certificate_archive']).read_bytes()).hexdigest()==check['seed_archive_sha256']
    for check in completion_checks:
        for field in ('source_sha256','raw_sha256'):
            for name,expected in check.get(field,{}).items():
                path=ROOT/name
                if not path.exists():path=ROOT/'experiments'/name
                assert digest(path)==expected,(field,name)
        for field in ('certificate_archive','tower_archive'):
            if field in check:
                assert digest(ROOT/check[field])==check[field+'_sha256'],check[field]
    assert regression['version']==wheel['version']==version=='0.22.0'
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
    # Restrict to the actual representation cases: a Levi fiber-transfer
    # audit is a different artifact from a specialized Springer-map fiber.
    fiber_paths=sorted({p for cartan_type,representation in degrees
        for p in (ROOT/'experiments').glob(
            cartan_type.lower()+'-'+representation+'-fiber-*.json')})
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
    paths=fiber_paths+runtime+completion_paths+[ROOT/name for name in other['source_sha256']]+[ROOT/'other-theories-verification.json']+[ROOT/name for name in (
        'EXCEPTIONAL_DEGREES.md','experiments/exceptional_degree_geometry.py',
        'experiments/exceptional-degree-geometry.json','experiments/exceptional_degree_fibers.py',
        'experiments/exceptional-degree-bounded-attempts.json',
        'experiments/e7_adjoint_module_fiber.py',
        'experiments/e7_adjoint_shifted_fiber.py',
        'experiments/e7_adjoint_tower_fiber.py',
        'experiments/e6_matrix_oracle.py','experiments/e7_matrix_oracle.py',
        'README.md','../README.md','../articles/README.md',
        'CLASSICAL_EQUIVARIANT.md','SYMPLECTIC_K_FITTING.md',
        'SYMPLECTIC_K_MINORS.md','symplectic-minors-verification.json',
        'tools/verify_symplectic_minors.py',
        'SYMPLECTIC_K_BALANCING.md','symplectic-balancing-verification.json',
        'tools/verify_symplectic_balancing.py',
        'equivariant-verification.json','tools/verify_equivariant.py',
        'tests/test_equivariant.py',
        '../articles/classical-equivariant-exposition/README.md',
        '../articles/classical-equivariant-exposition/main.tex',
        '../articles/classical-equivariant-exposition/verify_pdf.py',
        '../articles/classical-equivariant-exposition/verify_examples.py',
        '../articles/classical-equivariant-exposition/verification.json',
        '../articles/classical-equivariant-exposition/examples-verification.json',
        '../articles/classical-equivariant-exposition/output/pdf/Classical_Springer_Borel_equivariant_exposition.pdf',
        'E6_GEOMETRY.md','E7_GEOMETRY.md',
        'RESEARCH_LEDGER.md','AGENTS.md',
        'LOCAL_IMAGE_BOUND.md','EQUIVARIANT_IMAGE_PILOT.md','LOCAL_IMAGE_AUDIT.md',
        'E6_ROOT_WEIGHT_COMPARISON.md','E6_ROOT_WEIGHT_AUDIT.md',
        'E6_REMAINING_CASES.md','E6_REMAINING_AUDIT.md','E6_REMAINING_TABLES.md',
        'E7_ROOT_WEIGHT_COMPARISON.md','E7_ROOT_WEIGHT_AUDIT.md',
        'E7_REFINEMENT_MULTIPLICITIES.md',
        'E7_FULL_FLAG_COMPARISON.md','E7_D5_REFINEMENT.md','E7_D4_REFINEMENT.md','E7_REMAINING_PARABOLICS.md','E7_D5_A1_COMPARISON.md',
        'E7_LOCAL_SOURCES.md','E7_A5_FAMILY.md','E7_D4_A1_COMPARISON.md','E7_ALL_PARABOLICS.md',
        'experiments/e7_d5a1_sources.py','experiments/extend_e7_d5a1_sources.py',
        'experiments/audit_e7_d5a1_sources.py','experiments/e7_d5a1_norm.py','experiments/audit_e7_d5a1_norm.py',
        'experiments/e7_d5a1_images.py','experiments/audit_e7_d5a1_images.py','experiments/audit_e7_d5a1_moments.py',
        'experiments/e7_d5a1_source_algebra.py','experiments/e7_d5a1_free_source.py','experiments/summarize_e7_d5a1.py',
        'experiments/e7-d5a1-sources.json','experiments/e7-d5a1-sources-audit.json',
        'experiments/e7-d5a1-sources-22.json','experiments/e7-d5a1-sources-22-audit.json',
        'experiments/e7-d5a1-norm.json','experiments/e7-d5a1-norm-audit.json',
        'experiments/e7-d5a1-images.json','experiments/e7-d5a1-images-audit.json','experiments/e7-d5a1-moments-audit.json',
        'experiments/e7-d5a1-source-algebra.json','experiments/e7-d5a1-free-source.json','experiments/e7-d5a1-summary.json',
        'experiments/e7-d5a1-sources-certificates.json.gz','experiments/e7-d5a1-sources-22-certificates.json.gz',
        'experiments/e7-d5a1-norm-certificates.json.gz','experiments/e7-d5a1-images-certificates.json.gz',
        'experiments/e7-d5a1-free-source-certificates.json.gz',
        'experiments/e7_d6_sources.py','experiments/extend_e7_d6_sources.py','experiments/pack_e7_d6_sources.py','experiments/e7_d6_images.py',
        'experiments/e7_d5a1_degree_two.py','experiments/audit_e7_d5a1_degree_two.py',
        'experiments/e7-d5a1-degree-two.json','experiments/e7-d5a1-degree-two-audit.json',
        'experiments/audit_e7_d6_sources.py','experiments/audit_e7_d6_moments.py','experiments/audit_e7_d6_images.py',
        'experiments/e7_remaining_conjugacy.py','experiments/audit_e7_remaining_conjugacy.py','experiments/summarize_e7_d6.py',
        'experiments/e7-d6-sources.json','experiments/e7-d6-sources-audit.json','experiments/e7-d6-images.json',
        'experiments/e7-d6-images-audit.json','experiments/e7-d6-moments-audit.json','experiments/e7-d6-summary.json',
        'experiments/e7-remaining-conjugacy.json','experiments/e7-remaining-conjugacy-audit.json',
        'experiments/e7-d6-source-seed20.json.gz','experiments/e7-d6-sources-certificates.json.gz','experiments/e7-d6-images-certificates.json.gz',
        'experiments/e7_d4_source.py','experiments/e7_d4_source_algebra.py','experiments/e7_d4_images.py',
        'experiments/audit_e7_d4_source.py','experiments/audit_e7_d4_images.py','experiments/e7_d4_summary_check.py',
        'experiments/e7-d4-source.json','experiments/e7-d4-source-algebra.json','experiments/e7-d4-images.json',
        'experiments/e7-d4-source-audit.json','experiments/e7-d4-images-audit.json','experiments/e7-d4-summary.json',
        'experiments/e7-d4-tower-certificates.json.gz','experiments/e7-d4-index-certificates.jsonl.gz',
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
        'tools/monomial_staircase.py','E7_RELATIVE_IMAGE_API.md',
        'tools/package_e7_relative_images.py','tools/verify_relative_images.py',
        'tests/test_relative_images.py','tests/data/e7_relative_audited_bases.json',
        'relative-image-verification.json',
        'E6_RELATIVE_IMAGE_API.md','G2_FULL_FLAG.md','E7_QUADRATIC_DERIVATION.md','E7_DEGREE11_LIFTING.md','CENTRAL_PARITY_KERNEL.md',
        'experiments/central_parity_kernel.py','experiments/audit_central_parity_kernel.py',
        'experiments/central-parity-kernel.json','experiments/central-parity-kernel-audit.json',
        'SPIN_PARITY_EXTENSION.md',
        'experiments/spin_parity_extension.py','experiments/audit_spin_parity_extension.py',
        'experiments/spin-parity-extension.json','experiments/spin-parity-extension-audit.json',
        'SPIN_PARITY_RECURSION.md',
        'experiments/spin_parity_recursion.py','experiments/audit_spin_parity_recursion.py',
        'experiments/check_spin_parity_quadratic.py',
        'experiments/spin-parity-recursion.json','experiments/spin-parity-recursion-audit.json',
        'experiments/spin-parity-quadratic-check.json','experiments/spin-parity-recursion-certificates.json.gz',
        'SPIN_PARITY_BOUNDARY.md',
        'experiments/spin_parity_boundary.py','experiments/audit_spin_parity_boundary.py',
        'experiments/check_spin_parity_even_source.py',
        'experiments/spin-parity-boundary.json','experiments/spin-parity-boundary-audit.json',
        'experiments/spin-parity-boundary-certificates.json.gz',
        'experiments/spin-parity-boundary-rank15.json','experiments/spin-parity-boundary-rank15-audit.json',
        'experiments/spin-parity-boundary-rank15-certificates.json.gz','experiments/spin-parity-even-source-check.json',
        'experiments/e7_degree11_lifting.py','experiments/audit_e7_degree11_lifting.py',
        'experiments/e7-degree11-lifting.json','experiments/e7-degree11-lifting-audit.json',
        'experiments/e7-degree11-lifting-certificates.json.gz',
        'experiments/e7_quadratic_derivation.py','experiments/audit_e7_quadratic_derivation.py',
        'experiments/e7-quadratic-derivation.json','experiments/e7-quadratic-derivation-audit.json',
        'experiments/e7-quadratic-derivation-certificates.json.gz',
        'experiments/audit_type_a_quadratic_derivation.py','experiments/type-a-quadratic-derivation-audit.json',
        'tools/package_e6_relative_images.py','tests/test_e6_relative_images.py',
        'tests/test_g2_full_flag.py','tests/data/e6_relative_audited_base.json',
        'experiments/summarize_e6_all_parabolics.py','experiments/e6-all-parabolics-summary.json',
        'experiments/g2_full_flag.py','experiments/audit_g2_full_flag.py',
        'experiments/g2-full-flag.json','experiments/g2-full-flag-audit.json',
        'G2_CHARACTER_FORM.md','experiments/g2-character-form-audit.json',
        'tests/test_exceptional_degrees.py','character-verification.json',
        'documentation-verification.json','geometry-wheel-verification.json')]+[Path(__file__).resolve()]
    report=dict(result='passed',version=version,
        odd_spin_boundary=dict(adjacent_rank_obstruction_theorem=True,
            boundary_obstructions_elementary_two_groups=True,
            exact_degreewise_obstruction_rank_and_first_failure_formula=True,
            even_source_C='M[z^2]',
            entire_kernel_criterion='log2[N:C]=2^(n-1)',
            degreewise_obstruction_formula='a_(d-1)-e_d-e_(d-1)',
            new_rank15_first_failure=37,independent_original_derivative_replay=True,
            lattice_audit='exact inclusion, diagonal 2-primary bound, full row rank modulo two; no HNF software',
            separate_same_rank_even_source_check=True,
            closed_numerical_rank_classification_proved=False,runtime_unchanged=True),
        odd_spin_recursive_kernel=dict(arbitrary_rank_recursive_presentation=True,
            positive_all_degree_rank_residues_modulo_4=[1,2],
            every_rank_stable_through_polynomial_degree_n=True,
            finite_entire_kernel_decisions_through_rank=12,
            additional_stable_rank=11,additional_first_failure={'12':34},
            native_Hermite_and_bit_certificate_audit=True,
            separate_quadratic_straightening_check=True,
            arbitrary_rank_closed_classification_complete=False,runtime_unchanged=True),
        odd_spin_central_parity=dict(coefficient_ring='Z_(2)',parabolic='B',
            all_rank_symmetric_source_reduction=True,
            degree_four_failure_exactly_rank_three=True,
            stable_entire_kernel_ranks=[4,5,6],
            first_failure_degrees={'3':4,'7':11,'8':16},
            independent_finite_flatness_certificates=True,
            arbitrary_rank_classification_complete=False,runtime_unchanged=True),
        central_parity_kernel=dict(type_C_all_rank_all_parabolic_theorem=True,
            type_A_sufficient_family='SL_(2r), any r, at least one odd Levi block',
            Spin7_full_flag_degree_four_counterexample=True,
            Spin7_independent_polynomial_proof=True,
            finite_C3_parabolic_checks=3,root_data_classification_complete=False),
        structural_quadratic_derivation=dict(E7_full_flag_maximum_degree=13,
            E7_uniform_family_size=112,E7_initial_uniform_family_certificate_maximum_degree=10,
            E7_uniform_family_maximum_degree=11,E7_good_Levi_maximum_degree=13,
            E7_uniform_all_degree_proof=False,
            type_A_arbitrary_rank_theorem='SL_(2r), r odd, odd Levi block; 2-local',
            type_A_independent_Grassmannian_ranks=[6,20,10,120,252],
            runtime_unchanged=True),
        E7_degree_eleven_lifting=dict(degree=11,uniform_family_size=112,
            initial_kernel_rank_mod2=461,primitive_kernel_rank_mod2=462,
            coefficient_syzygy_dimension=1,saturation_quotient_order=2,
            exceptional_cubic_nonlifting_cases=4,good_Levi_count=108,
            certificate='experiments/e7-degree11-lifting.json',
            independent_audit='experiments/e7-degree11-lifting-audit.json',
            uniform_all_degree_proof=False),
        e6_relative_image_api=dict(proper_standard_parabolics=63,nonunit_factors=[3,9,27],
            certificate='experiments/e6-all-parabolics-summary.json',
            independent_verification=relative_api['E6'],abstract_relative_groups_only=True),
        g2_full_flag=dict(complete_absolute_images=True,representations=['minimal','adjoint'],
            integral_kernel=['X^2-3XY+3Y^2','Y^6'],free_image_rank=12,
            minimal_smith_factors=[r['smith_factors'] for r in g2_full['cases'][0]['degrees']],
            independent_chevalley_columns=72,shared_full_source_theorem=True),
        e7_relative_image_api=dict(proper_standard_parabolics=127,improper_G_over_G=True,
            compressed_smith_multiplicities=True,scaled_root_weight_comparisons=True,
            named_representations=['minimal','adjoint'],abstract_relative_groups_only=True,
            independent_formula_and_lattice_audit=relative_api,
            scope='Packages the previously proved finite E7 theorem and elementary scalar consequences; no new image theorem or absolute Schubert embedding.'),
        current_e7_relative_image_status='Complete formulas in all degrees for all 127 proper standard E7 parabolics. Earlier progress fields record the stated proof stages; the final completion field is authoritative for current scope.',
        other_theories_exploration=dict(reports_checked=other['reports_checked'],
            proof_notes=other['proof_notes_and_program'],ledger_entries=other['ledger_entries'],
            source_hashes_checked=True,symbolic_SL3_integral_reductions=True,
            historical_exploration_runtime_changed=False,
            scope='Finite K/oriented comparisons, representation-dependent arithmetic, equivariant determinants/local patterns, quantum projective and affine A/C families, controlled folding and pointwise loop functoriality. Proof notes and finite exact corroboration; no priority certificate, full exceptional K Smith tables or central/rotation-compatible Kac-Moody construction.'),
        classical_equivariant_API=classical_equivariant,
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
            independent_relation_constructions=2,direct_bgg_degrees=10,all_degree_unknown_base_polynomials=0,all_degree_refinements_settled=64,
            source_obstruction_degree=16,scope='Uniform degree-nine theorem extended to all degrees on all 64 P7 refinements by separate P7/A6/D5/D4 proofs. No novelty claim.'),
        e7_full_flag_comparison=dict(formula='q P_W/(1+q)',intermediary='E7/P2 with A6 Levi',
            intermediary_rank=576,intermediary_dimension=42,intermediary_all_degrees_audited=True,
            full_flag_total_factors_2=1451520,proper_good_levi_parabolics=108,additional_P7_basis_case=1,
            proper_parabolics_settled=109,P7_refinements_settled=61,
            scope='Complete computer-assisted A6 image proof and classical common-fiber transfer. A6-transfer scope only; the D5 pair and D4 are completed separately below. No novelty or ambient H/N classification claim.'),
        e7_d5_refinement=dict(source_generators=10,source_complete_all_degrees=True,quotient_rank=1512,dimension=43,
            independent_image_degrees=43,total_relative_factors_2=756,proper_parabolics_settled=111,P7_refinements_settled=63,
            scope='Complete computer-assisted D5 full-source comparison, with separate native-Python source, weight and image audits. E6 opposition gives the second D5 case; D4 is completed separately below.'),
        e7_d4_completion=dict(source_generators=10,source_complete_all_degrees=True,relative_module_generators=34,
            quotient_rank=15120,dimension=51,independently_audited_degrees=52,total_relative_factors_2=7560,
            proper_parabolics_settled=112,P7_refinements_settled=64,other_proper_parabolics_remaining=15,
            central_character_nodes=[2,5,7],remaining_parabolics_have_zero_degree_one_gap=True,
            scope='The 112 cases omitting at least one of nodes 2,5,7 have polynomial q P_W/((1+q)P_WI), including all 64 P7 refinements. D4 uses the full finite module and equal image-index corrections, independently audited by source/Buchberger and native Smith elimination. The D5 base is a shared audited input. The other 15 require a different formula. Structural explanation, ambient H/N and priority remain separate.'),
        e7_d6_historical_checkpoint=dict(new_parabolic='E7/P1, D6 Levi',quotient_rank=126,dimension=33,
            full_weight_source_degree_bound=33,full_root_source_degree_bound=27,
            all_33_image_degrees_independently_audited=True,
            multiplicity_polynomial=e7_d6[5]['multiplicity_polynomial'],total_relative_factors_2=e7_d6[5]['total_relative_factors_2'],
            proper_parabolics_settled=113,other_proper_parabolics_remaining=14,remaining_distinct_Weyl_conjugacy_classes=10,
            partial_D5_A1_first_gap_degree=2,partial_D5_A1_first_relative_group='Z/2',
            scope='Historical D6 step: full source components and a full-image sandwich prove the stated formula, with separate native audits. Four Weyl-conjugacy reductions are checked in a second basis. The then-remaining fourteen diagrams are now completed in the final field below.'),
        e7_d5a1_historical_bounds=dict(parabolic='E7/P6',levi_type='D5+A1',dimension=42,schubert_rank=756,
            all_degree_exponent_divides_2=True,full_source_multiplicity_degree_bound=22,
            independent_image_degrees=42,new_necessary_source_generator_degrees=[21,22],
            incomplete_generator_algebra_first_missing_degree=21,
            higher_multiplicity_unknown_degrees=e7_d5a1[9]['unknown_multiplicity_degrees'],
            proper_parabolics_with_complete_formulas=113,
            scope='Historical checkpoint: bounded reports remain unchanged. The new all-degree Hilbert source theorem and degree-23 witness below make the earlier image bounds exact. Earlier free-source trials remain incomplete and are not used as proof.'),
        e7_all_parabolics_completion=dict(proper_standard_parabolics=127,
            uniform_formula_parabolics=112,remaining_formula_parabolics=15,
            all_nonunit_relative_factors=2,first_gap_counts=complete_e7['first_gap_counts'],
            d5a1_full_root_source_generators=25,d5a1_total_factors_2=366,
            d4a1_full_root_source_generators=18,d4a1_weight_source_module_generators=34,
            d4a1_independently_audited_degrees=50,d4a1_total_factors_2=3618,
            a5_family_standard_diagrams=12,even_block_base_totals=[4911,29982],
            all_degree_source_hilbert_identities_independently_verified=True,
            full_embedded_image_calculations_independently_verified=True,
            printed_multiplicity_lists_match_audited_data=True,
            scope='Current complete relative full root/weight comparison in every degree for all proper E7 parabolics. Earlier counts above record historical proof stages. The even-block common formula is exact for the finite E7 family; geometric folding, arbitrary-rank generalizations, ambient H/N and literature priority remain separate. Version 0.22 packages both E6 and E7 relative classifications.'),
        hash_semantics='SHA256 after CRLF to LF normalization for text; raw bytes for .gz archives and .pdf documents',
        audit_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(set(paths))})
    (ROOT/'release-0.22-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='audit_sha256'},indent=2))


if __name__=='__main__':
    main()
