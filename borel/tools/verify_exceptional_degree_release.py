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
    review_package=read('review/adjoint-generic-fibers/package-verification.json')
    review_manifest=json.loads((ROOT/'review/adjoint-generic-fibers/package-manifest.json').read_text())
    assert review_package['archived_bytes_checked']
    assert review_package['replay_from_extracted_archive']
    assert review_package['fresh_pinned_engine_build_passed']
    assert review_package['proof_commit']==review_manifest['proof_commit']
    assert review_package['input_files']==len(review_manifest['files'])
    for entry in review_manifest['files']:
        assert digest(ROOT/entry['path'])==entry['sha256'],entry['path']
    theorem_audit=read('adjoint-theorem-audit.json')
    assert theorem_audit['proper_E_subsets_independently_enumerated']==445
    assert theorem_audit['cubic_standard_subset_transfers_replayed']==136
    assert theorem_audit['reports_bound']==8
    assert all(theorem_audit['exact_scheme_controls'].values())
    for name,expected in theorem_audit['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    e7_final=read('e7-final-families-verification.json')
    e7_final_audit=read('e7-final-families-audit.json')
    g2_complete=read('g2-adjoint-fibers-verification.json')
    f4_complete=read('f4-adjoint-fibers-verification.json')
    f4_audit=read('f4-adjoint-fibers-audit.json')
    e7_membership=read('e7-companion-membership-verification.json')
    e7_geometry=read('e7-companion-geometry-audit.json')
    e8_weyl=read('e8-weyl-transfers-verification.json')
    e8_weyl_audit=read('e8-weyl-transfers-audit.json')
    e8_triple=read('e8-a1-triple-verification.json')
    e8_triple_audit=read('e8-a1-triple-audit.json')
    e8_pair=read('e8-a2-a1-pair-verification.json')
    e8_pair_audit=read('e8-a2-a1-pair-audit.json')
    e8_final=read('e8-final-cases-verification.json')
    e8_final_audit=read('e8-final-cases-audit.json')
    e8_four=[read('e8-a1-four-%d-verification.json'%p) for p in (32003,32009)]
    # Earlier reports below remain dated checkpoints. This is the current
    # last-case certificate, which supersedes their open completion flags.
    e8_complete=[read('e8-last-%d-verification.json'%p) for p in (32003,32009)]
    e8_complete_audit=read('e8-last-case-audit.json')
    for record in e8_complete+[e8_complete_audit]:
        assert record['covered_proper_E8_standard_subsets']==255
        assert record['remaining_Weyl_classes']==0
        assert record['characteristic_zero_affine_singular_dimension_bound']==2
        for name,expected in record['source_sha256'].items():assert digest(ROOT/name)==expected,name
        for name,expected in record['raw_sha256'].items():assert digest(ROOT/name)==expected,name
    for record in e8_complete:
        assert record['original_ideal_memberships_replayed']
        assert record['homogeneous_degree_upper_bound']==record['accounted_infinity_degree']==21
    assert e8_final['covered_proper_E8_standard_subsets']==247
    assert e8_final['remaining_Weyl_classes']==1 and e8_final['remaining_standard_subsets']==8
    assert e8_final['unresolved_classes']==['A2+A1^3'] and not e8_final['all_E8_Levis_resolved']
    assert len(e8_final['newly_completed_classes'])==4
    assert e8_final_audit['native_full_root_characters']==5
    assert e8_final_audit['independent_actual_scheme_eliminations']==5
    assert e8_final_audit['multiplier_homogeneous_dimension']==3
    for record in [e8_final,e8_final_audit]+e8_four:
        for name,expected in record['source_sha256'].items():assert digest(ROOT/name)==expected,name
        for name,expected in record.get('raw_sha256',{}).items():assert digest(ROOT/name)==expected,name
    e8_last=json.loads((ROOT/'e8-last-case-exploration.json').read_text())
    assert not e8_last['characteristic_zero_bound_proved'] and not e8_last['all_E8_completion_claim']
    assert e8_last['affine_jacobian_deficiency_dimension']==1 and e8_last['prime']==32003
    for name,expected in e8_last['source_sha256'].items():assert digest(ROOT/name)==expected,name
    assert e8_pair['new_standard_subsets']==28 and e8_pair['regular_fiber_dimension']==10
    assert e8_pair['covered_proper_E8_standard_subsets']==218
    assert e8_pair['remaining_Weyl_classes']==5 and e8_pair['remaining_standard_subsets']==37
    assert e8_pair_audit['two_complex_boundary_components_checked']
    assert e8_pair_audit['k_zero_intersections_retained']
    assert e8_pair_audit['all_six_ramification_wall_factorizations']
    assert e8_pair_audit['actual_commuting_E8_roots']==6
    assert e8_pair_audit['central_schur_determinant']==48384
    assert e8_triple['new_standard_subsets']==21 and e8_triple['regular_fiber_dimension']==6
    assert e8_triple['covered_proper_E8_standard_subsets']==190
    assert e8_triple['remaining_Weyl_classes']==6 and e8_triple['remaining_standard_subsets']==65
    assert e8_triple_audit['entire_boundary_control'] and e8_triple_audit['exact_boundary_ideal_identity']
    assert e8_triple_audit['central_logarithmic_hessian_determinant']==-870912
    assert e8_weyl['cubic_Weyl_classes']==10 and e8_weyl['unresolved_cubic_classes']==7
    assert e8_weyl['covered_proper_standard_subsets']==169
    assert e8_weyl['unresolved_standard_subsets']==86
    assert e8_weyl_audit['explicit_words_checked']==93
    assert e8_weyl_audit['full_orbits_reenumerated']==10
    assert e8_weyl_audit['unresolved_classes']==7
    assert not e8_weyl['all_simple_adjoint_conjecture_proved']
    assert e7_membership['original_affine_ideal_memberships']==6
    assert e7_membership['independent_prime']!=e7_membership['generation_prime']
    assert e7_geometry['monomial_dimension_bound']==1
    for name,expected in e7_membership['raw_sha256'].items():assert digest(ROOT/name)==expected,name
    assert f4_complete['all_proper_F4_Levis_covered'] and f4_audit['independent_root_characters']==15
    assert g2_complete['all_proper_G2_Levis_covered']
    assert g2_complete['proper_standard_Levis']==3
    assert e7_final['covered_standard_E7_Levis']==127 and e7_final['all_E7_Levis_resolved']
    assert e7_final['a2a2a1']['regular_fiber_dimension']==14
    for record in (e7_final,e7_final_audit,g2_complete,f4_complete,f4_audit,e7_membership,e7_geometry,e8_weyl,e8_weyl_audit,e8_triple,e8_triple_audit,e8_pair,e8_pair_audit):
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    exposition_root=ROOT.parent/'articles/springer-cayley-exposition'
    exposition=read('../articles/springer-cayley-exposition/verification.json')
    assert exposition['pages']==11 and exposition['single_manuscript_page_per_pdf_page']
    assert exposition['page_layout']=='SinglePage'
    assert not exposition['unresolved_references'] and not exposition['overfull_boxes']
    assert exposition['E7_resolved']==123 and exposition['E7_unresolved']==4
    assert exposition['E7_unresolved_Weyl_families']==2
    assert not exposition['companion_memberships_verified']
    assert not exposition['all_simple_adjoint_conjecture_proved']
    for name,expected in exposition['source_sha256'].items():
        assert digest(exposition_root/name)==expected,name
    exposition_sources=json.loads((exposition_root/'source-manifest.json').read_text())
    exposition_live_source_drift={}
    for name,expected in exposition_sources['source_sha256'].items():
        actual=digest(ROOT.parent/name)
        # Preserve the dated PDF's source manifest. These two status files
        # evolve independently; record drift rather than rewriting provenance.
        if name in ('borel/PROOF_STATUS.md','borel/RESEARCH_LEDGER.md') and actual!=expected:
            exposition_live_source_drift[name]=dict(artifact_snapshot=expected,current=actual)
        else:
            assert actual==expected,name
    binary=read('binary-branching-verification.json')
    binary_audit=read('binary-branching-audit.json')
    for record in (binary,binary_audit):
        assert record['new_standard_embeddings']==22
        assert record['resolved_proper_standard_Levis']==123
        assert record['remaining_embeddings']==4
        assert not record['all_simple_adjoint_conjecture_proved']
        assert not record['full_map_degree_claim']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    assert binary['remaining_Weyl_classes']==2
    assert binary['regular_fiber_dimensions']=={'A1^3 cubic class':6,'A2+A1^2 cubic class':10}
    assert not binary['companion_membership_certificate_used']
    assert not binary['all_E7_Levis_resolved']
    assert not binary_audit['uses_Sage']
    assert not binary_audit['companion_memberships_verified']
    assert binary_audit['generic_regular_fiber_dimensions']==[6,10]
    assert binary_audit['roots']['seven_A1_root_character']
    assert binary_audit['roots']['cover_indices']==[8,12]
    assert binary_audit['symbolic']['original_characters_reconstructed']==2
    assert all(binary_audit['symbolic'][key] for key in (
        'original_target_derivatives','primitive_and_squarefree_certificates',
        'boundary_identities','direct_logarithmic_hessians','ramification_target_hyperplane'))
    assert binary_audit['Weyl']['Weyl_words_checked']==22
    assert binary_audit['Weyl']['entire_subsystems_compared']
    assert binary_audit['Weyl']['coverage_complete_within_two_cubic_classes']
    fiber_review=read('cayley-fiber-review-audit.json')
    assert fiber_review['control_count']==8 and all(fiber_review['controls'].values())
    assert not fiber_review['uses_Sage']
    assert fiber_review['E7_resolved']==101 and fiber_review['E7_unresolved']==26
    assert not fiber_review['selected_memberships_verified']
    assert not fiber_review['new_group_classification']
    assert not fiber_review['all_simple_adjoint_conjecture_proved']
    assert not fiber_review['literature_priority_claim']
    for name,expected in fiber_review['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    uniform_obstruction=read('uniform-adjoint-obstruction-audit.json')
    assert uniform_obstruction['proposed_stronger_lemma_false']
    assert not uniform_obstruction['generic_conjecture_disproved']
    assert not uniform_obstruction['all_simple_adjoint_conjecture_proved']
    assert not uniform_obstruction['new_Levi_coverage']
    assert uniform_obstruction['transverse_length']==3
    assert uniform_obstruction['adjoint_differential_rank']==7
    assert uniform_obstruction['source_regular_semisimple']
    assert uniform_obstruction['local_scheme_exact']
    assert uniform_obstruction['adjoint_Casimir_sandwich_formula_verified']
    assert uniform_obstruction['conjugation_Casimir_formula_verified']
    assert uniform_obstruction['full_equivariance_operator_identity_verified']
    assert uniform_obstruction['kernel_Lie_subalgebra_shortcut_false']
    assert uniform_obstruction['kernel_control_coroot_eigenvalue']=='-7/4'
    higher=uniform_obstruction['higher_order']
    assert higher['corrected_quartic_symbolic_check']
    assert higher['Sp4_raw_quartic']=='0'
    assert higher['Sp4_effective_quartic']=='3*x^4/2'
    assert higher['Sp4_implicit_expansion_matches_exact_local_ideal']
    assert higher['E6_excess_residual_constant_to_all_orders']
    assert higher['four_jet_vanishing_not_sufficient']
    assert not higher['generic_bad_image_bound_proved']
    generic=uniform_obstruction['generic_target']
    assert generic['Sp4_generic_fiber_jacobian_unit']
    assert generic['Sp4_independent_branch_separation']
    assert generic['E6_original_generic_ideal_verified']
    assert generic['E6_smooth_excess_dimension']==2
    assert generic['source_coefficient_localization_can_erase_dominant_fiber']
    assert generic['altered_Levi_central_polynomial_squarefree']
    assert generic['altered_Levi_components_generically_smooth']
    assert generic['altered_Levi_primary_decomposition_verified']
    assert not generic['generic_adjoint_smoothness_proved']
    assert not generic['new_Levi_coverage']
    reduced=uniform_obstruction['reducedness']
    assert reduced['complete_intersection_does_not_force_R0']
    assert reduced['existing_embedded_control_R0_without_S1']
    assert reduced['target_torsion_control_annihilator']=='c-1'
    assert reduced['target_torsion_control_generic_ring_zero']
    assert reduced['Sp4_relative_incidence_squarefree']
    assert reduced['E6_relative_incidence_radical']
    assert reduced['E6_displayed_full_gradient_not_regular_sequence']
    assert not reduced['all_adjoint_R0_proved']
    assert not reduced['all_adjoint_S1_proved']
    assert not reduced['uniform_generic_reducedness_proved']
    assert not reduced['new_Levi_coverage']
    for name,expected in uniform_obstruction['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    nilpotent_search=read('adjoint-nilpotent-search-audit.json')
    assert nilpotent_search['result']=='passed'
    assert not nilpotent_search['counterexample_found']
    assert not nilpotent_search['uniform_conjecture_proved']
    assert nilpotent_search['first_model']['entire_relative_critical_surface_smooth']
    assert nilpotent_search['first_model']['all_central_boundaries_retained']
    assert not nilpotent_search['second_model_source_screen_alone_resolves_boundary']
    assert nilpotent_search['independent_source_screen']['standard_monomials']==38
    assert nilpotent_search['independent_source_screen']['jacobian_ideal_unit']
    for name,expected in nilpotent_search['sha256'].items():
        assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    source_screen=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    assert source_screen['generic_relative_algebra_finite_etale_by_spreading']
    assert source_screen['monic_generic_basis_specializes_exactly']
    assert source_screen['source_boundary_not_tested']
    assert not source_screen['target_generic_claim']
    boundary_direct=json.loads((ROOT/'adjoint-simple-boundaries.json').read_text())
    assert boundary_direct['status']=='completed'
    assert len(boundary_direct['records'])==6
    assert all(r['dimension']==1 and r['cannot_dominate_two_dimensional_target']
               for r in boundary_direct['records'])
    boundary_audit=read('adjoint-boundary-exclusions-audit.json')
    assert boundary_audit['independently_replayed_count']==6
    assert not boundary_audit['partial_audit']
    assert boundary_audit['direct_membership_boundaries']==4
    assert boundary_audit['potential_involution_boundaries']==2
    assert boundary_audit['target_nondominance_proved_for_six_boundaries']
    assert not boundary_audit['entire_remaining_boundary_resolved']
    for name,expected in boundary_audit['source_sha256'].items():
        path=ROOT/'tools'/name if name.endswith('.py') else ROOT/name
        assert sha256(path.read_bytes()).hexdigest()==expected,name
    boundary_certificate=read('adjoint-boundary-exclusions-verification.json')
    assert boundary_certificate['boundary_count']==6
    assert boundary_certificate['direct_membership_boundaries']==4
    assert boundary_certificate['potential_involution_boundaries']==2
    for name,expected in boundary_certificate['source_sha256'].items():
        path=ROOT/'tools'/name if name.endswith('.py') else ROOT/name
        assert sha256(path.read_bytes()).hexdigest()==expected,name
    boundary_lift=json.loads((ROOT/'adjoint-central-boundary-lift-forward.json').read_text())
    assert boundary_lift['both_ideal_inclusions_checked']
    assert boundary_lift['monic_basis_generators']==28
    assert len(boundary_lift['denominator_factors'])==18
    assert sha256((ROOT/boundary_lift['certificate_file']).read_bytes()).hexdigest()==boundary_lift['certificate_sha256']
    reverse_lift=json.loads((ROOT/'adjoint-central-boundary-lift-reverse.json').read_text())
    assert reverse_lift['both_ideal_inclusions_checked']
    assert reverse_lift['monic_basis_generators']==28
    assert sha256((ROOT/reverse_lift['certificate_file']).read_bytes()).hexdigest()==reverse_lift['certificate_sha256']
    flatness_audit=read('adjoint-central-flatness-audit.json')
    assert flatness_audit['result']=='passed'
    assert [r['length'] for r in flatness_audit['checks']]==[37,37,38]
    assert flatness_audit['finite_over_entire_b_open_ruled_out']
    assert not flatness_audit['whole_open_flatness_proved']
    assert not flatness_audit['generic_target_reducedness_proved']
    for name,expected in flatness_audit['source_sha256'].items():
        path=ROOT/'tools'/name if name.endswith('.py') else ROOT/name
        assert sha256(path.read_bytes()).hexdigest()==expected,name
    for filename in ('adjoint-central-flatness-screen.json','adjoint-central-flatness-curve-screen.json'):
        for name,expected in json.loads((ROOT/filename).read_text())['source_sha256'].items():
            path=ROOT/'tools'/name if name.endswith('.py') else ROOT/name
            assert sha256(path.read_bytes()).hexdigest()==expected,name
    intersections=json.loads((ROOT/'adjoint-central-flatness-cover-intersections.json').read_text())
    assert intersections['residual_pairs']==22
    assert all(r['dimension']==0 for r in intersections['records'])
    # The earlier screens retain their historical limited scope. The following
    # structural proof now controls every fiber on the allowed source open.
    quasifinite=read('adjoint-quasifinite-verification.json')
    quasifinite_audit=read('adjoint-quasifinite-audit.json')
    assert quasifinite['all_allowed_critical_fibers_zero_dimensional']
    assert quasifinite['relative_algebra_flat_and_reduced_on_b_open']
    assert quasifinite['relative_scheme_pure_dimension']==2
    assert quasifinite['generic_target_finite_etale_or_zero_by_dimension_argument']
    assert quasifinite['generic_target_smoothness_proved']
    assert not quasifinite['uniform_adjoint_conjecture_proved']
    assert quasifinite['compact_center_identity_verified']
    assert quasifinite_audit['compact_center_identity_independently_verified']
    assert quasifinite_audit['centered_line_obstruction_membership_verified']
    assert len(quasifinite_audit['rank_zero_branches'])==8
    assert not quasifinite_audit['uses_Sage'] and not quasifinite_audit['uses_modular_specialization']
    for report in (quasifinite,quasifinite_audit):
        for name,expected in report['source_sha256'].items():
            assert sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    proof_paths=read('adjoint-proof-paths-audit.json')
    proof_paths_sage=read('adjoint-proof-paths-verification.json')
    assert proof_paths['carriers']['cubic_highest_roots']==386
    assert proof_paths['carriers']['existing_character_records']==136
    assert proof_paths['carriers']['counts']==proof_paths_sage['carrier_counts']
    assert proof_paths['smooth_relative_fiber_bad_image_exclusion_proved_in_note']
    assert not proof_paths['root_carriers_force_squarefree_schur_pencils']
    assert proof_paths['controls']['E6_original_generic_ideal_recovered']
    assert proof_paths['controls']['E6_entire_pencil_zero_on_dominant_source_stratum']
    assert proof_paths['controls']['reduced_singular_axes_control']
    assert not proof_paths['geometric_argument_verified_by_computer']
    assert proof_paths_sage['E6_actual_generic_ideal_equal']
    assert proof_paths_sage['independent_of_native_weight_coordinates']
    assert not proof_paths['invariant_operators']['flat_invariant_operators_force_equal_root_multiplicities']
    assert proof_paths_sage['flat_invariant_operator_control']==proof_paths['invariant_operators']['root_moment_eigenvalues']
    assert proof_paths['overlap']['E6_full_coefficient_alignment']
    assert proof_paths_sage['overlap']['E6_reflection_forces_k_equals_one']
    assert proof_paths['overlap']['SL4_proper_centralizer_overlap_cancellation']
    assert proof_paths_sage['overlap']['SL4_matrix_trace_and_exact_number_field_check']
    for key in ('root_string_identities','A2_directional_syzygies','aligned_quadratic_original_ideal_chart_checks'):
        assert proof_paths['overlap'][key]==proof_paths_sage['overlap'][key]
    for report in (proof_paths,proof_paths_sage):
        assert report['overlap']['zero_discriminant_nilpotent_control']
        assert report['overlap']['aligned_quadratic_singular_scheme_target_factorization']
        assert not report['overlap']['general_remaining_bad_image_bound_proved']
    for report in (proof_paths,proof_paths_sage):
        assert not report['uniform_adjoint_conjecture_proved']
        assert not report['new_Levi_coverage']
        for name,expected in report['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    assert sha256((ROOT/'unit-constraint-cases.json.gz').read_bytes()).hexdigest()==proof_paths['archive_sha256']
    radial=read('adjoint-radial-schur-audit.json')
    radial_sage=read('adjoint-radial-schur-verification.json')
    assert radial['radial_character_identity_is_classical']
    assert radial_sage['independent_of_native_root_strings']
    assert set(radial['root_identities'])==set(radial_sage['root_identities'])=={'A3','B3','C3','D4','G2','F4','E6','E7','E8'}
    for kind,row in radial['root_identities'].items():
        assert row['coefficientwise_radial_identity']
        for key in ('root_count','casimir_eigenvalue','root_length_squares'):
            assert row[key]==radial_sage['root_identities'][kind][key]
    assert radial['root_identities']['E6']['altered_Levi_radial_defect']==radial_sage['root_identities']['E6']['altered_Levi_radial_defect']
    for report in (radial,radial_sage):
        assert report['schur']['E8_R_D_wedge']=='515/6'
        assert report['schur']['E8_Q_D_wedge']=='-525/2'
        assert report['schur']['universal_Schur_envelope_derivative_verified']
        assert report['schur']['SL2_radial_derivation_regular_on_quotient_discriminant']
        assert report['schur']['E8_pole_point_quadratic_equations_inconsistent']
        assert not report['general_singular_image_bound_proved'] and not report['new_Levi_coverage']
        for name,expected in report['source_sha256'].items():assert digest(ROOT/name)==expected,name
    singular_radial=read('adjoint-singular-radial-audit.json')
    singular_radial_sage=read('adjoint-singular-radial-verification.json')
    assert singular_radial['differential']['unrestricted_full_matrix_identity']
    assert singular_radial['differential']['external_root_pair_logarithmic_derivative']
    assert singular_radial['control']['bad_image_dominant']
    assert singular_radial['control']['actual_incidence_ideal']==['q^2','q*y']
    assert not singular_radial['control']['actual_adjoint_character']
    assert singular_radial_sage['checks']['nonzero_generic_fiber_nilpotent']
    assert singular_radial_sage['checks']['localization_keeps_generic_fiber']
    assert singular_radial_sage['checks']['generic_bad_support_length']==2
    assert not singular_radial_sage['native_audit_imported']
    for report in (singular_radial,singular_radial_sage):
        assert not report['abstract_radial_identity_suffices']
        assert not report['uniform_adjoint_conjecture_proved']
        assert not report['adjoint_counterexample'] and not report['new_Levi_coverage']
        for name,expected in report['source_sha256'].items():assert digest(ROOT/name)==expected,name
    character_discriminant=read('adjoint-character-discriminant-audit.json')
    character_discriminant_sage=read('adjoint-character-discriminant-verification.json')
    assert character_discriminant['uniform_Jacobian_formula_is_Steinberg_corollary']
    assert character_discriminant['checks']['A2_Jacobian_square_equals_four_external_discriminant']
    assert character_discriminant['checks']['failure_on_generic_critical_incidence']
    assert character_discriminant_sage['checks']['full_and_external_matrix_determinants_verified']
    assert character_discriminant_sage['checks']['generic_critical_plus_second_Adams_gradient_is_unit_ideal']
    assert not character_discriminant_sage['native_audit_imported']
    for report in (character_discriminant,character_discriminant_sage):
        assert not report['checks']['criticality_propagates_to_Adams_characters']
        assert report['checks']['internal_discriminant_retained']
        assert not report['uniform_adjoint_conjecture_proved']
        assert not report['adjoint_counterexample'] and not report['new_Levi_coverage']
        for name,expected in report['source_sha256'].items():assert digest(ROOT/name)==expected,name
    type_a_nil=read('type-a-incidence-nilradical-audit.json')
    type_a_nil_sage=read('type-a-incidence-nilradical-verification.json')
    assert type_a_nil['checks']['SL4_actual_incidence_nonreduced']
    assert type_a_nil['checks']['SL4_target_annihilates_nilradical']
    assert type_a_nil_sage['checks']['SL4_original_Laurent_ideal_radical_verified']
    assert type_a_nil_sage['checks']['SL4_nilpotent_and_target_annihilator_verified']
    assert type_a_nil_sage['checks']['SL4_generic_target_ideal_is_unit']
    assert not type_a_nil_sage['native_audit_imported']
    for report in (type_a_nil,type_a_nil_sage):
        assert report['all_type_A_nilradical_theorem_written']
        assert report['generic_type_A_reducedness_already_followed_from_quadratic_theorem']
        assert not report['uniform_all_simple_reducedness_proved']
        assert not report['new_Levi_coverage'] and not report['novelty_claim']
        for name,expected in report['source_sha256'].items():assert digest(ROOT/name)==expected,name
    annihilator_routes=read('adjoint-annihilator-routes-audit.json')
    annihilator_routes_sage=read('adjoint-annihilator-routes-verification.json')
    assert annihilator_routes['checks']['E6_full_root_character_reconstructed']
    assert annihilator_routes['checks']['no_nonzero_target_multiple_repairs_radial_descent']
    assert annihilator_routes['checks']['E6_exact_defect_at_z2']=='1647/52'
    assert annihilator_routes_sage['checks']['E6_radial_defect_exact']=='1647/52'
    assert annihilator_routes_sage['checks']['E6_regular_semisimple_source_verified']
    assert annihilator_routes_sage['checks']['altered_Levi_etale_center_retains_embedded_nilpotent']
    assert not annihilator_routes_sage['native_audit_imported']
    for report in (annihilator_routes,annihilator_routes_sage):
        assert report['checks']['E6_radial_eigenidentity_retained']
        assert not report['exceptional_target_annihilator_derived']
        assert not report['uniform_generic_reducedness_proved']
        assert not report['adjoint_counterexample'] and not report['new_Levi_coverage']
        for name,expected in report['source_sha256'].items():assert digest(ROOT/name)==expected,name
    aligned_nil=read('aligned-critical-annihilator-audit.json')
    aligned_nil_sage=read('aligned-critical-annihilator-verification.json')
    assert aligned_nil['checks']['E6_explicit_original_ideal_membership_certificate']
    assert aligned_nil['checks']['entire_three_ideal_intersection_verified']
    assert aligned_nil['checks']['squarefree_support_polynomial_alone_fails']
    assert aligned_nil_sage['E6']['original_Laurent_gradient_radical_verified']
    assert aligned_nil_sage['E6']['original_Laurent_annihilator_colon_verified']
    assert aligned_nil_sage['E6']['central_target_annihilates_entire_nilradical']
    assert aligned_nil_sage['excluded_delta_zero_embedded_control_verified']
    assert aligned_nil_sage['E6_contact_positive_excess_retained']
    assert aligned_nil_sage['altered_Levi_embedded_nilpotent_retained']
    assert not aligned_nil_sage['native_audit_imported']
    for report in (aligned_nil,aligned_nil_sage):
        assert report['aligned_family_exact_annihilator_proved_in_note']
        assert report['generic_smoothness_coverage_already_known']
        assert not report['all_exceptional_annihilators_proved']
        assert not report['new_Levi_coverage'] and not report['novelty_claim']
        for name,expected in report['source_sha256'].items():assert digest(ROOT/name)==expected,name
    shared=read('shared-square-couplings-audit.json')
    shared_sage=read('shared-square-couplings-verification.json')
    assert shared['checks']['central_derivative_chain_rule']
    assert shared['checks']['E8_quadratic_factor_ranks']==[4,4]
    assert shared_sage['checks']['smooth_component_dimensions']==[0,1,1]
    assert shared_sage['checks']['original_coupled_generic_ideal_exact_decomposition']
    assert shared_sage['checks']['opposite_simple_zero_embedded_annihilator_verified']
    assert shared_sage['checks']['zero_total_constant_reduced_singular_control']
    assert not shared_sage['native_audit_imported']
    for record in (shared,shared_sage):
        assert record['shared_square_family_theorem_in_note']
        assert record['checks']['universal_affine_shared_variable_normal_form_disproved']
        assert not record['full_root_normal_form_for_remaining_cases_proved']
        assert not record['new_exceptional_group_coverage']
        assert not record['uniform_conjecture_proved']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():assert digest(ROOT/name)==expected,name
    jets=read('weyl-jet-obstruction-audit.json')
    jets_sage=read('weyl-jet-obstruction-verification.json')
    assert jets['checks']['full_Weyl_invariance_verified']
    assert jets['checks']['first_two_derived_jets_vanish']
    assert jets['checks']['third_derived_jet_nonzero']
    assert jets_sage['checks']['positive_irreducible_decomposition']
    assert jets_sage['checks']['original_localized_gradient_target_ideal_equal_double_structure']
    assert jets_sage['checks']['nilradical_square_zero_nonzero']
    assert jets_sage['checks']['generic_covered_chart_length']==48
    assert jets_sage['checks']['reduced_covered_chart_length']==24
    assert not jets_sage['native_audit_imported']
    for record in (jets,jets_sage):
        assert record['checks']['representation_dimension']==640
        assert record['checks']['trace_form_HH']==10368
        assert len(record['checks']['adjoint_long_root_jet_projectors'])==6
        assert record['actual_simple_group_representation'] and record['representation_is_reducible']
        assert record['genericity_is_within_proper_Levi_central_stratum']
        assert record['infinite_tensor_power_family_proved_in_note']
        assert not record['representation_is_adjoint']
        assert not record['uniform_adjoint_conjecture_disproved']
        assert not record['new_adjoint_coverage'] and not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():assert digest(ROOT/name)==expected,name
    rigidity=read('central-rigidity-audit.json')
    assert rigidity['control_count']==8 and all(rigidity['controls'].values())
    assert rigidity['generic_reduced_central_image_follows_from_trace_identity'] is False
    assert rigidity['new_control_is_adjoint_counterexample'] is False
    assert rigidity['covers_all_simple_types']
    assert rigidity['covers_nonregular_source_points']
    assert rigidity['constancy_on_reduced_components']
    assert rigidity['finite_scheme_image_and_actual_fiber_decomposition_proved']
    assert rigidity['generic_abelianization_image_finite']
    assert rigidity['finite_rational_hyperplane_existence_proved']
    assert not rigidity['generic_smoothness_proved']
    assert not rigidity['scheme_level_constancy_claim']
    assert not rigidity['new_Levi_coverage']
    assert not rigidity['literature_priority_claim']
    for name,expected in rigidity['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    cubic_roots=read('cubic-root-coefficients-audit.json')
    assert cubic_roots['finite_cubic_character_records_checked']==136
    assert cubic_roots['root_checks']['mixed_cubic_monomials']==230
    assert cubic_roots['root_checks']['D4_subsystems_checked']==78
    assert cubic_roots['root_checks']['orthogonal_charge_pairs']==48
    assert cubic_roots['independent_of_stored_character_polynomial_inversion']
    assert cubic_roots['full_character_unit_quadratic_coefficients_verified']
    assert cubic_roots['scheme_elimination_inverse_verified']
    assert cubic_roots['root_checks']['nonlinear_coefficient_weight_sets']==1036
    assert cubic_roots['all_nonlinear_coefficient_simplex_or_cross_polytope_structure_verified']
    assert cubic_roots['isolated_centralizer_factors_verified']
    assert cubic_roots['individual_coefficient_boundary_formulas_verified']
    assert cubic_roots['common_fundamental_character_coordinate_frames_verified']
    assert cubic_roots['squarefree_coordinate_ideal_decompositions_verified']
    assert cubic_roots['simultaneous_coefficient_zero_schemes_reduced_proved']
    resonance=cubic_roots['linear_resonance_control']
    assert resonance['same_all_nonlinear_coefficients'] and resonance['same_constant_term']
    assert resonance['actual_reductive_Levi_representation']
    assert resonance['representation_dimension']==82 and resonance['central_trace_norm']==180
    assert resonance['both_local_ideal_inclusions_verified']
    assert resonance['all_original_derivatives_replayed_in_nonreduced_ring']
    assert resonance['flat_root_laplacian_identity_preserved']
    assert not resonance['full_ambient_Weyl_invariance_preserved']
    assert resonance['invariant_derivations_not_ordinary_quotient_derivatives']
    assert resonance['sl2_adjoint_regular_unipotent_derivative']==4
    assert not resonance['ambient_simple_adjoint_representation']
    assert not resonance['generic_adjoint_conjecture_disproved']
    assert not resonance['coefficient_boundary_structure_alone_implies_rank_bound']
    assert not cubic_roots['simultaneous_boundary_transversality_proved']
    assert not cubic_roots['generic_adjoint_smoothness_proved']
    assert not cubic_roots['new_Levi_coverage']
    assert not cubic_roots['literature_priority_claim']
    for key in ('source_sha256','raw_sha256'):
        for name,expected in cubic_roots[key].items():
            assert digest(ROOT/name)==expected,name
    e7_pencil=read('e7-symmetric-pencil-verification.json')
    e7_pencil_audit=read('e7-symmetric-pencil-audit.json')
    e7_candidate=json.loads((ROOT/'e7-quartic-boundary-exploration.json').read_text())
    assert e7_pencil['resolved_proper_standard_Levis']==101
    assert e7_pencil['remaining_cubic_embeddings']==26 and e7_pencil['remaining_cubic_families']==4
    assert e7_pencil['symmetric_pencil']['rank_stratum_bounds_verified']
    assert e7_pencil['symmetric_pencil']['regular_fiber_dimension']==8
    assert e7_pencil['symmetric_pencil']['collision_polynomials_coprime']
    assert not e7_pencil_audit['uses_Sage']
    assert e7_pencil_audit['symmetric_pencil']['full_character_independently_reconstructed']
    assert e7_pencil_audit['symmetric_pencil']['rank_one_curves_and_sign_changes']
    assert e7_pencil_audit['symmetric_pencil']['ordinary_logarithmic_hessian_diagonal']==-24
    assert not e7_pencil['companion']['generic_smoothness_claim']
    assert not e7_pencil_audit['companion']['generic_smoothness_claim']
    assert e7_pencil_audit['candidate']['original_Jacobian_ideal_reconstructed']
    assert e7_pencil_audit['candidate']['selected_leading_support_dimension']==1
    assert not e7_pencil_audit['candidate']['selected_memberships_verified']
    assert e7_candidate['status']=='exploratory'
    assert e7_candidate['computed_singular_locus_dimension']==1
    assert not e7_candidate['independent_reverse_ideal_inclusion_complete']
    assert not e7_candidate['companion_generic_smoothness_promoted']
    assert digest(ROOT/e7_candidate['archive'])==e7_candidate['archive_sha256']
    for record in (e7_pencil,e7_pencil_audit):
        assert not record['all_E7_Levis_resolved']
        assert not record['full_map_degree_claim'] and not record['literature_priority_claim']
    for record in (e7_pencil,e7_pencil_audit,e7_candidate):
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    e7_boundary=read('e7-multiplier-boundary-verification.json')
    e7_boundary_audit=read('e7-multiplier-boundary-audit.json')
    assert e7_boundary['cubic_cases']==32 and e7_boundary['cubic_classes']==7
    assert e7_boundary['resolved_proper_standard_cases']==99
    assert e7_boundary['remaining_cubic_cases']==28 and e7_boundary['remaining_cubic_classes']==5
    assert e7_boundary['common_maximal_Levi_refinement']
    assert e7_boundary['application']['regular_fiber_dimension']==16
    assert e7_boundary['application']['generic_regular_fibers_nonempty_smooth']
    assert e7_boundary['application']['original_target_derivatives_verified']
    assert not e7_boundary_audit['uses_Sage']
    assert e7_boundary_audit['weyl']['full_Weyl_orbits_reenumerated']==7
    assert e7_boundary_audit['weyl']['explicit_words_verified']==32
    assert e7_boundary_audit['weyl']['all_32_characters_reconstructed']
    assert e7_boundary_audit['algebra']['original_derivatives_preserved']
    assert e7_boundary_audit['algebra']['full_logarithmic_hessian_schur_complement_verified']
    assert e7_boundary_audit['algebra']['boundary_dimension']==1
    assert e7_boundary_audit['algebra']['torus_hessian_determinant']=='-9125752/16875'
    assert digest(ROOT/'unit-constraint-cases.json.gz')==e7_boundary['case_archive_sha256']
    for record in (e7_boundary,e7_boundary_audit):
        assert not record['all_E7_Levis_resolved']
        assert not record['full_map_degree_claim'] and not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    unit_class=read('unit-constraint-classification-verification.json')
    unit_audit=read('unit-constraint-classification-audit.json')
    assert digest(ROOT/unit_class['case_archive'])==unit_class['case_archive_sha256']
    assert [unit_class['summary'][k]['counts']['cubic'] for k in ('E6','E7','E8')]==[11,32,93]
    assert len(unit_class['cases'])==136
    assert len(unit_class['E6_Weyl_transfers'])==11
    assert unit_class['all_E6_proper_Levi_generic_regular_fibers_smooth_or_empty']
    assert not unit_audit['uses_Sage']
    assert unit_audit['all_445_proper_standard_subsets_enumerated']
    assert unit_audit['cases_audited']==136
    assert unit_audit['maximum_remaining_variables']==4
    assert unit_audit['E6_explicit_Weyl_transfers_verified']==11
    assert unit_audit['all_original_target_derivatives']
    for record in (unit_class,unit_audit):
        assert not record['recursive_reduction_classified']
        assert not record['all_adjoint_generic_smoothness_proved']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    multipliers=read('multiplier-branching-verification.json')
    multiplier_audit=read('multiplier-branching-audit.json')
    assert multipliers['helper']['central_matrix_derivative_retained']
    assert len(multipliers['helper']['negative_cases'])==2
    assert multipliers['symbolic']['no_source_boundary_removed']
    assert multipliers['quotient']['degree']==6
    assert multipliers['quotient']['every_exceptional_fiber_one_double_four_simple']
    assert multipliers['quotient']['nonzero_exceptional_targets']==10
    assert multipliers['root_point']['trace_norm']==1800
    assert multipliers['root_point']['full_rank']==220
    assert not multiplier_audit['uses_Sage']
    assert multiplier_audit['independent_character']['cover_degree']==30
    assert multiplier_audit['symbolic']['universal_nonzero_target_ideal_equality']
    assert multiplier_audit['symbolic']['only_target_inverted']
    assert multiplier_audit['root_point']['kernel_dimension']==28
    for record in (multipliers,multiplier_audit):
        assert not record['all_exceptional_Levis_resolved']
        assert not record['full_map_generic_degree_claim']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    pencil_scope=read('adjoint-pencil-scope-verification.json')
    pencil_scope_audit=read('adjoint-pencil-scope-audit.json')
    assert [r['highest_weight_node'] for r in pencil_scope['highest_root_obstructions']]==[2,1,8]
    assert pencil_scope['E7']['actual_quotient']=='3*w^2-c*w-3'
    assert pencil_scope['E7']['universal_nonzero_target_ideal_equality']
    assert pencil_scope['independent_original_root_point']['full_rank']==113
    assert pencil_scope['independent_original_root_point']['trace_norm']==432
    assert not pencil_scope_audit['uses_Sage']
    assert [r['cover_degree'] for r in pencil_scope_audit['independently_reconstructed_characters']]==[6,12,30]
    assert pencil_scope_audit['symbolic']['mu12_action_on_original_character']
    assert pencil_scope_audit['symbolic']['generic_regular_orbits']==2
    assert pencil_scope_audit['symbolic']['regular_orbit_dimension']==20
    assert pencil_scope_audit['original_root_point']['full_rank']==113
    for record in (pencil_scope,pencil_scope_audit):
        assert not record['arbitrary_affine_changes_classified']
        assert not record['all_exceptional_Levis_resolved']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    spectral=read('spectral-branching-verification.json')
    spectral_audit=read('spectral-branching-audit.json')
    assert [r['accepted'] for r in spectral['character']['shape_tests']]==[False,False,True]
    assert spectral['character']['eigenvalue_multiplicities']==[1,1]
    assert spectral['character']['spectral_target_maps_zero']
    assert all(spectral['ideals'].values())
    assert not spectral['all_exceptional_classification']
    assert not spectral['new_exceptional_case_claim']
    assert not spectral['simple_group_counterexample_claim']
    realization=spectral['reductive_realization']
    assert realization['representation_dimension']==36
    assert realization['full_differential_rank']==2
    assert realization['regular_fiber_local_dimension']==9
    assert realization['regular_fiber_tangent_dimension']==10
    assert not realization['group_is_simple']
    assert not spectral_audit['uses_Sage']
    assert all(spectral_audit['checks'].values())
    assert spectral_audit['reductive_realization']['dimension']==36
    assert spectral_audit['reductive_realization']['differential_rank']==2
    assert spectral_audit['reductive_realization']['trace_form_nondegenerate']
    for record in (spectral,spectral_audit):
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    cubic=read('cubic-levi-refinements-verification.json')
    cubic_audit=read('cubic-levi-refinements-audit.json')
    assert cubic['remaining_two_mechanism_tests_resolved']
    assert cubic['Fricke_rank_five_is_not_Springer_degree']
    assert cubic['identities']['boundary_critical_value']==-2
    assert cubic['identities']['five_relative_hessian_determinants']==[8,-32,-32,-32,-32]
    assert cubic['identities']['target_dominance_hessian_diagonal']==[2,2,2,-18,-18,-18]
    assert [r['kernel_dimension'] for r in cubic['independent_root_points']]==[10,6]
    assert not cubic_audit['uses_Sage']
    assert cubic_audit['universal_fricke']['S_pairs_checked']==15
    assert cubic_audit['universal_fricke']['rank']==5
    assert cubic_audit['universal_fricke']['rank_is_over_central_parameters_not_target_parameters']
    assert cubic_audit['product_transfer']['product_transfer_on_both_nonzero_targets']
    assert [r['full_rank'] for r in cubic_audit['native_root_checks']]==[68,72]
    for record in (cubic,cubic_audit):
        assert not record['all_exceptional_Levis_resolved']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    norm_fibers=read('branching-norm-fibers-verification.json')
    norm_audit=read('branching-norm-fibers-audit.json')
    assert len(norm_fibers['adjoint_mechanism_tests'])==3
    assert all(row['max_degree']==3 for row in norm_fibers['adjoint_mechanism_tests'])
    assert norm_fibers['E6_minimal_quadratic_control']['max_degree']==2
    assert norm_fibers['E6_minimal_quadratic_control']['dimension']==27
    assert all(norm_fibers['exact_suspension_checks'][key] for key in (
        'E6_original_gradient_ideal_equals_scalar_model',
        'E6_formula_matches_root_reconstructed_character',
        'suspension_uniform_on_nonzero_target','E6_generic_quotient_finite_smooth',
        'two_pair_suspension_free_coordinate_preserved'))
    assert norm_fibers['E6_direct_root_differential']['full_rank']==64
    assert norm_fibers['E6_direct_root_differential']['kernel_dimension']==14
    assert not norm_audit['uses_Sage']
    assert len(norm_audit['independent_character_reconstructions'])==4
    assert norm_audit['symbolic_transfer']['generic_regular_K_orbits']==2
    assert norm_audit['independent_finite_root_control']['full_rank']==64
    for record in (norm_fibers,norm_audit):
        assert not record['all_exceptional_adjoint_Levis_resolved']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    contact=read('exceptional-adjoint-contact-verification.json')
    contact_audit=read('exceptional-adjoint-contact-audit.json')
    assert len(contact['contact_gradings'])==len(contact_audit['contact_gradings'])==5
    assert [row['V_dimension'] for row in contact['contact_gradings']]==[4,14,20,32,56]
    assert [row['trace_normalization'] for row in contact_audit['contact_gradings']]==[16,36,48,72,120]
    assert not contact_audit['uses_Sage']
    assert contact['E6_nonzero_contact_fiber_regular_dimension']==32
    assert contact['E6_actual_quotient_critical_scheme']=='A^2'
    assert contact['E6_root_trace_audit']['full_differential_rank']==46
    assert contact_audit['E6_independent_point']['full_differential_rank']==46
    assert contact_audit['E6_independent_point']['kernel_dimension']==32
    assert contact['exact_ideals']['G2_generic_fiber_nonempty_finite_smooth']
    assert not contact['classical_dimension_formula_extends_unchanged']
    for record in (contact,contact_audit):
        assert not record['all_exceptional_Levis_classified']
        assert not record['literature_priority_claim']
        for name,expected in record['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    classical_adjoint=read('classical-adjoint-fibers-verification.json')
    assert not classical_adjoint['uses_Sage']
    assert classical_adjoint['classical_all_rank_smoothness_written_proof']
    assert classical_adjoint['mixed_block_nonemptiness_classified']
    assert classical_adjoint['generic_regular_dimension']=='dim K-rank G'
    assert not classical_adjoint['public_API_changed']
    assert not classical_adjoint['literature_priority_claim']
    assert len(classical_adjoint['matrix_checks'])==4
    assert sum(r['trace_form_pairings'] for r in classical_adjoint['matrix_checks'])==1425
    assert sum(r['Springer_trace_pairings'] for r in classical_adjoint['matrix_checks'])==69
    assert sum(r['full_differential_pairings'] for r in classical_adjoint['matrix_checks'])==1425
    assert sum(r['regular_unipotent_control'] for r in classical_adjoint['matrix_checks'])==1
    assert len(classical_adjoint['deformation']['symbolic_deformation_branches'])==11
    assert not classical_adjoint['deformation']['rank_sweep']
    for name,expected in classical_adjoint['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    quadratic_branching=read('quadratic-branching-fibers-verification.json')
    adjoint_matrix_audit=read('type-a-adjoint-fibers-audit.json')
    assert len(quadratic_branching['scheme_controls'])==8
    assert len(quadratic_branching['type_A_block_partitions'])==58
    assert len(quadratic_branching['dominance_checks'])==24
    assert len(quadratic_branching['classical_root_restrictions'])==347
    assert sum(not r['quadratic'] for r in quadratic_branching['classical_root_restrictions'])==37
    assert not quadratic_branching['quadratic_theorem_requires_constant_rank']
    assert quadratic_branching['generic_smoothness_includes_rank_drop_loci']
    assert quadratic_branching['type_A_nonemptiness_proof_is_all_rank']
    assert not quadratic_branching['classical_transfer_nonemptiness_classified']
    assert not quadratic_branching['cubic_control_is_actual_representation']
    assert not quadratic_branching['literature_priority_claim']
    assert not quadratic_branching['public_API_changed']
    assert not adjoint_matrix_audit['uses_Sage']
    assert len(adjoint_matrix_audit['original_adjoint_matrix_checks'])==8
    assert adjoint_matrix_audit['trace_form_pairings']==865
    assert adjoint_matrix_audit['trace_pairings']==142
    assert adjoint_matrix_audit['differential_pairings']==2882
    assert sum(r['nonsemisimple_source'] for r in adjoint_matrix_audit['original_adjoint_matrix_checks'])==1
    for row in adjoint_matrix_audit['original_adjoint_matrix_checks']:
        assert row['fiber_tangent_dimension']==2*row['double_blocks']
        assert row['full_differential_rank']==row['N']**2-1-row['fiber_tangent_dimension']
        assert row['source_centralizer_dimension']==row['N']-1
    for report in (quadratic_branching,adjoint_matrix_audit):
        for name,expected in report['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    affine_classification=read('affine-branching-classification-verification.json')
    affine_classification_audit=read('affine-branching-classification-audit.json')
    assert len(affine_classification['classical_comparisons'])==2453
    assert affine_classification['exceptional_negative_pairs']==3240
    assert affine_classification['exceptional_affine_pairs']==133
    assert affine_classification['independent_full_character_comparisons']==1462
    assert len(affine_classification['built_in_branching_checks'])==3
    assert affine_classification['E7_all_fundamentals_classified']
    assert affine_classification['E8_all_fundamentals_classified']
    assert affine_classification['nonaffinity_is_not_a_fiber_singularity_proof']
    assert affine_classification_audit['exceptional_pairs']==3373
    assert affine_classification_audit['negative_pairs']==3240
    assert affine_classification_audit['affine_pairs']==133
    assert affine_classification_audit['reflection_word_certificates']==637
    assert not affine_classification_audit['uses_Sage']
    assert affine_classification_audit['cartan_matrices_independently_specified']
    assert affine_classification_audit['norms_and_Weyl_dimensions_recomputed_with_fractions']
    expected_affine_counts={'G2':[2,1],'F4':[2,1,1,7],
                           'E6':[25,4,1,1,1,25],'E7':[6,2,1,1,1,1,34],
                           'E8':[2,1,1,1,1,1,1,8]}
    summaries={(r['type'],r['node']):r for r in affine_classification['exceptional_summary']}
    certificates={(r['type'],r['node']):r for r in affine_classification['exceptional_certificates']}
    expected_pairs={(kind,k) for kind,counts in expected_affine_counts.items()
                    for k in range(1,len(counts)+1)}
    assert set(summaries)==set(certificates)==expected_pairs
    assert len(affine_classification['exceptional_summary'])==len(expected_pairs)
    assert len(affine_classification['exceptional_certificates'])==len(expected_pairs)
    for pair in expected_pairs:
        kind,k=pair; summary=summaries[pair]; record=certificates[pair]
        assert summary['affine_count']==expected_affine_counts[kind][k-1]
        assert summary['proper_standard_Levis']==2**int(kind[1:])-1
        assert sorted(summary['affine_node_sets'])==sorted(r['levi_nodes'] for r in record['positive_restrictions'])
        assert len(summary['affine_node_sets'])==summary['affine_count']
    for report in (affine_classification,affine_classification_audit):
        for name,expected in report['source_sha256'].items():
            assert digest(ROOT/name)==expected,name
    affine_branching=read('affine-branching-fibers-verification.json')
    assert len(affine_branching['type_A_coefficient_checks'])==158
    assert len(affine_branching['type_C_coefficient_checks'])==78
    assert len(affine_branching['type_A_weight_differentials'])==5
    assert len(affine_branching['type_C_Weyl_weight_differentials'])==12
    assert all(affine_branching['algebraic_controls'].values())
    assert not affine_branching['all_fundamental_branchings_affine']
    assert not affine_branching['nonempty_generic_fibers_proved_for_all_admissible_cases']
    assert not affine_branching['exceptional_target_values_classified']
    for row in affine_branching['type_A_weight_differentials']:
        assert row['fiber_tangent_dimension']==row['m']*(row['m']-1)
        assert row['full_differential_rank']==row['N']**2-1-row['fiber_tangent_dimension']
    for row in affine_branching['type_C_Weyl_weight_differentials']:
        assert row['fiber_tangent_dimension']==2*row['p']**2
        assert row['full_differential_rank']==row['n']*(2*row['n']+1)-row['fiber_tangent_dimension']
    for name,expected in affine_branching['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    spin_excess=read('spin-fiber-excess-verification.json')
    assert len(spin_excess['original_matrix_checks'])==8
    assert len(spin_excess['independent_weyl_characters'])==6
    assert spin_excess['total_trace_pairings']==244
    assert spin_excess['total_differential_pairings']==9724
    assert spin_excess['symbolic_critical_scheme']['normal_hessian_determinant']==-4
    assert spin_excess['generic_morse_candidate_refuted']
    assert not spin_excess['generic_reducedness_refuted']
    assert not spin_excess['nonregular_fibers_classified']
    for row in spin_excess['original_matrix_checks']:
        n=row['rank']
        assert row['differential_rank']==4*n
        assert row['actual_tangent_dimension']==n*(2*n-3)
        assert row['actual_tangent_dimension']-row['proposed_fiber_dimension']==n-2
    for name,expected in spin_excess['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    symplectic_strata=read('symplectic-target-strata-verification.json')
    assert len(symplectic_strata['first_block_ideal_checks'])==136
    assert len(symplectic_strata['even_pairing_checks'])==16
    assert len(symplectic_strata['odd_residual_checks'])==16
    assert len(symplectic_strata['signed_root_counts'])==18
    assert len(symplectic_strata['primitive_matrix_checks'])==14
    assert len(symplectic_strata['odd_primitive_matrix_checks'])==2
    assert len(symplectic_strata['weyl_characters'])==4
    assert symplectic_strata['total_primitive_trace_pairings']==494
    assert not symplectic_strata['odd_rank_multiplicities_classified']
    assert not symplectic_strata['nonregular_components_classified']
    for name,expected in symplectic_strata['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    simple_branching=read('simple-group-cayley-branching-verification.json')
    assert len(simple_branching['weyl_character_comparisons'])==20
    assert len(simple_branching['block_derivatives'])==113
    assert len(simple_branching['hook_elementary_coefficients'])==78
    assert sum(x['trace_pairings'] for x in simple_branching['middle_exterior_original_matrices'])==270
    assert len(simple_branching['middle_exterior_original_matrices'])==10
    assert len(simple_branching['middle_exterior_residuals'])==13
    assert all(x['root_assignments_independently_enumerated'] for x in simple_branching['middle_exterior_residuals'] if x['k']<=7)
    assert all(simple_branching['endpoint_controls'].values())
    assert not simple_branching['exceptional_Levi_classification_claimed']
    for name,expected in simple_branching['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    critical_fiber=read('general-cayley-fiber-critical-verification.json')
    assert len(critical_fiber['linear_eliminations'])==6
    assert len(critical_fiber['type_A_potentials'])==11
    assert len(critical_fiber['type_C_perfect_pairings'])==12
    assert len(critical_fiber['defining_SLn_baselines'])==9
    assert sum(x['trace_pairings'] for x in critical_fiber['actual_tensor_matrix_checks'])==12
    assert critical_fiber['symbolic_Sp4_corrected_potential_pairings']==4
    assert critical_fiber['nonregular_SL2_control'] and critical_fiber['ramified_PGL2_control']
    assert not critical_fiber['nonregular_product_theorem_claimed']
    for name,expected in critical_fiber['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    symplectic_fiber=read('symplectic-fiber-component-verification.json')
    assert len(symplectic_fiber['unimodular_pairings'])==20
    assert len(symplectic_fiber['cyclotomic_target_checks'])==40
    assert all(x['equal_external_pairs']==0 for x in symplectic_fiber['cyclotomic_target_checks'])
    assert len(symplectic_fiber['full_differential_checks'])==4
    assert sum(x['pairings'] for x in symplectic_fiber['actual_primitive_trace_checks'])==1396
    assert symplectic_fiber['rank_one_singular_model']['symbolic_primitive_trace_pairings']==10
    assert symplectic_fiber['rank_one_singular_model']['nilpotence_order']==3
    assert not symplectic_fiber['global_reduced_germ_conjecture_proved']
    for name,expected in symplectic_fiber['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    symplectic_realization=read('symplectic-centralizer-realization-verification.json')
    assert len(symplectic_realization['coefficient_isomorphisms'])==10
    assert len(symplectic_realization['actual_primitive_cases'])==7
    assert symplectic_realization['actual_primitive_tangent_pairings']==174
    assert symplectic_realization['contraction_kernels_and_invariance_checked']
    assert symplectic_realization['trace_indices_and_exact_centralizers_checked']
    assert not symplectic_realization['full_fiber_component_theorem']
    for name,expected in symplectic_realization['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    middle_fiber=read('middle-exterior-fiber-component-verification.json')
    assert sum(x['pairings'] for x in middle_fiber['independent_exterior_derivatives'])==5419
    assert len(middle_fiber['differential_checks'])==5
    assert len(middle_fiber['exact_transverse_curve_checks'])==7
    assert len(middle_fiber['cyclotomic_collision_checks'])==15
    for x in middle_fiber['exact_transverse_curve_checks']:
        assert x['order_at_one']==(2 if x['m']%2 else 1)
    for name,expected in middle_fiber['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    irreducible_realization=read('irreducible-centralizer-realization-verification.json')
    assert len(irreducible_realization['actual_irreducible_matrix_cases'])==12
    assert irreducible_realization['exact_ambient_tangent_pairings']==500
    assert irreducible_realization['symbolic_SL6_tangent_pairings']==35
    assert irreducible_realization['group_determinants_and_exact_centralizers_checked']
    assert not irreducible_realization['uses_character_or_borel_runtime']
    for name,expected in irreducible_realization['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    character_degeneration=read('type-a-character-degeneration-verification.json')
    assert character_degeneration['exterior_matrix_trace_moments']==140
    assert character_degeneration['exact_cyclotomic_central_checks']==203
    assert character_degeneration['full_symbolic_SL4_tangent_pairings']==45
    assert len(character_degeneration['all_order_witnesses'])==16
    assert character_degeneration['explicit_tensor_square_zero_map_checked']
    for name,expected in character_degeneration['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    fiber_assessment=read('cayley-fiber-assessment-verification.json')
    assert fiber_assessment['symbolic_trace_projection_checked']
    assert fiber_assessment['square_map_coordinate_change_and_inverse_checked']
    assert fiber_assessment['zero_linear_nonzero_quadratic_moment_checked']
    assert fiber_assessment['zero_fiber_components']==2
    assert fiber_assessment['component_dimensions']==[4,4]
    assert fiber_assessment['zero_fiber_reduced']
    for name,expected in fiber_assessment['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    g2_rank=read('g2-borel-rank-verification.json')
    assert g2_rank['full_symbolic_trace_representations']==[7,14]
    assert g2_rank['centralizer_curves_checked']==12
    assert g2_rank['rational_normal_form_rank_checks']==24
    assert g2_rank['determinantal_strata_verified']
    for name,expected in g2_rank['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    root_borel=read('root-borel-verification.json')
    assert root_borel['test_methods']==8
    assert root_borel['classical_standard_matrix_cases']==12
    assert root_borel['exceptional_adjoint_trace_cases']==15
    assert root_borel['exceptional_minimal_conjugation_cases']==4
    assert root_borel['E8_dense_identity']['positive_roots']==120
    assert root_borel['E8_dense_identity']['full_adjoint_trace_agrees']
    for name,expected in root_borel['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    borel_reconstruction=read('borel-reconstruction-verification.json')
    assert borel_reconstruction['test_methods']==8
    assert borel_reconstruction['symbolic_full_representation_projection_cases']==9
    assert borel_reconstruction['low_level_additional_defining_matrix_checks']==['SL2','SL4']
    assert not borel_reconstruction['representation_matrices_used_by_runtime']
    assert borel_reconstruction['regular_extension_checked']
    assert borel_reconstruction['changed_collapse_loci_checked']
    for name,expected in borel_reconstruction['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    borel_coordinates=read('borel-coordinates-verification.json')
    assert borel_coordinates['test_methods']==6
    assert borel_coordinates['independent_full_trace_gram_dimensions']==[8,10]
    assert borel_coordinates['symbolic_all_parameter_matrix_identities']
    assert borel_coordinates['symbolic_collapse_stratum_ranks']
    assert borel_coordinates['exact_empty_and_positive_dimensional_fibers']
    for name,expected in borel_coordinates['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    literature=read('cayley-literature-verification.json')
    assert literature['bibliographic_entries']==8
    assert literature['notebook']['saved_SInv_rank']==7
    assert literature['proposed_Part_I_Example_3_scalar_correction']['corrected_identity_derivative']==1
    for name,expected in literature['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    exterior=read('exterior-regularity-verification.json')
    assert exterior['regular_test_point_projection_checks']==66
    assert exterior['rational_interior_counterexamples']==45
    assert exterior['nonsemisimple_full_matrix_trace_checks']==8
    assert not exterior['novelty_claim'] and not exterior['runtime_API_changed']
    for name,expected in exterior['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    type_a_regularity=read('type-a-regularity-verification.json')
    assert type_a_regularity['character_cases']==149
    assert type_a_regularity['algebraic_witness_polynomials']==134
    assert type_a_regularity['symmetric_power_cases']==50
    assert type_a_regularity['symmetric_power_complete_factor_checks']==50
    assert type_a_regularity['nontrivial_central_quotient_descent_cases']==82
    assert not type_a_regularity['novelty_claim'] and not type_a_regularity['runtime_API_changed']
    for name,expected in type_a_regularity['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    type_c_regularity=read('type-c-fundamental-regularity-verification.json')
    assert type_c_regularity['primitive_kernel_cases']==36
    assert type_c_regularity['exact_contraction_weight_blocks']==19676
    assert type_c_regularity['quadratic_algebra_counterexamples']==28
    assert type_c_regularity['rational_defining_counterexamples']==7
    assert not type_c_regularity['novelty_claim'] and not type_c_regularity['runtime_API_changed']
    for name,expected in type_c_regularity['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    all_c_regularity=read('type-c-regularity-verification.json')
    assert all_c_regularity['character_cases']==79
    assert all_c_regularity['universal_Qi_trace_checks']==73
    assert all_c_regularity['nonfundamental_exact_algebraic_witnesses']==68
    assert all_c_regularity['same_central_parity_reducible_checks']==8
    assert not all_c_regularity['novelty_claim'] and not all_c_regularity['runtime_API_changed']
    for name,expected in all_c_regularity['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    type_b_regularity=read('type-b-regularity-verification.json')
    assert type_b_regularity['character_cases']==55
    assert type_b_regularity['positive_and_witness_alternant_checks']==110
    assert type_b_regularity['universal_Qi_trace_checks']==55
    assert type_b_regularity['integer_weight_cases']==34
    assert type_b_regularity['spinorial_cases']==21
    assert type_b_regularity['exterior_monomial_character_checks']==13
    assert type_b_regularity['independent_spin_sign_cube_checks']==7
    assert type_b_regularity['same_central_parity_reducible_checks']==8
    assert type_b_regularity['closed_spin_formula_checks']==8
    assert type_b_regularity['rank_one_spin_positive_control']
    assert not type_b_regularity['novelty_claim'] and not type_b_regularity['runtime_API_changed']
    for name,expected in type_b_regularity['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    mixed_regularity=read('mixed-regularity-verification.json')
    assert mixed_regularity['Sp_tensor_character_cases']==9
    assert mixed_regularity['Sp_regular_torus_trace_checks']==27
    assert mixed_regularity['independent_tensor_binomial_models']==9
    assert mixed_regularity['Sp_nonsemisimple_boundary_ranks']==5
    assert mixed_regularity['Sp_regular_nonsemisimple_good_open_checks']==15
    assert mixed_regularity['full_exterior_matrix_boundary_ranks']==2
    assert mixed_regularity['full_exterior_matrix_Lie_basis_checks']==31
    assert mixed_regularity['mixed_B_direct_finite_algebra_witnesses']==18
    assert mixed_regularity['mixed_B_boundary_linear_factor_exclusions']==372
    assert mixed_regularity['B2_torus_positive_control']['dimension']==16
    assert not mixed_regularity['novelty_claim'] and not mixed_regularity['runtime_API_changed']
    for name,expected in mixed_regularity['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    harmonic_localization=read('symplectic-harmonic-localization-verification.json')
    assert harmonic_localization['full_matrix_differential_checks']==48
    assert harmonic_localization['independent_Cayley_coordinate_Jacobian_checks']==24
    assert harmonic_localization['regular_critical_corank_one_checks']==12
    assert harmonic_localization['actual_exterior_tensor_differential_checks']==12
    assert harmonic_localization['low_degree_harmonic_checks']['dimensions']==[1,10,54,210]
    assert harmonic_localization['low_degree_harmonic_checks']['actual_exterior_trace_pullback_evaluations']==1650
    assert harmonic_localization['explicit_necessity_conjugacy_pairs']==15
    harmonic_audit=harmonic_localization['torus_boundary_audit']
    assert len(harmonic_audit['regular_torus_subspace_checks'])==48
    assert len(harmonic_audit['boundary_strata_checks'])==36
    assert len(harmonic_audit['actual_boundary_trace_checks'])==9
    assert len(harmonic_audit['rank_one_boundary_controls'])==3
    assert harmonic_localization['known_input_corollary']
    assert not harmonic_localization['novelty_claim'] and not harmonic_localization['runtime_API_changed']
    for name,expected in harmonic_localization['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    symplectic_degree=read('symplectic-cayley-degrees-verification.json')
    assert symplectic_degree['scalar_fiber_cases']==80
    assert symplectic_degree['independent_weight_fiber_cases']==24
    assert symplectic_degree['original_Laurent_Groebner_fibers']==10
    assert symplectic_degree['singular_target_degree_cases']==108
    assert symplectic_degree['full_exterior_matrix_fiber_cases']==9
    assert symplectic_degree['polynomial_character_extension']
    assert symplectic_degree['surjectivity_proved']
    assert symplectic_degree['finite_flat_over_invertible_targets']
    assert not symplectic_degree['novelty_claim']
    assert not symplectic_degree['externally_reviewed'] and not symplectic_degree['runtime_API_changed']
    for name,expected in symplectic_degree['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    spin_degree=read('spin-cayley-degrees-verification.json')
    assert spin_degree['scalar_fiber_cases']==48
    assert spin_degree['independent_weight_fiber_cases']==16
    assert spin_degree['original_Spin_Laurent_Groebner_fibers']==11
    assert spin_degree['original_SO_Laurent_Groebner_fibers']==5
    assert spin_degree['finite_flat_over_explicit_open']
    assert not spin_degree['novelty_claim'] and not spin_degree['externally_reviewed']
    assert not spin_degree['runtime_API_changed']
    for name,expected in spin_degree['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    spin_polynomial=read('spin-polynomial-cayley-degrees-verification.json')
    assert spin_polynomial['scalar_fiber_cases']==48
    assert spin_polynomial['independent_weight_fiber_cases']==24
    assert spin_polynomial['original_Spin_Laurent_Groebner_fibers']==14
    assert spin_polynomial['original_SO_Laurent_Groebner_fibers']==5
    assert spin_polynomial['collapsed_character_level_checks']==6
    assert spin_polynomial['finite_flat_over_explicit_open']
    assert spin_polynomial['SO_descent_iff_F_even']
    assert not spin_polynomial['novelty_claim'] and not spin_polynomial['externally_reviewed']
    assert not spin_polynomial['runtime_API_changed']
    for name,expected in spin_polynomial['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    spin_geometry=read('spin-cayley-geometry-verification.json')
    assert spin_geometry['maximal_finite_locus_uniform_written_proof']
    assert spin_geometry['full_global_Jacobian_written_proof']
    assert spin_geometry['weight_root_Jacobian_checks']==48
    assert spin_geometry['Clifford_boundary_Jacobian_checks']==33
    assert spin_geometry['finite_algebra_ramification_checks']==32
    assert spin_geometry['nonfinite_boundary_target_checks']==8
    assert spin_geometry['regular_semisimple_critical_invariant_checks']==42
    assert spin_geometry['Spin4_discriminant_control']['primitive_element_artifact']=='p-v'
    assert not spin_geometry['novelty_claim'] and not spin_geometry['externally_reviewed']
    assert not spin_geometry['runtime_API_changed']
    for name,expected in spin_geometry['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    half_spin=read('half-spin-cayley-verification.json')
    assert half_spin['uniform_Pfaffian_formula']
    assert half_spin['complete_chart_zero_fiber_written_proof']
    assert half_spin['chart_zero_fiber_smooth_written_proof']
    assert half_spin['global_boundary_dimension_lower_bound']
    assert half_spin['highest_weight_moment_checks']==12
    assert half_spin['full_chiral_Clifford_checks']==18
    assert half_spin['exact_zero_weight_checks']==42
    assert half_spin['full_zero_differential_rank_checks']==12
    assert half_spin['boundary_family_checks']==10
    assert half_spin['central_kernel_checks']==16
    assert half_spin['Spin10']==dict(chart_noncentral_classes_per_half_spin=4,
        chart_class_dimension=20,differential_rank=25,boundary_zero_family_dimension=33,
        generic_degree_determined=False,full_zero_fiber_classified=False)
    assert all(row['Spin8_degree']==32 for row in half_spin['triality_controls'])
    assert all(row['Spin6_original_torus_fiber_length']==4 for row in half_spin['Spin6_fibers'])
    assert not half_spin['novelty_claim'] and not half_spin['externally_reviewed']
    assert not half_spin['runtime_API_changed']
    for name,expected in half_spin['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    half_boundary=read('half-spin-boundary-verification.json')
    assert half_boundary['all_rank_underlying_zero_fiber_written_proof']
    assert half_boundary['all_rank_dimension_written_proof']
    assert half_boundary['Spin10_zero_fiber_dimension']==33
    assert half_boundary['semisimple_weight_checks']==322
    assert half_boundary['partition_sl2_checks']==3266
    assert half_boundary['direct_matrix_centralizer_checks']==34
    assert half_boundary['GL_partition_bounds']==683
    assert half_boundary['full_chiral_boundary_checks']==112
    assert half_boundary['finite_rank_dimension_checks']==17
    assert len(half_boundary['Spin10_Jordan_strata'])==7
    assert sorted(row['stratum_dimension'] for row in half_boundary['Spin10_Jordan_strata'])==[24,24,30,30,32,32,33]
    assert not half_boundary['boundary_scheme_reduced_proved']
    assert not half_boundary['boundary_irreducible_components_classified']
    assert not half_boundary['generic_degree_determined']
    assert not half_boundary['novelty_claim'] and not half_boundary['externally_reviewed']
    assert not half_boundary['runtime_API_changed']
    for name,expected in half_boundary['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    half_degree=read('half-spin-degree-verification.json')
    assert half_degree['Spin10_half_spin_generic_degree']==148
    assert half_degree['characteristic_zero_toric_written_proof']
    assert half_degree['toric_intersection_number']==208
    assert half_degree['toric_boundary_length']==60
    assert half_degree['actual_Spin_character_lattice']
    assert half_degree['group_torus_quotient_degrees_agree']
    assert len(half_degree['polytope_checks'])==6
    assert half_degree['boundary_algebra']['facet_critical_algebra_length']==6
    assert len(half_degree['boundary_algebra']['symbolic_local_checks'])==6
    assert half_degree['finite_field_fiber_checks']==12
    assert len(half_degree['finite_field_original_fibers'])==12
    assert {r['model'] for r in half_degree['finite_field_original_fibers']}=={
        'actual_P_D5_basis','original_Spin_cover'}
    assert all(r['length']==148 and r['squarefree_multiplication_characteristic_polynomial']
               for r in half_degree['finite_field_original_fibers'])
    assert len(half_degree['separate_exceptional_specializations'])==2
    assert all(r['length']==148 and r['geometric_distinct_point_count']==147
               and r['nonreduced_cyclic_algebra_certified']
               for r in half_degree['separate_exceptional_specializations'])
    assert not half_degree['all_rank_degree_formula_proved']
    assert not half_degree['novelty_claim'] and not half_degree['externally_reviewed']
    assert not half_degree['runtime_API_changed']
    for name,expected in half_degree['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    half_rank=read('half-spin-rank-degrees-verification.json')
    assert half_rank['all_rank_characteristic_zero_written_proof']
    assert half_rank['deeper_face_noncontribution_written_proof']
    assert half_rank['rank_minimum']==3
    assert half_rank['Spin12_degree']==536
    assert half_rank['Spin12_total_intersection']==1376
    assert half_rank['Spin12_boundary_correction']==840
    assert half_rank['even_rank_representation_kernel_order']==2
    assert half_rank['finite_fiber_checks']==4
    assert half_rank['squarefree_fiber_certificates']==4
    assert half_rank['additional_completed_rank7_fiber']['fiber_length']==1692
    assert all(r['length']==536 for r in half_rank['independent_fibers'])
    assert len(half_rank['combinatorial_checks'])==38
    assert half_rank['local_checks']['Spin12_curve_contribution']==360
    assert half_rank['local_checks']['Spin12_triple_contribution']==480
    assert not half_rank['novelty_claim'] and not half_rank['externally_reviewed']
    assert not half_rank['runtime_API_changed']
    for name,expected in half_rank['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    general_degree=read('general-cayley-degree-verification.json')
    assert general_degree['E6_minimal']['degree']==1080
    assert general_degree['E6_minimal']['boundary_correction']==0
    assert general_degree['E6_minimal']['weight_lattice_index']==1
    assert [r['generic_degree'] for r in general_degree['generic_character_examples']]==[8,6,8,12,14]
    assert all(r['base_field']=='QQ(b1,b2)' for r in general_degree['generic_character_examples'])
    assert not general_degree['finite_field_inference']
    assert not general_degree['novelty_claim'] and not general_degree['external_review']
    for name,expected in general_degree['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    e6_deformation=read('e6-adjoint-deformation-verification.json')
    assert e6_deformation['boundary_correction']==1890
    assert e6_deformation['root_torus_degree']==918
    assert e6_deformation['simply_connected_degree']==2754
    assert e6_deformation['geometry']['facet_cover_degree']==2
    assert e6_deformation['geometry']['root_layers']==[16,40,16]
    assert e6_deformation['local']['middle_layer_full_gradient_zero']
    assert e6_deformation['local']['generic_crossing_escape_count']==0
    assert not e6_deformation['Segre_formula_used_to_derive_result']
    assert not e6_deformation['novelty_claim'] and not e6_deformation['external_review']
    for name,expected in e6_deformation['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    f4_degree=read('f4-minimal-degree-verification.json')
    assert f4_degree['geometry']['group_torus_and_quotient_degree']==96
    assert f4_degree['geometry']['boundary_correction']==0
    assert f4_degree['geometry']['normalized_cone_volumes']==[4]*24
    assert f4_degree['geometry']['noncentral_zero_witness']['conjugacy_class_dimension']==16
    assert f4_degree['finite_flat_torus_and_quotient_written_corollary']
    assert f4_degree['full_group_map_not_finite']
    assert len(f4_degree['independent_fibers'])==4
    assert all(r['saturated_fiber_length']==96 and r['reduced_by_squarefree_multiplication'] for r in f4_degree['independent_fibers'])
    assert f4_degree['same_target_multiplication_polynomials_agree']
    for prime in [32003,32009]:
        rows=[r for r in f4_degree['independent_fibers'] if r['prime']==prime]
        assert len(rows)==2 and rows[0]['characteristic_polynomial_coefficients']==rows[1]['characteristic_polynomial_coefficients']
    assert not f4_degree['finite_fibers_used_as_proof_of_generic_degree']
    assert not f4_degree['novelty_claim'] and not f4_degree['externally_reviewed']
    for name,expected in f4_degree['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    f4_mixture=read('f4-mixture-degree-verification.json')
    assert f4_mixture['formula']=={'b=0':96,'b>0,a=b':240,'b>0,a=5b':336,'b>0,other':384}
    assert f4_mixture['geometry']['primitive_facet_height']==2
    assert f4_mixture['geometry']['facet_cone_volumes']==[16]*24
    assert f4_mixture['geometry']['nonnegative_exceptional_ratios']==['a=b','a=5*b']
    assert len(f4_mixture['independent_fibers'])==20
    for prime in [32003,32009]:
        for a,b,expected_degree in [(0,1,384),(1,1,240),(5,1,336),(2,1,384),(1,0,96)]:
            rows=[r for r in f4_mixture['independent_fibers'] if (r['a'],r['b'],r['prime'])==(a,b,prime)]
            assert len(rows)==2 and {r['model'] for r in rows}=={'weight_lattice','product_character'}
            assert all(r['saturated_fiber_length']==expected_degree and r['reduced_by_squarefree_multiplication'] for r in rows)
            assert rows[0]['characteristic_polynomial_coefficients']==rows[1]['characteristic_polynomial_coefficients']
    assert not f4_mixture['finite_fibers_used_as_proof_of_generic_degree']
    assert not f4_mixture['novelty_claim'] and not f4_mixture['externally_reviewed']
    for name,expected in f4_mixture['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    bc_mixture=read('bc-mixture-degrees-verification.json')
    assert len(bc_mixture['independent_fibers'])==68
    assert len(bc_mixture['geometry']['polytope_and_character_checks'])==10
    assert bc_mixture['spin_degree_multiplier']==2
    assert bc_mixture['pure_C_adjoint_form_degree_divisor']==2
    assert bc_mixture['exceptional_B_torus_and_quotient_nonfinite_by_zero_curves']
    assert all(r['exceptional_zero_curves_in_both_constructions'] for r in bc_mixture['geometry']['sign_checks'])
    grouped={}
    for row in bc_mixture['independent_fibers']:
        key=(row['kind'],row['rank'],row['a'],row['b'],row['prime'])
        grouped.setdefault(key,[]).append(row)
        assert row['reduced_by_squarefree_multiplication']
    assert len(grouped)==34
    for (kind,n,a,b,prime),rows in grouped.items():
        assert len(rows)==2 and {r['model'] for r in rows}=={'weight_lattice','trace_coordinates'}
        assert rows[0]['characteristic_polynomial_coefficients']==rows[1]['characteristic_polynomial_coefficients']
        if b==0: expected_degree=2**n
        elif kind=='C': expected_degree=4**n
        else:
            from math import comb
            expected_degree=2**n*(2**n-n)-sum(2*n*comb(n-1,k) for k in range(n) if a==(4*k-2*n+1)*b)
        assert all(r['saturated_fiber_length']==expected_degree for r in rows)
    assert [r['saturated_fiber_length'] for r in bc_mixture['isogeny_controls']]==[68,8]
    assert all(r['reduced_by_squarefree_multiplication'] for r in bc_mixture['isogeny_controls'])
    assert not bc_mixture['finite_fibers_used_as_proof_of_generic_degree']
    assert not bc_mixture['novelty_claim'] and not bc_mixture['externally_reviewed']
    for name,expected in bc_mixture['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    regularity_classification=read('regularity-classification-verification.json')
    assert regularity_classification['all_simply_connected_rows_classified']
    assert regularity_classification['central_quotients_corollary']
    assert regularity_classification['root_system_cases']==32
    assert regularity_classification['root_lattice_weight_lattice_checks']==1942
    assert regularity_classification['abstract_endpoint_patterns']==dict(root_string_paths=2,fork_orientations=8,passing_fork_orientations=0)
    assert regularity_classification['A_independent_GT_Freudenthal_cases']==48
    assert regularity_classification['A_positive_reducible_trace_controls']==8
    assert regularity_classification['A_mixed_exact_witnesses']==4
    assert regularity_classification['remaining_classical_exceptional_exact_witnesses']==30
    assert regularity_classification['positive_root_wall_trace_checks']==2146
    assert regularity_classification['rank_one_determinant_character_controls']==8
    assert not regularity_classification['novelty_claim'] and not regularity_classification['externally_reviewed']
    assert not regularity_classification['runtime_API_changed']
    for name,expected in regularity_classification['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
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
    assert classical_equivariant['ordinary_smith']['all_rank_closed_formula']
    assert classical_equivariant['ordinary_smith']['half_block_and_API_ranks']==128
    assert classical_equivariant['ordinary_smith']['original_inclusion_ranks']==64
    assert classical_equivariant['ordinary_smith']['exact_2_local_Gram_certificates']==31
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
        'SYMPLECTIC_K_ORDINARY_SMITH.md','symplectic-ordinary-verification.json',
        'tools/verify_symplectic_ordinary.py',
        'CAYLEY_LITERATURE_AUDIT.md','cayley-literature-verification.json',
        'tools/verify_cayley_literature.py','experiments/cayley-literature-source-record.json',
        'EXTERIOR_REGULARITY.md','exterior-regularity-verification.json',
        'tools/verify_exterior_regularity.py',
        'TYPE_A_REGULARITY.md','type-a-regularity-verification.json',
        'tools/verify_type_a_regularity.py',
        'TYPE_C_FUNDAMENTAL_REGULARITY.md','type-c-fundamental-regularity-verification.json',
        'tools/verify_type_c_fundamental_regularity.py',
        'TYPE_C_REGULARITY.md','type-c-regularity-verification.json',
        'tools/verify_type_c_regularity.py',
        'TYPE_B_REGULARITY.md','type-b-regularity-verification.json',
        'tools/verify_type_b_regularity.py',
        'MIXED_REGULARITY.md','mixed-regularity-verification.json',
        'tools/verify_mixed_regularity.py',
        'SYMPLECTIC_HARMONIC_LOCALIZATION.md','symplectic-harmonic-localization-verification.json',
        'tools/verify_symplectic_harmonic_localization.py',
        'SYMPLECTIC_CAYLEY_DEGREES.md','symplectic-cayley-degrees-verification.json',
        'tools/verify_symplectic_cayley_degree.py',
        'SPIN_CAYLEY_DEGREES.md','spin-cayley-degrees-verification.json',
        'tools/verify_spin_cayley_degree.py',
        'SPIN_POLYNOMIAL_CAYLEY_DEGREES.md','spin-polynomial-cayley-degrees-verification.json',
        'tools/verify_spin_polynomial_degree.py',
        'SPIN_CAYLEY_GEOMETRY.md','spin-cayley-geometry-verification.json',
        'tools/verify_spin_cayley_geometry.py',
        'HALF_SPIN_CAYLEY.md','half-spin-cayley-verification.json',
        'tools/verify_half_spin_cayley.py',
        'HALF_SPIN_BOUNDARY.md','half-spin-boundary-verification.json',
        'tools/verify_half_spin_boundary.py',
        'HALF_SPIN_DEGREE.md','half-spin-degree-verification.json',
        'tools/verify_half_spin_degree.py',
        'HALF_SPIN_RANK_DEGREES.md','half-spin-rank-degrees-verification.json',
        'tools/verify_half_spin_rank_degrees.py',
        'GENERAL_CAYLEY_DEGREE.md','general-cayley-degree-verification.json',
        'E6_ADJOINT_DEFORMATION.md','e6-adjoint-deformation-verification.json',
        'F4_MINIMAL_DEGREE.md','f4-minimal-degree-verification.json',
        'tools/verify_f4_minimal_degree.py',
        'F4_MIXTURE_DEGREES.md','f4-mixture-degree-verification.json',
        'tools/verify_f4_mixture_degree.py',
        'BC_MIXTURE_DEGREES.md','bc-mixture-degrees-verification.json',
        'tools/verify_bc_mixture_degrees.py',
        'tools/verify_e6_adjoint_deformation.py',
        'tools/verify_general_cayley_degree.py',
        'tools/explore_half_spin_rank6_product.py','experiments/half-spin-rank7-fiber.json',
        'REGULARITY_CLASSIFICATION.md','regularity-classification-verification.json',
        'tools/verify_regularity_classification.py',
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
    paths += [ROOT/'borel-coordinates-verification.json']
    paths += [ROOT/'borel-reconstruction-verification.json']
    paths += [ROOT/'root-borel-verification.json']
    paths += [ROOT/'g2-borel-rank-verification.json']
    paths += [ROOT/'cayley-fiber-assessment-verification.json']
    paths += [ROOT/'type-a-character-degeneration-verification.json']
    paths += [ROOT/'irreducible-centralizer-realization-verification.json']
    paths += [ROOT/'middle-exterior-fiber-component-verification.json']
    paths += [ROOT/'symplectic-centralizer-realization-verification.json']
    paths += [ROOT/'symplectic-fiber-component-verification.json']
    paths += [ROOT/'general-cayley-fiber-critical-verification.json']
    paths += [ROOT/'simple-group-cayley-branching-verification.json']
    paths += [ROOT/'symplectic-target-strata-verification.json']
    paths += [ROOT/'spin-fiber-excess-verification.json']
    paths += [ROOT/'affine-branching-fibers-verification.json']
    paths += [ROOT/name for name in ('AFFINE_BRANCHING_CLASSIFICATION.md',
        'affine-branching-classification-verification.json','affine-branching-classification-audit.json',
        'tools/verify_affine_branching_classification.py','tools/audit_affine_branching_certificates.py')]
    paths += [ROOT/name for name in ('QUADRATIC_BRANCHING_FIBERS.md',
        'quadratic-branching-fibers-verification.json','type-a-adjoint-fibers-audit.json',
        'tools/verify_quadratic_branching_fibers.py','tools/audit_type_a_adjoint_fibers.py')]
    paths += [ROOT/name for name in ('CLASSICAL_ADJOINT_FIBERS.md',
        'classical-adjoint-fibers-verification.json','tools/verify_classical_adjoint_fibers.py')]
    paths += [ROOT/name for name in ('EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md',
        'exceptional-adjoint-contact-verification.json','exceptional-adjoint-contact-audit.json',
        'tools/verify_exceptional_adjoint_contact.py','tools/audit_exceptional_adjoint_contact.py')]
    paths += [ROOT/name for name in ('BRANCHING_NORM_FIBERS.md',
        'branching-norm-fibers-verification.json','branching-norm-fibers-audit.json',
        'tools/verify_branching_norm_fibers.py','tools/audit_branching_norm_fibers.py')]
    paths += [ROOT/name for name in ('CUBIC_LEVI_REFINEMENTS.md',
        'cubic-levi-refinements-verification.json','cubic-levi-refinements-audit.json',
        'tools/verify_cubic_levi_refinements.py','tools/audit_cubic_levi_refinements.py')]
    paths += [ROOT/name for name in ('SPECTRAL_BRANCHING_FIBERS.md',
        'spectral-branching-verification.json','spectral-branching-audit.json',
        'tools/verify_spectral_branching_fibers.py','tools/audit_spectral_branching_fibers.py')]
    paths += [ROOT/name for name in ('ADJOINT_PENCIL_SCOPE.md',
        'adjoint-pencil-scope-verification.json','adjoint-pencil-scope-audit.json',
        'tools/verify_adjoint_pencil_scope.py','tools/audit_adjoint_pencil_scope.py')]
    paths += [ROOT/name for name in ('MULTIPLIER_BRANCHING_FIBERS.md',
        'multiplier-branching-verification.json','multiplier-branching-audit.json',
        'tools/linear_character_elimination.py','tools/verify_multiplier_branching_fibers.py',
        'tools/audit_multiplier_branching_fibers.py')]
    paths += [ROOT/name for name in ('EXCEPTIONAL_QUARTIC_REDUCTION.md',
        'unit-constraint-classification-verification.json','unit-constraint-classification-audit.json',
        'unit-constraint-cases.json.gz','tools/unit_constraint_systems.py',
        'tools/verify_unit_constraint_classification.py','tools/audit_unit_constraint_classification.py')]
    paths += [ROOT/name for name in ('E7_MULTIPLIER_BOUNDARY.md',
        'e7-multiplier-boundary-verification.json','e7-multiplier-boundary-audit.json',
        'tools/verify_e7_multiplier_boundary.py','tools/audit_e7_multiplier_boundary.py')]
    paths += [ROOT/name for name in ('E7_SYMMETRIC_PENCIL.md',
        'e7-symmetric-pencil-verification.json','e7-symmetric-pencil-audit.json',
        'e7-quartic-boundary-exploration.json','e7-quartic-boundary-candidate.json.gz',
        'tools/verify_e7_symmetric_pencil.py','tools/audit_e7_symmetric_pencil.py',
        'tools/explore_e7_cayley_quartic.py')]
    paths += [ROOT/name for name in ('CAYLEY_FIBER_REVIEW.md',
        'tools/audit_cayley_fiber_review.py','cayley-fiber-review-audit.json')]
    paths += [ROOT/name for name in ('BINARY_BRANCHING_FIBERS.md',
        'tools/verify_binary_branching_fibers.py','tools/audit_binary_branching_fibers.py',
        'binary-branching-verification.json','binary-branching-audit.json')]
    paths += [ROOT/'../articles/springer-cayley-exposition'/name for name in ('README.md','main.tex','verify_note.py',
        'source-manifest.json','verification.json',
        'output/pdf/Springer_Cayley_maps_progress_and_open_questions.pdf')]
    paths += [ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md',ROOT/'AGENTS.md']
    paths += [ROOT/'tools/audit_uniform_adjoint_obstruction.py',
              ROOT/'uniform-adjoint-obstruction-audit.json']
    paths += [ROOT/name for name in ('tools/search_adjoint_nilpotents.py',
        'tools/audit_adjoint_nilpotent_search.py','adjoint-nilpotent-search-audit.json',
        'adjoint-nilpotent-source-screen.json','adjoint-nilpotent-specialized-screen.json')]
    paths += [ROOT/name for name in ('tools/explore_adjoint_central_boundary.py',
        'tools/verify_adjoint_boundary_exclusions.py','tools/audit_adjoint_boundary_exclusions.py',
        'adjoint-simple-boundaries.json','adjoint-central-boundary-exploration.json',
        'adjoint-central-boundary-lift-forward.json','adjoint-central-membership-forward.json.gz',
        'adjoint-boundary-dimension-certificates.json.gz','adjoint-boundary-exclusions-audit.json',
        'adjoint-boundary-exclusions-verification.json')]
    paths += [ROOT/name for name in ('tools/test_adjoint_central_flatness.py',
        'tools/audit_adjoint_central_flatness.py','adjoint-central-flatness-audit.json',
        'adjoint-central-flatness-screen.json','adjoint-central-flatness-curve-screen.json',
        'adjoint-central-flatness-cover-intersections.json',
        'adjoint-central-boundary-lift-reverse.json','adjoint-central-membership-reverse.json.gz')]
    paths += [ROOT/name for name in ('tools/verify_adjoint_quasifinite.py',
        'tools/audit_adjoint_quasifinite.py','adjoint-quasifinite-verification.json',
        'adjoint-quasifinite-audit.json')]
    paths += [ROOT/name for name in ('tools/audit_adjoint_proof_paths.py',
        'tools/verify_adjoint_proof_paths.py','adjoint-proof-paths-audit.json',
        'adjoint-proof-paths-verification.json')]
    paths += [ROOT/name for name in ('tools/audit_adjoint_radial_schur.py',
        'tools/verify_adjoint_radial_schur.py','adjoint-radial-schur-audit.json',
        'adjoint-radial-schur-verification.json')]
    paths += [ROOT/name for name in ('tools/audit_adjoint_singular_radial.py',
        'tools/verify_adjoint_singular_radial.py','adjoint-singular-radial-audit.json',
        'adjoint-singular-radial-verification.json')]
    paths += [ROOT/name for name in ('tools/audit_adjoint_character_discriminant.py',
        'tools/verify_adjoint_character_discriminant.py','adjoint-character-discriminant-audit.json',
        'adjoint-character-discriminant-verification.json')]
    paths += [ROOT/name for name in ('tools/audit_type_a_incidence_nilradical.py',
        'tools/verify_type_a_incidence_nilradical.py','type-a-incidence-nilradical-audit.json',
        'type-a-incidence-nilradical-verification.json')]
    paths += [ROOT/name for name in ('tools/audit_adjoint_annihilator_routes.py',
        'tools/verify_adjoint_annihilator_routes.py','adjoint-annihilator-routes-audit.json',
        'adjoint-annihilator-routes-verification.json')]
    paths += [ROOT/name for name in ('tools/audit_aligned_critical_annihilator.py',
        'tools/verify_aligned_critical_annihilator.py','aligned-critical-annihilator-audit.json',
        'aligned-critical-annihilator-verification.json')]
    paths += [ROOT/name for name in ('SHARED_SQUARE_COUPLINGS.md',
        'tools/audit_shared_square_couplings.py','tools/verify_shared_square_couplings.py',
        'shared-square-couplings-audit.json','shared-square-couplings-verification.json')]
    paths += [ROOT/name for name in ('WEYL_JET_OBSTRUCTION.md',
        'tools/audit_weyl_jet_obstruction.py','tools/verify_weyl_jet_obstruction.py',
        'weyl-jet-obstruction-audit.json','weyl-jet-obstruction-verification.json')]
    paths += [ROOT/name for name in ('E7_FINAL_FAMILIES.md',
        'tools/audit_e7_final_families.py','tools/verify_e7_final_families.py',
        'e7-final-families-audit.json','e7-final-families-verification.json')]
    paths += [ROOT/name for name in ('E8_WEYL_TRANSFERS.md',
        'tools/verify_e8_weyl_transfers.py','tools/audit_e8_weyl_transfers.py',
        'e8-weyl-transfers-verification.json','e8-weyl-transfers-audit.json')]
    paths += [ROOT/name for name in ('E8_A1_TRIPLE.md',
        'tools/verify_e8_a1_triple.py','tools/audit_e8_a1_triple.py',
        'e8-a1-triple-verification.json','e8-a1-triple-audit.json')]
    paths += [ROOT/name for name in ('E8_A2_A1_PAIR.md',
        'tools/verify_e8_a2_a1_pair.py','tools/audit_e8_a2_a1_pair.py',
        'e8-a2-a1-pair-verification.json','e8-a2-a1-pair-audit.json')]
    paths += [ROOT/name for name in ('E8_FINAL_CASES.md','e8-final-character-models.json',
        'tools/verify_e8_final_cases.py','tools/audit_e8_final_cases.py',
        'tools/verify_e8_a1_four_certificate.py','e8-final-cases-verification.json',
        'e8-final-cases-audit.json','e8-a1-four-32003-certificate.json.gz',
        'e8-a1-four-32009-certificate.json.gz','e8-a1-four-32003-verification.json',
        'e8-a1-four-32009-verification.json','tools/explore_e8_last_affine.py',
        'e8-last-case-exploration.json')]
    paths += [ROOT/name for name in ('E8_LAST_ADJOINT_FIBER.md',
        'tools/verify_e8_last_case.py','tools/audit_e8_last_case.py',
        'tools/e8_last_f4_driver.c','tools/build_e8_last_certificate_engine.sh',
        'e8-last-32003-certificate.json.gz','e8-last-32009-certificate.json.gz',
        'e8-last-32003-verification.json','e8-last-32009-verification.json',
        'e8-last-case-audit.json')]
    paths += [ROOT/name for name in ('tools/verify_e7_companion_membership.py',
        'tools/audit_e7_companion_geometry.py','tools/explore_e7_homogeneous_certificate.py',
        'e7-homogeneous-membership-certificate.json.gz',
        'e7-companion-membership-verification.json','e7-companion-geometry-audit.json')]
    paths += [ROOT/name for name in ('G2_ADJOINT_FIBERS.md',
        'tools/verify_g2_adjoint_fibers.py','g2-adjoint-fibers-verification.json')]
    paths += [ROOT/name for name in ('F4_ADJOINT_FIBERS.md','tools/explore_f4_completion.py',
        'tools/verify_f4_adjoint_fibers.py','tools/audit_f4_adjoint_fibers.py',
        'f4-adjoint-fibers-verification.json','f4-adjoint-fibers-audit.json')]
    paths += [ROOT/name for name in ('CENTRAL_RIGIDITY_FIBERS.md',
        'tools/audit_central_rigidity.py','central-rigidity-audit.json')]
    paths += [ROOT/name for name in ('tools/audit_cubic_root_coefficients.py',
        'cubic-root-coefficients-audit.json')]
    paths += [ROOT/name for name in ('ADJOINT_GENERIC_FIBER_THEOREM.md',
        'tools/audit_adjoint_theorem.py','adjoint-theorem-audit.json')]
    paths += [ROOT/name for name in (
        'review/adjoint-generic-fibers/README.md',
        'review/adjoint-generic-fibers/CERTIFICATE_GUIDE.md',
        'review/adjoint-generic-fibers/run_review.py',
        'review/adjoint-generic-fibers/package-manifest.json',
        'review/adjoint-generic-fibers/package-verification.json',
        'tools/build_adjoint_review_package.py')]
    report=dict(result='passed',version=version,
        adjoint_review_package=dict(proof_commit=review_package['proof_commit'],
            input_files=review_package['input_files'],
            archived_bytes_checked=True, extracted_full_replay_passed=True,
            fresh_pinned_engine_build_passed=True, external_review_obtained=False,
            guide='review/adjoint-generic-fibers/README.md'),
        adjoint_theorem_consolidation=dict(standard_E_subsets=445,cubic_transfers=136,
            scheme_controls_checked=True,dependency_bindings_checked=True,
            external_peer_review=False,heavy_historical_memberships_rerun=False,
            note='ADJOINT_GENERIC_FIBER_THEOREM.md'),
        active_research_target=dict(uniform_simple_adjoint_conjecture=False,
            automatic_Levi_case_completion_paused=False,
            general_proof_complete=True,proof_method='classical and exceptional cases',
            uniform_root_theoretic_proof_complete=False,exceptional_case_completion_resumed=True,
            all_simple_adjoint_generic_regular_fiber_theorem_proved=True,
            current_E7_standard_Levis_covered=127,all_proper_E7_Levis_covered=True,
            current_E8_standard_Levis_covered=255,current_E8_unresolved_Weyl_classes=0,
            all_proper_E8_Levis_covered=True,E8_last_affine_singular_dimension_bound=2,
            E8_four_further_classes_generic_smooth_or_empty_proved=True,
            E8_A1_triple_generic_nonempty_smooth_regular_fibers_proved=True,
            E8_A2_A1_pair_generic_nonempty_smooth_regular_fibers_proved=True,
            all_proper_G2_Levis_covered=True,
            all_proper_F4_Levis_covered=True,
            authorized_E8_test_candidate_completed=True,
            larger_centralizer_shortcut_disproved=True,
            actual_adjoint_exceptional_triple_point_verified=True,
            uniform_central_rigidity_proved=True,
            rigidity_constancy_on_reductions=True,
            finite_central_scheme_image_proved=True,
            simply_laced_cubic_coefficient_structure_proved=True,
            central_trace_coordinate_scheme_elimination_proved=True,
            nonlinear_coefficient_centralizer_character_structure_proved=True,
            individual_coefficient_boundary_shape_proved=True,
            simultaneous_coefficient_zero_schemes_reduced_proved=True,
            coefficient_only_rank_argument_disproved=True,
            full_adjoint_nonresonance_still_unproved=True,
            cubic_root_carrier_classification_proved=True,
            smooth_relative_fiber_bad_image_exclusion_proved=True,
            overlapping_carrier_bad_image_bound_proved=False,
            all_flat_invariant_operator_shortcut_excluded=True,
            classical_radial_character_identity_verified=True,
            radial_identity_Schur_transport_written=True,
            determinant_only_scalar_reconstruction_for_fixed_E8_block_disproved=True,
            radial_identity_forces_singular_image_rank_defect=False,
            singular_block_radial_conormal_identity_written=True,
            scalar_radial_sufficiency_disproved_by_abstract_control=True,
            singular_block_control_is_adjoint_character=False,
            singular_block_control_generic_nilpotents_independently_verified=True,
            external_discriminant_Jacobian_square_Steinberg_corollary_written=True,
            criticality_propagation_to_higher_Adams_characters_disproved=True,
            full_character_discriminant_singular_image_bound_proved=False,
            all_type_A_relative_incidence_nilradical_described=True,
            type_A_central_target_coordinates_annihilate_nilradical=True,
            type_A_generic_reducedness_previously_quadratic_corollary=True,
            actual_adjoint_relative_incidence_always_reduced=False,
            uniform_all_simple_generic_target_reducedness_proved=True,
            generic_target_reducedness_scope='proper-Levi centralizer-exact targets, regular source',
            E6_radial_critical_ideal_descent_after_target_scaling_disproved=True,
            all_exceptional_nilradical_target_annihilators_derived=False,
            aligned_family_exact_nilradical_and_colon_proved=True,
            aligned_family_exact_target_annihilator_proved=True,
            aligned_E6_full_critical_ideal_target_annihilator_derived=True,
            aligned_generic_smoothness_coverage_previously_known=True,
            shared_square_couplings_generic_smoothness_proved=True,
            shared_square_couplings_new_group_coverage=False,
            universal_affine_shared_variable_normal_form_disproved=True,
            Weyl_symmetry_actual_simple_representation_sufficiency_disproved=True,
            adjoint_long_root_jet_projectors_verified=True,
            root_jet_projectors_imply_uniform_smoothness_proved=False,
            finite_central_image_alone_suffices_for_norm_certificate=False,
            simultaneous_boundary_transversality_proved=False,
            note='E8_LAST_ADJOINT_FIBER.md'),
        E8_adjoint_quasifinite_candidate=dict(retained_nodes=[1,2,3,5,6,8],
            derived_Levi_type='A2+A1+A2+A1',
            every_nonempty_allowed_central_fiber_zero_dimensional=True,
            all_parameter_intersections_retained=True,
            relative_scheme_flat_reduced_pure_dimension_two=True,
            generic_target_finite_etale_or_zero=True,
            regular_source_generic_smoothness_for_this_Levi=True,
            independent_exact_characteristic_zero_identity_audit=True,
            earlier_boundary_independent_replays_completed=6,
            earlier_boundary_independent_replays_total=6,
            boundary_direct_membership_certificates=4,boundary_original_potential_involutions=2,
            generic_target_degree_computed=False,nonemptiness_proved=False,
            all_special_targets_smooth=False,uniform_simple_adjoint_conjecture_proved=False,
            novelty_claim=False,certificate='adjoint-quasifinite-verification.json',
            independent_audit='adjoint-quasifinite-audit.json'),
        springer_cayley_exposition=dict(pages=11,E7_resolved=123,E7_unresolved=4,
            source_manifest_verified=True,single_manuscript_page_per_pdf_page=True,
            dated_source_manifest_preserved=True,current_status_source_drift=exposition_live_source_drift,
            new_mathematical_claim=False,literature_priority_claim=False,
            artifact_verification='../articles/springer-cayley-exposition/verification.json'),
        binary_branching=dict(common_critical_equation_criterion=True,
            generic_regular_fibers_nonempty_smooth=True,
            new_cubic_Weyl_classes=['A1^3','A2+A1^2'],
            new_standard_embeddings=22,regular_fiber_dimensions=[6,10],
            resolved_standard_Levis=123,remaining_embeddings=4,remaining_Weyl_classes=2,
            companion_memberships_verified=False,all_simple_adjoint_conjecture_proved=False,
            full_map_degree_claim=False,literature_priority_claim=False,
            certificate='binary-branching-verification.json',
            independent_audit='binary-branching-audit.json'),
        fiber_review_checkpoint=dict(standard_bad_image_criterion=True,
            independent_exact_controls=8,E7_resolved=101,E7_unresolved=26,
            selected_memberships_verified=False,new_group_classification=False,
            all_simple_adjoint_conjecture_proved=False,literature_priority_claim=False,
            next_milestone='Common root/central-charge criterion for multiple unresolved families',
            audit='cayley-fiber-review-audit.json'),
        E7_symmetric_pencil=dict(rank_stratum_criterion=True,invertible_leading_matrix_required=False,
            new_Levi_type='A1^4',new_standard_embeddings=2,regular_fiber_dimension=8,
            generic_regular_fibers_nonempty_smooth=True,resolved_standard_Levis=101,
            remaining_cubic_embeddings=26,remaining_cubic_families=4,
            companion_Levi_type='A2+A1^3',companion_exact_reduction_checked=True,
            companion_independent_dimension_certificate_complete=False,
            companion_generic_smoothness_promoted=False,all_E7_Levis_resolved=False,
            full_map_degree_claim=False,literature_priority_claim=False,
            certificate='e7-symmetric-pencil-verification.json',
            independent_audit='e7-symmetric-pencil-audit.json'),
        E7_multiplier_boundary=dict(cubic_Weyl_classes=7,cubic_standard_embeddings=32,
            common_maximal_Levi_refinement=True,nonunit_multiplier_boundary_criterion=True,
            new_Levi_type='A3+A1+A1',new_standard_embeddings=3,
            generic_regular_fibers_nonempty_smooth=True,regular_fiber_dimension=16,
            resolved_standard_Levis=99,total_proper_standard_Levis=127,
            remaining_cubic_classes=5,remaining_cubic_embeddings=28,
            all_E7_Levis_resolved=False,full_map_degree_claim=False,literature_priority_claim=False,
            certificate='e7-multiplier-boundary-verification.json',
            independent_audit='e7-multiplier-boundary-audit.json'),
        exceptional_quartic_reduction=dict(uniform_all_A_factor_elimination=True,
            proper_E_type_subsets=445,all_A_subsets=420,quadratic_subsets=284,cubic_subsets=136,
            original_coordinate_simultaneous_classification_complete=True,
            maximum_remaining_affine_variables=4,E6_all_proper_Levi_generic_smooth_or_empty=True,
            E6_cubic_embeddings_covered_by_explicit_Weyl_words=11,
            E7_E8_all_Levi_smoothness_proved=False,recursive_reduction_classified=False,
            literature_priority_claim=False,certificate='unit-constraint-classification-verification.json',
            independent_audit='unit-constraint-classification-audit.json'),
        coupled_multipliers=dict(unit_matrix_elimination_written_proof=True,
            unique_matching_sufficient_criterion=True,
            polynomial_residual_degree='max(2*deg(A)+2*deg(B),4*deg(B)+1)',
            E8_Levi='A2+A1+A4, omit node 4',E8_target_normalization='c=1800*t',
            E8_generic_regular_orbits=6,E8_regular_orbit_dimension=28,
            E8_nonzero_exceptional_targets=10,E8_exceptional_quotient='one double and four simple',
            direct_root_differential_rank=220,all_exceptional_Levis_resolved=False,
            full_map_generic_degree_claim=False,literature_priority_claim=False,
            certificate='multiplier-branching-verification.json',
            independent_audit='multiplier-branching-audit.json'),
        adjoint_pencil_scope=dict(original_coordinate_all_A_class='E6 omit node 4',
            arbitrary_affine_changes_classified=False,multiplier_elimination_written_proof=True,
            anchored_product_transfer_written_proof=True,
            E7_Levi='A1+A2+A3, omit node 4',E7_actual_quotient='3*W^2-c*W-3, c=432*t',
            E7_generic_regular_orbits=2,E7_regular_orbit_dimension=20,
            E7_nonzero_double_point_targets='c^2+36=0',
            direct_root_differential_rank=113,all_exceptional_Levis_resolved=False,
            literature_priority_claim=False,certificate='adjoint-pencil-scope-verification.json',
            independent_audit='adjoint-pencil-scope-audit.json'),
        spectral_branching=dict(arbitrary_matrix_size_written_theorem=True,
            semisimple_eigenvalue_block_excess='2*d-1',ordinary_branch_excess=0,
            noncollision_and_boundary_hypotheses_required=True,
            same_support_nonreduced_control=True,
            reductive_representation_counterexample_dimension=36,
            counterexample_group='SL2 x SL3 x Gm',
            counterexample_regular_fiber_local_dimension=9,
            counterexample_regular_fiber_tangent_dimension=10,
            simple_group_counterexample_claim=False,
            E6_suspension_recovered_with_boundary=True,
            new_exceptional_case_claim=False,literature_priority_claim=False,
            certificate='spectral-branching-verification.json',
            independent_audit='spectral-branching-audit.json',
            scope='Sufficient coefficient-sensitive branching criterion and an actual reductive representation obstruction. No all-exceptional classification or simple-group counterexample.'),
        cubic_levi_refinements=dict(product_transfer_written_proof=True,
            Fricke_controlled_boundary_written_proof=True,universal_relative_algebra_rank=5,
            rank_five_is_not_Springer_degree=True,
            specified_E6_refinements=['A1^2+A2','A1^3'],
            generic_regular_fibers_nonempty_smooth=True,generic_regular_dimensions=[10,6],
            original_three_cubic_mechanism_tests_resolved=True,
            all_exceptional_Levis_resolved=False,literature_priority_claim=False,
            certificate='cubic-levi-refinements-verification.json',
            independent_audit='cubic-levi-refinements-audit.json',
            scope='Exact product reduction and universal Fricke algebra with explicit boundary-image exclusion. Embedding-sensitive statements; no whole-fiber or generic-degree claim.'),
        branching_norm_fibers=dict(charge_sensitive_degree_bound_written_proof=True,
            simply_laced_low_radius_all_Levi_quadraticity=True,
            E6_minimal_all_Levi_generic_regular_smoothness=True,
            adjoint_cubic_shape_written_proof=True,
            exact_nonzero_target_suspension_written_proof=True,
            adjoint_mechanism_tests=3,then_unresolved_tests_now_resolved=['A1^3','A1^2+A2'],
            E6_A2A1A2_actual_quotient='3*w^2+2*c*w-3=0, c=144*t',
            E6_A2A1A2_generic_regular_orbits=2,E6_A2A1A2_orbit_dimension=14,
            all_exceptional_adjoint_Levis_resolved=False,literature_priority_claim=False,
            certificate='branching-norm-fibers-verification.json',
            independent_audit='branching-norm-fibers-audit.json',
            scope='Structural norm deductions and exact suspension; all-Levi minimal E6 smoothness allows emptiness and excess. No full adjoint classification or public API change.'),
        exceptional_adjoint_contact=dict(contact_gradings_checked=5,
            independent_native_module_and_root_reconstruction=True,
            E6_actual_quotient_critical_scheme='A^2',E6_regular_fiber_dimension=32,
            E6_excess_over_orbit_dimension=2,E6_full_differential_rank=46,
            root_length_obstruction_written_proof=True,
            linear_coupling_transfer_written_proof=True,
            empty_contact_regular_parts=['F4/C3','E7/D6','E8/E7'],
            G2_contact_generic_regular_dimension=2,
            all_exceptional_Levis_classified=False,literature_priority_claim=False,
            certificate='exceptional-adjoint-contact-verification.json',
            independent_audit='exceptional-adjoint-contact-audit.json',
            scope='Published contact branching; written fiber deductions and central descent. E6 refutes unchanged classical adjoint dimension, not smoothness. Whole fibers and other exceptional Levis remain separate; public API unchanged.'),
        classical_adjoint_fibers=dict(all_rank_generic_smoothness=True,
            all_rank_generic_nonemptiness_classified=True,
            generic_regular_dimension='dim K-rank G',
            BCD_criterion='GL blocks <=4; residual rank <=2/1/3, D_n with n>=4',
            proof_mechanisms=['controlled square substitution','quartic spectral lemma',
                              'zero-target deformation with nonzero trace derivative'],
            independent_matrix_cases=4,trace_form_pairings=1425,
            Springer_trace_pairings=69,full_differential_pairings=1425,
            symbolic_deformation_branches=11,rank_sweep=False,
            certificate='classical-adjoint-fibers-verification.json',
            scope='Written complete classical adjoint theorem for generic regular parts over proper Levi strata. Whole fibers, special targets, counts, other representations and priority remain separate. Public API unchanged.'),
        general_cayley_fiber_critical=critical_fiber,
        simple_group_branching=simple_branching,
        symplectic_target_strata=symplectic_strata,
        spin_fiber_excess=spin_excess,
        affine_branching_fibers=affine_branching,
        quadratic_branching_fibers=dict(generic_smoothness_without_constant_rank=True,
            dimension_formula='dim K-rank G+ell-k-q',
            type_A_adjoint_all_rank_nonemptiness='all blocks at most two and at least one singleton',
            type_A_adjoint_dimension='twice the number of double blocks',
            scheme_controls=8,type_A_partitions=58,dominance_checks=24,
            classical_root_restrictions=347,cubic_exclusions=37,
            original_adjoint_matrix_cases=8,original_trace_form_pairings=865,
            original_Springer_trace_pairings=142,original_full_differential_pairings=2882,
            independent_audit='type-a-adjoint-fibers-audit.json',
            certificate='quadratic-branching-fibers-verification.json',
            scope='Written quadratic theorem includes rank changes. Complete generic regular-fiber criterion for SL_N adjoint; classical B/C/D quadratic coverage, with B1/D2 residual cubic exceptions. Wholly nonregular fibers, classical transfer nonemptiness beyond A, counts, exceptional target equations and literature priority remain separate. Runtime unchanged.'),
        affine_branching_classification=dict(classical_all_rank_written_proof=True,
            classical_full_character_checks=2453,exceptional_pairs=3373,
            exceptional_affine_pairs=133,exceptional_nonaffine_pairs=3240,
            exceptional_full_character_checks=1462,native_branching_checks=3,
            all_exceptional_certificates_independently_audited_without_Sage=True,
            distinct_extremal_word_certificates=637,
            certificate='affine-branching-classification-verification.json',
            independent_audit='affine-branching-classification-audit.json',
            scope='Every fundamental representation and proper standard Levi of a simple complex simply connected group. All-rank classical proof and exhaustive finite exceptional norm-gap/extremal-weight certificates. Nonaffinity does not imply fiber singularity. Standard representation-theory inputs; novelty unverified; runtime unchanged.'),
        symplectic_fiber_component=symplectic_fiber,
        symplectic_centralizer_realization=symplectic_realization,
        middle_exterior_fiber_component=middle_fiber,
        irreducible_centralizer_realization=irreducible_realization,
        type_A_character_degeneration=character_degeneration,
        Cayley_fiber_assessment=fiber_assessment,
        G2_vertical_rank=g2_rank,
        borel_coordinate_API=borel_coordinates,
        borel_reconstruction_API=borel_reconstruction,
        automatic_root_borel_API=root_borel,
        general_degree_audit=general_degree,
        e6_adjoint_deformation=e6_deformation,
        f4_minimal_degree=f4_degree,
        BC_defining_adjoint_degrees=dict(formula_B=bc_mixture['formula_B'],
            formula_C=bc_mixture['formula_C'],independent_reduced_fibers=68,
            isogeny_controls=2,arbitrary_rank_written_proof=True,
            exceptional_B_nonfinite_zero_curves=True,
            spin_degree_multiplier=2,pure_C_adjoint_form_degree_divisor=2,
            novelty_claim=False,runtime_API_changed=False),
        f4_mixture_degree=dict(formula=f4_mixture['formula'],primitive_facet_height=2,
            exceptional_boundary_corrections={'a=b':144,'a=5b':48},
            independent_reduced_fibers=20,same_target_multiplication_polynomials_agree=True,
            characteristic_zero_written_proof=True,novelty_claim=False,runtime_API_changed=False),
        half_spin_rank_degree=dict(rank_minimum=3,
            formula='2^(r-1)*(r^2-3*r-2)+4*r',
            characteristic_zero_deformation_written_proof=True,
            deeper_faces_excluded=True,Spin12_degree=536,
            independent_reduced_Spin12_fibers=4,novelty_claim=False),
        Spin10_half_spin_degree=dict(generic_degree=148,
            characteristic_zero_toric_written_proof=True,
            toric_intersection_number=208,boundary_length=60,
            exact_reduced_fibers_in_two_models=12,
            nonreduced_specialization_controls=2,
            all_rank_degree_formula=False,novelty_claim=False),
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
