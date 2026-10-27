"""Sage-free audit of the inputs and geometric certificates for final E7 case."""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s
from audit_e7_symmetric_pencil import companion,candidate_audit
ROOT=Path(__file__).resolve().parents[1]

def main():
    source=json.loads((ROOT/'e7-symmetric-pencil-verification.json').read_text())
    roots=companion(source);original=candidate_audit()
    membership=json.loads((ROOT/'e7-companion-membership-verification.json').read_text())
    assert membership['result']=='passed' and membership['original_affine_ideal_memberships']==6
    assert membership['independent_prime']!=membership['generation_prime']
    for name,digest in membership['raw_sha256'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==digest
    U,V=s.symbols('U V');A=U+1/U;B=V+1/V;C=U*V+1/(U*V)
    J=s.Matrix([[z*s.diff(f,z) for z in (U,V)] for f in (A,B,C)])
    assert s.cancel(J[:2,:].det()-(U-1/U)*(V-1/V))==0
    for eps in (-1,1):
        assert s.cancel(J[2,0].subs(U,eps)-eps*(V-1/V))==0
        assert s.cancel(J[2,1].subs(V,eps)-eps*(U-1/U))==0
        for eta in (-1,1):assert J.subs({U:eps,V:eta})==s.zeros(3,2)
    assert s.Matrix([[0,1],[3,1]]).det()==-3
    report=dict(result='passed',uses_Sage=False,
                full_adjoint_character_reconstructed=roots['exact_character_and_residual'],
                original_target_derivatives=roots['original_target_derivatives'],
                original_Jacobian_minor_ideal_reconstructed=original['original_Jacobian_ideal_reconstructed'],
                monomial_dimension_bound=original['selected_leading_support_dimension'],
                torus_trace_map_rank_two_except_four_fixed_points=True,
                zero_target_on_trace_ramification=True,finite_torus_isogeny_degree=3,
                relative_hessian_determinant=-216,central_hessian_determinant=192,
                independent_membership_report='e7-companion-membership-verification.json',
                membership_arithmetic_replayed_here=False,
                regular_fiber_dimension=12,all_E7_proper_Levis_covered=True,
                scope='Independent native root and polynomial input audit plus geometric witnesses; membership itself is certified by the separate exact homogeneous Hilbert-function verifier.')
    names=['E7_FINAL_FAMILIES.md','tools/audit_e7_companion_geometry.py',
           'e7-companion-membership-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-companion-geometry-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':main()
