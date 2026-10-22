"""Check the actual local lattice isomorphism transporting the D4 module.

The old source is on the E6 root span. E6 opposition, followed by
projection to the first six E7 weight coordinates, intertwines the
relevant D5 and D4 groups. This establishes the full-source transfer
before imposing either Borel quotient.
"""
from pathlib import Path
from hashlib import sha256
import gzip,json
import sympy as S

ROOT=Path(__file__).resolve().parents[1]


def main():
    C=2*S.eye(7)
    for i,j in ((1,3),(3,4),(2,4),(4,5),(5,6),(6,7)):
        C[i-1,j-1]=C[j-1,i-1]=-1
    C6=C[:6,:6]
    assert C6.det()==3
    opposition=(6,2,5,4,3,1)
    F=S.zeros(7,6)
    for j,i in enumerate(opposition):F[i-1,j]=1
    assert abs(F[:6,:].det())==1
    old_weights=C[:,:6]*C6.inv()
    assert old_weights[:6,:]==S.eye(6)
    # On the old E6 root span, Omega is an R-basis since det(C6)=3.
    # F explicitly sends Omega_j to omega_opposition(j).
    old_D5=(2,3,4,5,6)
    for j in old_D5:
        k=opposition[j-1]
        assert F*C6.col(j-1)==C.col(k-1)
        old_reflection=S.eye(6)-C6.col(j-1)*S.eye(6).row(j-1)
        new_reflection=S.eye(7)-C.col(k-1)*S.eye(7).row(k-1)
        assert F*old_reflection==new_reflection*F
    assert sorted(opposition[j-1] for j in old_D5)==[1,2,3,4,5]
    assert sorted(opposition[j-1] for j in (2,3,4,5))==[2,3,4,5]
    E=S.eye(6);om=S.eye(7)
    old_us=[-E.col(1)+E.col(2),E.col(1)+E.col(2)-E.col(3),
            E.col(3)-E.col(4),E.col(4)-E.col(5),E.col(5)]
    new_us=[om.col(0),om.col(2)-om.col(0),om.col(3)-om.col(2),
            om.col(1)+om.col(4)-om.col(3),om.col(4)-om.col(1)]
    assert [F*u for u in old_us]==new_us[::-1]
    assert F*E.col(0)==om.col(5) and F*E.col(5)==om.col(0)
    s7=S.eye(7)-C.col(6)*S.eye(7).row(6)
    assert s7*F==F and s7*om.col(6)==om.col(5)-om.col(6)
    inputs=[ROOT/'experiments'/n for n in ('e7-d4-source-audit.json',
        'e7-d5-source-audit.json','e7-d4a1-generators-audit.json')]
    checks=[json.loads(p.read_text()) for p in inputs]
    assert all(c['result']=='passed' for c in checks)
    assert checks[0]['full_source_complete_all_degrees']
    assert checks[0]['relative_module_generators']==34
    assert checks[1]['full_source_complete_all_degrees']
    archive=ROOT/'experiments/e7-d4-tower-certificates.json.gz'
    cert=json.loads(gzip.decompress(archive.read_bytes()))
    assert len(cert['module_exponents'])==34 and len(cert['source_preimages'])==10
    report=dict(result='passed',full_weight_source_transfer_verified=True,
        old_E6_root_span_index=3,opposition=list(opposition),
        target_weight_submodule='Z_(2){omega_1,...,omega_6}',
        integral_basis_determinant=1,D5_reflection_intertwiners=5,
        D4_retained_nodes=[2,3,4,5],D5_retained_nodes=[1,2,3,4,5],
        old_u_coordinates_map_to_reversed_new_u=True,
        full_weight_source='A_D4[t], t=omega_7*(omega_6-omega_7)',
        complete_module_generators=34,
        scope='An exact local lattice and reflection isomorphism transports the earlier full D4/D5 source module before Borel. The A1 reflection fixes the six-dimensional coefficient lattice and acts by w -> h-w, giving the full norm-variable extension. The earlier all-degree source/module theorems are shared audited inputs.',
        engine='Native SymPy exact Cartan lattice and full reflection-matrix identities; no Sage or primary experiment imports.',
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in [Path(__file__).resolve()]+inputs},
        raw_sha256={archive.relative_to(ROOT).as_posix():sha256(archive.read_bytes()).hexdigest()})
    (ROOT/'experiments/e7-d4a1-weight-transfer-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED D4+A1 full weight-source transport, 34 module generators')


if __name__=='__main__':main()
