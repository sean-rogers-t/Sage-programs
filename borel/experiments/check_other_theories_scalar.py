"""Exact Laurent check of the global tensor-scalar equivariant K family in A2."""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import QQ,matrix,prod
from sage_borel import BorelMap,SpringerCharacterMap
from other_theories_k_pilot import demazure

ROOT=Path(__file__).resolve().parents[1]


def main():
    std=SpringerCharacterMap('A2',((1,0),))
    tensor=SpringerCharacterMap('A2',((2,0),(0,1)))
    t1,t2=std.torus_ring.gens();chi=t1+t2/t1+1/t2;a=chi/3
    assert all(x==a*y for x,y in zip(tensor.theta,std.theta))
    primary=json.loads((ROOT/'experiments/other-theories-k-pilot.json').read_text())
    cert=json.loads((ROOT/'experiments/other-theories-k-certificates.json').read_text())
    i=next(i for i,c in enumerate(primary['cases']) if c['cartan_type']=='A2' and c['highest_weights']==[[1,0]])
    case=cert['cases'][i];H=matrix(QQ,case['H_rows']);pivots=H.pivots()
    degrees=[sum(case['exponents'][j]) for j in pivots]
    B=matrix(QQ,primary['cases'][i]['linear_lattice_rows']);b=BorelMap('A2')
    inputs=[prod(sum(B[u,v]*b.weights[u] for u in range(2))**p for v,p in enumerate(case['exponents'][j])) for j in pivots]
    matrices=[]
    for ch in (std,tensor):
        columns=[]
        for f in inputs:
            pull=ch.pullback(ch.ring(f));terms={tuple(map(int,p)):v for p,v in pull.dict().items()};column=[]
            for word in case['words']:
                d=terms
                for node in reversed(word):d=demazure(d,node-1,case['roots'][node-1])
                column.append(ch.torus_ring(d))
            columns.append(column)
        matrices.append(matrix(std.torus_ring,columns).transpose())
    assert len(matrices[0].det().dict())==1  # a rational coefficient times a Laurent monomial: an R-unit
    assert matrices[1]==matrices[0]*matrix.diagonal([a**d for d in degrees])
    assert sorted(degrees)==[0,1,1,2,2,3]
    report=dict(result='passed',date='2026-10-22',group='SL3',parabolic='B',representation='V tensor V',
                global_Laurent_column_scaling_verified=True,reference_equivariant_matrix_is_unit=True,
                diagonal_character_exponents=sorted(degrees),
                scope='Independent exact actual-character/Demazure corroboration of the homogeneous-basis proof. The arbitrary-m,k,P theorem is the separate argument, not an extrapolation from A2.',
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in ('experiments/check_other_theories_scalar.py','experiments/other_theories_k_pilot.py')})
    report['raw_sha256']={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in ('experiments/other-theories-k-pilot.json','experiments/other-theories-k-certificates.json')}
    (ROOT/'experiments/other-theories-scalar-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED exact global A2 tensor-scalar matrix identity')


if __name__=='__main__':main()
