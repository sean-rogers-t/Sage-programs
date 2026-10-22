"""Exact checks supporting a claim-specific literature audit.

Only standard-library rational arithmetic is used. No downloaded notebook
cell is executed. Literature-search completeness and priority are not
mathematical conclusions certified by this script.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def rref(M):
    M=[[Q(x) for x in row] for row in M]
    row=0;pivots=[]
    for col in range(len(M[0])):
        pivot=next((i for i in range(row,len(M)) if M[i][col]),None)
        if pivot is None:continue
        M[row],M[pivot]=M[pivot],M[row]
        a=M[row][col];M[row]=[x/a for x in M[row]]
        for i in range(len(M)):
            if i!=row:
                a=M[i][col];M[i]=[x-a*y for x,y in zip(M[i],M[row])]
        pivots.append(col);row+=1
        if row==len(M):break
    return M,pivots


def multiply(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def main():
    record=json.loads((ROOT/'experiments/cayley-literature-source-record.json').read_text())
    assert len(record['sources'])==8
    assert {s['key'] for s in record['sources']}=={'KM','R','K','KR','KX','LPR','P','LY'}
    for s in record['sources']:
        if s['local_original']:
            # Resolve repository copies portably, rather than requiring
            # the original Windows absolute path on another host.
            name=s['local_original'].replace('\\','/').split('/')[-1]
            assert sha256((ROOT.parent/'articles'/name).read_bytes()).hexdigest()==s['sha256']
    b=record['notebook'];A=b['cell3_cartan'];B=b['cell3_SInv'];n=len(A)
    identity=[[int(i==j) for j in range(n)] for i in range(n)]
    augmented,pivots=rref([row+unit for row,unit in zip(A,identity)])
    assert pivots==list(range(8))
    inverse=[row[n:] for row in augmented]
    assert multiply(A,inverse)==identity
    assert inverse[-1]==[2,3,4,6,5,4,3,2]
    assert B[-1]==B[1]
    _,pivots_B=rref(B)
    assert len(pivots_B)==7 and multiply(A,B)!=identity
    X=b['cell20_X']
    assert all(len(row)==7 for row in X) and all(len(row)==8 for row in A)
    assert b['cells']==30 and b['execution_counts_all_null'] and b['not_executed']

    # Direct adjoint trace projection on SL2: D=diag(2,0,-2).
    # Coefficients of z^2 and z^-2 in tr(rho(t)D)/tr(D^2).
    weights=[2,0,-2];gram=sum(w*w for w in weights)
    corrected=[Q(2,gram),Q(-2,gram)]
    printed=[Q(1,2),Q(-1,2)]
    assert gram==8 and corrected==[Q(1,4),Q(-1,4)]
    assert 2*corrected[0]-2*corrected[1]==1
    assert 2*printed[0]-2*printed[1]==2
    # KM 6.1's general coefficient at highest weight n=2.
    assert Q(3,2**3+3*2**2+2*2)*2==corrected[0]

    # C2 defining torus witness, with its full positive root list.
    z,w=Q(2),Q(-1,2)
    roots=[z/w,z*w,z*z,w*w]
    images=[(t-1/t)/2 for t in (z,w)]
    derivative=[(t+1/t)/2 for t in (z,w)]
    assert all(t!=1 for t in roots) and images==[Q(3,4)]*2
    assert derivative[0]*derivative[1]==Q(-25,16)

    # SL4 on exterior^2: build the representation's diagonal weights
    # on simple coroots and solve the actual 3-by-3 trace equations.
    t=[Q(2),Q(-1,2),Q(1),Q(-1)]
    assert len(set(t))==4 and t[0]*t[1]*t[2]*t[3]==1
    basis=[[int(i==j)-int(i==j+1) for i in range(4)] for j in range(3)]
    pairs=list(combinations(range(4),2))
    weights=[[basis[j][a]+basis[j][b] for j in range(3)] for a,b in pairs]
    eigenvalues=[t[a]*t[b] for a,b in pairs]
    S=[[sum(mu[i]*mu[j] for mu in weights) for j in range(3)] for i in range(3)]
    rhs=[sum(mu[i]*ev for mu,ev in zip(weights,eigenvalues)) for i in range(3)]
    solution,pivots=rref([row+[v] for row,v in zip(S,rhs)])
    assert pivots==[0,1,2]
    c=[row[-1] for row in solution]
    output=[sum(c[j]*basis[j][i] for j in range(3)) for i in range(4)]
    assert output==[0,0,Q(3,4),Q(-3,4)]

    paths=['CAYLEY_LITERATURE_AUDIT.md','tools/verify_cayley_literature.py',
           'experiments/cayley-literature-source-record.json']
    report=dict(result='passed',date='2026-10-22',bibliographic_entries=8,
        notebook=dict(cells=30,static_inspection_only=True,cartan_rank=8,
            saved_SInv_rank=7,duplicate_rows_zero_based=[1,7],
            correct_inverse_last_row=[int(x) for x in inverse[-1]],
            cell20_dimensions_incompatible=True,external_data_file_required=True),
        proposed_Part_I_Example_3_scalar_correction=dict(weights=[2,0,-2],
            trace_gram=gram,corrected_z2_coefficient='1/4',printed_z2_coefficient='1/2',
            corrected_identity_derivative=1,printed_identity_derivative=2,
            agrees_with_KM_6_1_general_formula=True,author_approved_erratum=False),
        regularity_witnesses=dict(Sp4=dict(torus=[str(z),str(w)],
            positive_root_characters=[str(x) for x in roots],image=[str(x) for x in images],
            torus_derivative_determinant='-25/16'),
            SL4_exterior_square=dict(torus=[str(x) for x in t],
                image_diagonal=[str(x) for x in output],method='Six representation eigenvalues and the simple-coroot trace Gram matrix')),
        scope='Exact elementary normalization, matrix and witness checks; source-record integrity. No proof of literature-search completeness, novelty or 2026 openness, and no execution of external notebook cells.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'cayley-literature-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED: 8 source records, notebook matrix checks, Part I scalar discrepancy, Sp4 and SL4 regularity witnesses')


if __name__=='__main__':main()
