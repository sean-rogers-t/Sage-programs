"""Explicit finite full-image witnesses; Sage constructs, native audit checks."""
from pathlib import Path
from hashlib import sha256
from math import comb
import json
from sage.all import ZZ, PolynomialRing, matrix

ROOT=Path(__file__).resolve().parents[1]
def digest(p):
    return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def rows(M):
    return [[int(x) for x in row] for row in M.rows()]

def main():
    cases=[];bindings={}
    for m in range(2,7):
        name='experiments/type-a-k-specialization-%d.json'%m
        data=json.loads((ROOT/name).read_text());bindings[name]=digest(ROOT/name)
        S=PolynomialRing(ZZ,names=['c%d'%i for i in range(1,m)])
        values={S.gen(i-1):ZZ(comb(m,i)) for i in range(1,m)}
        ev=lambda x:ZZ(S(x).subs(values))
        B=matrix(ZZ,[[ev(x) for x in row] for row in data['presentation']])
        R=matrix(ZZ,[[ev(x) for x in row] for row in data['syzygies']])
        K=B.right_kernel_matrix().transpose()
        D,U,V=R.smith_form()
        target=U*(m*K);Z=matrix(ZZ,R.ncols(),K.ncols())
        for i in range(R.nrows()):
            for j in range(K.ncols()):
                if i<min(D.nrows(),D.ncols()) and D[i,i]:
                    assert target[i,j]%D[i,i]==0
                    Z[i,j]=target[i,j]//D[i,i]
                else:assert target[i,j]==0
        W=V*Z
        assert B*K==0 and R*W==m*K
        cases.append(dict(m=m,generators=B.ncols(),relations=R.ncols(),
                          kernel_rank=K.ncols(),K=rows(K),W=rows(W)))
    bindings['tools/build_type_a_k_annihilator.py']=digest(Path(__file__))
    report=dict(result='passed',cases=cases,source_sha256=bindings,
                scope='Explicit R0 W=m K witnesses for existing complete full-image presentations, m=2,...,6; no arbitrary-rank annihilator proof.')
    (ROOT/'type-a-k-annihilator-witnesses.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'result':'passed','full_image_witnesses':len(cases)}))
if __name__=='__main__':main()
