"""One independent full-flag control for the all-rank free-source transfer.

Starts from both actual SL3 root differences in the rank-six flag algebra,
not the projective presentation or its predicted tensor decomposition.
This finite control is not the proof of the general transfer theorem.
"""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import ZZ, PolynomialRing, matrix, vector, identity_matrix, zero_matrix, block_matrix, singular

ROOT=Path(__file__).resolve().parents[1]


def main():
    S=PolynomialRing(ZZ,names=['a','b']);a,b=S.gens()
    I=identity_matrix(S,3);Z=zero_matrix(S,3)
    X=matrix(S,[[0,0,1],[1,0,-b],[0,1,a]])
    XX=block_matrix([[X,Z],[Z,X]])
    YY=block_matrix([[Z,-b*I+a*X-X*X],[I,a*I-X]])
    assert XX*YY==YY*XX
    ZZZ=a*identity_matrix(S,6)-XX-YY
    assert XX+YY+ZZZ==a*identity_matrix(S,6)
    assert XX*YY+XX*ZZZ+YY*ZZZ==b*identity_matrix(S,6)
    assert XX*YY*ZZZ==identity_matrix(S,6)
    roots=[XX-YY,XX-ZZZ]
    one=vector(S,[1,0,0,0,0,0])
    initial=matrix(S,[one]+[M*one for M in roots]).transpose()
    def literal(M):return ','.join(str(x).replace('**','^') for x in M.list())
    singular.eval('ring r=integer,(a,b),dp;')
    singular.eval('matrix M[6][3]='+literal(initial)+'; module G=std(module(M));')
    for j,M in enumerate(roots):singular.eval('matrix T%d[6][6]=%s;'%(j,literal(M)))
    rounds=0
    while True:
        rounds+=1;singular.eval('module extra;')
        for j in range(2):singular.eval('extra=extra,reduce(module(T%d*matrix(G)),G);'%j)
        if int(singular.eval('size(extra);'))==0:break
        singular.eval('G=std(G,extra);kill extra;')
        assert rounds<20
    for j in range(2):
        assert int(singular.eval('size(reduce(module(T%d*matrix(G)),G));'%j))==0
    B=singular('matrix(G)').sage()
    singular.eval('module sy=syz(G);')
    R=singular('matrix(sy)').sage();assert B*R==0
    evaluate=lambda x:ZZ(S(x)(3,3))
    B0=matrix(ZZ,B.nrows(),B.ncols(),[evaluate(x) for x in B.list()])
    R0=matrix(ZZ,R.nrows(),R.ncols(),[evaluate(x) for x in R.list()])
    torsion=[int(x) for x in R0.elementary_divisors() if abs(x)>1]
    ordinary=[int(x) for x in B0.elementary_divisors() if x]
    assert B0.rank()==6 and R0.rank()==B.ncols()-6
    assert torsion==[3]*4 and ordinary==[1,1,3,3,3,3]
    names=['TYPE_A_K_INTEGRAL_RECURRENCE.md','tools/verify_type_a_k_flag_transfer.py']
    report=dict(result='passed',scope='One actual full-flag control; the arbitrary-rank proof is free source base change.',
        m=3,ambient_rank=6,full_source_generators='x-y,x-z',
        determinant_one=True,closure_rounds=rounds,full_integral_closure=True,
        specialization_kernel_factors=torsion,ordinary_smith=ordinary,
        projective_presentation_imported=False,
        presentation=[[str(x) for x in row] for row in B.rows()],
        syzygies=[[str(x) for x in row] for row in R.rows()],
        source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in names})
    (ROOT/'type-a-k-flag-transfer-verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('presentation','syzygies','source_sha256')},indent=2))


if __name__=='__main__':main()
