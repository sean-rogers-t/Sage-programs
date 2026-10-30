"""Test a structural generator hypothesis, not a new multiplicity census.

Candidate: the saturated B square ring is generated over its original
image by all high monomial classes m_nu/2^length, parts>=ceil(n/2).
These divisions have a separate complete-flag/quadric explanation.
"""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import Partitions, Partition, ZZ, QQ, GF, matrix, diagonal_matrix
from verify_type_b_adapted_basis import basis,S,SF

ROOT=Path(__file__).resolve().parents[1]

def main():
    result=[];M=SF.monomial()
    for A,l in [(3,2),(4,2),(3,3),(4,3),(5,3),(4,4),(5,4),(5,5),(7,3)]:
        n=A+l;c=(n+1)//2
        generators=[()]
        for m in range(c,A*l+1):
            generators.extend(tuple(p) for p in Partitions(m,min_part=c,max_part=n-1,max_length=l))
        rows=[]
        for m in range(A*l+1):
            ps,es,V,F,fs=basis(A,l,m);columns=[]
            for nu in generators:
                if sum(nu)>m:continue
                for lam in Partitions(m-sum(nu),max_part=A,max_length=l):
                    product=S(lam)*S(M(nu))
                    z=matrix(ZZ,len(ps),1,[product[Partition(p)] for p in ps])
                    coord=diagonal_matrix(ZZ,[2**e for e in es])*V*z/2**len(nu)
                    coord=coord.change_ring(ZZ)
                    columns.append(coord.column(0))
            C=matrix(ZZ,len(ps),len(columns),lambda i,j:columns[j][i])
            rank=C.change_ring(GF(2)).rank()
            rows.append(dict(m=m,rank=len(ps),candidate_mod2_rank=int(rank),spans=rank==len(ps)))
        case=dict(A=A,ell=l,generator_count=len(generators),degrees=rows,
                  all_degrees_span=all(r['spans'] for r in rows))
        result.append(case)
        print(json.dumps(dict(A=A,ell=l,spans=case['all_degrees_span'],
                             failures=[r for r in rows if not r['spans']])),flush=True)
    report=dict(status='completed',scope='Finite test of a candidate generator theorem; not a proof.',
                cases=result,source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'classical-divided-generator-test.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
