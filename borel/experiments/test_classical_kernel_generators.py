"""Test an integral ideal presentation for every primitive-prefix kernel.

Candidate generators: m_nu with >=e parts, all >=ceil(n/2).
This is a theorem-search control; agreement alone is not a proof.
"""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import Partitions, Partition, ZZ, matrix
from verify_type_b_adapted_basis import basis,S,SF

ROOT=Path(__file__).resolve().parents[1]

def main():
    M=SF.monomial();cases=[]
    for A,l in [(3,2),(4,2),(3,3),(4,3),(5,3),(4,4),(5,4),(5,5),(7,3)]:
        n=A+l;c=(n+1)//2;gs=[()]
        for m in range(c,A*l+1):
            gs.extend(tuple(p) for p in Partitions(m,min_part=c,max_part=n-1,max_length=l))
        rows=[]
        for m in range(A*l+1):
            ps,es,V,F,fs=basis(A,l,m)
            columns=[]
            for nu in gs:
                if not nu or sum(nu)>m:continue
                for lam in Partitions(m-sum(nu),max_part=A,max_length=l):
                    value=S(lam)*S(M(nu))
                    col=V*matrix(ZZ,len(ps),1,[value[Partition(p)] for p in ps])
                    assert all(col[i,0]==0 for i,e in enumerate(es) if e<len(nu))
                    columns.append((len(nu),col.column(0)))
            for e in range(1,l+1):
                ix=[i for i,v in enumerate(es) if v>=e]
                selected=[col for length,col in columns if length>=e]
                C=matrix(ZZ,len(selected),len(ix),lambda i,j:selected[i][ix[j]])
                exact=C.row_module()==ZZ**len(ix)
                rows.append(dict(m=m,e=e,kernel_rank=len(ix),generated_integrally=exact))
        case=dict(A=A,ell=l,all_equal=all(r['generated_integrally'] for r in rows),components=rows)
        cases.append(case)
        print(json.dumps(dict(A=A,ell=l,all_equal=case['all_equal'],
                             failures=[r for r in rows if not r['generated_integrally']])),flush=True)
    report=dict(status='completed',scope='Finite ideal test, not an arbitrary-rank proof.',cases=cases,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'classical-kernel-generator-test.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
