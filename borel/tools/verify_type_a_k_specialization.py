"""Character/Schubert input checks and stored complete module certificates."""
from pathlib import Path
from hashlib import sha256
from math import comb
from itertools import combinations
import json
from sage.all import ZZ, GF, PolynomialRing, matrix, prod
from sage_borel import type_a_equivariant

ROOT=Path(__file__).resolve().parents[1]
def digest(name):return sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()

def main():
    cases=[];restrictions=0;expansions=0
    expected={2:[],3:[3,3],4:[2,4,4],5:[5]*20,6:[3]*10+[6]*6}
    bound=['TYPE_A_K_SPECIALIZATION.md','tools/verify_type_a_k_specialization.py']
    for n in range(2,7):
        k=type_a_equivariant(n,omitted_nodes=1,theory='K')
        ts=k.diagonal_parameters
        # Sage 9.5's pre-existing one-variable Laurent coercion in the API
        # fails for SL2; its direct two-by-two control is checked natively.
        for j,g in enumerate(k.source_generators() if n>2 else (),1):
            direct=k.restrictions(g)
            for w,value in direct.items():
                selected=ts[w[0]-1]
                differences=[selected-t for i,t in enumerate(ts) if i!=w[0]-1]
                want=sum(prod(c) for c in combinations(differences,j))
                assert value==want
                restrictions+=1
            assert k.expand(g).restrictions()==direct
            expansions+=1
        name='experiments/type-a-k-specialization-%d.json'%n
        row=json.loads((ROOT/name).read_text());bound.append(name)
        assert row['result']=='passed' and row['full_module_closure']
        assert row['alternate_order_check']['same_kernel_and_ordinary_factors']
        for key,value in row['source_sha256'].items():assert digest(key)==value
        S=PolynomialRing(ZZ,names=['c'+str(i) for i in range(1,n)])
        B=matrix(S,[[S(x) for x in a] for a in row['presentation']])
        Q=matrix(S,[[S(x) for x in a] for a in row['syzygies']])
        assert B*Q==0
        ev={S.gen(i-1):comb(n,i) for i in range(1,n)}
        Q0=matrix(ZZ,Q.nrows(),Q.ncols(),[ZZ(x.subs(ev)) for x in Q.list()])
        factors=[int(abs(x)) for x in Q0.elementary_divisors() if abs(x)>1]
        assert factors==row['specialization_kernel_factors']==expected[n]
        local=[]
        for p in ZZ(n).prime_divisors():
            min_generators=B.ncols()-Q0.change_ring(GF(p)).rank()
            assert min_generators>=n
            assert (min_generators==n)==(n==2)
            local.append({'p':int(p),'minimal_local_generators':int(min_generators),'generic_rank':n})
        cases.append({'n':n,'kernel_factors':factors,'local_checks':local})
    report=dict(result='passed',fixed_point_checks=restrictions,schubert_expansions=expansions,
                API_rank_range=[3,6],SL2_API_limitation='Existing Sage 9.5 univariate Laurent coercion; SL2 checked by direct presentation, API unchanged.',
                cases=cases,scope='Independent character/Schubert inputs and exact certificate checks; syzygy completeness relies on the recorded strong integral standard-basis calculation.',
                source_sha256={name:digest(name) for name in bound})
    (ROOT/'type-a-k-specialization-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
