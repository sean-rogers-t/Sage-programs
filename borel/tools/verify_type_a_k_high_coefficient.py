"""Independent Sage replay of the native recursive coefficient witnesses."""
from hashlib import sha256
from math import comb
from pathlib import Path
import gzip
import json
from sage.all import QQ, ZZ, PolynomialRing

ROOT=Path(__file__).resolve().parents[1]


def main():
    archive=ROOT/'experiments/type-a-k-high-coefficient-witnesses.json.gz'
    data=json.loads(gzip.decompress(archive.read_bytes()))
    count=0
    for case in data['cases']:
        e,k,d=case['excess'],case['coefficient_degree'],case['weight']
        assert d==k+e and ((e==1 and k>=3) or (e>=2 and k>=e+1))
        R=PolynomialRing(QQ,names=['M','L']+['a%d'%i for i in range(1,d+1)])
        M,L,*aa=R.gens()
        qs=[M*L**j+sum(comb(j,i)*aa[i-1]*L**(j-i) for i in range(1,j+1)) for j in range(1,e+3)]
        images=[M]+list(aa)+qs
        for witness in case['witnesses']:
            total=R(0)
            for numerator,denominator,powers in witness['terms']:
                assert sum(powers)==k
                assert all(p<=max(2,e) for p in ZZ(denominator).prime_divisors())
                if witness['coefficient_ideal_witness']:
                    assert any(powers[1:d+1])
                total+=QQ(numerator)/denominator*R.prod(g**a for g,a in zip(images,powers))
            expected=M**witness['mass_power']*L**witness['L_power']*R.prod(
                aa[i-1] for i in witness['coefficient_indices'])
            assert total==expected
            count+=1
    names=['TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md','tools/verify_type_a_k_high_coefficient.py']
    report=dict(result='passed',exact_symbolic_witnesses=count,native_builder_imported=False,
                method='Independent Sage polynomial substitution of saved rational coefficient lists',
                all_rank_claim_from_finite_checks=False,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in names},
                raw_sha256={str(archive.relative_to(ROOT)).replace('\\','/'):sha256(archive.read_bytes()).hexdigest()})
    (ROOT/'type-a-k-high-coefficient-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',exact_symbolic_witnesses=count)))


if __name__=='__main__':
    main()
