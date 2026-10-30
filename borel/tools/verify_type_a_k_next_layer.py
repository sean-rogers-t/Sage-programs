"""Fixed tensor obstruction and universal moment-lattice corroboration."""
from hashlib import sha256
from math import comb
from pathlib import Path
import json
from sage.all import ZZ, QQ, GF, PolynomialRing, matrix, lcm
from verify_type_a_k_stable_multiplicities import exponents
from verify_type_a_k_stable_multiplicities import expected
from verify_type_a_k_first_correction import build_image, primary
from binary_tensor_certificates import binary_tensor_certificate, cubic_offset_obstruction

ROOT = Path(__file__).resolve().parents[1]


def sector(p, k, raw=False):
    d = p+3
    R = PolynomialRing(ZZ,names=['L']+['a%d'%i for i in range(1,d+1)])
    L,*aa = R.gens()
    q = [L**j+sum(comb(j,i)*aa[i-1]*L**(j-i) for i in range(1,j+1)) for j in range(1,d+1)]
    # M=1 encodes homogenization to coefficient degree k, not mass specialization.
    basis = [e for e in exponents([1]+list(range(1,d+1)),d) if sum(e[1:]) <= k]
    index = {e:i for i,e in enumerate(basis)}
    source, ideal = [],[]
    for e in exponents(list(range(1,d+1))*2,d):
        if sum(e) > k:
            continue
        polynomial = R.prod(g**t for g,t in zip(list(aa)+q,e))
        row = [ZZ(0)]*len(basis)
        for power,value in polynomial.dict().items():
            row[index[tuple(power)]] = value
        source.append(row)
        if any(e[:d]):
            ideal.append(row)
    if raw:
        return basis,matrix(ZZ,source),matrix(ZZ,ideal)
    J = matrix(ZZ,source).hermite_form(include_zero_rows=False)
    # Reduce the ideal in source coordinates to avoid ambient coefficient swell.
    coordinates = matrix(ZZ,J.transpose().solve_right(matrix(ZZ,ideal).transpose()).transpose())
    C = coordinates.hermite_form(include_zero_rows=False)*J
    return basis,J,C


def p_factors(A,p):
    return [int(p**ZZ(x).valuation(p)) for x in A.elementary_divisors()
            if x and ZZ(x).valuation(p)>0]


def verify_moment_block(p):
    basis,J,C = sector(p,3)
    basis4,J4,C4 = sector(p,4)
    relative = matrix(ZZ,J.transpose().solve_right(C.transpose()).transpose())
    assert p_factors(J,p) == p_factors(C,p) == [p]
    assert p_factors(relative,p) == []
    assert p_factors(J4,p) == p_factors(C4,p) == []
    saturated = C.saturation()
    inC = C.transpose().solve_right(saturated.transpose())
    inJ = J.transpose().solve_right(saturated.transpose())
    valuations = []
    for j in range(saturated.nrows()):
        a = int(ZZ(lcm(x.denominator() for x in inC.column(j))).valuation(p))
        b = int(ZZ(lcm(x.denominator() for x in inJ.column(j))).valuation(p))
        assert a == b
        valuations.append(a)
    assert 1 in valuations
    # Replay multiplication by M as the same coefficient polynomial in the next sector.
    idx4 = {e:i for i,e in enumerate(basis4)}
    for old,new in [(J,J4),(C,C4)]:
        sat = old.saturation()
        coords = old.transpose().solve_right(sat.transpose())
        selected = next(j for j in range(sat.nrows())
                        if ZZ(lcm(x.denominator() for x in coords.column(j))).valuation(p)>0)
        order = ZZ(lcm(x.denominator() for x in coords.column(selected)))
        representative = (order//p)*sat[selected]
        row = [ZZ(0)]*len(basis4)
        for e,value in zip(basis,representative):
            row[idx4[e]] = value
        image = new.transpose().solve_right(matrix(ZZ,[row]).transpose())
        assert all(x.denominator().valuation(p)==0 for x in image.list())
    return dict(p=p, degree=p+3, coefficient_degree_three_ambient_rank=len(basis),
                source_primary_torsion=[p], ideal_primary_torsion=[p],
                relative_primary_torsion=[], natural_torsion_map_isomorphism=True,
                source_and_ideal_mass_action_zero_on_primary_torsion=True,
                coefficient_degree_four_primary_torsion=[])


def main():
    tensors = [binary_tensor_certificate(k,3,proper) for k in [3,4] for proper in [False,True]]
    assert [t['gluing_obstruction_dimension'] for t in tensors] == [1,1,0,0]
    assert all(t['coefficient_kernel_dimension']==1 for t in tensors)
    witness = cubic_offset_obstruction()
    moments = []
    for p in [5,7]:
        print('Checking universal moment block at p=%d' % p,flush=True)
        moments.append(verify_moment_block(p))
        print(json.dumps(moments[-1]),flush=True)
    higher = []
    for k in [3,4]:
        basis,J,C=sector(13,k,raw=True)
        rational=[int(A.change_ring(QQ).rank()) for A in [J,C]]
        modular=[int(A.change_ring(GF(13)).rank()) for A in [J,C]]
        assert [a-b for a,b in zip(rational,modular)]==([1,1] if k==3 else [0,0])
        higher.append(dict(p=13,degree=16,coefficient_degree=k,
                           rational_ranks=rational,modular_ranks=modular,
                           smith_factors_computed=False))
    masses=[]
    for m,p in [(5,5),(25,5),(14,7)]:
        _,basis,J,C,_=build_image(m,p+3,'moment')
        factors=primary(J,C,p)
        assert factors==[int(p**ZZ(m).valuation(p))]*expected(p+3)
        masses.append(dict(m=m,p=p,degree=p+3,ambient_rank=len(basis),
                           actual_unrooted_moment_factors=factors))
    names = ['TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md','tools/verify_type_a_k_next_layer.py',
             'tools/binary_tensor_certificates.py','tools/verify_type_a_k_stable_multiplicities.py',
             'tools/verify_type_a_k_first_correction.py']
    report = dict(result='passed',tensor_certificates=tensors,explicit_obstruction=witness,
                  moment_blocks=moments, higher_weight_rank_control=higher,
                  actual_mass_moment_images=masses,
                  actual_elementary_degree_p_plus_three_formula=False,
                  all_rank_claim_from_finite_checks=False,
                  source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                                 for name in names})
    (ROOT/'type-a-k-next-layer-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',moment_blocks=len(moments))))


if __name__=='__main__':
    main()
