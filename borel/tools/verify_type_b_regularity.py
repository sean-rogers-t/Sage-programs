"""Exact corroboration of the uniform type B central--Weyl obstruction.

Run with Sage Python. Freudenthal weights are checked against an independent
Weyl determinant and its derivatives, exterior monomials and spin sign cubes.
All witness calculations take place in Q(i), without floating point roots.
"""
from collections import defaultdict
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

from sage.all import QQ, QuadraticField, WeylCharacterRing, matrix, prod
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]
K = QuadraticField(-1, 'ii')
ii = K.gen()


def labels(n, bound):
    for a in product(range(bound+1), repeat=n):
        if 0 < sum(a) <= bound:
            yield a


def oracle(n, a, rings):
    if n not in rings:
        rings[n] = WeylCharacterRing('B'+str(n), style='coroots')
    character = rings[n](tuple(a))
    weights = {}
    for w, m in character.weight_multiplicities().items():
        doubled = tuple(2*v for v in w.to_vector())
        assert all(v == int(v) for v in doubled)
        weights[tuple(int(v) for v in doubled)] = int(m)
    p = tuple(2*sum(a[j:n-1])+a[-1] for j in range(n))
    assert weights[p] == 1 and sum(weights.values()) == character.degree()
    assert all(all(v % 2 == a[-1] % 2 for v in w) for w in weights)
    return p, weights


def trace(weights, point):
    n = len(point)
    gram = matrix(QQ, n, n, lambda i,j: sum(
        QQ(m)*w[i]*w[j]/4 for w,m in weights.items()))
    S = gram[0,0]
    assert S > 0 and gram == S*matrix.identity(QQ,n)
    values = [(w,m,prod(a**v for a,v in zip(point,w))) for w,m in weights.items()]
    character = sum(m*v for w,m,v in values)
    gradient = [sum(QQ(m)*w[j]*v/2 for w,m,v in values) for j in range(n)]
    return character, gradient, S, [v/S for v in gradient]


def alternant(p, point):
    """Weyl character and a_j/2 times its a_j derivative, no weight oracle."""
    n = len(p)
    def determinant_and_derivatives(exponents):
        M = matrix(K, [[a**k-a**(-k) for a in point] for k in exponents])
        derivatives = []
        for j,a in enumerate(point):
            N = matrix(M)
            for i,k in enumerate(exponents):
                N[i,j] = QQ(k)/2*(a**k+a**(-k))
            derivatives.append(N.det())
        return M.det(), derivatives
    numerator, dn = determinant_and_derivatives([p[j]+2*n-2*j-1 for j in range(n)])
    denominator, dd = determinant_and_derivatives([2*n-2*j-1 for j in range(n)])
    assert denominator
    return numerator/denominator, [(a*denominator-numerator*b)/denominator**2
                                  for a,b in zip(dn,dd)]


def regular(point):
    z = [a*a for a in point]
    assert all(a != 1 for a in z)
    assert all(z[i] != z[j] and z[i]*z[j] != 1
               for i in range(len(z)) for j in range(i+1,len(z)))


def check(n, a, p, weights):
    witness = [ii]+[K(j) for j in range(2,n+1)]
    regular(witness)
    row = {}
    for name,point in [('positive',[K(j) for j in range(2,n+2)]),('witness',witness)]:
        chi, gradient, S, image = trace(weights,point)
        chi2, gradient2 = alternant(p,point)
        assert chi == chi2 and gradient == gradient2
        if name == 'witness':
            parity = (-1)**a[-1]
            inverted = [1/point[0]]+point[1:]
            for w in weights:
                ev = prod(v**d for v,d in zip(point,w))
                ev2 = prod(v**d for v,d in zip(inverted,w))
                assert ev2 == parity*ev
            if parity == 1:
                assert image[0] == 0
            else:
                assert all(y == 0 for y in image[1:])
            row = dict(rank=n, dynkin_labels=list(a), dimension=sum(weights.values()),
                       central_scalar=parity, trace_Gram_scalar=str(S),
                       image=[str(y) for y in image],
                       vanishing_root='epsilon_1' if parity == 1 else 'epsilon_2')
    return row


def exterior_weights(n,k):
    vector = [(0,)*n]
    for j in range(n):
        for sign in [-2,2]:
            vector.append(tuple(sign if i == j else 0 for i in range(n)))
    result = defaultdict(int)
    for subset in combinations(vector,k):
        result[tuple(sum(w[i] for w in subset) for i in range(n))] += 1
    return dict(result)


def main():
    rings, cache, rows = {}, {}, []
    exterior_count = spin_count = 0
    for n,bound in [(2,4),(3,3),(4,2),(5,1)]:
        for a in labels(n,bound):
            p, weights = oracle(n,a,rings)
            cache[n,a] = p,weights
            rows.append(check(n,a,p,weights))
        for k in range(1,n):
            a = tuple(1 if j == k-1 else 0 for j in range(n))
            assert exterior_weights(n,k) == cache[n,a][1]
            exterior_count += 1
        # Top vector exterior power has highest weight 2 omega_n, not omega_n.
        if bound >= 2:
            a = (0,)*(n-1)+(2,)
            assert exterior_weights(n,n) == cache[n,a][1]
            exterior_count += 1
        a = (0,)*(n-1)+(1,)
        assert {w:1 for w in product([-1,1],repeat=n)} == cache[n,a][1]
        spin_count += 1
    # Higher-rank basic spin cases: an explicit third construction of weights.
    for n in range(6,9):
        a = (0,)*(n-1)+(1,)
        p,weights = oracle(n,a,rings)
        assert {w:1 for w in product([-1,1],repeat=n)} == weights
        rows.append(check(n,a,p,weights))
        spin_count += 1
    sums = []
    for n in range(2,6):
        for last in [0,1]:
            a = (1,)+(0,)*(n-2)+(last,)
            b = (0,)*(n-1)+(2+last,)
            _,wa = oracle(n,a,rings)
            _,wb = oracle(n,b,rings)
            weights = defaultdict(int,wa)
            for w,m in wb.items():
                weights[w] += m
            if last == 0:
                weights[(0,)*n] += 3
            point = [ii]+[K(j) for j in range(2,n+1)]
            _,_,S,image = trace(weights,point)
            if last == 0:
                assert image[0] == 0
            else:
                assert all(y == 0 for y in image[1:])
            # Independent additivity uses each actual trace form, not an
            # unweighted sum of the separately normalized Springer maps.
            pa = tuple(2*sum(a[j:n-1])+a[-1] for j in range(n))
            pb = tuple(2*sum(b[j:n-1])+b[-1] for j in range(n))
            da, db = alternant(pa,point)[1],alternant(pb,point)[1]
            assert image == [(x+y)/S for x,y in zip(da,db)]
            sums.append(dict(rank=n, central_scalar=(-1)**last,
                             summands=[list(a),list(b)], trivial_multiplicity=3 if last == 0 else 0,
                             image=[str(y) for y in image]))
    # Rank-one counter-control, and the explicit spin formula at every rank.
    spin_formula_checks = []
    for n in range(1,9):
        weights = {w:1 for w in product([-1,1],repeat=n)}
        point = [ii]+[K(j) for j in range(2,n+1)]
        _,_,S,image = trace(weights,point)
        assert S == QQ(2)**(n-2)
        expected = [QQ(2)**(1-n)*(point[j]-1/point[j])*
                    prod(point[k]+1/point[k] for k in range(n) if k != j)
                    for j in range(n)]
        assert image == expected and image[0] != 0
        assert all(y == 0 for y in image[1:])
        spin_formula_checks.append(dict(rank=n,image=[str(y) for y in image]))
    assert spin_formula_checks[0]['image'] == ['2*ii']
    assert spin_formula_checks[1]['image'] == ['5/2*ii','0']
    paths = ['TYPE_B_REGULARITY.md','tools/verify_type_b_regularity.py']
    report = dict(result='passed', date='2026-10-22',sage_version=SAGE_VERSION,
        arbitrary_rank_proof='TYPE_B_REGULARITY.md, Sections 2-4',
        group='Spin_(2n+1)(C), n>=2; and all nontrivial SO_(2n+1) representations',
        character_cases=len(rows),positive_and_witness_alternant_checks=2*len(rows),
        universal_Qi_trace_checks=len(rows),
        integer_weight_cases=sum(r['central_scalar']==1 for r in rows),
        spinorial_cases=sum(r['central_scalar']==-1 for r in rows),
        exterior_monomial_character_checks=exterior_count,
        independent_spin_sign_cube_checks=spin_count,
        same_central_parity_reducible_checks=len(sums),
        closed_spin_formula_checks=len(spin_formula_checks),
        finite_Dynkin_height_bounds={'2':4,'3':3,'4':2,'5':1},
        additional_basic_spin_ranks=[6,7,8],rank_one_spin_positive_control=True,
        witnesses=rows,reducible_witnesses=sums,spin_formula_controls=spin_formula_checks,
        runtime_API_changed=False,novelty_claim=False,
        scope='Exact finite corroboration; not extrapolation, external review or a mixed-parity classification.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'type-b-regularity-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('witnesses','reducible_witnesses','spin_formula_controls','source_sha256')},indent=2))


if __name__ == '__main__':
    main()
