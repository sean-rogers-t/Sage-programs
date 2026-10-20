"""Finite corroboration of ORTHOGONAL_JOINT_IMAGES.md, not an all-rank proof.

Run from borel/: PYTHONPATH=src sage -python experiments/orthogonal_joint_structure.py
Checks integral kernels, image comparisons, Gaussian ranks and Cartesian BGG.
"""

from functools import lru_cache
from itertools import combinations
import hashlib
import json
from pathlib import Path
import sys
import time

from sage.all import PolynomialRing, ZZ, matrix, prod
from sage.env import SAGE_VERSION
from sage_borel import BorelMap, GeneratorAlgebra

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.oracle import OrthonormalBGG

POLYS = PolynomialRing(ZZ, 't')
t = POLYS.gen()


@lru_cache(None)
def gaussian(n, r):
    if r < 0 or r > n:
        return POLYS.zero()
    if r == 0 or r == n:
        return POLYS.one()
    return gaussian(n-1, r-1) + t**r * gaussian(n-1, r)


def ranks(kind, a, b, k):
    n = a+b+k
    G = lambda N, r: gaussian(N, r)(t**2)
    Q = G(n, b)*G(n-b, k)
    if kind == 'B':
        F = prod(1+t**j for j in range(1, b+1))
        return F*G(n, b), G(n, k), F*Q
    F = prod(1+t**j for j in range(1, b))
    middle = F*(G(n, b)+t**b*G(n-1, b))
    tail = G(n, k)+t**k*G(n-1, k)
    joint = F*(Q + t**b*G(n-1, b)*G(n-b, k)
               + t**k*G(n-1, k)*G(n-k, b)
               + t**(b+k)*G(n-1, b)*G(n-1-b, k))
    return POLYS(middle), POLYS(tail), POLYS(joint)


def elementary(values, j):
    return sum(prod(term) for term in combinations(values, j))


def algebras(kind, a, b, k):
    flag = BorelMap('%s%d' % (kind, a+b+k), omitted_nodes=(a, a+b))
    x = flag.classical_coordinates()
    middle = {'c%d' % j: elementary(x[a:a+b], j) for j in range(1, b+1)}
    tail = {'e%d' % j: elementary(tuple(y*y for y in x[-k:]), j)
            for j in range(1, k+1 if kind == 'B' else k)}
    if kind == 'D':
        tail['p'] = prod(x[-k:])
    return flag, tuple(GeneratorAlgebra(flag, inputs)
                       for inputs in (middle, tail, dict(middle, **tail)))


def relations(algebra, kind, a, b, k):
    ring = algebra.domain
    c = [ring.one()] + list(ring.gens()[:b])
    middle_squared = [ring.one()] + [
        sum((-1)**(i+r)*c[r]*c[2*i-r]
            for r in range(max(0, 2*i-b), min(b, 2*i)+1))
        for i in range(1, b+1)]
    tail = [ring.one()] + list(ring.gens()[b:])
    if kind == 'D':
        tail[-1] = tail[-1]**2
    total = [sum(middle_squared[i]*tail[j] for i in range(b+1)
                 for j in range(k+1) if i+j == r) for r in range(b+k+1)]
    # h_r is the coefficient of 1/E(-z); all coefficients remain integral.
    h = [ring.one()]
    for r in range(1, a+b+k+1):
        h.append(sum((-1)**(i+1)*total[i]*h[r-i]
                     for i in range(1, min(r, b+k)+1)))
    result = [(2*r, h[r]) for r in range(a+1, a+b+k+1)]
    if kind == 'D':
        result.append((2*a+b+k, c[b]*ring.gen(b+k-1)*h[a]))
    return result


def check_kernel(algebra, report, generators):
    columns = []
    for degree, relation in generators:
        if degree <= report.degree:
            for monomial in algebra.monomials(report.degree-degree):
                polynomial = relation*monomial
                columns.append([ZZ(polynomial.monomial_coefficient(m)) for m in report.monomials])
    relations_matrix = matrix(ZZ, len(report.monomials), len(columns),
                              lambda i, j: columns[j][i])
    assert report.matrix*relations_matrix == 0
    rank = relations_matrix.rank()
    assert rank == report.kernel_dimension
    # A primitive subgroup of the correct rank inside the kernel is the
    # whole integer kernel, not just a rationally spanning relation list.
    if rank:
        diagonal = relations_matrix.smith_form()[0]
        assert all(abs(diagonal[i, i]) == 1 for i in range(rank))
    return int(rank)


def cartesian_check(algebra, report, oracle):
    inputs = [oracle.from_fundamental(f) for f in algebra.generator_polynomials]
    for monomial in report.monomials:
        powers = monomial.exponents()[0]
        polynomial = prod(f**e for f, e in zip(inputs, powers))
        independent = oracle.expansion(polynomial, report.degree)
        actual = algebra.evaluate(monomial)
        assert {tuple(w.reduced_word()): c for w, c in independent.items()} == {
            tuple(w.reduced_word()): c for w, c in actual.items()}
    return len(report.monomials)


def main():
    start = time.time()
    historical = json.loads((ROOT/'two-step-integral-results.json').read_text())
    historical_rows = 0
    for case in historical['cases']:
        if case['type'] not in ('B', 'D'):
            continue
        predicted = ranks(case['type'], case['a'], case['b'], case['k'])
        for row in case['degrees']:
            for key, polynomial in zip(('middle', 'tail', 'joint'), predicted):
                assert row[key]['rank'] == polynomial[row['degree']]
            historical_rows += 1
    cases = [('B',1,1,1,9), ('B',1,2,1,9), ('B',2,1,1,9),
             ('B',1,2,2,8), ('B',2,2,2,8),
             ('D',1,1,2,12), ('D',1,2,2,10), ('D',2,1,2,10),
             ('D',1,1,3,10), ('D',2,2,2,8), ('D',1,3,2,8), ('D',1,2,3,8)]
    records = []
    independent_columns = kernel_rows = 0
    oracles = {}
    for kind, a, b, k, limit in cases:
        flag, (middle, tail, joint) = algebras(kind, a, b, k)
        predicted = ranks(kind, a, b, k)
        generators = relations(joint, kind, a, b, k)
        rows = []
        for d in range(limit+1):
            m, e, j = [algebra.integral_degree(d) for algebra in (middle, tail, joint)]
            for report, polynomial in zip((m, e, j), predicted):
                assert report.rank == polynomial[d], (kind,a,b,k,d,report.rank,polynomial[d])
                assert all(int(f) & (int(f)-1) == 0 for f in report.smith_factors)
            comparison = m.compare(e)
            assert comparison.intersection.rank == (1 if d == 0 else 0)
            quotient = j.lattice().quotient_by(comparison.sum)
            assert quotient.torsion_factors == ()
            assert quotient.free_rank == j.rank-m.rank-e.rank+(1 if d == 0 else 0)
            relation_rank = check_kernel(joint, j, generators)
            kernel_rows += relation_rank
            independent = ((a+b+k <= 4 and d <= 6) or
                           (kind == 'B' and a+b+k == 3 and d == 7) or
                           (kind == 'D' and a+b+k == 4 and d == 8) or
                           (kind,a,b,k,d) == ('D',1,2,2,6))
            if independent:
                label = '%s%d' % (kind, a+b+k)
                if label not in oracles:
                    oracles[label] = OrthonormalBGG(label)
                independent_columns += cartesian_check(joint, j, oracles[label])
            rows.append({'degree': d, 'middle_rank': int(m.rank), 'tail_rank': int(e.rank),
                         'joint_rank': int(j.rank), 'joint_smith_factors': list(map(int,j.smith_factors)),
                         'joint_mod_sum_free_rank': int(quotient.free_rank),
                         'joint_mod_sum_torsion': [], 'integral_kernel_rank': relation_rank})
        records.append({'type': kind, 'a': a, 'b': b, 'k': k, 'max_degree': limit, 'degrees': rows})
        print('Verified %s (a,b,k)=%s through d=%d' % (kind,(a,b,k),limit), flush=True)
    sources = [Path(__file__).resolve(), ROOT/'ORTHOGONAL_JOINT_IMAGES.md',
               ROOT/'TYPE_C_JOINT_PRIMITIVITY.md', ROOT/'tests/oracle.py',
               ROOT/'two-step-integral-results.json'] + sorted((ROOT/'src/sage_borel').glob('*.py'))
    report = {'result': 'passed', 'sage_version': SAGE_VERSION,
              'seconds': round(time.time()-start,3),
              'command': 'PYTHONPATH=src sage -python experiments/orthogonal_joint_structure.py',
              'scope': 'Finite checks of the separate all-rank argument, not a proof by enumeration.',
              'historical_degree_rows': historical_rows,
              'fresh_degree_rows': sum(len(case['degrees']) for case in records),
              'independent_cartesian_columns': independent_columns,
              'integral_kernel_ranks_total': kernel_rows,
              'cases': records,
              'source_sha256': {str(p.relative_to(ROOT)).replace('\\','/'):
                                hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    (ROOT/'experiments/orthogonal_joint_structure_verification.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print('All checks passed in %.3fs' % report['seconds'], flush=True)


if __name__ == '__main__':
    main()
