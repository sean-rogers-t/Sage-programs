"""Finite checks of the type C joint-image proof and its simultaneous basis.

Run from borel/: PYTHONPATH=src sage -python experiments/type_c_joint_primitivity.py
The proof is in TYPE_C_JOINT_PRIMITIVITY.md; this script only corroborates it.
"""

from functools import lru_cache
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path
import sys
import time

from sage.all import ZZ, SchubertPolynomialRing, matrix, prod
from sage.env import SAGE_VERSION
from sage_borel import two_step_symplectic

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tests.oracle import OrthonormalBGG


def gaussian(n, k):
    """Gaussian coefficients from partition counting, not a rational product."""
    @lru_cache(None)
    def coefficients(rows, width):
        if rows == 0 or width == 0:
            return (1,)
        first, second = coefficients(rows, width-1), coefficients(rows-1, width)
        out = [0]*(rows*width+1)
        for i, x in enumerate(first):
            out[i] += x
        for i, x in enumerate(second):
            out[i+width] += x
        return tuple(out)
    return coefficients(k, n-k)


def multiply(left, right):
    out = [0]*(len(left)+len(right)-1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i+j] += x*y
    return out


def doubled(poly):
    out = [0]*(2*len(poly)-1)
    out[::2] = poly
    return out


def ranks(n, b, k):
    fiber = [1]
    for i in range(1, b+1):
        fiber = multiply(fiber, [1]+[0]*(i-1)+[1])
    middle = multiply(fiber, doubled(gaussian(n, b)))
    tail = doubled(gaussian(n, k))
    joint = multiply(middle, doubled(gaussian(n-b, k)))
    return middle, tail, joint


def qtilde(parts, chern):
    """Integral Pfaffian representative for a Lagrangian Schubert class."""
    def c(i):
        return chern[i] if 0 <= i < len(chern) else 0
    parts = tuple(parts)
    if len(parts) % 2:
        parts += (0,)

    @lru_cache(None)
    def pfaffian(indices):
        if not indices:
            return chern[0]
        i = indices[0]
        answer = 0
        for position in range(1, len(indices)):
            j = indices[position]
            entry = c(i)*c(j) + 2*sum((-1)**p*c(i+p)*c(j-p) for p in range(1, j+1))
            answer += (-1)**(position+1)*entry*pfaffian(indices[1:position]+indices[position+1:])
        return answer
    return pfaffian(parts)


def basis_polynomials(flag, max_degree):
    x = flag.classical_coordinates()
    # Reorder M,A,T, so BOTH selected quaternionic factors are standard
    # Grassmannian projections from Fl(b,n-k;n).
    z = tuple(t*t for t in x[flag.a:flag.a+flag.b]+x[:flag.a]+x[-flag.k:])
    chern = [flag.ring.one()]+list(flag.middle_generators().values())
    strict = [(tuple(reversed(parts)), sum(parts))
              for size in range(flag.b+1) for parts in combinations(range(1, flag.b+1), size)]
    lagrangian = [(parts, degree, qtilde(parts, chern)) for parts, degree in strict]
    schubert = SchubertPolynomialRing(ZZ)
    result = {d: [] for d in range(max_degree+1)}
    for word in permutations(range(1, flag.n+1)):
        descents = {i+1 for i in range(flag.n-1) if word[i] > word[i+1]}
        if not descents <= {flag.b, flag.n-flag.k}:
            continue
        length = sum(word[i] > word[j] for i in range(flag.n) for j in range(i+1, flag.n))
        if 2*length > max_degree:
            continue
        polynomial = schubert(list(word)).expand()
        base = polynomial(*z[:polynomial.parent().ngens()])
        for parts, degree, u in lagrangian:
            d = 2*length+degree
            if d <= max_degree:
                result[d].append((parts, word, descents <= {flag.b},
                                  not parts and descents <= {flag.n-flag.k}, u*base))
    return result


def main():
    start = time.time()
    historical = json.loads((ROOT/'two-step-integral-results.json').read_text())
    historical_rows = 0
    for case in historical['cases']:
        if case['type'] != 'C':
            continue
        coefficients = ranks(case['n'], case['b'], case['k'])
        for row in case['degrees']:
            d = row['degree']
            for kind, poly in zip(('middle','tail','joint'), coefficients):
                assert row[kind]['rank'] == (poly[d] if d < len(poly) else 0)
            assert row['joint_mod_sum']['torsion_factors'] == []
            historical_rows += 1
    # Several new block sizes, including middle size three and a balanced C6.
    cases = [(1,1,1,9), (1,2,1,8), (2,2,1,8), (1,3,1,8),
             (2,1,2,8), (1,1,3,8), (2,2,2,8), (1,3,2,6), (2,3,2,6)]
    records = []
    basis_columns = oracle_columns = 0
    oracle = OrthonormalBGG('C3')
    for a,b,k,limit in cases:
        flag = two_step_symplectic(a+b+k,b,k)
        middle,tail,joint = flag.middle_algebra(), flag.tail_algebra(), flag.joint_algebra()
        expected = ranks(flag.n,b,k)
        bases = basis_polynomials(flag,limit)
        rows = []
        for d in range(limit+1):
            reports = [algebra.integral_degree(d) for algebra in (middle,tail,joint)]
            m,t,j = reports
            for report, poly in zip(reports,expected):
                assert report.rank == (poly[d] if d < len(poly) else 0)
                assert report.smith_factors == (1,)*report.rank
            separate = m.compare(t)
            assert separate.intersection.rank == (1 if d == 0 else 0)
            assert separate.sum.is_saturated
            quotient = j.lattice().quotient_by(separate.sum)
            assert quotient.torsion_factors == ()
            assert quotient.free_rank == j.rank-m.rank-t.rank+(1 if d==0 else 0)
            images = []
            for parts, word, in_middle, in_tail, polynomial in bases[d]:
                image = flag.expand(polynomial)
                assert j.contains(image)
                assert m.contains(image) == in_middle
                assert t.contains(image) == in_tail
                assert separate.sum.contains(image) == (in_middle or in_tail)
                if flag.n == 3:
                    independent = oracle.expansion(oracle.from_fundamental(polynomial),d)
                    actual = {tuple(w.reduced_word()):c for w,c in image.items()}
                    expected_image = {tuple(w.reduced_word()):c for w,c in independent.items()}
                    assert actual == expected_image
                    oracle_columns += 1
                images.append(image)
                basis_columns += 1
            assert len(images) == j.rank
            support = tuple(sorted(set().union(*(set(v) for v in images)),
                                   key=lambda w: tuple(w.reduced_word())))
            basis = matrix(ZZ,len(support),len(images),
                           {(i,h):value[w] for i,w in enumerate(support) for h,value in enumerate(images)})
            assert basis.rank() == len(images)
            if images:
                diagonal = basis.smith_form()[0]
                assert all(abs(diagonal[i,i])==1 for i in range(len(images)))
            rows.append({'degree':d,'middle_rank':int(m.rank),'tail_rank':int(t.rank),
                         'joint_rank':int(j.rank),'mixed_quotient_rank':int(quotient.free_rank),
                         'basis_columns':len(images),'joint_and_sum_saturated':True})
        records.append({'a':a,'b':b,'k':k,'n':a+b+k,'degrees':rows})
        print('Verified (a,b,k)=%s through d=%d: %d simultaneous basis columns' %
              ((a,b,k),limit,sum(row['basis_columns'] for row in rows)),flush=True)
    source_files = [Path(__file__).resolve(), ROOT/'TYPE_C_JOINT_PRIMITIVITY.md',ROOT/'tests/oracle.py']
    source_files += sorted((ROOT/'src/sage_borel').glob('*.py'))
    report = {'result':'passed','sage_version':SAGE_VERSION,'seconds':round(time.time()-start,3),
              'command':'PYTHONPATH=src sage -python experiments/type_c_joint_primitivity.py',
              'scope':'Finite corroboration of the separate geometric proof; not a proof by enumeration.',
              'historical_degree_rows_checked':historical_rows,
              'fresh_degree_rows':sum(len(c['degrees']) for c in records),
              'simultaneous_basis_columns':basis_columns,'independent_cartesian_basis_columns':oracle_columns,
              'cases':records,'source_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):
                                            hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}}
    (ROOT/'experiments/type_c_joint_primitivity_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('All checks passed in %.3fs'%report['seconds'],flush=True)


if __name__=='__main__':
    main()
