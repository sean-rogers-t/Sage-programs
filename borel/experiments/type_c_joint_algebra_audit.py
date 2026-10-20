"""Independent finite check of the integral base-change step, without sage_borel.

Run from borel/: sage -python experiments/type_c_joint_algebra_audit.py
Uses formal Chern classes and integer determinants, not BGG or Schubert images.
The general argument is in TYPE_C_JOINT_PRIMITIVITY.md, not this enumeration.
"""

from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

from sage.all import PolynomialRing, ZZ, matrix, prod
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def exponents(weights, degree):
    if not weights:
        if degree == 0:
            yield ()
        return
    for power in range(degree // weights[0] + 1):
        for rest in exponents(weights[1:], degree - power * weights[0]):
            yield (power,) + rest


def check_block(b, limit):
    ring = PolynomialRing(ZZ, b, names=['c%s' % i for i in range(1, b + 1)])
    chern = (ring.one(),) + ring.gens()

    def c(i):
        return chern[i] if 0 <= i <= b else ring.zero()

    # Extract squared elementary classes from C(v)C(-v), independently
    # of the Pfaffian implementation and the library's root coordinates.
    p = tuple(sum((-1)**(i + r) * c(r) * c(2*i-r)
                  for r in range(2*i+1)) for i in range(1, b+1))

    @lru_cache(None)
    def q(parts):
        if not parts:
            return ring.one()
        if len(parts) == 1:
            return c(parts[0])
        if len(parts) % 2:
            return sum((-1)**i * c(part) * q(parts[:i] + parts[i+1:])
                       for i, part in enumerate(parts))
        left = parts[0]
        return sum((-1)**(j-1) * (
            c(left)*c(right) + 2*sum((-1)**r*c(left+r)*c(right-r)
                                    for r in range(1, right+1)))
            * q(parts[1:j]+parts[j+1:])
            for j, right in enumerate(parts[1:], 1))

    records = []
    for degree in range(limit+1):
        monomials = tuple(exponents(tuple(range(1, b+1)), degree))
        columns = []
        for powers in monomials:
            parts = tuple(i for i in range(b, 0, -1) for _ in range(powers[i-1]))
            strict = tuple(i for i in range(b, 0, -1) if powers[i-1] % 2)
            factored = prod(p[i]**(powers[i]//2) for i in range(b)) * q(strict)
            assert q(parts) == factored, (b, degree, parts)
            # Sage ETuple keys do not have Python tuple lookup semantics.
            columns.append({tuple(key): value for key, value in ring(factored).dict().items()})
        transition = matrix(ZZ, len(monomials), len(monomials),
                            lambda i, j: columns[j].get(monomials[i], 0))
        determinant = transition.det()
        assert abs(determinant) == 1, (b, degree, determinant)
        records.append({'degree': degree, 'rank': len(monomials),
                        'determinant': int(determinant)})
    return {'b': b, 'max_degree': limit, 'degrees': records}


def main():
    start = time.time()
    cases = []
    for b in range(1, 6):
        cases.append(check_block(b, 14))
        print('Verified b=%d through polynomial degree 14' % b, flush=True)
    sources = [Path(__file__).resolve(), ROOT/'TYPE_C_JOINT_PRIMITIVITY.md',
               ROOT/'TYPE_C_JOINT_AUDIT.md']
    report = {
        'result': 'passed', 'sage_version': SAGE_VERSION,
        'seconds': round(time.time()-start, 3),
        'command': 'sage -python experiments/type_c_joint_algebra_audit.py',
        'scope': 'Finite formal-polynomial audit; no sage_borel imports. Does not prove arbitrary rank.',
        'unimodular_transition_matrices': sum(len(case['degrees']) for case in cases),
        'pair_removal_identities': sum(row['rank'] for case in cases for row in case['degrees']),
        'cases': cases,
        'source_sha256': {str(p.relative_to(ROOT)).replace('\\', '/'):
                          hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    }
    (ROOT/'experiments/type_c_joint_algebra_audit_verification.json').write_text(
        json.dumps(report, indent=2)+'\n')
    print('All checks passed: %d matrices, %d identities in %.3fs' % (
        report['unimodular_transition_matrices'], report['pair_removal_identities'],
        report['seconds']), flush=True)


if __name__ == '__main__':
    main()
