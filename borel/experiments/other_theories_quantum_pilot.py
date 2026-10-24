"""Integral quantum projective-space images, independently checked by gcd DP.

This is a research pilot, not a production quantum Schubert backend.
The ordinary closed projective image formula is an input proved elsewhere.
"""
from hashlib import sha256
import json
from math import comb, gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def coefficients(m, through):
    # All monomials in the actual Levi generators E_j, without using the
    # claimed closed image. A quantum power q can also be supplied freely.
    exact = [1]
    quantum = []
    for d in range(through + 1):
        if d:
            value = 0
            for j in range(1, min(m - 1, d) + 1):
                value = gcd(value, comb(m, j) * exact[d-j])
            exact.append(value)
        value = 0
        for power in range(d // m + 1):
            value = gcd(value, exact[d - power * m])
        quantum.append(value)
    return exact, quantum


def main():
    reports = []; degrees = 0; products = 0
    for m in range(2, 81):
        through = 8*m
        exact, quantum = coefficients(m, through)
        factors = [m // gcd(m, r) for r in range(m)]
        assert factors[0] == 1
        assert exact[:m] == factors
        assert quantum == [factors[d % m] for d in range(through+1)]
        for a in range(m):
            for b in range(m):
                assert (factors[a]*factors[b]) % factors[(a+b) % m] == 0
                products += 1
        degrees += through+1
        reports.append(dict(m=m,degree_through=through,
                            image_basis_coefficients=factors,
                            faithful_q_zero_specialization=True))
    report = dict(result='passed',date='2026-10-22',
                  uniform_family='SL_m standard representation, literal root-lattice Levi source, P1',
                  quantum_ring='Z[q,h]/(h^m-q)',
                  generator_images='E_j -> binomial(m,j)*h^j, 1<=j<m',
                  independent_degree_gcd_cases=degrees,
                  multiplication_integrality_checks=products,cases=reports,
                  scope='Direct actual-generator recurrence checks the written arbitrary-m image formula; not a computation of a new ambient quantum ring.',
                  source_sha256={'experiments/other_theories_quantum_pilot.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-quantum-pilot.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED',degrees,'degree gcd checks and',products,'multiplication checks')


if __name__ == '__main__':
    main()
