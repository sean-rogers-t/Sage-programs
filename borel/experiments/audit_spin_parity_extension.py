"""Native independent audit: rational membership + finite flatness, no Sage.

The supplied integer polynomials are certificates, not trusted ideal generators.
We prove they generate the whole local kernel from their rational membership,
their reduced colength 2**n, and the known rational Borel rank. Thus we do not
repeat or trust the primary integer saturation algorithm.
"""
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    primary_path = ROOT/'experiments/spin-parity-extension.json'
    primary = json.loads(primary_path.read_text())
    cases = []
    for case in primary['cases']:
        n = case['rank']
        xs = sp.symbols(' '.join(case['source_variables']))
        z = xs[0]; cs = [sp.Integer(1), 2*z]+list(xs[1:])+[sp.Integer(0)]*n
        rels = [sp.expand(sum((-1)**k*cs[k]*cs[2*j-k]
                             for k in range(max(0, 2*j-n), min(n, 2*j)+1)))
                for j in range(1, n+1)]
        lifts = [sp.sympify(s.replace('^', '**'), locals=dict(zip(case['source_variables'], xs)))
                 for s in case['integral_lifts']]
        print('rational membership', n, flush=True)
        rational = sp.groebner(rels, *reversed(xs), order='grevlex', domain=sp.QQ)
        assert all(rational.reduce(f)[1] == 0 for f in lifts)
        modular = sp.groebner(lifts, *xs, order='grevlex', modulus=2)
        leading = [p.LM(order=modular.order).exponents for p in modular.polys]
        bounds = [min(e[j] for e in leading if e[j] and sum(e) == e[j]) for j in range(n)]
        standard = [e for e in product(*(range(b) for b in bounds))
                    if not any(all(a >= b for a, b in zip(e, lm)) for lm in leading)]
        assert len(standard) == 2**n
        weights = case['source_weights']
        hilbert = [0]*(n*(n+1)//2+1)
        for e in standard:
            hilbert[sum(a*b for a, b in zip(e, weights))] += 1
        expected = [1]
        for j in range(1, n+1):
            old = expected[:]; expected += [0]*j
            for k, coefficient in enumerate(old): expected[k+j] += coefficient
        assert hilbert == expected == case['symmetric_hilbert_coefficients']
        failures = []
        for f in lifts:
            remainder = modular.reduce(sp.diff(f, z))[1]
            if remainder != 0:
                terms = sp.Poly(f, *xs).monoms()
                degree = max(sum(a*b for a, b in zip(e, weights)) for e in terms)
                failures.append(dict(degree=degree, derivative_remainder=str(remainder)))
        first = min((f['degree'] for f in failures), default=None)
        assert (not failures) == case['entire_kernel_preserved']
        assert first == case['first_failure_degree']
        # Audited deduction: each lift lies in K; reduced quotient is finite
        # of dimension 2**n. Graded Nakayama gives finiteness over Z_(2).
        # Its rational rank is >= 2**n, by the Borel quotient. The structure
        # theorem over Z_(2) forces that rank exactly and zero torsion. Thus
        # the supplied ideal equals K, in every degree, not just sampled ones.
        cases.append(dict(rank=n, entire_kernel_preserved=not failures,
                          first_failure_degree=first, rational_memberships=len(lifts),
                          reduced_quotient_length=len(standard),
                          finite_flatness_certificate=True, failures=failures))
        print('certified', n, first, flush=True)
    # Closed Spin9 presentation, independently specified rather than read.
    z, c2, c3, c4 = sp.symbols('z c2 c3 c4')
    closed = [2*z**2-c2, z**2*c2-2*z*c3+c4, 2*c2*c4-c3**2,
              z**8+z*c3*c4-c2**2*c4]
    assert sp.expand(2*closed[1]-(c2**2-4*z*c3+2*c4)-c2*closed[0]) == 0
    assert sp.expand(closed[2]+c3**2-2*c2*c4) == 0
    lower = sp.groebner(closed[:3], c2, c4, c3, z, domain=sp.QQ)
    assert lower.reduce(closed[3]-c4**2/4)[1] == 0
    G = sp.groebner(closed, z, c2, c3, c4, modulus=2)
    assert all(G.reduce(f)[1] == 0 for f in [c2, c4, c3**2, z**8])
    assert all(G.reduce(sp.diff(f, z))[1] == 0 for f in closed)
    report = dict(result='passed', date='2026-10-22',
                  method='SymPy rational Groebner membership, independent modular Groebner and staircase, finite-flatness proof',
                  shared_inputs='Root datum, rational Borel theorem and integer certificate polynomials; no Sage, integer saturation or runtime Smith data.',
                  all_rank_quartic_statement='The written symmetric-source proof establishes failure in degree four iff n=3; the finite certificates do not prove the arbitrary-rank assertion.',
                  Spin9_closed_presentation_checked=True,
                  source_sha256={
                      'experiments/audit_spin_parity_extension.py': digest(Path(__file__)),
                      'experiments/spin-parity-extension.json': digest(primary_path)}, cases=cases)
    (ROOT/'experiments/spin-parity-extension-audit.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
