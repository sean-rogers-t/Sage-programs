"""Exact entire-kernel certificates for odd spin full flags, ranks 3--8.

The source is the actual simply connected weight lattice, localized at 2.
Eliminate a variable satisfying 2*t=1 over ZZ, then retain a small list of
integral lifts generating the entire reduced ideal. A separate native audit
uses rational membership and a finite-flatness argument instead of saturation.
"""
from hashlib import sha256
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, ZZ, prod

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def degree(poly, weights):
    return max((sum(a*b for a, b in zip(e, weights)) for e in poly.dict()), default=-1)


def staircase(gb, ring):
    leading = [tuple(g.lm().exponents()[0]) for g in gb]
    bounds = []
    for j in range(ring.ngens()):
        bounds.append(min(e[j] for e in leading if e[j] and sum(e) == e[j]))
    exponents = [tuple([0]*ring.ngens())]
    for j, bound in enumerate(bounds):
        exponents = [e[:j]+(k,)+e[j+1:] for e in exponents for k in range(bound)]
    return [e for e in exponents if not any(all(a >= b for a, b in zip(e, lm)) for lm in leading)]


def main():
    cases = []
    for n in range(3, 9):
        names = ['z']+['c'+str(i) for i in range(2, n+1)]
        weights = [1]+list(range(2, n+1))
        P = PolynomialRing(ZZ, n+1, names=['t']+names, order='lex')
        t, z = P.gens()[:2]
        cs = [P(1), 2*z]+list(P.gens()[2:])+[P(0)]*n
        relations = [sum((-1)**k*cs[k]*cs[2*j-k]
                         for k in range(max(0, 2*j-n), min(n, 2*j)+1))
                     for j in range(1, n+1)]
        saturated = [g for g in P.ideal(relations+[2*t-1]).groebner_basis() if g.degree(t) == 0]
        F = PolynomialRing(GF(2), n, names=names, order='degrevlex')
        reduced = F.ideal([F(str(g)) for g in saturated])
        chosen = []
        current = F.ideal([])
        for g in sorted(saturated, key=lambda g: (degree(g, [0]+weights), len(g.dict()), str(g))):
            bar = F(str(g))
            if current.reduce(bar):
                chosen.append(g)
                current = F.ideal([F(str(h)) for h in chosen])
        assert current == reduced
        basis = list(current.groebner_basis())
        standard = staircase(basis, F)
        assert len(standard) == 2**n
        hilbert = [0]*(n*(n+1)//2+1)
        for e in standard:
            hilbert[sum(a*b for a, b in zip(e, weights))] += 1
        Q = PolynomialRing(ZZ, 'q'); q = Q.gen()
        expected = prod(1+q**j for j in range(1, n+1))
        assert hilbert == expected.list()
        failures = []
        for g in chosen:
            bar = F(str(g)); rem = current.reduce(bar.derivative(F.gen(0)))
            if rem:
                failures.append(dict(degree=degree(g, [0]+weights), integral_relation=str(g),
                                     derivative_remainder=str(rem)))
        cases.append(dict(rank=n, group='Spin('+str(2*n+1)+')', parabolic='B',
                          entire_kernel_preserved=not failures,
                          first_failure_degree=min((f['degree'] for f in failures), default=None),
                          source_variables=names, source_weights=weights,
                          integral_lifts=[str(g) for g in chosen],
                          reduced_groebner_basis=[str(g) for g in basis],
                          symmetric_quotient_rank=len(standard),
                          symmetric_hilbert_coefficients=hilbert,
                          failures=failures))
        print(n, 'stable' if not failures else 'first failure '+str(cases[-1]['first_failure_degree']), flush=True)
    result = dict(result='passed', date='2026-10-22', coefficient_ring='Z_(2)',
                  source='Sym of the full Spin weight lattice; E=R[z,c2,...,cn] and S free over E',
                  method='integer elimination (r_j,2*t-1), entire reduced ideal, Artin source reduction',
                  scope='Entire kernels for the six specified ranks, not an arbitrary-rank classification.',
                  source_sha256={'experiments/spin_parity_extension.py': digest(Path(__file__))}, cases=cases)
    (ROOT/'experiments/spin-parity-extension.json').write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
