"""Independent native verification of the Pieri lattice certificates.

SymPy column Hermite forms independently establish the image lattices.
Exact triangular integer division and bit arithmetic check the derivative
factorizations, or a relation with nonzero derivative. No Sage is used.
The all-rank lowering proof is checked separately by coefficient identities
and sparse raising/lowering commutators, not inferred from the finite table.
"""
from hashlib import sha256
import gzip
from itertools import combinations
import json
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    raw = path.read_bytes()
    return sha256(raw if path.suffix == '.gz' else raw.replace(b'\r\n', b'\n')).hexdigest()


def levels(n):
    result = [[] for _ in range(n*(n+1)//2+1)]
    for length in range(n+1):
        for p in combinations(range(1, n+1), length): result[sum(p)].append(p)
    return result


def covers(parts, n):
    if not parts or parts[0] != 1: yield (1,)+parts
    for j, part in enumerate(parts):
        if part != n and part+1 not in parts:
            yield parts[:j]+(part+1,)+parts[j+1:]


def up_apply(source, target, H, n, columns):
    out = [[0]*columns for _ in target]; ix = {p: i for i, p in enumerate(target)}
    for i, p in enumerate(source):
        for q in covers(p, n):
            for j, c in enumerate(H[i]): out[ix[q]][j] += c
    return out


def coords(H, A):
    """Solve a lower triangular integral lattice basis, without QQ solving."""
    if not H:
        assert not A
        return []
    m = len(H); width = len(A[0]); result = []
    for i in range(m):
        assert H[i][i] > 0 and not any(H[i][i+1:])
        row = A[i][:]
        for k in range(i):
            if H[i][k]:
                row = [a-H[i][k]*b for a, b in zip(row, result[k])]
        assert all(c % H[i][i] == 0 for c in row)
        result.append([c//H[i][i] for c in row])
    return result


def decode(rows, width):
    return [[(int(row, 16) >> j) & 1 for j in range(width)] for row in rows]


def multiply(A, B, width):
    out = [[0]*width for _ in A]
    for i, row in enumerate(A):
        for k, c in enumerate(row):
            if c:
                for j, b in enumerate(B[k]): out[i][j] += c*b
    return out


def lowering_checks():
    n, i = sp.symbols('n i', integer=True)
    top = n*(n+1)/2
    b = n*(n+1)-i*(i+1)
    assert sp.expand(b.subs(i, 1)-2*top+2) == 0
    assert sp.expand(b-b.subs(i, i-1)+2*i) == 0
    assert sp.expand(b.subs(i, n-1)-2*n) == 0
    checked = 0
    for rank in range(1, 10):
        D = rank*(rank+1)//2
        def down(p):
            if p and p[0] == 1: yield p[1:], D
            for j, part in enumerate(p):
                if part > 1 and part-1 not in p:
                    yield p[:j]+(part-1,)+p[j+1:], (rank-part+1)*(rank+part)
        for level in levels(rank):
            for p in level:
                commutator = {}
                for q in covers(p, rank):
                    for r, c in down(q): commutator[r] = commutator.get(r, 0)+c
                for q, c in down(p):
                    for r in covers(q, rank): commutator[r] = commutator.get(r, 0)-c
                commutator = {r: c for r, c in commutator.items() if c}
                expected = D-2*sum(p)
                assert commutator == ({p: expected} if expected else {})
                for q, c in down(p):
                    assert (2**len(p)*c) % (2*2**len(q)) == 0
                checked += 1
    # Infinite Schur-P derivation: delta(p1)=1, delta(pi)=2*p_(i-1).
    for j in range(1, 10):
        p = (sp.Integer(1),)+sp.symbols('p1:'+str(2*j+1))
        relation = p[j]**2+2*sum((-1)**k*p[j+k]*p[j-k] for k in range(1, j))+(-1)**j*p[2*j]
        derivative = sp.diff(relation, p[1])+sum(2*p[k-1]*sp.diff(relation, p[k]) for k in range(2, 2*j+1))
        assert sp.expand(derivative) == 0
    return dict(symbolic_commutator_coefficients=True, strict_partitions_checked=checked,
                infinite_quadratic_derivation_identities_checked=9,
                all_rank_argument='The symbolic coefficients and the written operator proof establish all n; finite commutators are corroboration.')


def main():
    primary_path = ROOT/'experiments/spin-parity-recursion.json'
    primary = json.loads(primary_path.read_text())
    archive = ROOT/primary['certificate_archive']
    assert digest(archive) == primary['certificate_archive_sha256']
    raw = json.loads(gzip.decompress(archive.read_bytes()))
    cases = []
    for source, expected in zip(raw['cases'], primary['cases']):
        n = source['rank']; lev = levels(n); top = len(lev)-1
        hs = [[[1]]]; ds = [[]]; first = None
        for record in source['degrees'][1:]:
            d = record['degree']; target = lev[d] if d <= top else []
            m = len(target); prev = len(lev[d-1])
            ZH = up_apply(lev[d-1], target, hs[-1], n, prev)
            A = [[(2**len(p) if i == j else 0) for j in range(m)]+ZH[i]
                 for i, p in enumerate(target)]
            H = record['lattice_rows']; T = coords(H, A)
            if m:
                # M is a diagonal sublattice of the generated image, so its
                # determinant is a proved multiple of the unknown Hermite
                # determinant. This bound enables native modular Hermite
                # without assuming the supplied basis is correct.
                bound = 2**sum(len(p) for p in target)
                native = hermite_normal_form(sp.Matrix(A), D=bound, check_rank=True)
                assert abs(int(native.det())) == abs(sp.prod(H[i][i] for i in range(m)))
            # A is contained in H by exact coordinates; equal native Hermite
            # determinant proves equality of the two full-rank lattices.
            jet = [r[:] for r in hs[-1]]
            if d >= 2:
                product = multiply(hs[-2], ds[-1], prev)
                additional = up_apply(lev[d-2], lev[d-1], product, n, prev)
                jet = [[a+b for a, b in zip(r, s)] for r, s in zip(jet, additional)]
            B = [[0]*m+r for r in jet]
            C = coords(hs[-1], B)
            if record['derivative_rows'] is None:
                x = int(record['relation_bits'], 16)
                assert all(sum(c*((x >> j) & 1) for j, c in enumerate(row)) % 2 == 0 for row in T)
                remainder = sum((sum(c*((x >> j) & 1) for j, c in enumerate(row)) % 2) << i
                                for i, row in enumerate(C))
                assert remainder == int(record['derivative_remainder_bits'], 16) != 0
                first = d
                break
            delta = decode(record['derivative_rows'], m)
            reconstructed = multiply(delta, T, m+prev)
            assert all((a-b) % 2 == 0 for r, s in zip(reconstructed, C) for a, b in zip(r, s))
            hs.append(H); ds.append(delta)
        assert first == expected['first_failure_degree']
        assert first is not None or source['degrees'][-1]['degree'] == top+1
        cases.append(dict(rank=n, entire_kernel_preserved=first is None, first_failure_degree=first,
                          degrees_checked=source['degrees'][-1]['degree'], native_lattices_verified=True))
        print(n, first or 'stable in every degree', flush=True)
    report = dict(result='passed', date='2026-10-22', cases=cases,
                  lowering_proof=lowering_checks(),
                  method='Native SymPy Hermite determinants, exact triangular coordinates and bit relation identities',
                  shared_inputs='Classical Schur-Q lattice and one-box Pieri rule; supplied lattice/derivative certificates independently verified, not trusted. No Sage or polynomial saturation.',
                  arbitrary_rank_closed_classification_complete=False,
                  source_sha256={'experiments/audit_spin_parity_recursion.py': digest(Path(__file__)),
                                 'experiments/spin-parity-recursion.json': digest(primary_path)})
    (ROOT/'experiments/spin-parity-recursion-audit.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
