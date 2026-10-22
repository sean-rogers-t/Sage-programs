"""Native G2/B audit: root/coroot Chevalley multiplication and minor gcds.

Uses Python's standard library only. No Sage, runtime or primary experiment
imports, and no divided differences. The rank-two Cartan datum and the
previously proved full torus source are shared mathematical inputs.
"""
from collections import deque
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import gcd
from pathlib import Path
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
ONE = (1,0,0,1)
SIMPLE = ((-1,0,1,1), (1,3,0,-1))
ROOTS = ((2,-1), (-3,2))


def mul(a, b):
    return (a[0]*b[0]+a[1]*b[2], a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2], a[2]*b[1]+a[3]*b[3])


def inverse(a):
    det = a[0]*a[3]-a[1]*a[2]
    assert abs(det) == 1
    return (a[3]//det, -a[1]//det, -a[2]//det, a[0]//det)


def action(a, x):
    return (a[0]*x[0]+a[1]*x[1], a[2]*x[0]+a[3]*x[1])


def word_matrix(word):
    out = ONE
    for i in word:
        out = mul(out, SIMPLE[i-1])
    return out


def smith(rows):
    """All G2/B target degrees have rank <=2; use determinantal divisors."""
    first = 0
    for row in rows:
        for value in row:
            first = gcd(first, value)
    if not first:
        return []
    if len(rows) == 1:
        return [first]
    assert len(rows) == 2
    second = 0
    for i, j in combinations(range(len(rows[0])), 2):
        second = gcd(second, rows[0][i]*rows[1][j]-rows[0][j]*rows[1][i])
    assert second and second % first == 0
    return [first, second//first]


def main():
    start = perf_counter()
    source = ROOT/'experiments/g2-full-flag.json'
    data = json.loads(source.read_text())
    assert data['result'] == 'passed' and data['complete']
    lengths, queue = {ONE:0}, deque([ONE])
    while queue:
        w = queue.popleft()
        for s in SIMPLE:
            v = mul(w, s)
            if v not in lengths:
                lengths[v] = lengths[w]+1
                queue.append(v)
    assert len(lengths) == 12
    assert [sum(v == d for v in lengths.values()) for d in range(7)] == [1,2,2,2,2,2,1]
    roots = {}
    for w in lengths:
        inv = inverse(w)
        for i, alpha in enumerate(ROOTS):
            a = action(w, alpha)
            # Positive roots have nonnegative simple-root coordinates.
            if min(2*a[0]+3*a[1], a[0]+2*a[1]) < 0:
                continue
            reflection = mul(mul(w, SIMPLE[i]), inv)
            coroot = (inv[i*2], inv[i*2+1])  # (w^-1)^T e_i
            if a in roots:
                assert roots[a] == (reflection, coroot)
            roots[a] = (reflection, coroot)
    assert len(roots) == 6

    @lru_cache(None)
    def monomial(a, b):
        if a+b == 0:
            return {ONE:1}
        node = 0 if a else 1
        lower = monomial(a-1, b) if a else monomial(a, b-1)
        out = {}
        for w, coefficient in lower.items():
            for reflection, coroot in roots.values():
                v = mul(w, reflection)
                if lengths[v] == lengths[w]+1 and coroot[node]:
                    out[v] = out.get(v, 0)+coefficient*coroot[node]
        return out

    checked = 0
    assert [case['representation'] for case in data['cases']] == ['minimal','adjoint']
    for case in data['cases']:
        sx, sy = case['linear_generators']
        assert (sx,sy) == ((6,2) if case['representation'] == 'minimal' else (24,8))
        assert [row['degree'] for row in case['degrees']] == list(range(8))
        for row in case['degrees']:
            d = row['degree']
            assert row['cohomological_degree'] == 2*d
            assert row['source_exponents'] == [[a,d-a] for a in range(d+1)]
            normal_forms = {(i,j) for i in (0,1) for j in range(6) if i+j == d}
            assert set(map(tuple,row['quotient_basis_exponents'])) == normal_forms
            labels = [word_matrix(word) for word in row['schubert_words']]
            assert len(set(labels)) == len(labels)
            assert set(labels) == {w for w,d in lengths.items() if d == row['degree']}
            assert row['rank'] == row['ambient_rank'] == len(labels)
            columns = [{w: c*sx**a*sy**b for w,c in monomial(a,b).items()}
                       for a,b in row['source_exponents']]
            rebuilt = [[column.get(w,0) for column in columns] for w in labels]
            assert rebuilt == row['coefficient_matrix'], (case['representation'],row['degree'])
            assert smith(rebuilt) == row['smith_factors']
            checked += len(columns)
        # Check the two exact kernel relations directly in this different ring model.
        def scaled(a,b):
            return {w:c*sx**a*sy**b for w,c in monomial(a,b).items()}
        quadratic = {}
        for a,b,c in ((2,0,1),(1,1,-3),(0,2,3)):
            for w,v in scaled(a,b).items():
                quadratic[w] = quadratic.get(w,0)+c*v
        assert not any(quadratic.values()) and not scaled(0,6)
        # The twelve integer normal forms have no rational dependency.
        for row in case['degrees'][:7]:
            labels = [word_matrix(word) for word in row['schubert_words']]
            columns = [scaled(a,b) for a,b in row['quotient_basis_exponents']]
            normal = [[column.get(w,0) for column in columns] for w in labels]
            assert smith(normal) == row['smith_factors']
    assert checked == 72
    report = dict(result='passed', complete=True, sage_imported=False, runtime_imported=False,
                  algorithm='Independent rank-two Weyl matrices, positive-root/coroot Chevalley products and gcds of all rank-one/two minors',
                  group_elements=12, positive_roots=6, checked_monomial_columns=checked,
                  representations=['minimal','adjoint'], all_degrees=True,
                  integral_kernel_verified=True,
                  scope='All arithmetic reconstructed independently. Shared Cartan datum, characteristic-map theorem and previously proved full torus source; no independent reproof of that source theorem or novelty claim.',
                  seconds=round(perf_counter()-start,3),
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in (Path(__file__).resolve(),source)})
    (ROOT/'experiments/g2-full-flag-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED independent G2/B Chevalley audit:',checked,'columns')


if __name__ == '__main__':
    main()
