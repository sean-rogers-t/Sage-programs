"""Independent exact BGG checks for the corrected type A statements.

Uses the independently implemented BGG calculator from the full table audit,
never the dissertation's historical code. Run with Sage 9.5 or later.
"""
from pathlib import Path
import sys, json
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / 'papers' / 'table_audit'))
from recompute import Calculator
from sage.all import SymmetricFunctions, QQ, Partitions, prod
from itertools import combinations

schur = SymmetricFunctions(QQ).schur()


class TypeACalculator(Calculator):
    # Sage constructs a univariate ring for A1; support that boundary too.
    def expansion(self,f,parabolic):
        degree = f.total_degree() if hasattr(f,'total_degree') else f.degree()
        cache = {():f}
        out = {}
        for w in self.representatives(parabolic,degree):
            prefix = ()
            for i in reversed(w.reduced_word()):
                nxt = prefix+(int(i),)
                if nxt not in cache:
                    cache[nxt] = self.divided(cache[prefix],i)
                prefix = nxt
            if cache[prefix]:
                out[w] = QQ(cache[prefix])
        return out


def elementary(values, degree):
    return sum(prod(c) for c in combinations(values, degree))


def permutation_word(values):
    values = list(values)
    sorting = []
    while values != sorted(values):
        for a in range(len(values)-1):
            if values[a] > values[a+1]:
                values[a], values[a+1] = values[a+1], values[a]
                sorting.append(a+1)
    return list(reversed(sorting))


def grassmannian_word(mu, r, n):
    mu = list(mu) + [0] * (r-len(mu))
    first = [a + mu[r-a] for a in range(1,r+1)]
    tail = [a for a in range(1,n+1) if a not in first]
    return permutation_word(first+tail)


counts = {'elementary_rows': 0, 'schur_rows': 0, 'outside_rectangle_rows': 0}
for n in range(2,9):
    calc = TypeACalculator('A'+str(n-1))
    x = [calc.x[0]]
    x += [calc.x[j]-calc.x[j-1] for j in range(1,n-1)]
    x += [-calc.x[-1]]
    for r in range(1,n):
        d = n-r
        A,B = x[:r],x[r:]
        for i in range(1,r+1):
            actual = calc.expansion(calc.R(elementary(A,i)),r)
            expected = {calc.W.from_reduced_word(list(range(r-i+1,r+1))): QQ(1)}
            assert actual == expected, ('head',n,r,i,actual,expected)
            counts['elementary_rows'] += 1
        for i in range(1,d+1):
            actual = calc.expansion(calc.R(elementary(B,i)),r)
            expected = {calc.W.from_reduced_word(list(range(r+i-1,r-1,-1))): QQ((-1)**i)}
            assert actual == expected, ('tail',n,r,i,actual,expected)
            counts['elementary_rows'] += 1
        if n <= 6:
            for size in range(0,7):
                for lam in Partitions(size,max_length=d):
                    raw = schur(lam).expand(d)
                    f = calc.R(raw(*B)) if raw.parent().ngens()>1 else calc.R(raw(B[0]))
                    actual = calc.expansion(f,r)
                    if lam and lam[0]>r:
                        expected = {}
                        counts['outside_rectangle_rows'] += 1
                    else:
                        mu = lam.conjugate()
                        w = calc.W.from_reduced_word(grassmannian_word(mu,r,n))
                        expected = {w:QQ((-1)**size)}
                    assert actual == expected, ('schur',n,r,list(lam),actual,expected)
                    counts['schur_rows'] += 1
print(json.dumps({'status':'PASS', 'checks':counts},indent=2))
(ROOT/'type_a_checks.json').write_text(json.dumps({'status':'PASS','checks':counts},indent=2),encoding='utf-8')
