"""Proof-boundary checks, separate from the library and its Sage survey.

Uses Python exact arithmetic and the previously documented independent root
utilities. Exhausts nonempty support masks, distinct-generator word orders,
and Kostant path coefficients for all five fixed exceptional diagrams.
"""
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd, prod
from pathlib import Path
from time import perf_counter
import json

from audit_trace_source_survey import (root_data,weights_from_strings,
                                      determinant,inverse,matmul)
from audit_irreducible_lattices import kostant_near_highest, same_lattice

ROOT = Path(__file__).resolve().parents[1]


def shortest_path(C,highest,target):
    pending = deque([(target,)])
    while pending:
        path = pending.popleft()
        if highest[path[-1]]:
            return tuple(reversed(path))
        pending.extend(path+(j,) for j in range(len(C))
                       if j not in path and C[path[-1]][j])
    raise AssertionError('disconnected support')


def nonzero_word(C,highest,word):
    # Each generator is used once. Its raising operator therefore commutes
    # with all previously used lowering operators and kills that vector.
    # The finite-dimensional sl2 criterion is exactly a positive h_i-weight.
    mu = list(highest)
    for j in word:
        if mu[j] <= 0:
            return False
        mu = [mu[i]-C[i][j] for i in range(len(C))]
    return True


def main():
    started = perf_counter()
    masks = paths = prefixes = words = 0
    profiles = {}
    for kind in ('G2','F4','E6','E7','E8'):
        C,d,positive,invC = root_data(kind)
        n = len(C)
        local = dict(support_masks=0,target_paths=0,path_prefixes=0,word_orders=0)
        for mask in range(1,1<<n):
            highest = tuple((mask>>i)&1 for i in range(n))
            local['support_masks'] += 1
            for target in range(n):
                path = shortest_path(C,highest,target)
                assert path[-1] == target and highest[path[0]] == 1
                assert all(highest[j]==0 for j in path[1:])
                local['target_paths'] += 1
                delta = [0]*n
                for length in range(1,len(path)+1):
                    prefix = path[:length]
                    delta[prefix[-1]] = 1
                    assert kostant_near_highest(kind,highest,tuple(delta)) == 1
                    valid = []
                    for word in permutations(prefix):
                        local['word_orders'] += 1
                        if nonzero_word(C,highest,word):
                            valid.append(word)
                    assert valid == [prefix]
                    local['path_prefixes'] += 1
        masks += local['support_masks']; paths += local['target_paths']
        prefixes += local['path_prefixes']; words += local['word_orders']
        profiles[kind] = local
        print(kind,local,flush=True)
    # A path through another supported node need not have multiplicity one.
    assert kostant_near_highest('A2',(1,1),(1,1)) == 2
    assert kostant_near_highest('A3',(1,0,1),(1,1,1)) == 3
    # Check ordinary coordinate duality versus normalized metric duality in A1.
    for a in range(1,33):
        weights = list(range(a,-a-1,-2))
        content = gcd(*weights)
        M = sum(w*w for w in weights)
        index = F(a*(a+1)*(a+2),6)
        assert index*2/content == F(M,content)
        assert content == gcd(a,2)
    # The extremal weights alone can miss simple roots: V(4*omega) in A1.
    assert gcd(4,-4) == 4 and gcd(4,2,0,-2,-4) == 2
    # Central-charge tests follow from integer matrices, not sampled weights.
    for kind,prime,charge in [('E6',3,[1,0,-1,0,1,-1]),
                              ('E7',2,[0,1,0,0,1,0,1])]:
        C = root_data(kind)[0]
        assert abs(determinant(C)) == prime and gcd(*charge) == 1
        assert all(sum(charge[i]*C[i][j] for i in range(len(C)))%prime == 0
                   for j in range(len(C)))
    # The embedded coroot lattice equals the LONG-root lattice with this scale.
    long_lattices = {}
    for kind,expected_index in [('G2',3),('F4',4)]:
        C,d,positive,invC = root_data(kind)
        n = len(C)
        lengths = [F(x,max(d)) for x in d]
        H = [[F(C[i][j])/lengths[j] for j in range(n)] for i in range(n)]
        long_roots = [w for a,w,norm in positive if F(norm,max(d))==2]
        invH = inverse(H)
        assert abs(determinant(H)) == expected_index
        for w in long_roots:
            assert all(x.denominator==1 for row in matmul(invH,[[x] for x in w]) for x in row)
        basis = next(list(map(list,zip(*selected))) for selected in combinations(long_roots,n)
                     if abs(determinant(list(zip(*selected))))==expected_index)
        same_lattice(H,basis)
        long_lattices[kind] = dict(index_in_P=expected_index,positive_long_roots=len(long_roots))
    # Independent construction of the F4 same-index direct-sum counterexample.
    minimal = weights_from_strings('F4',(0,0,0,1))
    adjoint = weights_from_strings('F4',(1,0,0,0))
    candidates = [w for w,m in minimal.items() if m==1]
    D = next(list(map(list,zip(*selected))) for selected in combinations(candidates,4)
             if abs(determinant(list(zip(*selected))))==1)
    # These actual multiplicity-one weights give D_min=P; all other generators
    # are integral. D_triple=3D_min exactly. The adjoint contains the simple roots.
    assert abs(determinant(root_data('F4')[0])) == 1
    moment = lambda W:[[sum(m*w[i]*w[j] for w,m in W.items()) for j in range(4)] for i in range(4)]
    Mmin = moment(minimal)
    Mtriple = moment({w:3*m for w,m in minimal.items()})
    Madj = moment(adjoint)
    assert Mtriple == Madj == [[3*x for x in row] for row in Mmin]
    Bmin = matmul(Mmin,list(map(list,zip(*inverse(D)))))
    Btriple = matmul(Mtriple,list(map(list,zip(*inverse([[3*x for x in row] for row in D])))))
    Badj = matmul(Madj,list(map(list,zip(*inverse(root_data('F4')[0])))))
    same_lattice(Bmin,Btriple)
    same_lattice(Badj,[[3*x for x in row] for row in Btriple])
    assert abs(determinant(Badj)) == 3**4*abs(determinant(Btriple))
    files = [Path(__file__).resolve(),ROOT/'experiments/audit_irreducible_lattices.py',
             ROOT/'experiments/audit_trace_source_survey.py',ROOT/'IRREDUCIBLE_LATTICE_AUDIT.md',
             ROOT/'IRREDUCIBLE_CHARACTER_LATTICES.md']
    report = dict(result='passed',engine='Python standard library only',
        exceptional_support_masks=masks,target_paths=paths,kostant_prefix_checks=prefixes,
        distinct_generator_orders_checked=words,profiles=profiles,rank_one_cases=32,
        nonshortest_path_controls=[dict(type='A2',multiplicity=2),dict(type='A3',multiplicity=3)],
        long_root_lattices=long_lattices,F4_reducible_source_index_ratio=3**4,
        scope='All nonempty support masks on the five fixed exceptional diagrams; word-order and Kostant checks. Adversarial path, extremal-only, metric-dual, central-charge, long-root-scale and direct-sum controls. Not an external or formal proof audit.',
        seconds=round(perf_counter()-started,3),source_hash_semantics='SHA256 with CRLF normalized to LF.',
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/irreducible-lattice-edge-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print({k:v for k,v in report.items() if k not in ('source_sha256','profiles')},flush=True)


if __name__ == '__main__':
    main()
