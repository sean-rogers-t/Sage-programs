"""Separate exact-arithmetic audit; no Sage or runtime-library imports.

Reuses the earlier independent auditor's root/reflection and Freudenthal
utilities. Adds truncated Kostant character coefficients to check the path
weights by a different multiplicity formula, even for huge highest weights.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import lcm, prod
from pathlib import Path
from time import perf_counter
import json

from audit_trace_source_survey import (root_data,weights_from_strings,determinant,
                                      inverse,matmul)

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(None)
def partitions(kind,delta):
    """Kostant partition function on a small rectangular truncation."""
    positive = [a for a,w,norm in root_data(kind)[2] if all(x<=y for x,y in zip(a,delta))]
    zero = (0,)*len(delta)
    coefficients = {zero:1}
    for alpha in positive:
        new = dict(coefficients)
        for beta,c in coefficients.items():
            gamma = tuple(x+y for x,y in zip(beta,alpha))
            while all(x<=y for x,y in zip(gamma,delta)):
                new[gamma] = new.get(gamma,0)+c
                gamma = tuple(x+y for x,y in zip(gamma,alpha))
        coefficients = new
    return coefficients.get(delta,0)


@lru_cache(None)
def kostant_near_highest(kind,highest,delta):
    """Weyl numerator in the box 0<=lambda+rho-w(lambda+rho)<=delta."""
    C = root_data(kind)[0]
    n = len(highest)
    zero = (0,)*n
    signs, pending = {zero:1}, deque([zero])
    while pending:
        a = pending.popleft()
        for i in range(n):
            pairing = highest[i]+1-sum(C[i][j]*a[j] for j in range(n))
            if pairing <= 0:
                continue
            b = list(a)
            b[i] += pairing
            b = tuple(b)
            if not all(x<=y for x,y in zip(b,delta)):
                continue
            if b in signs:
                assert signs[b] == -signs[a]
            else:
                signs[b] = -signs[a]
                pending.append(b)
    return sum(sign*partitions(kind,tuple(x-y for x,y in zip(delta,a)))
               for a,sign in signs.items())


def same_lattice(A,B):
    T = matmul(inverse(A),B)
    assert all(x.denominator==1 for r in T for x in r)
    assert abs(determinant(T)) == 1


def main():
    started = perf_counter()
    path = ROOT/'experiments/irreducible-lattice-survey.json'
    data = json.loads(path.read_text())
    assert data['result'] == 'passed'
    full_count = path_count = 0
    summary = {}
    for row in data['cases']:
        kind,highest = row['cartan_type'],tuple(row['highest_weight'])
        C,d,positive,invC = root_data(kind)
        n = len(highest)
        lengths = [F(x,max(d)) for x in d]
        H = [[F(C[i][j],1)/lengths[j] for j in range(n)] for i in range(n)]
        dimension = prod(F(sum(lengths[i]*a[i]*(highest[i]+1) for i in range(n)),
                           sum(lengths[i]*a[i] for i in range(n))) for a,w,norm in positive)
        invH = inverse(H)
        casimir = sum(highest[i]*invH[i][j]*(highest[j]+2) for i in range(n) for j in range(n))
        index = dimension*casimir/(n+2*len(positive))
        central = lcm(*(x.denominator for r in matmul(invC,[[x] for x in highest]) for x in r))
        assert dimension == int(row['dimension']) and index == F(row['trace_index'])
        assert central == row['central_character_order']
        M = [[index*x for x in r] for r in H]
        assert M == [[F(x) for x in r] for r in row['trace_gram']]
        D,B = [[[F(x) for x in r] for r in row[key]] for key in ('moment_basis','linear_basis')]
        assert abs(determinant(D)) == abs(determinant(C))/central
        for v in list(zip(*C))+[highest]:
            assert all(x.denominator==1 for r in matmul(inverse(D),[[x] for x in v]) for x in r)
        same_lattice(B,matmul(M,list(map(list,zip(*inverse(D))))))
        paths = row['multiplicity_one_paths']
        assert len(paths) == n
        for target,path_nodes in enumerate(paths):
            nodes = [j-1 for j in path_nodes]
            assert nodes[-1] == target and highest[nodes[0]] > 0
            assert all(highest[j]==0 for j in nodes[1:]) and len(set(nodes)) == len(nodes)
            assert all(C[a][b] for a,b in zip(nodes,nodes[1:]))
            assert all(not C[a][b] for k,a in enumerate(nodes) for b in nodes[k+2:])
            delta = [0]*n
            for j in nodes:
                delta[j] += 1
                assert kostant_near_highest(kind,highest,tuple(delta)) == 1
                path_count += 1
        if kind in ('G2','F4','E6','E7','E8'):
            expected = 'Q_long' if kind in ('G2','F4') else ('P' if central==1 else 'Q')
            assert expected == row['lattice_type']
            base = H if expected=='Q_long' else (C if expected=='Q' else
                    [[int(i==j) for j in range(n)] for i in range(n)])
            same_lattice(B,[[index*x for x in r] for r in base])
            if kind=='E6':
                assert central == (1 if (highest[0]-highest[2]+highest[4]-highest[5])%3==0 else 3)
            if kind=='E7':
                assert central == (1 if (highest[1]+highest[4]+highest[6])%2==0 else 2)
        if row['full_weight_check']:
            weights = weights_from_strings(kind,highest)
            assert sum(weights.values()) == dimension and len(weights)==row['distinct_weights']
            for path_nodes in paths:
                mu = highest
                assert weights[mu] == 1
                for j in path_nodes:
                    mu = tuple(mu[i]-C[i][j-1] for i in range(n))
                    assert weights[mu] == 1
            actual = [[sum(m*w[i]*w[j] for w,m in weights.items()) for j in range(n)] for i in range(n)]
            assert actual == M
            # Integrality from every weighted coefficient, with the determinant
            # and multiplicity-one paths above certifying completeness.
            K = matmul(list(map(list,zip(*B))),inverse(M))
            for w,m in weights.items():
                assert all(x.denominator==1 for r in matmul(K,[[m*x] for x in w]) for x in r)
            same_lattice(B,[[F(x) for x in r] for r in row['full_linear_basis']])
            full_count += 1
        summary[kind] = summary.get(kind,0)+1
        if sum(summary.values())%100==0:
            print('Audited',sum(summary.values()),'elapsed',round(perf_counter()-started,2),flush=True)
    # Same index and central character do not classify reducible representations.
    minimal = weights_from_strings('F4',(0,0,0,1))
    adjoint = weights_from_strings('F4',(1,0,0,0))
    gram = lambda weights:[[sum(m*w[i]*w[j] for w,m in weights.items()) for j in range(4)] for i in range(4)]
    assert gram({w:3*m for w,m in minimal.items()}) == gram(adjoint)
    assert data['reducible_control']['trace_index'] == data['reducible_control']['adjoint_trace_index'] == 18
    assert data['reducible_control']['linear_scale'] == 6 and data['reducible_control']['adjoint_linear_scale'] == 18
    assert full_count == data['full_weight_cases'] and sum(summary.values()) == data['total_cases']
    files = [path,Path(__file__).resolve(),ROOT/'experiments/audit_trace_source_survey.py',
             ROOT/'IRREDUCIBLE_CHARACTER_LATTICES.md']
    result = dict(result='passed',engine='Python standard library only',cases=summary,
        total_cases=data['total_cases'],full_weight_reconstructions=full_count,
        kostant_path_prefix_checks=path_count,reducible_same_index_control='passed',
        seconds=round(perf_counter()-started,3),
        scope='All recorded root/index/central-class and lattice data; Kostant path coefficients for every case; separate root-string/Freudenthal reconstruction for the dimension-bounded cases. Large cases do not enumerate full weights. No Borel-image calculations.',
        source_hash_semantics='SHA256 with CRLF normalized to LF.',
        source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in files})
    (ROOT/'experiments/irreducible-lattice-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k!='source_sha256'},flush=True)


if __name__ == '__main__':
    main()
