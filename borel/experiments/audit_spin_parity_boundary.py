"""Independent SymPy/Python audit of adjacent-rank boundary certificates.

Rebuilds generating columns and proves lattice equality using exact inclusion
and mod-two full row rank. Then independently runs derivative factorization
test on the projected lattices. Also checks Schur-P lowering by Pfaffians,
without using the generating-series proof or Pieri matrices for that check.
"""
import argparse
from functools import lru_cache
import gzip
import json
from pathlib import Path

import sympy as sp
from audit_spin_parity_recursion import digest, coords, up_apply

ROOT = Path(__file__).resolve().parents[1]


def rank_bits(rows):
    piv = {}
    for x in rows:
        while x:
            k = x.bit_length()-1
            if k in piv: x ^= piv[k]
            else:
                piv[k] = x
                break
    return len(piv)


def factor(T, C):
    piv = {}
    for i, row in enumerate(T):
        x = sum((c & 1) << j for j,c in enumerate(row)); coefficient = 1 << i
        while x:
            k = x.bit_length()-1
            if k in piv:
                x ^= piv[k][0]; coefficient ^= piv[k][1]
            else:
                piv[k] = (x, coefficient)
                break
        assert x, 'T must present a full image basis modulo two'
    remainders = []; answers = []
    for x in C:
        coefficient = 0; remainder = 0
        while x:
            k = x.bit_length()-1
            if k not in piv:
                remainder |= 1 << k
                x ^= 1 << k
            else:
                x ^= piv[k][0]; coefficient ^= piv[k][1]
        remainders.append(remainder); answers.append(coefficient)
    return rank_bits(remainders), answers


def pfaffian_checks():
    max_part = 8
    p = (sp.Integer(1),)+sp.symbols('p1:'+str(2*max_part+1))
    # Eliminate even special classes using the monic integral rho_j.
    sub = {}
    for j in range(1,max_part+1):
        value = (-1)**(j+1)*(p[j]**2+2*sum((-1)**k*p[j+k]*p[j-k] for k in range(1,j)))
        sub[p[2*j]] = sp.expand(value.subs(sub, simultaneous=True))
    def normal(f): return sp.expand(f.subs(sub, simultaneous=True))
    @lru_cache(None)
    def pair(a,b):
        if b == 0: return p[a]
        return p[a]*p[b]+2*sum((-1)**k*p[a+k]*p[b-k] for k in range(1,b))+(-1)**b*p[a+b]
    @lru_cache(None)
    def P(lam):
        if not lam: return sp.Integer(1)
        if len(set(lam)) != len(lam): return sp.Integer(0)
        if len(lam) & 1: lam = lam+(0,)
        # P_pfaffian uses P_ab in even size; entries with zero are P_a.
        return sp.expand(sum((-1)**(j+1)*pair(lam[0],lam[j])*P(lam[1:j]+lam[j+1:])
                             for j in range(1,len(lam))))
    from itertools import combinations
    count = 0
    for k in range(1,max_part+1):
        for increasing in combinations(range(1,max_part+1),k):
            if sum(increasing) > 15: continue
            lam = tuple(reversed(increasing)); f = P(lam)
            df = sp.diff(f,p[1])+sum(2*p[i-1]*sp.diff(f,p[i]) for i in range(2,2*max_part+1))
            expected = sp.Integer(0)
            for j,v in enumerate(lam):
                q = lam[:j]+((v-1,) if v > 1 else ())+lam[j+1:]
                expected += (2 if v > 1 else 1)*P(q)
            assert normal(df-expected) == 0, lam
            count += 1
    return count


def audit(report):
    archive = ROOT/report['certificate_archive']
    assert digest(archive) == report['certificate_archive_sha256']
    raw = json.loads(gzip.decompress(archive.read_bytes()))
    results = []
    for case, cert in zip(report['cases'], raw['cases']):
        n = case['rank']; r = n+1
        from itertools import combinations
        lev = [[] for _ in range(r*(r+1)//2+1)]
        for k in range(r+1):
            for p in combinations(range(1,r+1),k): lev[sum(p)].append(p)
        for ps in lev: ps.sort(key=lambda p:(r in p,p))
        previous = [[1]]; projected = [[1]]
        prev_z_coordinates = [[]]; derivative = []
        degree_rows = []; witness = None
        for rec, expected in zip(cert['degrees'][1:],case['degree_rows']):
            d = rec['degree']; ps = lev[d]; old = lev[d-1]; m = len(ps)
            a = sum(r not in p for p in ps); ap = sum(r not in p for p in old)
            UH = up_apply(old,ps,previous,r,len(old))
            A = [[(2**len(p) if i == j else 0) for j in range(m)]+UH[i] for i,p in enumerate(ps)]
            H = rec['lattice_rows']; inclusion = coords(H,A)
            # A contains diag(2**length), so ambient/A and H/A are 2-primary.
            # Full row rank modulo two of H^-1 A makes H/A have odd order;
            # therefore H/A=0. This proves equality without any HNF routine.
            assert rank_bits([sum((c & 1)<<j for j,c in enumerate(row))
                              for row in inclusion]) == m
            assert not any(H[i][j] for i in range(a) for j in range(a,m))
            # Rebuild the boundary map without reading a saved obstruction matrix.
            twiceB = [[0]*(m-a) for _ in range(ap)]
            ix = {p:i for i,p in enumerate(old[:ap])}
            for i,p in enumerate(ps[a:],a):
                if n not in p:
                    for j in range(a,m): twiceB[ix[p[:-1]+(n,)]][j-a] = 2*H[i][j]
            J = coords(projected,twiceB)
            boundary_rank = rank_bits([sum((c & 1)<<j for j,c in enumerate(row)) for row in J])
            # Independent original recurrence, using the projected rank-n image.
            F = [row[:a] for row in H[:a]]
            UF = up_apply(old[:ap],ps[:a],projected,n,ap)
            An = [[(2**len(p) if i == j else 0) for j in range(a)]+UF[i] for i,p in enumerate(ps[:a])]
            T = coords(F,An)
            C = []
            for i in range(ap):
                x = 1 << i
                for j,c in enumerate(prev_z_coordinates[i]):
                    if c & 1: x ^= derivative[j]
                C.append(x << a)
            derivative_rank, dd = factor(T,C)
            assert boundary_rank == derivative_rank == expected['obstruction_rank']
            degree_rows.append(dict(degree=d,obstruction_rank=boundary_rank,
                                    independent_derivative_defect=derivative_rank))
            previous = H; projected = F
            if derivative_rank:
                j = next(j for j in range(m-a) if any(row[j] & 1 for row in J))
                witness = dict(boundary_schubert_terms=[dict(partition=list(ps[i]),coefficient=H[i][a+j])
                                                       for i in range(a,m) if H[i][a+j]],
                               twice_lowered_boundary_image_coordinates_mod2=[i for i,row in enumerate(J) if row[j] & 1],
                               coordinate_basis='Top-left block of the preceding adjacent-rank lattice certificate')
                break
            derivative = dd
            prev_z_coordinates = [row[a:] for row in T]
        assert degree_rows[-1]['degree'] == case['tested_through_degree']
        prior = json.loads((ROOT/'experiments/spin-parity-recursion.json').read_text())
        matches = [p for p in prior['cases'] if p['rank']==n]
        if matches: assert matches[0]['first_failure_degree'] == case['first_failure_degree']
        results.append({k:case[k] for k in ('rank','first_failure_degree','tested_through_degree','entire_kernel_decision','entire_kernel_preserved')}
                       | dict(degree_rows=degree_rows,first_failure_witness=witness))
        print(n,'independently audited',flush=True)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stem',default='spin-parity-boundary')
    args = parser.parse_args()
    path = ROOT/'experiments'/(args.stem+'.json'); source = json.loads(path.read_text())
    results = audit(source)
    count = pfaffian_checks()
    report = dict(result='passed',date='2026-10-22',coefficient_ring='Z_(2)',parabolic='B',
                  independent_arithmetic='Python exact triangular inclusion and full row rank modulo two prove lattice equality, using the diagonal power-of-two sublattice; no HNF routine, FLINT or Sage. SymPy is used only for Pfaffian polynomial identities.',
                  independent_derivative_factorization=True,schur_P_pfaffian_lowering_checks=count,
                  closed_numerical_rank_classification_proved=False,cases=results,
                  source_sha256={'experiments/audit_spin_parity_boundary.py':digest(Path(__file__)),
                                 'experiments/audit_spin_parity_recursion.py':digest(ROOT/'experiments/audit_spin_parity_recursion.py')},
                  raw_sha256={path.relative_to(ROOT).as_posix():digest(path),source['certificate_archive']:digest(ROOT/source['certificate_archive'])})
    (ROOT/'experiments'/(args.stem+'-audit.json')).write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__': main()
