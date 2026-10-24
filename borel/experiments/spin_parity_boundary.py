"""Adjacent-rank boundary criterion; no derivative matrices are propagated.

Requires native python-flint (tested 0.9.0). Default exact cases corroborate
the prior derivative recurrence. --ranks and --through allow bounded runs;
an unfinished range is never reported as an entire-kernel decision.
"""
import argparse
from hashlib import sha256
import gzip
from itertools import combinations
import json
from pathlib import Path
from time import monotonic

from flint import fmpz_mat

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    raw = path.read_bytes()
    return sha256(raw if path.suffix == '.gz' else raw.replace(b'\r\n', b'\n')).hexdigest()


def levels(r):
    result = [[] for _ in range(r*(r+1)//2+1)]
    for k in range(r+1):
        for p in combinations(range(1, r+1), k): result[sum(p)].append(p)
    for row in result: row.sort(key=lambda p: (r in p, p))
    return result


def coordinates(H, B):
    """Exact lower triangular division; reject an unproved integrality step."""
    result = []
    for i in range(len(H)):
        assert H[i][i] > 0 and not any(H[i][i+1:])
        row = B[i][:]
        for k in range(i):
            if H[i][k]:
                row = [a-H[i][k]*b for a, b in zip(row, result[k])]
        assert all(a % H[i][i] == 0 for a in row)
        result.append([a//H[i][i] for a in row])
    return result


def bit_rank(rows):
    pivots = {}
    for row in rows:
        x = sum((int(c) & 1) << j for j, c in enumerate(row))
        while x:
            i = x.bit_length()-1
            if i in pivots: x ^= pivots[i]
            else:
                pivots[i] = x
                break
    return len(pivots)


def case(n, through):
    start = monotonic(); r = n+1; lev = levels(r)
    prev = [[1]]; records = [dict(degree=0, lattice_rows=prev)]
    rows = []; failure = None
    last = min(n*(n+1)//2+1, through)
    for d in range(1, last+1):
        ps = lev[d]; old = lev[d-1]; size = len(ps)
        a = sum(r not in p for p in ps); ap = sum(r not in p for p in old)
        ix = {p: i for i, p in enumerate(ps)}
        U = fmpz_mat(size, len(old))
        for j, p in enumerate(old):
            if not p or p[0] > 1: U[ix[(1,)+p], j] = 1
            for k, v in enumerate(p):
                if v < r and v+1 not in p:
                    U[ix[p[:k]+(v+1,)+p[k+1:]], j] = 1
        UH = U*fmpz_mat(prev)
        A = fmpz_mat([[(2**len(ps[i]) if i == j else 0) for j in range(size)]
                      +[int(UH[i,j]) for j in range(UH.ncols())] for i in range(size)])
        G = A.transpose().hnf()
        H = [[int(G[j,i]) for j in range(size)] for i in range(size)]
        assert not any(H[i][j] for i in range(a) for j in range(a, size))
        B = [[0]*(size-a) for _ in range(ap)]
        target_index = {p: i for i, p in enumerate(old[:ap])}
        for i, p in enumerate(ps[a:], a):
            if n not in p:
                q = p[:-1]+(n,)
                for j in range(a, size): B[target_index[q]][j-a] = H[i][j]
        F = [row[:ap] for row in prev[:ap]]
        J = coordinates(F, [[2*c for c in row] for row in B])
        defect = bit_rank(J)
        rows.append(dict(degree=d, adjacent_schubert_rank=size,
                         boundary_lattice_rank=size-a, obstruction_rank=defect))
        records.append(dict(degree=d, lattice_rows=H))
        print(n, d, 'obstruction', defect, 'seconds', round(monotonic()-start,2), flush=True)
        prev = H
        if defect:
            failure = d
            break
    return dict(rank=n, adjacent_rank=r, tested_through_degree=d,
                entire_kernel_decision=failure is not None or d >= n*(n+1)//2+1,
                entire_kernel_preserved=False if failure else (True if d >= n*(n+1)//2+1 else None),
                first_failure_degree=failure, degree_rows=rows), dict(rank=n, degrees=records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ranks', nargs='+', type=int, default=[3,4,7,8,11,12])
    parser.add_argument('--through', type=int, default=1000)
    parser.add_argument('--stem', default='spin-parity-boundary')
    args = parser.parse_args()
    cases = []; certs = []
    for n in args.ranks:
        assert n > 0 and args.through > 0
        result, cert = case(n, args.through)
        cases.append(result); certs.append(cert)
    archive = ROOT/'experiments'/(args.stem+'-certificates.json.gz')
    archive.write_bytes(gzip.compress(json.dumps(dict(cases=certs), separators=(',', ':')).encode(), mtime=0))
    report = dict(result='passed', date='2026-10-22', coefficient_ring='Z_(2)',
                  parabolic='B', source='Full Spin weight lattice, after free Artin symmetric-source reduction',
                  method='Adjacent-rank image Hermite blocks, boundary lowering, rank of twice-boundary coordinates modulo two',
                  closed_numerical_rank_classification_proved=False,
                  certificate_archive=archive.relative_to(ROOT).as_posix(),
                  certificate_archive_sha256=digest(archive),
                  source_sha256={'experiments/spin_parity_boundary.py': digest(Path(__file__))}, cases=cases)
    (ROOT/'experiments'/(args.stem+'.json')).write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__': main()
