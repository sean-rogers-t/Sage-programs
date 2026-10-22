"""Check the same-rank even-source index formula against boundary ranks.

Uses separate one-step and two-step image recurrences in rank n. Does not
read the adjacent-rank lattice certificates or use the boundary map.
Requires python-flint; shares its arithmetic with the boundary calculation.
"""
from itertools import combinations
import json
from pathlib import Path

from flint import fmpz_mat
from spin_parity_boundary import digest

ROOT = Path(__file__).resolve().parents[1]


def valuation(x):
    x = abs(int(x)); assert x and x & (x-1) == 0
    return x.bit_length()-1


def lattice(M, Z):
    rows = M.nrows()
    A = fmpz_mat([[int(M[i,j]) for j in range(rows)]+[int(Z[i,j]) for j in range(Z.ncols())]
                 for i in range(rows)])
    G = A.transpose().hnf()
    return fmpz_mat([[int(G[j,i]) for j in range(rows)] for i in range(rows)])


def case(n, through=None):
    top = n*(n+1)//2; last = min(top+1,through or top+1)
    lev = [[] for _ in range(top+1)]
    for k in range(n+1):
        for p in combinations(range(1,n+1),k): lev[sum(p)].append(p)
    H = [fmpz_mat([[1]])]; C = [fmpz_mat([[1]])]; ups = []; es = [0]
    rows = [dict(degree=0,rank=1,index_exponent=0,obstruction_rank=0)]
    for d in range(1,last+1):
        ps = lev[d] if d <= top else []; prev = lev[d-1]
        ix = {p:i for i,p in enumerate(ps)}; m = len(ps)
        U = fmpz_mat(m,len(prev))
        for j,p in enumerate(prev):
            if (not p or p[0]>1) and (1,)+p in ix: U[ix[(1,)+p],j]=1
            for k,v in enumerate(p):
                q = p[:k]+(v+1,)+p[k+1:]
                if v<n and v+1 not in p and q in ix: U[ix[q],j]=1
        ups.append(U)
        if m:
            M = fmpz_mat(m,m)
            for i,p in enumerate(ps): M[i,i]=2**len(p)
            h = lattice(M,U*H[-1])
            c = lattice(M,U*ups[-2]*C[-2] if d>=2 else fmpz_mat(m,0))
            eh = sum(valuation(h[i,i]) for i in range(m))
            ec = sum(valuation(c[i,i]) for i in range(m))
            e = ec-eh
        else:
            h=c=fmpz_mat(0,0);e=0
        # N/C has exponent two; verify both containments exactly using FLINT.
        if m:
            X=h.solve(c);Y=c.solve(2*h)
            assert all(x.q==1 for x in X.entries()) and all(x.q==1 for x in Y.entries())
        obstruction = len(prev)-e-es[-1]
        assert obstruction>=0
        rows.append(dict(degree=d,rank=m,index_exponent=e,obstruction_rank=obstruction))
        es.append(e);H.append(h);C.append(c)
        print(n,d,obstruction,flush=True)
    first = next((x['degree'] for x in rows if x['obstruction_rank']),None)
    complete = last>=top+1
    if complete:
        assert sum(x['obstruction_rank'] for x in rows)==2**n-2*sum(es)
        assert (first is None)==(sum(es)==2**(n-1))
    return dict(rank=n,tested_through_degree=last,complete_index_polynomial=complete,
                first_failure_degree=first,total_index_exponent=sum(es) if complete else None,
                degree_rows=rows)


def main():
    expected=[]
    paths=['experiments/spin-parity-boundary.json','experiments/spin-parity-boundary-rank15.json']
    for path in paths: expected.extend(json.loads((ROOT/path).read_text())['cases'])
    results=[]
    for prior in expected:
        n=prior['rank'];result=case(n,37 if n==15 else None)
        assert result['first_failure_degree']==prior['first_failure_degree']
        for p in prior['degree_rows']:
            assert result['degree_rows'][p['degree']]['obstruction_rank']==p['obstruction_rank']
        results.append(result)
    report=dict(result='passed',date='2026-10-22',cases=results,
                coefficient_ring='Z_(2)',parabolic='B',
                method='Separate same-rank image recurrences N=M+p1*N and C=M+p1^2*C; exact image indices and formula a_(d-1)-e_d-e_(d-1)',
                independent_of_boundary_lattices=True,shared_FLINT_arithmetic=True,
                closed_numerical_rank_classification_proved=False,
                source_sha256={'experiments/check_spin_parity_even_source.py':digest(Path(__file__)),
                               'experiments/spin_parity_boundary.py':digest(ROOT/'experiments/spin_parity_boundary.py')},
                raw_sha256={p:digest(ROOT/p) for p in paths})
    (ROOT/'experiments/spin-parity-even-source-check.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
