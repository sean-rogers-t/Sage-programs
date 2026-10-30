"""Construct a single adapted basis, then test its actual Borel embedding.

This is research verification, not a new public API. Sage symmetric
functions construct rows; stored divided-difference columns are independent
of that construction. New small columns are also computed directly.
"""
import json
import sys
from collections import Counter
from math import comb
from pathlib import Path
from hashlib import sha256
from functools import lru_cache
from sage.all import ZZ, QQ, Partition, Partitions, SymmetricFunctions, matrix, diagonal_matrix

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments'))
from verify_type_b_general_smith import bead_statistic
from verify_type_b_filtration import schur_product_matrix

SF = SymmetricFunctions(ZZ)
S, H = SF.schur(), SF.homogeneous()

def corner(A,ell,t,g):
    b=ell-t
    n=A+ell
    r,c=n//2,(n+1)//2
    for j in range(max(0,b-r),min(b,c,(b+n%2)//2)+1):
        bounds=(c-j,)*j+(r-b+j,)*(b-j)
        if all(x<=y for x,y in zip(g[t:],bounds)):
            return j
    return None

@lru_cache(None)
def basis(A, ell, m):
    parts = tuple(sorted((tuple(p) for p in Partitions(m, max_part=A, max_length=ell)),
                         key=lambda p: (sum(x*x for x in p), p)))
    exps, rows, hrows = [], [], []
    for part in parts:
        g = part + (0,)*(ell-len(part))
        e = bead_statistic(A, ell, g)
        j = corner(A, ell, e, g)
        assert j is not None
        mu = tuple(x for x in g[:e] if x)
        alpha = tuple(x for x in g[e:e+j] if x)
        beta = tuple(x for x in g[e+j:] if x)
        value = S(H(mu))*S(alpha)*S(beta)
        rows.append([value[Partition(p)] for p in parts])
        hv = H(value)
        hrows.append([hv[Partition(p)] for p in parts])
        exps.append(e)
    V = matrix(ZZ, rows)
    U = matrix(ZZ, hrows)
    assert all(U[i,j]==0 for i in range(len(parts)) for j in range(i))
    assert list(U.diagonal()) == [1]*len(parts)
    assert abs(V.det()) == 1
    F = V.inverse().change_ring(ZZ)
    polynomials = tuple(sum((F[i,j]*S(p) for i,p in enumerate(parts)), S.zero())
                        for j in range(len(parts)))
    return parts, tuple(exps), V, F, polynomials

def decode(data):
    return matrix(ZZ, data['rows'], data['columns'],
                  {(e['row'],e['column']):ZZ(e['value']) for e in data['entries']})

def check_embedding(A, ell, m, B, source_parts):
    parts, exps, V, F, polys = basis(A,ell,m)
    # Stored columns are tail Schur classes. Inverse total Chern series
    # gives s_lambda(tail^2)=(-1)^m s_lambda'(head^2).
    indices = {tuple(Partition(p).conjugate()):i for i,p in enumerate(source_parts)}
    head = (-1)**m * B.matrix_from_columns([indices[p] for p in parts])
    primitive = head * F * diagonal_matrix(QQ,[QQ(1)/2**e for e in exps])
    primitive = primitive.change_ring(ZZ)
    assert all(abs(x)==1 for x in primitive.elementary_divisors() if x)
    assert primitive*diagonal_matrix(ZZ,[2**e for e in exps])*V == head
    return dict(A=A,ell=ell,m=m,rank=len(parts),integral=True,primitive=True,
                reconstructed=True)

def main():
    checks, products, negative = [], [], []
    # Compare every prefix, rather than only the final rank.
    for A,ell in [(3,2),(2,3),(4,2),(2,4),(3,3),(4,3),(4,4),(5,3)]:
        count = 0
        for m in range(A*ell+1):
            parts, es, V, F, fs = basis(A,ell,m)
            for t in range(ell+1):
                actual, _ = schur_product_matrix(A,ell,m,t,parts)
                selected=V.matrix_from_rows([i for i,e in enumerate(es) if e<=t])
                assert selected.row_module()==actual.row_module()
                count += 1
        # The source products, and scaled saturation products, are computed
        # in the rectangular Schur quotient, without a Borel matrix.
        terms = 0
        for u in range(A*ell+1):
            ps,es,V,F,fs = basis(A,ell,u)
            for v in range(u,A*ell-u+1):
                qs,ds,W,G,gs = basis(A,ell,v)
                rs,ts,X,J,hs = basis(A,ell,u+v)
                for i,f in enumerate(fs):
                    for j,g in enumerate(gs):
                        z=f*g
                        col=matrix(ZZ,len(rs),1,[z[Partition(p)] for p in rs])
                        coeff=X*col
                        for k,c in enumerate(coeff.column(0)):
                            value=QQ(c)*QQ(2)**(ts[k]-es[i]-ds[j])
                            assert value.denominator()==1,(A,ell,u,v,i,j,k,c)
                            terms+=bool(c)
        checks.append(dict(A=A,ell=ell,prefixes=count))
        products.append(dict(A=A,ell=ell,nonzero_structure_constants=terms,integral=True))
        print(json.dumps(products[-1]),flush=True)
    recorded=json.loads((ROOT/'experiments/type_b_multiplicity_search_results.json').read_text())
    embeddings=[]
    for case in recorded['cases']:
        for row in case['degrees']:
            if row['status']!='computed':continue
            embeddings.append(check_embedding(case['k'],case['a'],row['half_polynomial_degree'],
                                               decode(row['B_matrix']),row['partitions']))
    corner_checks=[]
    M=SF.monomial()
    for A,ell in [(2,2),(3,2),(4,2),(3,3),(4,3),(5,3),(4,4),(5,4),(5,5),(7,3)]:
        c=(A+ell+1)//2
        family=[()]
        for length in range(1,ell+1):
            family.extend(tuple(p) for m in range(length*c,length*A+1)
                          for p in Partitions(m,min_part=c,max_part=A,length=length))
        for nu in family:
            ps,es,V,F,fs=basis(A,ell,sum(nu))
            value=S(M(nu))
            column=matrix(ZZ,len(ps),1,[value[Partition(p)] for p in ps])
            expected=matrix(ZZ,len(ps),1,[int(p==nu) for p in ps])
            assert V*column==expected
            assert es[ps.index(nu)]==len(nu)
        pc=0
        for i,lam in enumerate(family):
            for mu in family[i:]:
                union=tuple(sorted(lam+mu,reverse=True))
                if sum(union)>A*ell:continue
                value=S(M(lam)*M(mu))
                coefficient=1
                for j in range(c,A+1):
                    coefficient*=comb(lam.count(j)+mu.count(j),lam.count(j))
                expected=S(M(union))*coefficient if len(union)<=ell else S.zero()
                for p in Partitions(sum(union),max_part=A,max_length=ell):
                    assert value[p]==expected[p]
                pc+=1
        corner_checks.append(dict(A=A,ell=ell,classes=len(family),products=pc))
    report=dict(status='passed',scope='Finite checks of a uniform constructive deduction; no external review.',
                prefixes=checks,multiplication=products,actual_embeddings=embeddings,
                divided_power_corners=corner_checks,
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'type-b-adapted-basis-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status='passed',embedding_cases=len(embeddings))),flush=True)

if __name__=='__main__':main()
