"""Finite original-coordinate constraint classifier; no generic inversion.

Sparse Laurent determinants classify every choice of equally sized,
disjoint multiplier/paired coordinate sets. Coordinates are not changed
before this test. A successful determinant is exactly one Laurent term.
"""
from itertools import combinations, permutations


def add(P,Q,scale=1):
    out=dict(P)
    for key,value in Q.items():
        out[key]=out.get(key,0)+scale*value
        if not out[key]:del out[key]
    return out


def multiply(P,Q):
    out={}
    for a,c in P.items():
        for b,d in Q.items():
            key=tuple(x+y for x,y in zip(a,b))
            out[key]=out.get(key,0)+c*d
    return {k:v for k,v in out.items() if v}


def determinant(M,central_rank):
    out={};n=len(M)
    for perm in permutations(range(n)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term={(0,)*central_rank:sign}
        for i,j in enumerate(perm):term=multiply(term,M[i][j])
        out=add(out,term)
    return out


def classify(row):
    n=len(row['nodes']);s=len(row['central_cocharacters']);terms=row['terms']
    linear=[i for i in range(n) if all(t['exponent'][i]<=1 for t in terms)]
    # D_ij(phi) lies in the Laurent base ring precisely when no term of
    # affine degree >=3 contains both i,j. Distinct charges are retained.
    coefficients={};blocked=set()
    for i,j in combinations(range(n),2):
        C={}
        for term in terms:
            e=term['exponent']
            if not e[i] or not e[j]:continue
            if sum(e)!=2:blocked.add((i,j));break
            C=add(C,{tuple(term['charge']):term['coefficient']})
        coefficients[i,j]=C
    results=[]
    for k in range(1,n//2+1):
        for tau in combinations(linear,k):
            if any(sum(t['exponent'][i] for i in tau)>1 for t in terms):continue
            for paired in combinations([j for j in range(n) if j not in tau],k):
                if any(tuple(sorted((i,j))) in blocked for i in tau for j in paired):continue
                M=[[coefficients[tuple(sorted((i,j)))] for j in paired] for i in tau]
                det=determinant(M,s)
                if len(det)!=1:continue
                charge,coeff=next(iter(det.items()))
                results.append(dict(multipliers=list(tau),paired=list(paired),
                                    determinant_charge=list(charge),determinant_coefficient=coeff))
    maximum=max((len(r['paired']) for r in results),default=0)
    best=[r for r in results if len(r['paired'])==maximum]
    return dict(accepted_systems=len(results),maximum_pairs=maximum,
                maximum_systems=best,
                chosen=min(best,key=lambda r:(r['multipliers'],r['paired'])) if best else None)
