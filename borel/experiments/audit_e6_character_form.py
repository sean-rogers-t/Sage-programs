"""Pure-Python audit of the E6 integral-form certificate and low-degree lattices.

Uses fractions, integer arithmetic and dictionaries only; no Sage or library
imports. Checks the all-degree leading-monomial argument and independently
reconstructs Laurent coefficients through degree four.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
from time import perf_counter
from math import gcd
from itertools import combinations

ROOT=Path(__file__).resolve().parents[1]
ZERO=(0,)*6
C=[[2,0,-1,0,0,0],[0,2,0,-1,0,0],[-1,0,2,-1,0,0],
   [0,-1,-1,2,-1,0],[0,0,0,-1,2,-1],[0,0,0,0,-1,2]]
Y=[[3,0,0,0,0,0],[-1,-2,2,0,0,0],[-1,2,2,-2,0,0],
   [-1,0,0,2,-2,0],[-1,0,0,0,2,-2],[-1,0,0,0,0,2]]


def determinant(A):
    A=[list(map(F,row)) for row in A]
    answer=F(1)
    for i in range(len(A)):
        j=next((j for j in range(i,len(A)) if A[j][i]),None)
        if j is None:return F(0)
        if j!=i:A[i],A[j]=A[j],A[i];answer=-answer
        pivot=A[i][i];answer*=pivot
        for j in range(i+1,len(A)):
            q=A[j][i]/pivot
            for k in range(i+1,len(A)):A[j][k]-=q*A[i][k]
    return answer


def inverse(A):
    n=len(A)
    B=[[F(c) for c in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for i in range(n):
        j=next(j for j in range(i,n) if B[j][i]);B[i],B[j]=B[j],B[i]
        q=B[i][i];B[i]=[v/q for v in B[i]]
        for j in range(n):
            if j!=i:
                q=B[j][i];B[j]=[a-q*b for a,b in zip(B[j],B[i])]
    return [r[n:] for r in B]


def combine(polys,coeffs):
    out={}
    for p,c in zip(polys,coeffs):
        if c:
            for e,a in p.items():out[e]=out.get(e,0)+c*a
    return {e:c for e,c in out.items() if c}


def multiply(a,b):
    out={}
    for u,c in a.items():
        for v,d in b.items():
            w=tuple(x+y for x,y in zip(u,v));out[w]=out.get(w,0)+c*d
    return {e:c for e,c in out.items() if c}


def modular_rank(rows,p):
    pivots={}
    width=len(rows[0])
    for row in rows:
        row=[int(x)%p for x in row]
        for k in sorted(pivots):
            if row[k]:
                q=row[k];row=[(x-q*y)%p for x,y in zip(row,pivots[k])]
        k=next((i for i,v in enumerate(row) if v),None)
        if k is not None:
            inv=pow(row[k],-1,p);pivots[k]=[(x*inv)%p for x in row]
        if len(pivots)==width:return width
    return len(pivots)


def matmul(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def smith_square(H):
    n=len(H);previous=1;factors=[]
    for k in range(1,n+1):
        divisor=0
        for rows in combinations(range(n),k):
            for cols in combinations(range(n),k):
                divisor=gcd(divisor,int(determinant([[H[i][j] for j in cols] for i in rows])))
        factors.append(divisor//previous);previous=divisor
    return factors


def check_image_basis(M,H):
    assert all(F(c).denominator==1 for row in H for c in row)
    coordinates=matmul(inverse(H),M)
    assert all(F(c).denominator==1 for row in coordinates for c in row)
    divisor=0
    for cols in combinations(range(len(M[0])),len(M)):
        divisor=gcd(divisor,int(determinant([[row[j] for j in cols] for row in coordinates])))
        if divisor==1:break
    assert divisor==1


def main():
    start=perf_counter()
    path=ROOT/'experiments/e6-character-form-results.json'
    data=json.loads(path.read_text())
    assert data['result']=='passed'
    weights={(1,0,0,0,0,0)};pending=list(weights)
    while pending:
        w=pending.pop()
        for i in range(6):
            v=tuple(w[j]-w[i]*C[i][j] for j in range(6))
            if v not in weights:weights.add(v);pending.append(v)
    assert len(weights)==27
    gram=[[sum(w[i]*w[j] for w in weights) for j in range(6)] for i in range(6)]
    assert gram==[[6*c for c in row] for row in C]
    moments=[{w:w[i] for w in weights if w[i]} for i in range(6)]
    theta=[combine(moments,row) for row in inverse(gram)]
    certificate=data['torus_certificate'];U=certificate['unimodular_change']
    assert determinant(U)==-1
    transformed=[combine(moments,row) for row in U]
    ordered=sorted(weights,reverse=True)
    assert [list(w) for w in ordered]==certificate['ordered_weights']
    assert [[p.get(w,0) for w in ordered] for p in transformed]==certificate['coefficient_rows']
    for cert in certificate['leading']:
        prime=cert['characteristic']
        leading=[max(w for w,c in p.items() if c and (not prime or c%prime)) for p in transformed]
        assert [list(w) for w in leading]==cert['exponents']
        assert [p[w] for p,w in zip(transformed,leading)]==cert['coefficients']
        assert determinant(leading)==cert['exponent_determinant']!=0
    # Weyl equivariance of the moment map, checked coefficientwise.
    checks=0
    for i in range(6):
        for j in range(6):
            reflected={tuple(w[k]-w[i]*C[i][k] for k in range(6)):c for w,c in moments[j].items()}
            expected=combine([moments[j],moments[i]],[1,-C[j][i]])
            assert reflected==expected;checks+=1
    ys=[combine(theta,row) for row in Y]
    squares=[multiply(p,p) for p in ys[1:]]
    inputs=[ys[0],combine(squares,[1]*5),
            combine([multiply(squares[i],squares[j]) for i in range(5) for j in range(i+1,5)],[1]*10)]
    @lru_cache(None)
    def monomial(p):
        if not any(p):return {ZERO:F(1)}
        i=next(i for i,e in enumerate(p) if e);q=list(p);q[i]-=1
        return multiply(inputs[i],monomial(tuple(q)))
    last=next(r for r in data['direct_laurent_rows'] if r['degree']==4 and r['mode']=='filtered')
    allpowers=[tuple(p) for p in last['monomial_exponents']]
    support=sorted(set().union(*(set(monomial(p)) for p in allpowers)))
    coefficient_checks=0
    for row in data['direct_laurent_rows']:
        polys=[monomial(tuple(p)) for p in row['monomial_exponents']]
        B=[[F(x) for x in r] for r in row['source_basis_matrix']]
        assert determinant(B)
        integral=[combine(polys,[B[i][j] for i in range(len(B))]) for j in range(len(B))]
        assert all(c.denominator==1 for p in integral for c in p.values())
        N=[[int(p.get(w,0)) for p in integral] for w in support]
        cert=row['completeness_certificate']
        divisor=0
        for minor in cert['maximal_minors']:
            det=abs(determinant([N[i] for i in minor['rows']]))
            assert det==int(minor['determinant'])
            divisor=gcd(divisor,int(det))
            assert modular_rank(N,minor['selection_prime'])==len(B)
        assert divisor==cert['minor_gcd']==1
        coefficient_checks+=len(N)*len(B)
    oldpath=ROOT/'experiments/exceptional-minuscule-results.json'
    old=json.loads(oldpath.read_text())
    image_checks=0
    for row in data['degrees']:
        H=row['image_basis']
        check_image_basis(row['image_matrix'],H)
        assert smith_square(H)==row['smith_factors']
        image_checks+=1
        for comparison in row['comparisons']:
            containment=comparison['containment_matrix']
            assert smith_square(containment)==comparison['quotient_factors']
            oldH=matmul(H,inverse(containment))
            case=next(c for c in old['cases'] if c['kind']=='E6' and c['form']==comparison['form'])
            prior=case['degrees'][row['degree']]
            M=[prior['matrix'][prior['schubert_words'].index(word)] for word in row['schubert_words']]
            check_image_basis(M,oldH)
            assert smith_square(oldH)==comparison['old_smith_factors']
            assert matmul(oldH,containment)==H
            image_checks+=1
    result=dict(result='passed',engine='Python standard library only',weights=27,
      reflection_coefficient_checks=checks,leading_characteristic_cases=2,
      complete_laurent_lattices_checked=len(data['direct_laurent_rows']),
      laurent_coefficient_positions_checked=coefficient_checks,
      image_lattices_and_smith_lists_checked=image_checks,embedded_containments_checked=2*len(data['degrees']),
      seconds=round(perf_counter()-start,3),
      source_hash_semantics='SHA256 with CRLF normalized to LF.',
      source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                     for p in (path,oldpath,Path(__file__).resolve())})
    (ROOT/'experiments/e6-character-form-audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
