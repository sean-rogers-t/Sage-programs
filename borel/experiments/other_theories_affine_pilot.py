"""Integral divisor powers on Gr_SL(n), and independent affine A1 Chevalley.

The published symmetric-function Hopf model gives the first calculation.
The prime-power formula is checked against a direct multinomial gcd DP.
An unrelated infinite-dihedral Chevalley calculation checks full affine
A1 divisor-image Smith factors ((d-1)!,d!). No infinite trace is used.
"""
from functools import lru_cache
from hashlib import sha256
import json
from math import comb,factorial,gcd
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def primes(n):
    sieve=[True]*(n+1)
    for p in range(2,n+1):
        if sieve[p]:
            yield p
            for j in range(p*p,n+1,p):sieve[j]=False


def vp_factorial(n,p):
    answer=0
    while n:n//=p;answer+=n
    return answer


def closed(n,d):
    result=1
    for p in primes(d):
        power=1
        while power*p<n:power*=p
        result*=p**vp_factorial(d//power,p)
    return result


def multinomial_gcds(n,through):
    gs=[1]
    for d in range(1,through+1):
        value=0
        for j in range(1,min(n-1,d)+1):value=gcd(value,comb(d,j)*gs[d-j])
        gs.append(value)
    return gs


def symplectic_closed(n,d):
    result=2**vp_factorial(d,2)
    for p in primes(d):
        if p==2:continue
        power=1
        while power*p<2*n:power*=p
        result*=p**vp_factorial(d//power,p)
    return result


def symplectic_gcds(n,through):
    gs=[1]
    for d in range(1,through+1):
        value=0
        for j in range(1,min(2*n-1,d)+1,2):
            value=gcd(value,comb(d,j)*2**(j-1)*gs[d-j])
        gs.append(value)
    return gs


def reduced(word):
    result=[]
    for i in word:
        if result and result[-1]==i:result.pop()
        else:result.append(i)
    return tuple(result)


def words(d):return [()] if d==0 else [tuple((start+j)%2 for j in range(d)) for start in (0,1)]


def reflection_coroot(word):
    assert len(word)%2==1 and word==word[::-1]
    middle=len(word)//2;root=[0,0];root[word[middle]]=1
    for i in reversed(word[:middle]):
        pairing=2*root[i]-2*root[1-i]
        root[i]-=pairing
    if max(root)<=0:root=[-c for c in root]
    assert min(root)>=0
    return root


def up(d,node):
    source=words(d);target=words(d+1);result=[[0]*len(source) for _ in target]
    for i,v in enumerate(target):
        for j,w in enumerate(source):
            reflection=reduced(w[::-1]+v)
            result[i][j]=reflection_coroot(reflection)[node]
    return result


def apply(A,v):return [sum(a*b for a,b in zip(row,v)) for row in A]


def affine_A1(through=100):
    # Every degree-d divisor monomial eps0^a eps1^(d-a), evaluated via
    # the affine Chevalley formula and the actual real coroot coefficients.
    columns=[[1]];reports=[]
    for d in range(1,through+1):
        U0=up(d-1,0);U1=up(d-1,1)
        columns=[apply(U1,v) for v in columns]+[apply(U0,columns[-1])]
        g=0;detg=0
        for v in columns:
            for x in v:g=gcd(g,x)
        for i,v in enumerate(columns):
            for w in columns[i+1:]:detg=gcd(detg,v[0]*w[1]-v[1]*w[0])
        actual=[g,detg//g];expected=[factorial(d-1),factorial(d)]
        assert actual==expected,(d,actual,expected)
        reports.append(dict(degree=d,smith_factors=[str(x) for x in actual]))
    return reports


def main():
    cases=[];checked=0
    for n in range(2,33):
        values=multinomial_gcds(n,200)
        assert values==[closed(n,d) for d in range(201)]
        checked+=len(values)
        first=next(d for d,g in enumerate(values) if g>1)
        threshold=min(p**a for p in primes(n+20) for a in range(1,10) if p**a>=n)
        assert first==threshold
        cases.append(dict(n=n,checked_through=200,first_nonsaturated_divisor_power=first,
                          initial_contents=[str(g) for g in values[:min(31,len(values))]]))
    a1=affine_A1();symplectic=[];symplectic_checked=0
    for n in range(1,25):
        values=symplectic_gcds(n,200)
        assert values==[symplectic_closed(n,d) for d in range(201)]
        assert all(c%closed(2*n,d)==0 for d,c in enumerate(values))
        symplectic_checked+=len(values)
        symplectic.append(dict(rank=n,degree_through=200,
                               initial_contents=[str(c) for c in values[:31]]))
    report=dict(result='passed',date='2026-10-22',type_A_content_cases=cases,
                multinomial_gcd_vs_prime_power_cases=checked,independent_A1_Chevalley_degrees=a1,
                type_C_content_cases=symplectic,type_C_weighted_gcd_vs_prime_power_cases=symplectic_checked,
                scope='Finite arithmetic checks corroborate the written arbitrary-n prime-power proof and evaluation/level construction; no affine Springer trace map is asserted.',
                source_sha256={'experiments/other_theories_affine_pilot.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-affine-pilot.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED',checked,'A and',symplectic_checked,'C content cases and',len(a1),'independent affine Chevalley degrees')


if __name__=='__main__':main()
