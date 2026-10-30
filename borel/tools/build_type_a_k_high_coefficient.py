"""Construct Hasse-moment source witnesses by the uniform deficit recurrence.

No lattice solver is used. Denominators are supported at primes <= max(2,e).
The finite archive permits an independent polynomial replay in Sage.
"""
from functools import lru_cache
from hashlib import sha256
from math import comb, prod
from pathlib import Path
import gzip
import json
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(None)
def exponents(weights,n):
    if not weights:
        return ((),) if n==0 else ()
    return tuple((a,)+tail for a in range(n//weights[0]+1)
                 for tail in exponents(weights[1:],n-a*weights[0]))


def build(e,k):
    assert (e==1 and k>=3) or (e>=2 and k>=e+1)
    d=k+e
    M,L=s.symbols('M L')
    aa=s.symbols('a1:'+str(d+1))
    Q=s.symbols('Q1:'+str(e+3))
    variables=(M,)+aa+Q
    ML=Q[0]-aa[0]

    @lru_cache(None)
    def lift(u,lam,t):
        h=t-u
        coefficient=prod(aa[i-1] for i in lam)
        if h<=0:
            return coefficient*M**(u-t)*ML**t
        assert h<=e
        if u>=1 and 1 in lam and len(lam)>=2:
            smaller=list(lam)
            smaller.remove(1)
            smaller=tuple(smaller)
            leading=prod(aa[i-1] for i in smaller)*ML**(u-1)*(M*Q[h+1]-ML*Q[h])
            return leading-sum(comb(h+1,i-1)*lift(u,tuple(sorted(smaller+(i,))),t+1-i)
                               for i in range(2,h+3))
        if u>=2:
            # This also handles the pure ordinary monomial. Its first
            # positive-coefficient successor has u-1>=2 when lam is empty.
            return coefficient*ML**(u-1)*Q[h]-sum(
                comb(h+1,i)*lift(u-1,tuple(sorted(lam+(i,))),t-i) for i in range(1,h+2))
        assert u==0 and 1 in lam and len(lam)>=2
        if e>=2 and h==e:
            # The sole top-deficit target has all coefficient indices one.
            # q_2*q_e supplies coefficient 2e instead of the unavailable e+1.
            assert lam==(1,)*k and t==e
            leading=aa[0]**(k-2)*Q[1]*Q[e-1]
            correction=0
            for i in range(3):
                for j in range(e+1):
                    if i==j==1:continue
                    new_u=int(i==0)+int(j==0)
                    new_lam=tuple(sorted((1,)*(k-2)+tuple(x for x in (i,j) if x)))
                    correction+=comb(2,i)*comb(e,j)*lift(new_u,new_lam,e+2-i-j)
            return (leading-correction)/(2*e)
        smaller=list(lam)
        smaller.remove(1)
        smaller=tuple(smaller)
        leading=prod(aa[i-1] for i in smaller)*Q[h]
        return (leading-lift(1,smaller,t+1)-sum(
            comb(h+1,i)*lift(0,tuple(sorted(smaller+(i,))),t+1-i) for i in range(2,h+2)))/(h+1)

    substitution={Q[j-1]:M*L**j+sum(comb(j,i)*aa[i-1]*L**(j-i) for i in range(1,j+1))
                  for j in range(1,e+3)}
    records=[]
    for exp in exponents(tuple([1]+list(range(1,d+1))),d):
        t=exp[0]
        lam=tuple(i for i,n in enumerate(exp[1:],1) for _ in range(n))
        u=k-len(lam)
        if u<0 or t<=u:
            continue
        polynomial=s.Poly(s.expand(lift(u,lam,t)),*variables)
        target=M**u*prod(aa[i-1] for i in lam)*L**t
        assert s.expand(polynomial.as_expr().subs(substitution)-target)==0
        terms=[]
        for powers,c in polynomial.terms():
            assert sum(powers)==k
            assert sum(i*powers[i] for i in range(1,d+1))+sum(
                j*powers[d+j] for j in range(1,e+3))==d
            assert not lam or any(powers[1:d+1])
            assert all(int(p)<=max(2,e) for p in s.factorint(c.q))
            terms.append([int(c.p),int(c.q),list(powers)])
        records.append(dict(mass_power=u,coefficient_indices=list(lam),L_power=t,
                            coefficient_ideal_witness=bool(lam),terms=terms))
    return dict(excess=e,coefficient_degree=k,weight=d,variable_names=[str(x) for x in variables],
                witnesses=records)


def main():
    cases=[]
    for e,k in [(1,3),(2,4),(3,5),(4,6),(3,7),(5,7),(2,3),(3,4),(4,5),(5,6)]:
        cases.append(build(e,k))
        print(json.dumps(dict(excess=e,coefficient_degree=k,witnesses=len(cases[-1]['witnesses']))),flush=True)
    archive=ROOT/'experiments/type-a-k-high-coefficient-witnesses.json.gz'
    archive.write_bytes(gzip.compress(json.dumps(dict(cases=cases),separators=(',',':')).encode(),mtime=0))
    names=['TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md','tools/build_type_a_k_high_coefficient.py']
    report=dict(result='passed',method='Recursive symbolic source witnesses; no matrix membership solver',
                cases=[dict(excess=c['excess'],coefficient_degree=c['coefficient_degree'],weight=c['weight'],
                            witness_count=len(c['witnesses'])) for c in cases],
                total_witnesses=sum(len(c['witnesses']) for c in cases),
                all_positive_coefficient_targets_have_literal_ideal_witnesses=True,
                denominator_prime_bound='max(2,e); no division by M or representation rank',
                boundary_cases_k_equal_e_plus_one=4,
                all_rank_claim_from_finite_checks=False,
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in names},
                raw_sha256={str(archive.relative_to(ROOT)).replace('\\','/'):sha256(archive.read_bytes()).hexdigest()})
    (ROOT/'type-a-k-high-coefficient-build.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
