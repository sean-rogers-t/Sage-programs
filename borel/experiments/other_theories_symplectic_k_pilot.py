"""Independent projective K-image matrices and a recursive quotient model.

Uses only native exact integer/rational arithmetic and SymPy. No Sage,
BGG, trace-weight implementation or saved exceptional table is imported.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
try:
    from flint import fmpz_mat
except ImportError:
    fmpz_mat=None

ROOT=Path(__file__).resolve().parents[1]


def smith_factors(M):
    if fmpz_mat is not None:
        S=fmpz_mat([[int(x) for x in row] for row in M.tolist()]).snf()
        return [abs(int(S[i,i])) for i in range(M.rows)]
    S=smith_normal_form(M,domain=sp.ZZ)
    return [abs(int(S[i,i])) for i in range(M.rows)]


def choose(n,k):
    return (-1)**k*comb(k-n-1,k) if n<0 else (comb(n,k) if n>=k else 0)


def multiply(a,b,m):
    return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(m)]


def run(m,recursive=False):
    b=[0]+[2 if j==1 else (-1)**(j+1) for j in range(1,m)]
    columns=[];p=[1]+[0]*(m-1)
    for j in range(m):columns.append(p);p=multiply(p,b,m)
    M=sp.Matrix(columns).T
    Euler=sp.Matrix(m,m,lambda i,j:sum((-1)**k*comb(j,k)*choose(j-2*k+i,i) for k in range(j+1)))
    Pascal=sp.Matrix(m,m,lambda i,j:comb(i,j) if i>=j else 0)
    assert Euler==Pascal*M
    factors=smith_factors(M)
    assert all(x&(x-1)==0 for x in factors)
    exponents=[x.bit_length()-1 for x in factors]
    assert sum(exponents)==m*(m-1)//2
    assert exponents.count(0)==(m+1)//2
    j=(m-1)//2
    predicted=1 if j==0 else 4*j-j.bit_count()
    assert max(exponents)==predicted
    assert exponents.count(predicted)==(2 if m%2==0 and m>=4 else 1)
    if recursive:
        # Coker M = Z[b]/(c_m,b*c_(m-1)), c_0=0,c_1=1.
        x=sp.Symbol('b');previous=sp.Poly(0,x);current=sp.Poly(1,x)
        for _ in range(1,m):previous,current=current,sp.Poly((x-2)*current.as_expr()+x*previous.as_expr(),x)
        g=sp.Poly(x*previous.as_expr(),x);rank=m-1
        R=sp.Matrix(rank,rank,lambda i,j:sp.rem(g*sp.Poly(x**j,x),current).nth(i))
        assert [1]+smith_factors(R)==factors
    return dict(projective_dimension=m-1,symplectic_rank=(m//2 if m%2==0 else None),
                smith_2_exponents=exponents,Euler_matrix_agrees=True,
                recursive_quotient_agrees=recursive,exponent_formula_verified=True)


def main():
    cases=[]
    for m in range(2,65):
        cases.append(run(m,recursive=(m<=32)))
        if m%8==0:print('m',m,'checked',flush=True)
    # An independent rational recurrence for the inverse coordinate
    # v(b)=b/2+sqrt(1+b^2/4)-1 verifies its exact denominator formula.
    c=Fraction(1);denominator_checks=0
    for j in range(1,1001):
        c*=Fraction(3-2*j,2*j)
        coefficient=c/(4**j)
        assert coefficient.denominator==2**(4*j-j.bit_count())
        denominator_checks+=1
    report=dict(result='passed',date='2026-10-22',cases=cases,
                inverse_coefficient_denominator_checks=denominator_checks,
                scope='Complete finite K-image Smith factors checked by v-powers, independent projective Euler characteristics, and a two-polynomial quotient. The all-rank exponent and largest-factor multiplicity are proved separately; no closed full Smith formula is claimed.',
                source_sha256={'experiments/other_theories_symplectic_k_pilot.py':sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()})
    (ROOT/'experiments/other-theories-symplectic-k-pilot.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED',len(cases),'Smith cases, 31 recursive quotients and',denominator_checks,'inverse denominators')


if __name__=='__main__':main()
