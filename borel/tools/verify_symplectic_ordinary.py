"""Independent integer Smith, Euler-basis and 2-local Gram certificates.

No API formula supplies the expected Smith exponents: factorial products
are compared to both independently constructed matrices and the API's
binary-digit formula. The finite ranges corroborate the written proof.
"""
from hashlib import sha256
import json
from math import comb, factorial
from pathlib import Path
from time import perf_counter
from sage.all import ZZ, QQ, PolynomialRing, matrix, diagonal_matrix
from sage_borel.equivariant import SymplecticProjectiveKImage

ROOT=Path(__file__).resolve().parents[1]


def valuation(x):
    return int(QQ(x).valuation(2))


def predicted(n):
    r=n//2
    exponents=[(3+3*(n%2)+6*j+
                valuation(ZZ(factorial(r+j+n%2))//factorial(r-j-1)))
               for j in range(r)]
    return [0]*n+([1] if n%2 else [])+[e for e in exponents for _ in range(2)]


def half_data(n):
    S=PolynomialRing(ZZ,'u');u=S.gen();r=n//2
    E=S(sum(comb(n,2*j)*(-2)**(n-2*j)*u**j for j in range(r+1)))
    O=S(sum(comb(n,2*j+1)*(-2)**(n-2*j-1)*u**j for j in range((n+1)//2)))
    p,q=(E,O) if n%2==0 else (O,E)
    M=matrix(ZZ,r,r,lambda i,j:((2*q*u**j)%p)[i])
    return p,q,M


def original(n,euler=False):
    m=2*n
    # Direct truncated v-power multiplication, no even/odd quotient.
    b=[0]+[2 if j==1 else (-1)**(j+1) for j in range(1,m)]
    columns=[];power=[1]+[0]*(m-1)
    for _ in range(m):
        columns.append(power)
        power=[sum(power[i]*b[k-i] for i in range(k+1)) for k in range(m)]
    M=matrix(ZZ,columns).transpose()
    if euler:
        def choose(a,k):
            return (-1)**k*comb(k-a-1,k) if a<0 else (comb(a,k) if a>=k else 0)
        # chi(L^i*b^j), i=0,...,m-1, in projective dimension m-1.
        Euler=matrix(ZZ,m,m,lambda i,j:
            sum((-1)**k*comb(j,k)*choose(j-2*k+i,i) for k in range(j+1)))
        Pascal=matrix(ZZ,m,m,lambda i,j:comb(i,j) if i>=j else 0)
        assert Pascal.det()==1 and Euler==Pascal*M
    return M


def continued_fraction(n):
    S=PolynomialRing(QQ,'t');t=S.gen();F=S.fraction_field();value=F(1)
    for k in range(n-1,0,-1):
        c=QQ(n*n-k*k)/((2*k-1)*(2*k+1))
        value=1/(1+c*t*value)
    A=sum(comb(n,2*k)*t**k for k in range(n//2+1))
    B=sum(comb(n,2*k+1)*t**k for k in range((n+1)//2))
    assert value==F(B)/(n*A)


def gram(n,p,q):
    r=n//2;S=PolynomialRing(QQ,'u');u=S.gen()
    moments=[((q*p.parent().gen()**k)%p)[r-1] for k in range(2*r-1)]
    Hankel=matrix(QQ,r,r,lambda i,j:2*moments[i+j])
    def alpha(k):
        return -4*QQ(n*n-k*k)/((2*k-1)*(2*k+1)) if k<n else QQ(0)
    def beta(j):
        if n%2:
            return alpha(2*j+1)+alpha(2*j+2)
        return alpha(1) if j==0 else alpha(2*j)+alpha(2*j+1)
    def lam(j):
        return alpha(2*j-1+n%2)*alpha(2*j+n%2)
    previous,current=S(0),S(1);polynomials=[];diagonal=[];norm=2*moments[0]
    for j in range(r):
        polynomials.append(current);diagonal.append(norm)
        previous,current=current,(u-beta(j))*current-(lam(j)*previous if j else 0)
        norm*=lam(j+1)
    U=matrix(QQ,r,r,lambda i,j:polynomials[i][j])
    assert U.det()==1
    assert all(x.denominator()%2 for x in U.list())
    assert U*Hankel*U.transpose()==diagonal_matrix(QQ,diagonal)
    assert [valuation(x) for x in diagonal]==predicted(n)[n+n%2::2]
    return dict(rank=n,dimension=r,unimodular_over_Z_local_2=True,
                exact_Gram_identity=True,diagonal_2_exponents=[valuation(x) for x in diagonal])


def main():
    start=perf_counter();half_cases=[];original_cases=[];gram_cases=[]
    for n in range(1,129):
        expected=predicted(n)
        # Bypass the irrelevant spectral-polynomial constructor: the
        # public method depends only on the already-validated rank.
        api=object.__new__(SymplecticProjectiveKImage);api.rank=n
        actual=[valuation(x) for x in api.ordinary_smith_factors()]
        assert actual==expected
        assert sum(expected)==n*(2*n-1) and expected.count(0)==n
        if n>=2:
            assert expected[-1]==4*(n-1)-int(n-1).bit_count()
            p,q,M=half_data(n)
            half=[valuation(x) for x in M.elementary_divisors()]
            assert [0]*n+([1] if n%2 else [])+[e for e in half for _ in range(2)]==expected
            if n<=32:gram_cases.append(gram(n,p,q))
        half_cases.append(dict(rank=n,full_smith_2_exponents=expected))
        if n<=64:
            assert [valuation(x) for x in original(n,n<=32).elementary_divisors()]==expected
            original_cases.append(n)
        if n<=40:continued_fraction(n)
        if n%16==0:print('ordinary ranks through',n,'passed',flush=True)
    paths=['tools/verify_symplectic_ordinary.py','SYMPLECTIC_K_ORDINARY_SMITH.md',
           'src/sage_borel/equivariant.py']
    report=dict(result='passed',seconds=round(perf_counter()-start,3),
        half_block_and_API_cases=half_cases,original_projective_inclusion_ranks=original_cases,
        independent_Euler_basis_ranks=list(range(1,33)),
        exact_finite_continued_fraction_ranks=list(range(1,41)),
        exact_2_local_Gram_certificates=gram_cases,
        proof_scope='All-rank formula is deduced in the proof note from Euler E750 VI-VII and orthogonal Hankel diagonalization (Miller-Stanton Theorem 1), with every recurrence denominator odd. These finite calculations are independent corroboration, not a proof of the arbitrary-rank quantifier or of priority.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'symplectic-ordinary-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASSED',len(half_cases),'half/API,',len(original_cases),'original inclusions,',len(gram_cases),'2-local Gram certificates',flush=True)


if __name__=='__main__':main()
