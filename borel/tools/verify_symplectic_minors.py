"""Independent integral Schur certificates for the monic multiplication minor theorem.

This corroborates the written arbitrary-ring argument; finite checks do not
prove its quantifiers or establish priority. No production minor routine is used.
"""
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path

from sage.all import ZZ, GF, PolynomialRing, matrix, vector, prod, gcd

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b"\r\n",b"\n")).hexdigest()


def compound(M,j):
    indices = list(combinations(range(M.nrows()),j))
    return matrix(M.base_ring(),len(indices),len(indices),
                  lambda a,b:M.matrix_from_rows_and_columns(indices[a],indices[b]).det(algorithm="df"))


def schur_certificate(r,j):
    R = PolynomialRing(ZZ,2*r,names=["p"+str(i) for i in range(r)]+["q"+str(i) for i in range(r)])
    U = PolynomialRing(R,"u")
    u = U.gen()
    p = u**r+sum(R.gen(i)*u**i for i in range(r))
    q = sum(R.gen(r+i)*u**i for i in range(r))
    C = matrix(R,r,r,lambda a,b:((u**(b+1))%p)[a])
    H = matrix(R,r,r,lambda a,b:((q*u**b)%p)[a])
    assert C*H==H*C
    T = PolynomialRing(R,"t")
    t = T.gen()
    exterior = compound(matrix(T,r,r,lambda a,b:int(a==b)+t*C[a,b]),j)
    indices = list(combinations(range(r),j))
    d = len(indices)
    E = [matrix(R,d,d,lambda a,b:exterior[a,b][k]) for k in range(j+1)]
    identity,zero = matrix.identity(R,d),matrix(R,d,d)
    h = [identity]
    for k in range(1,r):
        h.append(sum(((-1)**(ell-1)*E[ell]*h[k-ell] for ell in range(1,min(j,k)+1)),zero))
    W = compound(H,j)
    seed = vector(R,[1]+[0]*(d-1))
    checked = 0
    for column,I in enumerate(indices):
        partition = tuple(I[j-1-a]-(j-1-a) for a in range(j))
        operator = matrix(R,d,d)
        for perm in permutations(range(j)):
            degrees = [partition[a]-a+perm[a] for a in range(j)]
            if min(degrees,default=0)<0:
                continue
            sign = (-1)**sum(perm[a]>perm[b] for a in range(j) for b in range(a+1,j))
            term = identity
            for degree in degrees:
                term = term*h[degree]
            operator += sign*term
        expected = vector(R,[int(a==column) for a in range(d)])
        assert operator*seed==expected
        # Explicit integral ideal membership for every noninitial column:
        # each full j-minor is an R-linear combination of initial-column minors.
        assert operator*W.column(0)==W.column(column)
        checked += d
    return dict(rank=r,minor_size=j,full_minors=checked,initial_minors=d,
                integral_schur_identities=True)


def low_rank_formulas():
    R = PolynomialRing(ZZ,6,names=("P","Q","R0","a","b","c"))
    P,Q,R0,a,b,c = R.gens()
    U = PolynomialRing(R,"u")
    u = U.gen()
    p,q = u**3+P*u**2+Q*u+R0,a*u**2+b*u+c
    H = matrix(R,3,3,lambda i,j:((q*u**j)%p)[i])
    K = (c*c-Q*a*c+R0*a*b,b*c-P*a*c+R0*a*a,
         b*b-P*a*b-a*c+Q*a*a)
    actual = tuple(H.matrix_from_rows_and_columns(rows,[0,1]).det(algorithm="df")
                   for rows in combinations(range(3),2))
    assert K==actual
    adj = H.adjugate()
    # Gives I_1(H)*I_3(H) contained in I_2(H)^2 without denominators.
    assert adj.adjugate()==H.det(algorithm="df")*H
    # Rank five has the same universal 2-by-2 block as rank four.
    S = PolynomialRing(ZZ,5,names=["s"+str(i) for i in range(1,6)])
    s = S.gens()
    e = [S.one()]+[sum(prod(x) for x in combinations(s,k)) for k in range(1,6)]
    A,B = e[1]*e[4]-e[5],e[1]*e[2]-e[3]
    F = prod(s[i]+s[j] for i in range(5) for j in range(i+1,5))
    assert A*A-e[2]*A*B+e[4]*B*B==-F
    return dict(rank_five_resultant_identity=True,cubic_initial_minors=list(map(str,K)),
                cubic_adjugate_identity=True)


def specialization_checks():
    counts = dict(integer=0,finite_field=0)
    for r in range(1,7):
        S = PolynomialRing(ZZ,"u")
        u = S.gen()
        for p,q in ((u**r,sum((i+1)*u**i for i in range(r))),
                    (prod(u-i for i in range(r)),prod(u-i for i in range(r//2))),
                    (u**r+2*u+2,2*u**max(r-1,0)+4)):
            H = matrix(ZZ,r,r,lambda a,b:((q*u**b)%p)[a])
            smith = H.elementary_divisors()
            for j in range(1,r+1):
                initial = [H.matrix_from_rows_and_columns(rows,list(range(j))).det()
                           for rows in combinations(range(r),j)]
                assert abs(gcd(initial))==abs(prod(smith[:j]))
                counts["integer"]+=1
            for prime in (2,3,5):
                M = H.change_ring(GF(prime))
                for j in range(1,r+1):
                    assert M[:,:j].rank()==min(M.rank(),j)
                    counts["finite_field"]+=1
    return counts


def original_inclusion_checks():
    """Every closed rank 5/6/7 ideal after integer base change, from 2n columns."""
    checked = 0
    for n in (5,6,7):
        for values in ([2]*n,list(range(1,n+1)),[1,-1]+[2]*(n-2),[0]*n):
            U = PolynomialRing(ZZ,"L")
            L = U.gen()
            relation = prod(L*L-s*L+1 for s in values)
            inverse = -sum(relation[k]*L**(k-1) for k in range(1,2*n+1))
            b = (L-inverse)%relation
            columns,power = [],U.one()
            for j in range(2*n):
                columns.append([power[i] for i in range(2*n)])
                power = power*b%relation
            smith = matrix(ZZ,columns).transpose().elementary_divisors()
            e = [ZZ.one()]+[sum(prod(c) for c in combinations(values,j)) for j in range(1,n+1)]
            F = prod(values[i]+values[j] for i in range(n) for j in range(i+1,n))
            if n==5:
                A,B = e[1]*e[4]-e[5],e[1]*e[2]-e[3]
                g = abs(gcd(A,B))
                expected = [32*F*F,16*abs(F)*g,8*g*g,4*g,2]
            else:
                P,Q,R0 = e[2],e[4],e[6]
                a,b,c = ((e[1],e[3],e[5]) if n==6 else
                         (e[1]*e[2]-e[3],e[1]*e[4]-e[5],e[1]*e[6]-e[7]))
                g = abs(gcd([a,b,c]))
                h = abs(gcd([c*c-Q*a*c+R0*a*b,b*c-P*a*c+R0*a*a,
                             b*b-P*a*b-a*c+Q*a*a]))
                expected = [2**n*F*F,2**(n-1)*abs(F)*h,2**(n-2)*h*h,
                            2**(n-3)*g*h,2**(n-4)*g*g,2**(n-5)*g]
                if n==7:
                    expected.append(2)
            for i,value in enumerate(expected):
                assert abs(prod(smith[:2*n-i]))==value,(n,values,i)
                checked += 1
    return checked


def main():
    certificates = [schur_certificate(r,j) for r in range(1,5) for j in range(1,r+1)]
    low_rank = low_rank_formulas()
    fibers = specialization_checks()
    original = original_inclusion_checks()
    paths = ["tools/verify_symplectic_minors.py","SYMPLECTIC_K_MINORS.md"]
    report = dict(result="passed",symbolic_certificates=certificates,
                  full_minor_identities=sum(c["full_minors"] for c in certificates),
                  low_rank=low_rank,specializations=fibers,
                  original_inclusion_closed_ideal_checks=original,
                  scope="Universal coefficient matrix identities through degree four and finite arithmetic fibers; the arbitrary-ring proof is written separately. Not external review or a priority claim.",
                  source_sha256={name:digest(ROOT/name) for name in paths})
    (ROOT/"symplectic-minors-verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
