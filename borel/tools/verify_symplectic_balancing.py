"""Independent integral certificates and arithmetic checks for Hankel balancing.

The all-rank proof uses the cited Hankel product theorems and a monic
Groebner lift. The universal quartic identities below prove that finite
polynomial case over Z; finite specializations only corroborate the tables.
"""
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from sage.all import ZZ, GF, PolynomialRing, matrix, vector, prod, gcd, IntegerVectors

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b"\r\n",b"\n")).hexdigest()


def hankelization(r):
    R = PolynomialRing(ZZ,2*r,names=["p"+str(i) for i in range(r)]+["q"+str(i) for i in range(r)])
    U = PolynomialRing(R,"u")
    u = U.gen()
    p = u**r+sum(R.gen(i)*u**i for i in range(r))
    q = sum(R.gen(r+i)*u**i for i in range(r))
    H = matrix(R,r,r,lambda i,j:((q*u**j)%p)[i])
    T = matrix(R,r,r,lambda i,j:((u**(i+j))%p)[r-1])
    reverse = T.matrix_from_rows(list(reversed(range(r))))
    assert all(reverse[i,i]==1 for i in range(r))
    assert all(reverse[i,j]==0 for i in range(r) for j in range(i))
    moments = [((q*u**i)%p)[r-1] for i in range(2*r-1)]
    assert T*H==matrix(R,r,r,lambda i,j:moments[i+j])
    return r


def quartic_universal_identities():
    R = PolynomialRing(ZZ,7,names=["h"+str(i) for i in range(7)])
    h = R.gens()
    H = matrix(R,4,4,lambda i,j:h[i+j])
    def distinct_minors(j):
        answer,seen = [],set()
        for rows in combinations(range(4),j):
            for cols in combinations(range(4),j):
                value = H.matrix_from_rows_and_columns(rows,cols).det(algorithm="df")
                if value and value not in seen and -value not in seen:
                    answer.append((rows,cols,value))
                    seen.add(value)
        return answer
    K,L = distinct_minors(2),distinct_minors(3)
    monomials = list(IntegerVectors(4,7))
    lookup = {tuple(e):i for i,e in enumerate(monomials)}
    def coefficients(f):
        answer = vector(ZZ,len(monomials))
        for e,c in f.dict().items():
            answer[lookup[tuple(e)]] = c
        return answer
    pairs = [(i,j) for i in range(len(K)) for j in range(i,len(K))]
    M = matrix(ZZ,[coefficients(K[i][2]*K[j][2]) for i,j in pairs]).transpose()
    basis,transform = M.transpose().hermite_form(include_zero_rows=False,transformation=True)
    assert transform*M.transpose()==basis
    targets = [(i,j,h[i]*L[j][2]) for i in range(7) for j in range(len(L))]
    target_matrix = matrix(ZZ,[coefficients(t[2]) for t in targets]).transpose()
    coordinates = basis.transpose().solve_right(target_matrix)
    assert coordinates.denominator()==1
    identities = transform.transpose()*matrix(ZZ,coordinates)
    assert M*identities==target_matrix
    records = []
    for column,(i,j,target) in enumerate(targets):
        sparse = [(pairs[k][0],pairs[k][1],int(c)) for k,c in enumerate(identities.column(column)) if c]
        # Check as literal polynomial identities independently of the HNF coordinates.
        assert sum(c*K[a][2]*K[b][2] for a,b,c in sparse)==target
        records.append(dict(moment=i,three_minor=j,products=sparse))
    return dict(identity_count=len(records),integer_coefficients=True,
                maximum_products=max(len(r["products"]) for r in records),
                two_minors=[dict(rows=list(rows),columns=list(cols)) for rows,cols,_ in K],
                three_minors=[dict(rows=list(rows),columns=list(cols)) for rows,cols,_ in L],
                identities=records)


def quartic_generators():
    R = PolynomialRing(ZZ,8,names=("P","Q","R0","S","a","b","c","d"))
    P,Q,R0,S,a,b,c,d = R.gens()
    U = PolynomialRing(R,"u")
    u = U.gen()
    p,q = u**4+P*u**3+Q*u*u+R0*u+S,a*u**3+b*u*u+c*u+d
    H = matrix(R,4,4,lambda i,j:((q*u**j)%p)[i])
    B = (d*d-a*R0*d+a*c*S,c*d-a*Q*d+a*b*S,b*d-a*P*d+a*a*S,
         c*c-a*Q*c-b*d+a*b*R0,b*c-a*P*c-a*d+a*a*R0,b*b-a*P*b-a*c+a*a*Q)
    assert B==tuple(H.matrix_from_rows_and_columns(rows,[0,1]).det(algorithm="df")
                    for rows in combinations(range(4),2))
    v = (-S*(b-a*P),-a*S-R0*(b-a*P),d-a*R0-Q*(b-a*P),c-a*Q-P*(b-a*P))
    assert tuple(H.column(2))==v
    D = (v[0]*B[3]-v[1]*B[1]+v[2]*B[0],v[0]*B[4]-v[1]*B[2]+v[3]*B[0],
         v[0]*B[5]-v[2]*B[2]+v[3]*B[1],v[1]*B[5]-v[2]*B[4]+v[3]*B[3])
    assert D==tuple(H.matrix_from_rows_and_columns(rows,[0,1,2]).det(algorithm="df")
                    for rows in combinations(range(4),3))
    return dict(two_generators=6,three_generators=4,universal_integer_identities=True)


def original_inclusion_checks():
    checks = 0
    for n in (8,9):
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
            factors = matrix(ZZ,columns).transpose().elementary_divisors()
            e = [ZZ.one()]+[sum(prod(c) for c in combinations(values,j)) for j in range(1,n+1)]
            V = PolynomialRing(ZZ,"u")
            u = V.gen()
            p = u**4+e[2]*u**3+e[4]*u*u+e[6]*u+e[8]
            coefficients = ((e[1],e[3],e[5],e[7]) if n==8 else
                (e[1]*e[2]-e[3],e[1]*e[4]-e[5],e[1]*e[6]-e[7],e[1]*e[8]-e[9]))
            q = sum(coefficients[j]*u**(3-j) for j in range(4))
            H = matrix(ZZ,4,4,lambda i,j:((q*u**j)%p)[i])
            ideals = [ZZ.one()]+[abs(gcd([H.matrix_from_rows_and_columns(rows,list(range(j))).det()
                         for rows in combinations(range(4),j)])) for j in range(1,5)]
            for i in range(n):
                k = n-i
                ell = k-n%2
                expected = 2**k*ideals[ell//2]*ideals[(ell+1)//2]
                assert abs(prod(factors[:2*n-i]))==expected,(n,values,i)
                checks += 1
    return checks


def polynomial_fiber_checks():
    checks = 0
    for prime in (2,3):
        for r in (4,5):
            R = PolynomialRing(GF(prime),r,names=["q"+str(i) for i in range(r)])
            U = PolynomialRing(R,"u")
            u = U.gen()
            q = sum(R.gen(i)*u**i for i in range(r))
            for p in (u**r,u**r+u+1):
                H = matrix(R,r,r,lambda i,j:((q*u**j)%p)[i])
                ideals = {j:[H.matrix_from_rows_and_columns(rows,list(range(j))).det(algorithm="df")
                             for rows in combinations(range(r),j)] for j in range(1,r+1)}
                for a in range(1,r):
                    for b in range(a+2,r+1):
                        c,d = (a+b)//2,(a+b+1)//2
                        monomials = list(IntegerVectors(a+b,r))
                        lookup = {tuple(e):i for i,e in enumerate(monomials)}
                        def coefficients(f):
                            answer = vector(GF(prime),len(monomials))
                            for e,v in f.dict().items():
                                answer[lookup[tuple(e)]]=v
                            return answer
                        products = [x*y for x in ideals[c] for y in ideals[d]]
                        span = matrix(GF(prime),[coefficients(f) for f in products]).row_space()
                        assert all(coefficients(x*y) in span for x in ideals[a] for y in ideals[b])
                        checks += 1
    return checks


def main():
    ranks = [hankelization(r) for r in range(1,7)]
    quartic = quartic_universal_identities()
    generators = quartic_generators()
    original = original_inclusion_checks()
    fibers = polynomial_fiber_checks()
    paths = ["tools/verify_symplectic_balancing.py","SYMPLECTIC_K_BALANCING.md"]
    report = dict(result="passed",universal_hankelization_ranks=ranks,
        quartic=quartic,quartic_generators=generators,original_inclusion_checks=original,
        polynomial_fiber_balancing_checks=fibers,
        proof_scope="All-rank deduction from published generic Hankel product theorems and a monic integral lift. Universal quartic Z identities and finite fibers are checked separately; no primary decomposition is asserted after specialization.",
        source_sha256={name:digest(ROOT/name) for name in paths})
    (ROOT/"symplectic-balancing-verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k not in ("quartic","source_sha256")},indent=2))
    print("Universal quartic integer polynomial identities:",quartic["identity_count"])


if __name__=="__main__":
    main()
