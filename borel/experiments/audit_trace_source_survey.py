"""Independent standard-library audit: root strings, Freudenthal and minors.

No Sage or sage_borel imports, and no imports from other repository auditors.
The proof covers all degrees/primes; the finite data check its hypotheses
and independently reconstruct every reported representation's weights.
"""
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path
from time import perf_counter
import json

ROOT = Path(__file__).resolve().parents[1]


def determinant(A):
    A = [list(map(int, r)) for r in A]
    if not A:
        return 1
    sign, previous = 1, 1
    for k in range(len(A)-1):
        i = next((i for i in range(k,len(A)) if A[i][k]), None)
        if i is None:
            return 0
        if i != k:
            A[i], A[k] = A[k], A[i]
            sign = -sign
        pivot = A[k][k]
        for i in range(k+1,len(A)):
            lower = A[i][k]
            for j in range(k+1,len(A)):
                q, remainder = divmod(A[i][j]*pivot-lower*A[k][j], previous)
                assert remainder == 0
                A[i][j] = q
            A[i][k] = 0
        previous = pivot
    return sign*A[-1][-1]


def inverse(A):
    n = len(A)
    B = [[F(c) for c in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        i = next(i for i in range(k,n) if B[i][k])
        B[k], B[i] = B[i], B[k]
        pivot = B[k][k]
        B[k] = [a/pivot for a in B[k]]
        for i in range(n):
            if i != k:
                q = B[i][k]
                B[i] = [a-q*b for a,b in zip(B[i],B[k])]
    return [r[n:] for r in B]


def matmul(A,B):
    return [[sum(a*b for a,b in zip(r,c)) for c in zip(*B)] for r in A]


def cartan(kind):
    n = int(kind[1:])
    if kind == "G2":
        return [[2,-3],[-1,2]]
    if kind == "F4":
        return [[2,-1,0,0],[-1,2,-1,0],[0,-2,2,-1],[0,0,-1,2]]
    if kind == "B2":
        return [[2,-1],[-2,2]]
    C = [[2*int(i==j) for j in range(n)] for i in range(n)]
    edges = [(i,i+1) for i in range(n-1)] if kind[0] == "A" else [
        (0,2),(2,3),(1,3)]+[(i,i+1) for i in range(3,n-1)]
    for i,j in edges:
        C[i][j] = C[j][i] = -1
    return C


@lru_cache(None)
def root_data(kind):
    C = cartan(kind)
    n = len(C)
    d = [None]*n
    d[0] = F(1)
    pending = [0]
    while pending:
        i = pending.pop()
        for j in range(n):
            if i != j and C[i][j] and d[j] is None:
                d[j] = d[i]*F(C[i][j],C[j][i])
                pending.append(j)
    scale = lcm(*(a.denominator for a in d))
    lengths = [int(a*scale) for a in d]
    assert all(lengths[i]*C[i][j] == lengths[j]*C[j][i] for i in range(n) for j in range(n))
    roots = {tuple(int(i==j) for i in range(n)) for j in range(n)}
    pending = list(roots)
    while pending:
        alpha = pending.pop()
        for i in range(n):
            beta = list(alpha)
            beta[i] -= sum(C[i][j]*alpha[j] for j in range(n))
            beta = tuple(beta)
            if beta not in roots:
                roots.add(beta)
                pending.append(beta)
    positive = []
    for a in sorted(roots):
        if all(v>=0 for v in a):
            w = tuple(sum(C[i][j]*a[j] for j in range(n)) for i in range(n))
            norm = sum(a[i]*lengths[i]*w[i] for i in range(n))
            assert norm > 0
            positive.append((a,w,norm))
    return C, lengths, positive, inverse(C)


@lru_cache(None)
def weights_from_strings(kind, highest):
    C,d,positive,invC = root_data(kind)
    n = len(C)
    # Complete simple-root strings through known weights; finite norm bound.
    support, pending = {highest}, deque([highest])
    while pending:
        mu = pending.popleft()
        for i in range(n):
            for k in range(1,mu[i]+1):
                nu = tuple(mu[j]-k*C[j][i] for j in range(n))
                if nu not in support:
                    support.add(nu)
                    pending.append(nu)
    def root_coordinates(mu):
        return [sum(invC[i][j]*(highest[j]-mu[j]) for j in range(n)) for i in range(n)]
    height = {mu:sum(root_coordinates(mu)) for mu in support}
    assert all(h.denominator==1 and h>=0 for h in height.values())
    multiplicities = {highest:1}
    for mu in sorted(support,key=lambda mu:height[mu]):
        if mu == highest:
            continue
        a = root_coordinates(mu)
        assert all(v.denominator==1 and v>=0 for v in a)
        denominator = sum(a[i]*d[i]*(highest[i]+mu[i]+2) for i in range(n))
        total = 0
        for alpha,w,norm in positive:
            pairing = sum(d[i]*alpha[i]*mu[i] for i in range(n))
            nu, k = tuple(mu[i]+w[i] for i in range(n)), 1
            while nu in support:
                assert nu in multiplicities
                total += (pairing+k*norm)*multiplicities[nu]
                k += 1
                nu = tuple(nu[i]+w[i] for i in range(n))
        value = F(2*total,denominator)
        assert value.denominator==1 and value>0
        multiplicities[mu] = int(value)
    # A second classical identity checks the total dimension independently.
    dimension = F(1)
    for alpha,w,norm in positive:
        dimension *= F(sum(d[i]*alpha[i]*(highest[i]+1) for i in range(n)),
                       sum(d[i]*alpha[i] for i in range(n)))
    assert dimension.denominator==1 and sum(multiplicities.values())==dimension
    for mu,m in multiplicities.items():
        for i in range(n):
            reflected = tuple(mu[j]-mu[i]*C[j][i] for j in range(n))
            assert multiplicities.get(reflected) == m
    return multiplicities


def row_reduce(A,p):
    A = [[int(c)%p for c in r] for r in A]
    pivot, columns = 0, []
    for j in range(len(A[0])):
        i = next((i for i in range(pivot,len(A)) if A[i][j]),None)
        if i is None:
            continue
        A[pivot], A[i] = A[i], A[pivot]
        q = pow(A[pivot][j],-1,p)
        A[pivot] = [q*c%p for c in A[pivot]]
        for i in range(len(A)):
            if i != pivot and A[i][j]:
                q = A[i][j]
                A[i] = [(a-q*b)%p for a,b in zip(A[i],A[pivot])]
        columns.append(j)
        pivot += 1
        if pivot == len(A):
            break
    return A,columns


def multiply(a,b):
    out = {}
    for u,c in a.items():
        for v,d in b.items():
            e = tuple(x+y for x,y in zip(u,v))
            out[e] = out.get(e,0)+c*d
    return {e:c for e,c in out.items() if c}


def bounded_polynomials(polys,degree):
    zero = (0,)*len(next(iter(polys[0])))
    def powers(n,d):
        if n==1:
            yield (d,)
        else:
            for a in range(d+1):
                for tail in powers(n-1,d-a):
                    yield (a,)+tail
    result = []
    for d in range(degree+1):
        for e in powers(len(polys),d):
            f = {zero:1}
            for g,k in zip(polys,e):
                for _ in range(k):
                    f = multiply(f,g)
            result.append(f)
    return result


def polynomial_column_rank(polys,p):
    basis = {}
    support = sorted(set().union(*(set(f) for f in polys)))
    for e in support:
        row = [f.get(e,0)%p for f in polys]
        for j,r in sorted(basis.items()):
            if row[j]:
                q = row[j]
                row = [(a-q*b)%p for a,b in zip(row,r)]
        j = next((j for j,c in enumerate(row) if c),None)
        if j is not None:
            q = pow(row[j],-1,p)
            basis[j] = [q*c%p for c in row]
            if len(basis) == len(polys):
                return len(basis)
    return len(basis)


def main():
    started = perf_counter()
    path = ROOT/"experiments/trace-source-survey-results.json"
    data = json.loads(path.read_text())
    assert data["result"] == "passed" and len(data["cases"]) == 59
    prime_count, dimensions, nonzero_weights = 0, 0, 0
    for number,row in enumerate(data["cases"]):
        kind = row["cartan_type"]
        C,d,positive,invC = root_data(kind)
        n = len(C)
        assert C == [[int(v) for v in r] for r in row["cartan_matrix"]]
        weights = {}
        for highest in row["highest_weights"]:
            for mu,m in weights_from_strings(kind,tuple(highest)).items():
                weights[mu] = weights.get(mu,0)+m
        assert weights == {tuple(mu):m for mu,m in row["weights"]}
        assert sum(weights.values()) == row["dimension"]
        M = [[sum(m*mu[i]*mu[j] for mu,m in weights.items()) for j in range(n)] for i in range(n)]
        assert M == [[int(v) for v in r] for r in row["trace_gram"]]
        B = [[int(v) for v in r] for r in row["linear_basis"]]
        assert determinant(M) and determinant(B)
        K = matmul(list(map(list,zip(*B))),inverse(M))
        support = sorted((mu for mu in weights if any(mu)),reverse=True)
        assert [list(mu) for mu in support] == row["support"]
        coefficients = [[weights[mu]*sum(K[j][i]*mu[i] for i in range(n)) for j in range(n)] for mu in support]
        assert all(v.denominator==1 for r in coefficients for v in r)
        coefficients = [[int(v) for v in r] for r in coefficients]
        assert coefficients == [[int(v) for v in r] for r in row["coefficients"]]
        divisor = 0
        for minor in row["saturation_certificate"]["maximal_minors"]:
            det = abs(determinant([coefficients[i] for i in minor["rows"]]))
            assert det == int(minor["determinant"])
            divisor = gcd(divisor,det)
        assert divisor == 1 == row["saturation_certificate"]["minor_gcd"]
        cert = row["certificate"]
        assert cert == dict(certified=True,rank=n,support_size=len(support),
                            primitive_coefficients=True,trace_factorization=True)
        # Both verify the Sage basis and independently choose a different order.
        reverse_order = sorted(range(len(support)),key=lambda i:tuple(reversed(support[i])),reverse=True)
        for check in row["prime_checks"]:
            p = check["prime"]
            U = [[int(v) for v in r] for r in check["change"]]
            assert determinant(U)%p
            transformed = matmul(U,list(map(list,zip(*coefficients))))
            leading = [next(j for j,c in enumerate(r) if c%p) for r in transformed]
            assert leading == check["pivots"]
            assert [list(support[i]) for i in leading] == check["leading_exponents"]
            det = determinant([support[i] for i in leading])
            assert det == int(check["exponent_determinant"]) and det
            assert determinant([coefficients[i] for i in leading])%p
            A = [[coefficients[i][j] for i in reverse_order] for j in range(n)]
            reduced,pivots = row_reduce(A,p)
            assert len(pivots)==n
            E = [support[reverse_order[i]] for i in pivots]
            assert determinant(E)
            prime_count += 1
        dimensions += row["dimension"]
        nonzero_weights += len(support)
        print("Audited",number+1,kind,row["highest_weights"],"dimension",row["dimension"],
              "elapsed",round(perf_counter()-started,2),flush=True)
    expansions = {}
    for check in data["bounded_searches"]:
        case = check["case"]
        if case not in expansions:
            row = data["cases"][case]
            fs = [{tuple(mu):int(r[j]) for mu,r in zip(row["support"],row["coefficients"]) if int(r[j])}
                  for j in range(len(row["linear_basis"]))]
            expansions[case] = bounded_polynomials(fs,check["max_degree"])
        fs = expansions[case]
        assert polynomial_column_rank(fs,check["prime"]) == len(fs)
        assert check["relation"] is None
    assert prime_count == 354 and len(data["bounded_searches"]) == 18
    # Adversarial control: primitive coefficients without trace factorization.
    assert determinant([[0,1],[1,0]]) and determinant([[2,0],[1,0]]) == 0
    f = {(1,0):1}
    g = {(2,0):1,(0,1):2}
    assert polynomial_column_rank(bounded_polynomials([f,g],2),2) < 6
    # A zero Jacobian in characteristic two does not imply a relation.
    assert determinant([[2,0],[0,2]]) == 4
    assert polynomial_column_rank(bounded_polynomials([{(2,0):1},{(0,2):1}],2),2) == 6
    files = [path,Path(__file__).resolve(),ROOT/"TRACE_SOURCE_THEOREM.md"]
    result = dict(result="passed",engine="Python standard library only",
                  independently_reconstructed_representations=59,
                  prime_adapted_checks=prime_count, bounded_searches=18,
                  total_dimensions=dimensions, nonzero_weight_rows=nonzero_weights,
                  controls="Non-trace primitive Laurent counterexample and inseparable independent polynomials.",
                  scope="Independent weight strings/Freudenthal, Weyl dimensions, trace and linear lattices, small-minor saturation, two orders of leading bases, bounded expansions. No Borel-image computations are asserted.",
                  seconds=round(perf_counter()-started,3),
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n",b"\n")).hexdigest() for p in files})
    (ROOT/"experiments/trace-source-survey-audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="source_sha256"},indent=2),flush=True)


if __name__ == "__main__":
    main()
