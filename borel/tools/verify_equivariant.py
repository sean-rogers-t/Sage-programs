"""Bind classical equivariant code/proofs to regression and rank-nine checks."""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from itertools import combinations
from sage.all import ZZ, PolynomialRing, matrix, block_diagonal_matrix, prod
from sage_borel import type_a_equivariant, symplectic_projective_k_image

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b"\r\n",b"\n")).hexdigest()


def verify_doubling(n):
    """Independent universal-coefficient matrix certificates, not API matrices."""
    R = PolynomialRing(ZZ,n,names=["e"+str(i) for i in range(1,n+1)])
    e = (R.one(),)+R.gens()
    S = PolynomialRing(R,"z")
    z = S.gen()
    even = S(sum((-1)**(n-2*j)*e[n-2*j]*(z*z+4)**j for j in range(n//2+1)))
    q = S(sum((-1)**(n-2*j-1)*e[n-2*j-1]*(z*z+4)**j for j in range((n+1)//2)))
    G = even+z*q
    def multiplication(multiplier,modulus):
        d = modulus.degree()
        return matrix(R,d,d,lambda i,j:((multiplier*z**j)%modulus)[i])
    N = multiplication(q,G)
    if n%2==0:
        T = matrix(R,n,n,lambda i,j:(z**j+z*(q*z**j).quo_rem(even)[0])[i])
        assert T.det()==1 and N*T==multiplication(q,even)
        parity_matrix,modulus,multiplier = multiplication(q,even),even,q
    else:
        B = matrix(R,n,n,lambda i,j:(q if j==0 else z**(j-1))[i])
        assert B.det() in (R.one(),-R.one())
        A = matrix(R,B.inverse()*N)
        assert list(A.column(0))==[1]+[0]*(n-1)
        K = matrix(R,n-1,n-1,lambda i,j:(q*z**(j+1)).quo_rem(G)[0][i])
        assert K.det()==1
        parity_matrix = multiplication(even,q)
        assert A[1:,1:]==-parity_matrix*K
        modulus,multiplier = q,even
    r = n//2
    order = list(range(0,2*r,2))+list(range(1,2*r,2))
    parity_matrix = parity_matrix.matrix_from_rows_and_columns(order,order)
    # Even polynomials descend to v=z^2. The translation u=v+4 is unimodular.
    V = PolynomialRing(R,"v")
    v = V.gen()
    p = V([modulus[2*j] for j in range(r+1)])
    f = V([multiplier[2*j] for j in range(n//2+1)])
    Hshift = matrix(R,r,r,lambda i,j:((f*v**j)%p)[i])
    assert parity_matrix==block_diagonal_matrix([Hshift,Hshift])
    p0,f0 = p(v-4),f(v-4)
    H = matrix(R,r,r,lambda i,j:((f0*v**j)%p0)[i])
    U = matrix(R,r,r,lambda i,j:((v+4)**j)[i])
    assert U.det()==1 and Hshift*U==U*H
    return n


def verify_rank_four_minor_ideals():
    R = PolynomialRing(ZZ,4,names=("a","y","v","w"))
    a,y,v,w = R.gens()
    H = matrix(R,[[y,-a*w],[a,y-a*v]])
    D = y*y-a*v*y+a*a*w
    assert H.det()==D
    M = block_diagonal_matrix([H,H])
    minors = {}
    for k in range(1,5):
        minors[k] = {M.matrix_from_rows_and_columns(rows,cols).det()
                     for rows in combinations(range(4),k) for cols in combinations(range(4),k)}
    # Entry and size-two ideals: all entries have (a,y)-degree one;
    # displayed a,y and their three quadratic products occur among the minors.
    assert {a,y} <= minors[1]
    assert {a*a,a*y,y*y} <= minors[2]
    for k in (1,2):
        assert all(exponent[0]+exponent[1]>=k for value in minors[k] for exponent in value.dict())
    assert all(value==0 or any(value==sign*D*entry for entry in H.list() for sign in (-1,1))
               for value in minors[3])
    assert all(any(sign*value in minors[3] for sign in (-1,1)) for value in (D*a,D*y))
    assert minors[4]=={D*D}
    return True


def main():
    regression = json.loads((ROOT/"character-verification.json").read_text())
    assert regression["result"]=="passed"
    for name in ("src/sage_borel/equivariant.py","src/sage_borel/__init__.py","tests/test_equivariant.py"):
        assert regression["source_sha256"][name]==digest(ROOT/name)
    minors = json.loads((ROOT/"symplectic-minors-verification.json").read_text())
    assert minors["result"]=="passed" and minors["full_minor_identities"]==94
    assert minors["original_inclusion_closed_ideal_checks"]==72
    for name,expected in minors["source_sha256"].items():
        assert digest(ROOT/name)==expected,name
    balancing = json.loads((ROOT/"symplectic-balancing-verification.json").read_text())
    assert balancing["result"]=="passed" and balancing["quartic"]["identity_count"]==70
    assert balancing["original_inclusion_checks"]==68
    for name,expected in balancing["source_sha256"].items():
        assert digest(ROOT/name)==expected,name
    doubling_ranks = [verify_doubling(n) for n in range(1,9)]
    rank_four_ideals = verify_rank_four_minor_ideals()

    # The displayed rank-three unimodular reduction, over independent
    # elementary parameters rather than a substituted numerical point.
    ring = PolynomialRing(ZZ,3,names=("e1","e2","e3"))
    e1,e2,e3 = ring.gens()
    d,c,q = e2+4,4*e1+e3,e3-e1*e2
    N = matrix(ring,[[d,c,e1*c],[0,0,q],[1,e1,e1**2]])
    N.swap_rows(0,2)
    N.add_multiple_of_column(1,0,-e1)
    N.add_multiple_of_column(2,0,-e1**2)
    N.add_multiple_of_row(2,0,-d)
    N.swap_columns(1,2)
    N.add_multiple_of_row(2,1,-e1)
    assert N==matrix.diagonal(ring,[1,q,q])

    start = perf_counter()
    a = type_a_equivariant(10,1,theory="K")
    f = 10*a.weights[0]
    answer = a.expand(f)
    assert answer.restrictions()==a.restrictions(f)
    type_a_seconds = perf_counter()-start
    assert len(a.basis)==10

    start = perf_counter()
    c = symplectic_projective_k_image(9,spectral=True)
    factors = c.ordinary_smith_factors()
    assert c._matrix is None and len(factors)==18
    assert sum(ZZ(s).valuation(2) for s in factors)==9*17
    symplectic_seconds = perf_counter()-start
    assert factors[-1]==2**(4*8-ZZ(8).popcount())

    paths = [
        "src/sage_borel/equivariant.py","src/sage_borel/__init__.py",
        "tests/test_equivariant.py","CLASSICAL_EQUIVARIANT.md",
        "SYMPLECTIC_K_FITTING.md","tools/verify_equivariant.py",
        "SYMPLECTIC_K_MINORS.md","tools/verify_symplectic_minors.py",
        "symplectic-minors-verification.json",
        "SYMPLECTIC_K_BALANCING.md","tools/verify_symplectic_balancing.py",
        "symplectic-balancing-verification.json",
        "character-verification.json",
    ]
    report = dict(
        result="passed",new_test_methods=16,
        full_regression_methods=regression["test_methods_run"],
        type_A=dict(group="SL10",parabolic="P1",root_rank=9,fixed_points=10,
                    theory="K",exact_restriction_roundtrip=True,seconds=round(type_a_seconds,4)),
        type_C=dict(group="Sp18",parabolic="P1",rank=9,
                    ordinary_smith=list(factors),index_exponent=153,
                    symbolic_matrix_avoided=True,seconds=round(symplectic_seconds,4)),
        symbolic_rank_three_unimodular_reduction=True,
        symbolic_doubling_unimodular_certificates=doubling_ranks,
        symbolic_rank_four_minor_ideal_certificates=rank_four_ideals,
        initial_column_minor_theorem=dict(symbolic_identities=minors["full_minor_identities"],
            original_inclusion_checks=minors["original_inclusion_closed_ideal_checks"],
            literature_input="Integral Schur cyclicity of exterior powers; Laksov-Thorup and Gatto-Santiago"),
        hankel_balancing=dict(universal_quartic_integer_identities=70,
            original_inclusion_checks=68,all_rank_integral_proof=True,
            literature_input="Published generic Hankel product/monic Groebner theorems; coefficient-pairing equivalence and integral lift"),
        proof_scope="All-rank doubled half-size presentation, initial-column/Schur generators and a single balanced product for every Fitting ideal; explicit lists through rank nine. Global cyclic decompositions only in ranks 1-3; specialized primary structure and full ordinary Smith multiplicities remain.",
        API_scope="Type A all-parabolic input evaluation and literal integrality; automatic complete source generators only on full flags and P1. Type C projective integral image/cokernel.",
        source_sha256={name:digest(ROOT/name) for name in paths},
    )
    (ROOT/"equivariant-verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:v for k,v in report.items() if k!="source_sha256"},indent=2))


if __name__=="__main__":
    main()
