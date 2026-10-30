"""Sage-free checks using the original adjoint matrices on sl_N."""
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


@lru_cache(None)
def algebra(N):
    basis=[]
    for i in range(N):
        for j in range(N):
            if i!=j:
                E=sp.zeros(N); E[i,j]=1; basis.append(E)
    for i in range(N-1):
        H=sp.zeros(N); H[i,i]=1; H[N-1,N-1]=-1; basis.append(H)
    def coordinates(M):
        assert sp.trace(M)==0
        return sp.Matrix([M[i,j] for i in range(N) for j in range(N) if i!=j]+
                         [M[i,i] for i in range(N-1)])
    ad=[sp.Matrix.hstack(*(coordinates(X*Y-Y*X) for Y in basis)) for X in basis]
    entries=[[(i,j,A[i,j]) for i in range(A.rows) for j in range(A.cols) if A[i,j]] for A in ad]
    gram=sp.Matrix(len(basis),len(basis),lambda i,j:2*N*sp.trace(basis[i]*basis[j]))
    assert gram.det()!=0
    for i,A in enumerate(ad):
        for j,terms in enumerate(entries):
            assert sum(A[b,a]*c for a,b,c in terms)==gram[j,i]
    return basis,coordinates,ad,entries,gram


def theta(g):
    N=g.rows; inverse=g.inv()
    return (sp.trace(inverse)*g-sp.trace(g)*inverse)/(2*N)


def two_block_point(t,p,unipotent=False):
    N=2*t+p; w=sp.Rational(2)
    vs=[sp.Rational(j+2) for j in range(p-1)]
    vs.append(1/(w**N*sp.prod(vs)))
    ss=[sp.Rational(i+3) for i in range(t-1)]
    if unipotent: assert t>=2; ss[0]=sp.Rational(2)
    ss.append(-sum(v+1/v for v in vs)/2-sum(ss))
    blocks=[w*sp.Matrix([[0,-1],[1,s]]) for s in ss]
    blocks += [sp.Matrix([[w*v]]) for v in vs]
    g=sp.diag(*blocks); assert g.det()==1
    H=sum(1/v-v for v in vs)/(4*N)
    lambdas=[H*s for s in ss]; nus=[H*(v+1/v) for v in vs]
    assert H!=0 and len(set(lambdas+nus))==t+p
    target=sp.diag(*[value for lam in lambdas for value in (lam,lam)],*nus)
    assert theta(g)==target and sp.trace(target)==0
    return g,dict(double_blocks=t,singletons=p,nonsemisimple_source=unipotent,
                  target_eigenvalues=[str(a) for a in lambdas+nus],expected_tangent_dimension=2*t)


def torus_point(N):
    values=[sp.Rational(i+2) for i in range(N-1)]
    values.append(1/sp.prod(values)); g=sp.diag(*values); h=theta(g)
    assert len(set(h.diagonal()))==N
    return g,dict(double_blocks=0,singletons=N,nonsemisimple_source=False,
                  target_eigenvalues=[str(a) for a in h.diagonal()],expected_tangent_dimension=0)


def check(g,description):
    N=g.rows; basis,coordinates,ad,entries,gram=algebra(N); dim=len(basis)
    inverse=g.inv(); A=sp.trace(g); B=sp.trace(inverse)
    Ad=sp.Matrix.hstack(*(coordinates(g*Y*inverse) for Y in basis))
    assert sp.trace(Ad)==A*B-1
    point=theta(g)
    for Y,terms in zip(basis,entries):
        assert sum(Ad[b,a]*c for a,b,c in terms)==2*N*sp.trace(point*Y)
    direct_pairings=sp.zeros(dim); derivative_columns=[]
    for j,(X,adX) in enumerate(zip(basis,ad)):
        derivative=(-sp.trace(X*inverse)*g+B*g*X-sp.trace(g*X)*inverse+A*X*inverse)/(2*N)
        derivative_columns.append(coordinates(derivative))
        product=Ad*adX
        for i,terms in enumerate(entries): direct_pairings[i,j]=sum(product[b,a]*c for a,b,c in terms)
    D=sp.Matrix.hstack(*derivative_columns)
    assert direct_pairings==gram*D
    rank=D.rank(); tangent=dim-rank
    assert tangent==description['expected_tangent_dimension'],(N,description,tangent)
    commutators=sp.Matrix.hstack(*(coordinates(g*X-X*g) for X in basis))
    assert dim-commutators.rank()==N-1
    return dict(N=N,representation_dimension=dim,**description,
                original_matrix_trace_pairings=dim,original_matrix_differential_pairings=dim*dim,
                full_differential_rank=rank,fiber_tangent_dimension=tangent,
                source_centralizer_dimension=N-1,determinant_one=True,
                independent_original_representation=True)


def main():
    points=[torus_point(N) for N in (3,4,5)]
    points += [two_block_point(t,p) for t,p in [(1,1),(1,2),(1,3),(2,1)]]
    points += [two_block_point(2,1,True)]
    rows=[]
    for g,description in points:
        row=check(g,description); rows.append(row)
        print('Original adjoint matrices passed: SL'+str(g.rows),description,flush=True)
    # A regular zero-fiber control: it must not be assigned a proper target Levi.
    g=sp.diag(sp.Matrix([[0,-1],[1,3]]),sp.Matrix([[0,-1],[1,-3]]))
    assert theta(g)==sp.zeros(4)
    basis,coordinates,_,_,_=algebra(4)
    centralizer=len(basis)-sp.Matrix.hstack(*(coordinates(g*X-X*g) for X in basis)).rank()
    assert centralizer==3
    report=dict(result='passed',date='2026-10-23',uses_Sage=False,
                original_adjoint_matrix_checks=rows,
                trace_form_pairings=sum((N*N-1)**2 for N in (3,4,5)),
                trace_pairings=sum(r['original_matrix_trace_pairings'] for r in rows),
                differential_pairings=sum(r['original_matrix_differential_pairings'] for r in rows),
                regular_zero_fiber_control=dict(N=4,centralizer_dimension=centralizer,
                    target_centralizer='SL4; not a proper Levi',adjoint_not_fundamental=True),
                scope='Exact original Ad/ad matrices and full group differentials at eight rational points, including one regular nonsemisimple source. Corroborates the separate all-rank proofs, not a numerical extrapolation.',
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for name in ['QUADRATIC_BRANCHING_FIBERS.md',
                                            'tools/audit_type_a_adjoint_fibers.py']})
    (ROOT/'type-a-adjoint-fibers-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','original_adjoint_matrix_checks')},indent=2))


if __name__=='__main__': main()
