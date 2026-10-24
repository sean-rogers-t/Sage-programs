"""Exact checks of the type C component proof, separate from its all-rank argument."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import (QQ, CyclotomicField, PolynomialRing, matrix,
                      identity_matrix, block_diagonal_matrix, binomial)
from sage.env import SAGE_VERSION
from verify_symplectic_centralizer_realization import (
    spectrum, invariant_form, contraction, symplectic_basis)
from verify_type_a_character_degeneration import wedge_lie, compound

ROOT=Path(__file__).resolve().parents[1]


def elementary(g):
    p=g.charpoly(); n=g.nrows()
    return [(-1)**j*p[n-j] for j in range(n+1)]


def raw(g,k,e=None):
    if e is None: e=elementary(g)
    return sum(((-1)**(j-1)*e[k-j]*g**j for j in range(1,k+1)),
               matrix(g.base_ring(),g.nrows()))


def index(m):
    return QQ(binomial(4*m-2,2*m-1)-binomial(4*m-2,2*m-3))


def theta(g,J,m):
    F=raw(g,2*m)-raw(g,2*m-2)
    return (F-J.inverse()*F.transpose()*J)/(2*index(m))


def jacobian(g,J,m,basis):
    n=g.nrows(); e=elementary(g); inv=J.inverse()
    powers=[identity_matrix(QQ,n)]+[g**j for j in range(1,2*m+1)]
    gradients=[raw(g,j,e) for j in range(2*m)]
    result=[]
    for X in basis:
        dg=g*X; dF=matrix(QQ,n)
        for k,sign in ((2*m,1),(2*m-2,-1)):
            for j in range(1,k+1):
                de=(gradients[k-j]*X).trace()
                dp=sum((powers[r]*dg*powers[j-1-r] for r in range(j)),matrix(QQ,n))
                dF+=sign*(-1)**(j-1)*(de*powers[j]+e[k-j]*dp)
        result.append(((dF-inv*dF.transpose()*J)/(2*index(m))).list())
    return matrix(QQ,result).transpose()


def main():
    pairings=[]
    R=PolynomialRing(QQ,'z'); z=R.gen()
    for m in range(1,21):
        v=[z**k+z**(2*m-k) for k in range(1,m)]+[z**m]
        L=matrix(QQ,[[((1-z*z)*u*w)[2*m] for w in v] for u in v])
        E=identity_matrix(QQ,m)
        for j in range(m-2): E[j,j+2]=-1
        assert L==E*E.transpose() and L.det()==1
        Q0=sum(z**(2*j) for j in range(m+1))
        assert (1-z*z)*Q0==1-z**(2*m+2)
        assert all(((1-z*z)*Q0*u)[2*m]==0 for u in v)
        pairings.append(dict(m=m,determinant=int(L.det()),factorization_verified=True))

    collisions=[]
    for m in range(1,41):
        F=CyclotomicField(2*m+2); w=F.gen()
        beta=[w**j for j in range(1,2*m+2) if j!=m+1]
        values=[(b*b-1)*(1-b)**(2*m)/(2*index(m)) for b in beta]
        assert all(values)
        equal=[(j,k) for j in range(2*m) for k in range(j+1,2*m) if values[j]==values[k]]
        # Check the closed target formula directly from the torus character.
        RF=PolynomialRing(F,'z'); zz=RF.gen()
        D=(1+zz)**(2*m)*sum(zz**(2*j) for j in range(m+1))
        for b,h in zip(beta,values):
            left,rem=D.quo_rem(1+b*zz); assert rem==0
            right,rem=D.quo_rem(1+zz/b); assert rem==0
            direct=((1-zz*zz)*zz*(b*left-right/b))[2*m]/(2*index(m))
            assert direct==h
        collisions.append(dict(m=m,equal_external_pairs=len(equal)))

    ranks=[]; traces=[]
    for m in range(1,5):
        _,Q,B=spectrum(m,[QQ(0)]*m)
        J0=matrix(QQ,2*m)
        for j in range(m): J0[j,m+j]=1; J0[m+j,j]=-1
        J=block_diagonal_matrix(J0,invariant_form(B))
        Y=matrix(QQ,2*m)
        for j in range(m-1): Y[j,j+1]=1; Y[m+j+1,m+j]=-1
        Y[m-1,2*m-1]=1
        I=identity_matrix(QQ,2*m)
        g=block_diagonal_matrix((I+Y)*(I-Y).inverse(),B)
        t=block_diagonal_matrix(I,B)
        h=theta(t,J,m); assert theta(g,J,m)==h
        basis=symplectic_basis(J)
        kdim=len(basis)-matrix(QQ,[(h*X-X*h).list() for X in basis]).rank()
        dimension=kdim-2*m
        D=jacobian(g,J,m,basis)
        assert len(basis)-D.rank()==dimension==2*m*m
        ranks.append(dict(m=m,target_centralizer_dimension=kdim,
                          component_dimension=dimension,differential_rank=D.rank(),
                          ambient_dimension=len(basis)))
        if m<=2:
            C=contraction(J,2*m); V=C.right_kernel().basis_matrix().transpose()
            rows=list(V.transpose().pivots()); left=V.matrix_from_rows(rows).inverse()
            def restrict(A):
                result=left*(A*V).matrix_from_rows(rows)
                assert V*result==A*V
                return result
            rep=restrict(compound(g,2*m))
            lie=[restrict(wedge_lie(X,2*m)) for X in basis]
            trace=matrix(QQ,[[(rep*dX*dZ).trace() for dX in lie] for dZ in lie])
            polynomial=matrix(QQ,[[index(m)*(matrix(QQ,4*m,4*m,D.column(j))*Z).trace()
                                   for j in range(len(basis))] for Z in basis])
            assert trace==polynomial
            traces.append(dict(m=m,primitive_dimension=V.ncols(),
                               pairings=len(basis)**2,rank=trace.rank()))
        print('Full ambient derivative passed:',ranks[-1],flush=True)

    # Symbolic original trace construction on the entire m=1 centralizer chart.
    P=PolynomialRing(QQ,names=('p','r','s','b')); FF=P.fraction_field()
    p,r,s,b=map(FF,P.gens())
    A=matrix(FF,[[p,r],[s,(1+r*s)/p]])
    B=matrix(FF,[[b,0],[0,1/b]])
    J2=matrix(QQ,[[0,1],[-1,0]]); J=block_diagonal_matrix(J2,J2)
    V=contraction(J,2).right_kernel().basis_matrix().transpose()
    rows=list(V.transpose().pivots()); left=V.matrix_from_rows(rows).inverse()
    def restrict_symbolic(M):
        answer=left*(M*V).matrix_from_rows(rows)
        assert V*answer==M*V
        return answer
    rep=restrict_symbolic(compound(block_diagonal_matrix(A,B),2))
    expected=block_diagonal_matrix((b+1/b)*(A-A.inverse())/4,
                                  A.trace()*(B-B.inverse())/4)
    assert rep.trace()==1+A.trace()*(b+1/b)
    for X in symplectic_basis(J):
        assert (rep*restrict_symbolic(wedge_lie(X,2))).trace()==2*(expected*X).trace()

    # The exact m=1 singular local model, after an invertible formal coordinate change.
    S=PolynomialRing(QQ,names=('x','y','z','w')); x,y,z,w=S.gens()
    ideal=S.ideal([w*x,w*y,w*z,x*x+y*z-w*w])
    radical=S.ideal([w,x*x+y*z])
    assert ideal.radical()==radical
    assert ideal.reduce(w*w)!=0 and ideal.reduce(w**3)==0
    assert all(ideal.reduce(v*w*w)==0 for v in S.gens())
    assert ideal.dimension()==2
    paths=['SYMPLECTIC_FIBER_COMPONENT.md','tools/verify_symplectic_fiber_component.py',
           'tools/verify_symplectic_centralizer_realization.py',
           'tools/verify_type_a_character_degeneration.py']
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                unimodular_pairings=pairings,cyclotomic_target_checks=collisions,
                full_differential_checks=ranks,actual_primitive_trace_checks=traces,
                rank_one_singular_model=dict(radical_verified=True,nilpotence_order=3,
                    nonzero_socle_verified=True,dimension=2,
                    symbolic_primitive_trace_pairings=10),
                arbitrary_rank_component_proof='Written separately; finite checks are corroboration.',
                global_reduced_germ_conjecture_proved=False,
                target_noncollision_all_rank_proved=False,
                scope='Exact rational primitive traces and polynomial derivatives; cyclotomic target comparisons; exact singular ideal. Shared Sage and wedge helpers, no external audit or novelty claim.',
                source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'symplectic-fiber-component-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
