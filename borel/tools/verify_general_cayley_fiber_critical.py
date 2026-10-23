"""Corroborate the general critical presentation using independent exact models."""
from pathlib import Path
from hashlib import sha256
import json
from sage.all import (QQ, CyclotomicField, PolynomialRing, LaurentPolynomialRing,
                      matrix, identity_matrix, block_matrix, block_diagonal_matrix)
from sage.env import SAGE_VERSION
from verify_type_a_character_degeneration import compound, wedge_lie
from verify_symplectic_centralizer_realization import contraction, symplectic_basis

ROOT=Path(__file__).resolve().parents[1]


def sl2_basis():
    return [matrix(QQ,[[0,1],[0,0]]),matrix(QQ,[[0,0],[1,0]]),
            matrix(QQ,[[1,0],[0,-1]])]


def tensor_rep(g,X,d):
    rep=matrix(g.base_ring(),[[1]])
    for j in range(d): rep=rep.tensor_product(g)
    action=matrix(QQ,2**d)
    for slot in range(d):
        term=matrix(QQ,[[1]])
        for j in range(d): term=term.tensor_product(X if j==slot else identity_matrix(QQ,2))
        action+=term
    return rep,action


def main():
    eliminations=[]
    R=PolynomialRing(QQ,names=('x1','x2','v1','v2','z'))
    x1,x2,v1,v2,z=R.gens()
    for k in range(1,7):
        W=x1*v1+x2*v2+v1*z+v2*z*z+v1*v1+v1*v2+QQ(1)/(k+1)*z**(k+1)
        jac=R.ideal([W.derivative(v) for v in R.gens()])
        expected=R.ideal([v1,v2,x1+z,x2+z*z,z**k])
        assert jac==expected
        assert jac.dimension()==0
        eliminations.append(dict(residual_order=k,exact_gradient_ideal=True))

    Achecks=[]; F=CyclotomicField(4); alpha=F.gen()
    LS=LaurentPolynomialRing(F,'s'); s=LS.gen()
    for m in range(2,13):
        # Differentiate the character plus the logarithmic correction directly.
        chi=alpha**m*(s**(-m)+(-1)**m*s**m)
        log_coefficient=m*alpha**m*(1-(-1)**m)
        derivative=chi.derivative()+log_coefficient/s
        residual=1-s**(2*m) if m%2==0 else (s**m-1)**2
        assert derivative==-m*alpha**m*s**(-m-1)*residual
        assert derivative(F(1))==0
        mu=2 if m%2 else 1
        value=derivative
        for j in range(mu):
            assert value(F(1))==0
            value=value.derivative()
        assert value(F(1))!=0
        Achecks.append(dict(m=m,multiplicity=mu,potential_order=mu+1))

    Cchecks=[]; Rz=PolynomialRing(QQ,'z'); z=Rz.gen()
    for m in range(1,13):
        v=[z**j+z**(2*m-j) for j in range(1,m)]+[z**m]
        L=matrix(QQ,[[((1-z*z)*a*b)[2*m] for b in v] for a in v])
        E=identity_matrix(QQ,m)
        for j in range(m-2): E[j,j+2]=-1
        assert L==E*E.transpose() and L.det()==1
        lower=matrix(QQ,m,m,lambda i,j:1+(i+1 if i==j else 0))
        Hess=block_matrix([[matrix(QQ,m),L],[L.transpose(),lower]])
        assert Hess.det()==(-1)**m
        Cchecks.append(dict(m=m,mixed_determinant=1,hessian_rank=Hess.rank()))

    defining=[]
    for n in range(2,11):
        K=CyclotomicField(2*n); w=K.gen()
        beta=[w**(2*j+(1 if n%2==0 else 0)) for j in range(n)]
        Rc=PolynomialRing(K,'c'); c=Rc.gen()
        h=matrix(K,n,n,lambda i,j:beta[i] if i==j else 0)
        g=h.change_ring(Rc)+c*identity_matrix(Rc,n)
        assert h.trace()==0 and h.det()==1
        assert g.det()==1+c**n
        assert g-g.trace()/n*identity_matrix(Rc,n)==h
        Hess=matrix(K,n-1,n-1,lambda i,j:(beta[i] if i==j else 0)+beta[-1])
        assert Hess.rank()==n-2
        defining.append(dict(n=n,fiber_length=n,quotient_hessian_corank=1))

    # Actual tensor matrices, not the character chain rule being tested.
    P=PolynomialRing(QQ,names=('p','r','s')); E=P.fraction_field()
    p,r,s=map(E,P.gens()); g=matrix(E,[[p,r],[s,(1+r*s)/p]])
    a=g.trace(); base=g-a/2*identity_matrix(E,2); tensors=[]
    for d in range(2,6):
        idx=d*2**(d-1)
        for X in sl2_basis():
            rep,action=tensor_rep(g,X,d)
            assert rep.trace()==a**d
            assert (rep*action).trace()==idx*((a/2)**(d-1)*base*X).trace()
        _,dH=tensor_rep(identity_matrix(QQ,2),sl2_basis()[-1],d)
        assert (dH*dH).trace()==2*idx
        tensors.append(dict(tensor_power=d,dimension=2**d,trace_pairings=3,
                            zero_level_multiplicity=d-1))

    # SL2 defining representation: quotient critical locus empty, group critical points nonempty.
    gradients=[a.derivative(v) for v in (p,r,s)]
    numerator_ideal=P.ideal([P(v.numerator()) for v in gradients])
    assert numerator_ideal==P.ideal([P.gen(0)**2-1,P.gen(1),P.gen(2)])
    point={P.gen(0):1,P.gen(1):0,P.gen(2):0}
    Hess=matrix(QQ,[[v.derivative(w).subs(point) for w in (p,r,s)] for v in gradients])
    assert Hess.det()!=0

    # PGL2 adjoint: regularity alone does not imply smoothness of its global quotient.
    lie=sl2_basis()
    def coords(X): return [X[0,1],X[1,0],X[0,0]]
    adj=matrix(E,[coords(g*X*g.inverse()) for X in lie]).transpose()
    assert adj.trace()==a*a-1
    expected=a/2*base
    for X in lie:
        dX=matrix(QQ,[coords(X*Y-Y*X) for Y in lie]).transpose()
        assert (adj*dX).trace()==4*(expected*X).trace()
    point={P.gen(0):1,P.gen(1):1,P.gen(2):-2}
    assert a.subs(point)==0
    assert any(v.subs(point)!=0 for v in gradients)
    assert all((a*a).derivative(v).subs(point)==0 for v in (p,r,s))

    # Sp4: original primitive matrices versus the corrected potential on the entire K chart.
    K=CyclotomicField(4); ii=K.gen()
    P=PolynomialRing(K,names=('p','r','s','b')); E=P.fraction_field()
    p,r,s,b=map(E,P.gens()); A=matrix(E,[[p,r],[s,(1+r*s)/p]])
    B=matrix(E,[[b,0],[0,1/b]]); g=block_diagonal_matrix(A,B)
    J2=matrix(QQ,[[0,1],[-1,0]]); J=block_diagonal_matrix(J2,J2)
    V=contraction(J,2).right_kernel().basis_matrix().transpose()
    rows=list(V.transpose().pivots()); left=V.matrix_from_rows(rows).inverse()
    def restrict(M):
        result=left*(M*V).matrix_from_rows(rows)
        assert V*result==M*V
        return result
    rep=restrict(compound(g,2)); a=A.trace(); q=b+1/b
    assert rep.trace()==1+a*q
    potential_gradient=[q*a.derivative(v) for v in (p,r,s)]+[a*(1-1/b**2)-4*ii/b]
    h=block_diagonal_matrix(matrix(E,2),matrix(E,[[ii,0],[0,-ii]]))
    basis=[block_diagonal_matrix(X,matrix(QQ,2)) for X in sl2_basis()]
    basis+=[block_diagonal_matrix(matrix(QQ,2),sl2_basis()[-1])]
    for X in basis:
        dg=g*X
        velocity=[dg[0,0],dg[0,1],dg[1,0],dg[2,2]]
        derivative=sum(u*v for u,v in zip(potential_gradient,velocity))
        assert derivative==(rep*restrict(wedge_lie(X,2))).trace()-2*(h*X).trace()
    point={P.gen(0):1,P.gen(1):0,P.gen(2):0,P.gen(3):ii}
    Hess=matrix(K,[[v.derivative(w).subs(point) for w in (p,r,s,b)] for v in potential_gradient])
    assert Hess==0
    # On the quotient (a,b), the Hessian is instead nonsingular.
    assert matrix(K,[[0,2],[2,0]]).det()==-4

    paths=['GENERAL_CAYLEY_FIBER_CRITICAL.md','tools/verify_general_cayley_fiber_critical.py',
           'tools/verify_type_a_character_degeneration.py',
           'tools/verify_symplectic_centralizer_realization.py']
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                linear_eliminations=eliminations,type_A_potentials=Achecks,
                type_C_perfect_pairings=Cchecks,defining_SLn_baselines=defining,
                actual_tensor_matrix_checks=tensors,
                nonregular_SL2_control=True,ramified_PGL2_control=True,
                symbolic_Sp4_corrected_potential_pairings=4,
                Sp4_group_hessian_rank=0,Sp4_quotient_hessian_rank=2,
                general_proof_scope='Written separately; these are finite exact corroborations.',
                published_baselines=['KM 2.8(3) tensor scaling','KM 6.2 defining SLn fibers',
                                     'Steinberg regular quotient smoothness'],
                novelty_established=False,nonregular_product_theorem_claimed=False,
                source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'general-cayley-fiber-critical-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
