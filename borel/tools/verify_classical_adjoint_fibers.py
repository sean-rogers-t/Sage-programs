"""Bounded, Sage-free audits of the all-rank classical adjoint proofs.

No rank sweep. Matrix representations are built on tensor-square bases,
independently of the character API and of the quotient equations.
"""
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import json

import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def zero(value):
    assert s.cancel(value) == 0, value


def symbolic_checks():
    tau, L, u, v, a = s.symbols('tau L u v a')
    p = u**4-tau*u**3+L*u**2+tau*u-1
    # Independently multiply the selected and complementary quadratics.
    factored = s.expand((u**2-a*u+v)*(u**2-(tau-a)*u-1/v))
    cubic_relation = a*(v**2+1)-tau*v*(v+1)
    zero((factored.coeff(u, 1)-tau)*v-cubic_relation)
    chosen_sum = tau*v*(v+1)/(v**2+1)
    zero(chosen_sum*(1+1/v)-tau*(v+1)**2/(v**2+1))
    assert s.diff(tau*(v+1)**2/(v**2+1), v) != 0
    zero(s.expand(-u**4*p.subs(u, 1/u))-p.subs(L, -L))
    disc = s.discriminant(p, u)
    assert s.Poly(disc, L).LC() != 0
    # Reconstruct the quartic from both normalized eigenvalue formulas.
    T, d, lam = s.symbols('T d lam', nonzero=True)
    orth = (T*(u-1/u)-(u**2-1/u**2))/(2*d)
    symp = (T*(u-1/u)+(u**2-1/u**2))/(2*d)
    zero(-2*d*u**2*(orth-lam)-p.subs({tau:T,L:2*d*lam}))
    zero(2*d*u**2*(symp-lam)-p.subs({tau:-T,L:-2*d*lam}))
    zero(p.subs(L, 0)-(u**2-1)*(u**2-tau*u+1))

    x, y, z, q, c = s.symbols('x y z q c')
    P = x*x/2+z
    A = x
    phi = P+A*y*y
    psi = P+A*q
    zero(s.diff(phi,y)-2*y*s.diff(psi,q).subs(q,y*y))
    zero(s.diff(phi,x)-s.diff(psi,x).subs(q,y*y))
    zero(z*s.diff(phi,z)-(z*s.diff(psi,z)).subs(q,y*y))
    actual = s.groebner([s.diff(phi,x),s.diff(phi,y),z*s.diff(phi,z)-c],x,y,z)
    expected = s.groebner([x+y*y,y**3,z-c],x,y,z)
    assert actual == expected
    assert actual.reduce(y**3)[1] == 0 and actual.reduce(y*y)[1] != 0
    # This prevents an incorrect assertion that quadratic descent alone suffices.

    # Minimal B1 and D2 residual models, with the GL2 central character retained.
    w = s.symbols('w')
    H = (z+w)*x
    P0 = x*x+z*z+w*w
    # Laurent differentiation before imposing z*w=1.
    D = lambda f: z*s.diff(f,z)-w*s.diff(f,w)
    Pb = P0-H-1
    Ab = H+1
    field = s.QQ.frac_field(c)
    bad_b = s.groebner([s.diff(Pb,x),Ab,D(Pb)-c,z*w-1],x,z,w,domain=field)
    assert list(bad_b) == [1]
    for sign in (-1,1):
        bad_d = s.groebner([s.diff(P0,x),H-2*sign,D(P0)-c,z*w-1],x,z,w,domain=field)
        assert list(bad_d) == [1]

    # Weight construction independent of the fundamental-character expression.
    r, b = s.symbols('r b', nonzero=True)
    yr, yb = r+1/r, b+1/b
    residual_b = [r*r,s.Integer(1),r**-2]
    residual_d = [r*b,r/b,b/r,1/(r*b)]
    wedge_char = lambda eig: sum(i*j for i,j in combinations(eig,2))
    zero(sum(residual_b)-(yr*yr-1))
    zero(wedge_char(residual_b)-(yr*yr-1))
    zero(sum(residual_d)-yr*yb)
    zero(wedge_char(residual_d)-(yr*yr+yb*yb-2))
    e = s.symbols('e', nonzero=True)
    h_eig = [z*e,z/e,1/(z*e),e/z]
    x_eig = e+1/e
    h = (z+1/z)*x_eig
    pp = x_eig*x_eig+z*z+z**-2
    zero(wedge_char(h_eig+residual_b)-(pp-h-1+(h+1)*yr*yr))
    zero(wedge_char(h_eig+residual_d)-(pp+yr*yr+yb*yb-2+h*yr*yb))
    zero((yr*yr+yb*yb+h*yr*yb)-
         ((2+h)*(yr+yb)**2+(2-h)*(yr-yb)**2)/4)

    t = s.symbols('t', integer=True, positive=True)
    positive_rows = [('B',0,s.Integer(1),1/(1-2*t)),
                     ('B',1,s.Integer(-1),1/(2*t-1)),
                     ('B',2,T-1,1/(2*t)),
                     ('C',0,s.Integer(0),s.Integer(0)),
                     ('C',1,-T,s.Integer(0)),
                     ('D',0,s.Integer(0),s.Integer(0)),
                     ('D',2,s.Integer(0),s.Integer(0)),
                     ('D',3,T,s.Integer(0))]
    rows=[]
    for family,m,R,tvalue in positive_rows:
        block = (-2 if family=='C' else 2)*t*T
        zero((T-block-R).subs(T,tvalue))
        rows.append(dict(family=family,residual_rank=m,T=str(tvalue),
                         symbolic_all_t_trace_identity=True))
    return dict(vieta_and_reciprocal_identities=True,
                quartic_discriminant=s.factor(disc).__str__(),
                square_chain_rule_checked=True,
                persistent_cubic_control_basis=[str(f) for f in actual],
                cubic_control_is_representation=False,
                generic_bad_ideal_controls=['B1 with GL2','D2 plus with GL2','D2 minus with GL2'],
                residual_character_weight_identities=6,
                all_rank_four_block_trace_identities=rows)


def deformation_checks():
    """Symbolic mechanism checks, not a sweep of ranks or block partitions."""
    T,L,u=s.symbols('T L u')
    for sign in (1,-1):
        p=u**4-sign*T*u**3+L*u**2+sign*T*u-1
        derivative=-(1-u**-2)*s.diff(p,T)/s.diff(p,u)
        numerator=s.fraction(s.cancel((derivative-sign).subs(L,0)))[0]
        assert s.rem(numerator,u*u-sign*T*u+1,u)==0
        for root in (1,-1):
            zero(derivative.subs({L:0,u:root}))
    Q=s.symbols('Q',integer=True,nonnegative=True)
    A=s.symbols('a_C',integer=True,nonnegative=True)
    # Each row is an entire symbolic branch, with restrictions in the note.
    branches=[('C all blocks',(1+Q+A)*T,s.Integer(0)),
              ('B0 balanced',(1-Q)*T-1,1/(1-Q)),
              ('B1 balanced',(1-Q)*T+1,-1/(1-Q)),
              ('B2 balanced',-Q*T+1,1/Q),
              ('D0/D2 balanced',(1-Q)*T,s.Integer(0)),
              ('D3 balanced',-Q*T,s.Integer(0)),
              ('B0 single GL1',T-3,s.Integer(3)),
              ('B0 single GL3',-T-3,s.Integer(-3)),
              ('B1 single GL1',T-1,s.Integer(1)),
              ('B1 single GL3',-T-1,s.Integer(-1)),
              ('D2 single GL3',-2*T,s.Integer(0))]
    rows=[]
    for label,F,point in branches:
        zero(F.subs(T,point))
        derivative=s.diff(F,T)
        assert derivative!=0
        rows.append(dict(branch=label,T0=str(point),derivative=str(derivative)))
    # All local block selections used by the proof, not ambient-rank examples.
    choices=[(1,0,0,1),(2,1,1,0),(2,0,0,2),(3,1,1,1),(4,1,1,2),
             (1,1,0,0),(3,1,0,2)]
    assert all(b==plus+minus+q and plus in (0,1) and minus in (0,1)
               and q in (0,1,2) for b,plus,minus,q in choices)
    assert all(plus==minus for b,plus,minus,q in choices[:5])
    return dict(implicit_root_derivatives_checked_for_both_signs=True,
                symbolic_deformation_branches=rows,
                local_root_selection_rules=[list(row) for row in choices],
                rank_sweep=False,
                excluded_D_low_rank_patterns=[dict(residual_rank=0,blocks=[1]),
                                             dict(residual_rank=0,blocks=[3]),
                                             dict(residual_rank=2,blocks=[1])])


def tensor_square(g, symmetric=False, differential=False):
    """Matrices on e_i wedge e_j or the monomials e_i e_j."""
    n=g.rows
    pairs=list(combinations_with_replacement(range(n),2) if symmetric
               else combinations(range(n),2))
    def product(i,j,k,l):
        if differential:
            return g[i,k]*int(j==l)+int(i==k)*g[j,l]
        return g[i,k]*g[j,l]
    def entry(row,col):
        i,j=pairs[row]; k,l=pairs[col]
        if symmetric and i==j:
            return product(i,j,k,l)
        return product(i,j,k,l)+(1 if symmetric else -1)*product(j,i,k,l)
    return s.MutableSparseMatrix(len(pairs),len(pairs),entry)


def tr_product(A,B):
    return sum(value*B[j,i] for (i,j),value in A.todok().items())


def so_case(n):
    basis=[]
    for i,j in combinations(range(n),2):
        X=s.zeros(n);X[i,j]=1;X[j,i]=-1;basis.append(X)
    g=s.eye(n)
    for j in range(n//2):
        t=s.Integer(j+2); cc=(1-t*t)/(1+t*t);ss=2*t/(1+t*t)
        g[2*j:2*j+2,2*j:2*j+2]=s.Matrix([[cc,-ss],[ss,cc]])
    assert g.T*g==s.eye(n) and g.det()==1
    return basis,g,s.eye(n),False


def sp_cases():
    n=2; basis=[]
    for i in range(n):
        for j in range(n):
            X=s.zeros(4);X[i,j]=1;X[n+j,n+i]=-1;basis.append(X)
    for i,j in combinations_with_replacement(range(n),2):
        for upper in (True,False):
            X=s.zeros(4)
            if upper:
                X[i,n+j]=1;X[j,n+i]=1
            else:
                X[n+i,j]=1;X[n+j,i]=1
            basis.append(X)
    J=s.zeros(4);J[:2,2:]=s.eye(2);J[2:,:2]=-s.eye(2)
    unipot=s.eye(4);unipot[:2,2:]=s.Matrix([[1,2],[2,3]])
    g=s.diag(2,3,s.Rational(1,2),s.Rational(1,3))*unipot
    N=s.zeros(4);N[0,1]=1;N[3,2]=-1;N[1,3]=1
    assert N**4==s.zeros(4) and N**3!=s.zeros(4)
    regular_unipot=s.eye(4)+N+N**2/2+N**3/6
    return [(basis,g,J,True),(basis,regular_unipot,J,True)]


def matrix_checks():
    rows=[]
    for label,case in [('SO7 exterior square',so_case(7)),
                       ('SO8 exterior square',so_case(8))]+[
                           ('Sp4 symmetric square '+str(i),case) for i,case in enumerate(sp_cases())]:
        basis,g,J,symmetric=case; n=g.rows; dim=len(basis)
        assert g.T*J*g==J
        assert all(X.T*J+J*X==s.zeros(n) for X in basis)
        R=tensor_square(g,symmetric)
        dr=[tensor_square(X,symmetric,True) for X in basis]
        assert R.rows==dim
        d=n+(2 if symmetric else -2);sgn=1 if symmetric else -1
        inv=g.inv();T=s.trace(g)
        theta=(T*(g-inv)+sgn*(g*g-inv*inv))/(2*d)
        zero(s.trace(R)-(T*T+sgn*s.trace(g*g))/2)
        for i,X in enumerate(basis):
            zero(tr_product(R,dr[i])-d*s.trace(theta*X))
            for j,Y in enumerate(basis):
                zero(tr_product(dr[i],dr[j])-d*s.trace(X*Y))
        for j,Y in enumerate(basis):
            derivative=(s.trace(g*Y)*(g-inv)+T*(g*Y+Y*inv)
                        +sgn*(g*Y*g+g*g*Y+Y*inv*inv+inv*Y*inv))/(2*d)
            direct=R*dr[j]
            for i,X in enumerate(basis):
                zero(tr_product(direct,dr[i])-d*s.trace(derivative*X))
        rows.append(dict(case=label,defining_dimension=n,adjoint_dimension=dim,
                         trace_form_pairings=dim*dim,Springer_trace_pairings=dim,
                         full_differential_pairings=dim*dim,
                         regular_unipotent_control='square 1' in label))
        print('Independent matrix check:',label,flush=True)
    return rows


def main():
    symbolic=symbolic_checks();deformation=deformation_checks()
    print('Symbolic proof and deformation identities passed',flush=True)
    matrices=matrix_checks()
    paths=['CLASSICAL_ADJOINT_FIBERS.md','tools/verify_classical_adjoint_fibers.py',
           'QUADRATIC_BRANCHING_FIBERS.md','GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report=dict(result='passed',date='2026-10-23',uses_Sage=False,sympy_version=s.__version__,
                symbolic=symbolic,deformation=deformation,matrix_checks=matrices,
                classical_all_rank_smoothness_written_proof=True,
                generic_regular_dimension='dim K-rank G',
                mixed_block_nonemptiness_classified=True,
                nonemptiness_criterion='all GL blocks <=4; residual m<=2 in B, <=1 in C, <=3 in D (n>=4)',
                public_API_changed=False,literature_priority_claim=False,
                scope='Finite checks audit identities and normalization; all-rank proofs are in the note. No new rank sweep or external proof referee.',
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in paths})
    (ROOT/'classical-adjoint-fibers-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',matrix_cases=len(matrices),
                         trace_form_pairings=sum(r['trace_form_pairings'] for r in matrices),
                         Springer_pairings=sum(r['Springer_trace_pairings'] for r in matrices),
                         differential_pairings=sum(r['full_differential_pairings'] for r in matrices)),indent=2))


if __name__=='__main__':
    main()
