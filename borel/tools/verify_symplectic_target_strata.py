"""Exact corroboration of target-stratum theorems for last type C modules."""
from hashlib import sha256
from itertools import permutations, product
from math import factorial, prod, comb
from pathlib import Path
import json
from sage.all import (QQ, PolynomialRing, matrix, vector, identity_matrix,
                      diagonal_matrix, block_diagonal_matrix, companion_matrix,
                      CyclotomicField, Partitions, WeylCharacterRing)
from sage.env import SAGE_VERSION
from verify_symplectic_centralizer_realization import contraction, symplectic_basis
from verify_type_a_character_degeneration import compound, wedge_lie

ROOT=Path(__file__).resolve().parents[1]


def elementary(values,k):
    if k<0: return QQ(0)
    a=[QQ(1)]+[QQ(0)]*k
    for x in values:
        for j in range(k,0,-1): a[j]+=x*a[j-1]
    return a[k]


def reciprocal(R,z,rank,coeff):
    if rank==0: return R(1)
    return 1+z**(2*rank)+sum(coeff[j-1]*(z**j+z**(2*rank-j)) for j in range(1,rank))+coeff[-1]*z**rank


def first_block_checks():
    rows=[]; odd=[]; even=[]
    for n in range(1,17):
        for p in range(1,n+1):
            q=n-p
            R=PolynomialRing(QQ,names=tuple('b'+str(j) for j in range(1,max(q,1)+1)))
            Z=PolynomialRing(R,'z'); z=Z.gen()
            Q=reciprocal(Z,z,q,list(R.gens())[:q]); C=(1-z*z)*Q
            v=[z**j+z**(2*p-j) for j in range(1,p)]+[z**p]
            direct=[(C*f)[n] for f in v]
            def c(j): return C[j] if 0<=j<=C.degree() else R(0)
            triangular=[c(q-j) for j in range(p)]
            assert R.ideal(direct)==R.ideal(triangular)
            if p>q:
                assert R.ideal(direct)==R.ideal([1])
            else:
                M=matrix(QQ,[[f.derivative(b) for b in R.gens()] for f in direct])
                assert M.rank()==p
            rows.append(dict(n=n,p=p,q=q,empty=(p>q)))
    R=PolynomialRing(QQ,names=('t','z')); t,z=R.gens()
    for m in range(1,17):
        Q=sum(z**(2*j) for j in range(m+2))+t*sum(z**(2*j+1) for j in range(m+1))
        C=(1-z*z)*Q
        assert C==1+t*z-t*z**(2*m+3)-z**(2*m+4)
        # Extract in z by substitution into a univariate ring over QQ[t].
        T=PolynomialRing(QQ,'t'); U=PolynomialRing(T,'z'); zz=U.gen(); tt=T.gen()
        Qu=sum(zz**(2*j) for j in range(m+2))+tt*sum(zz**(2*j+1) for j in range(m+1))
        assert ((1-zz*zz)*Qu*(1+zz**(2*m)))[2*m+1]==tt
        odd.append(dict(m=m,residual_character='t'))
        Q0=sum(zz**(2*j) for j in range(m+1))
        assert Q0.gcd(Q0.derivative()).degree()==0 and Q0(1)!=0 and Q0(-1)!=0
        vv=[zz**j+zz**(2*m-j) for j in range(1,m)]+[zz**m]
        L=matrix(QQ,[[((1-zz*zz)*a*b)[2*m] for b in vv] for a in vv])
        E=identity_matrix(QQ,m)
        for j in range(m-2): E[j,j+2]=-1
        assert L==E*E.transpose() and L.det()==1
        even.append(dict(m=m,mixed_determinant=1,spectrum_squarefree=True))
    return rows,odd,even


def signed_assignments(m,partition):
    result=set()
    for perm in permutations(range(1,m+1)):
        for signs in product((-1,1),repeat=m):
            values=[a*b for a,b in zip(perm,signs)]; start=0; blocks=[]
            for size in partition:
                blocks.append(tuple(sorted(values[start:start+size]))); start+=size
            result.add(tuple(blocks))
    return sorted(result)


def counts():
    rows=[]
    for m in range(1,6):
        for part in Partitions(m):
            actual=len(signed_assignments(m,part))
            expected=2**m*factorial(m)//prod(factorial(int(r)) for r in part)
            assert actual==expected
            rows.append(dict(m=m,multiplicities=list(map(int,part)),components=actual))
    return rows


def standard_form(m):
    J=matrix(QQ,2*m)
    for j in range(m): J[j,m+j]=1; J[m+j,j]=-1
    return J


def symplectic_companion(poly):
    """Construct a cyclic symplectic matrix, checking both form and polynomial.

    Solve for an invariant alternating form, then change to the fixed
    standard form. This does not invoke the character-gradient formula.
    """
    K=poly.base_ring(); A=companion_matrix(poly,format='right'); d=A.nrows(); m=d//2
    skew=[]
    for i in range(d):
        for j in range(i+1,d):
            S=matrix(K,d); S[i,j]=1; S[j,i]=-1; skew.append(S)
    equations=matrix(K,[(A.transpose()*S*A-S).list() for S in skew]).transpose()
    forms=[sum((x*S for x,S in zip(v,skew)),matrix(K,d)) for v in equations.right_kernel().basis()]
    candidates=forms+[sum((K(c)**j*S for j,S in enumerate(forms)),matrix(K,d)) for c in range(1,12)]
    J=next(S for S in candidates if S.det()!=0)
    remaining=list(identity_matrix(K,d).columns()); vs=[]; ws=[]
    pairing=lambda v,w:(v*J).dot_product(w)
    while remaining:
        v=remaining.pop(0); j=next(j for j,w in enumerate(remaining) if pairing(v,w)!=0)
        w=remaining.pop(j); w=w/pairing(v,w); vs.append(v); ws.append(w)
        remaining=[x-pairing(v,x)*w+pairing(w,x)*v for x in remaining]
    S=matrix(K,vs+ws).transpose()
    assert S.transpose()*J*S==standard_form(m)
    result=S.inverse()*A*S
    assert result.transpose()*standard_form(m)*result==standard_form(m)
    assert result.charpoly().list()==poly.list()
    return result


def actual_primitive_checks():
    records=[]
    for m in (1,2):
        K=CyclotomicField(2*m+2); w=K.gen(); roots=[w**j for j in range(1,m+1)]
        J=block_diagonal_matrix(standard_form(m),standard_form(m)); n=4*m
        C=contraction(J,2*m); V=C.right_kernel().basis_matrix().transpose()
        pivots=list(V.transpose().pivots()); left=V.matrix_from_rows(pivots).inverse()
        def restrict(A):
            AV=A*V; result=left*AV.matrix_from_rows(pivots)
            assert V*result==AV
            return result
        basis=symplectic_basis(J); actions=[restrict(wedge_lie(X,2*m)) for X in basis]
        # Sparse trace contractions avoid multiplying matrices only to take a trace.
        trace_terms=[[(j,i,act[i,j]) for (i,j) in act.nonzero_positions()] for act in actions]
        I=QQ(comb(4*m-2,2*m-1)-(comb(4*m-2,2*m-3) if 2*m>=3 else 0))
        for X,act in zip(basis,actions):
            assert (act*act).trace()==I*(X*X).trace()
        Z=PolynomialRing(K,'z'); z=Z.gen()
        vv=[z**j+z**(2*m-j) for j in range(1,m)]+[z**m]
        L=matrix(K,[[((1-z*z)*a*b)[2*m] for b in vv] for a in vv])
        d0=vector(K,[((1-z*z)*(1+z**(2*m))*v)[2*m] for v in vv])
        for part in Partitions(m):
            targets=[K(2*j+1) for j,size in enumerate(part) for unused in range(size)]
            for assignment in signed_assignments(m,part):
                labels=[x for block in assignment for x in block]
                beta=[roots[abs(x)-1]**(1 if x>0 else -1) for x in labels]
                alphabet=beta+[1/b for b in beta]
                jac=matrix(K,m,m,lambda i,j:
                    beta[i]*elementary(alphabet[:i]+alphabet[i+1:],j)
                    -1/beta[i]*elementary(alphabet[:m+i]+alphabet[m+i+1:],j))
                assert jac.det()!=0
                a=L.transpose().solve_right(jac.solve_right(vector(K,[2*I*x for x in targets]))-d0)
                coefficients=[K(1)]+list(a)+list(reversed(list(a)[:-1]))+[K(1)]
                poly=sum((-1)**j*coefficients[j]*z**(2*m-j) for j in range(2*m+1))
                A=symplectic_companion(poly)
                B=diagonal_matrix(K,alphabet); g=block_diagonal_matrix(A,B)
                h=diagonal_matrix(K,[K(0)]*(2*m)+targets+[-x for x in targets])
                assert g.transpose()*J*g==J and g.det()==1
                rep=restrict(compound(g,2*m))
                for X,terms in zip(basis,trace_terms):
                    observed=sum(rep[i,j]*c for i,j,c in terms)
                    assert observed==I*(h*X).trace(),(m,part,assignment)
                # A is cyclic; no first-block eigenvalue can equal an external one.
                assert all(poly(b)!=0 for b in alphabet)
                records.append(dict(m=m,multiplicities=list(map(int,part)),assignment=list(labels),
                                    primitive_dimension=V.ncols(),trace_pairings=len(basis),
                                    cyclic_first_block=True,no_cross_block_spectrum=True))
            print('Primitive matrices passed:',m,list(part),flush=True)
    return records


def weyl_characters():
    rows=[]
    for n in range(2,6):
        W=WeylCharacterRing('C'+str(n),style='coroots'); weights=W(*([0]*(n-1)+[1])).weight_multiplicities()
        for seed in (2,5):
            b=[QQ(seed+j) for j in range(n)]; alphabet=b+[1/x for x in b]
            actual=sum(mult*prod(b[j]**int(w[j]) for j in range(n)) for w,mult in weights.items())
            assert actual==elementary(alphabet,n)-elementary(alphabet,n-2)
        rows.append(dict(n=n,dimension=int(sum(weights.values())),points=2))
    return rows


def odd_primitive_checks():
    # The primitive Sp6 module has dimension 14. Obtain the target from
    # original matrices, then independently check the residual derivative.
    m=1; K=CyclotomicField(6); w=K.gen(); beta=[w,w*w]
    alphabet=beta+[1/b for b in beta]
    J=block_diagonal_matrix(standard_form(1),standard_form(2))
    C=contraction(J,3); V=C.right_kernel().basis_matrix().transpose()
    pivots=list(V.transpose().pivots()); left=V.matrix_from_rows(pivots).inverse()
    def restrict(M):
        MV=M*V; result=left*MV.matrix_from_rows(pivots)
        assert V*result==MV
        return result
    basis=symplectic_basis(J); actions=[restrict(wedge_lie(X,3)) for X in basis]
    jac=matrix(K,2,2,lambda i,j:beta[i]*elementary(alphabet[:i]+alphabet[i+1:],j)
               -1/beta[i]*elementary(alphabet[:2+i]+alphabet[2+i+1:],j))
    log_derivative=jac.transpose().solve_right(vector(K,[1,0]))
    # Independently differentiate Q_t(-beta)=0 at t=0 to check signs.
    assert list(log_derivative)==[(1+b*b)/(2*b+4*b**3) for b in beta]
    records=[]; I=QQ(5)
    for a in (3,4):
        A=matrix(K,[[0,-1],[1,a]])
        g=block_diagonal_matrix(A,diagonal_matrix(K,alphabet)); rep=restrict(compound(g,3))
        targets=[]
        for j in range(2):
            H=matrix(QQ,6); H[2+j,2+j]=1; H[4+j,4+j]=-1
            dH=restrict(wedge_lie(H,3))
            assert (dH*dH).trace()==I*(H*H).trace()
            targets.append((rep*dH).trace()/(2*I))
        assert all(x!=0 for x in targets) and targets[0]!=targets[1] and targets[0]!=-targets[1]
        h=diagonal_matrix(K,[0,0]+targets+[-x for x in targets])
        for X,act in zip(basis,actions):
            assert (rep*act).trace()==I*(h*X).trace()
        assert 1-2*I*vector(K,targets).dot_product(log_derivative)==0
        records.append(dict(first_block_trace=a,primitive_dimension=V.ncols(),trace_pairings=len(basis),
                            residual_derivative_zero=True,source_root_sign_checked=True))
    return records


def main():
    blocks,odd,even=first_block_checks()
    original=actual_primitive_checks(); odd_original=odd_primitive_checks()
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                first_block_ideal_checks=blocks,odd_residual_checks=odd,even_pairing_checks=even,
                signed_root_counts=counts(),primitive_matrix_checks=original,weyl_characters=weyl_characters(),
                odd_primitive_matrix_checks=odd_original,
                total_primitive_trace_pairings=sum(x['trace_pairings'] for x in original+odd_original),
                all_rank_claims_have_written_proofs=True,
                odd_rank_multiplicities_classified=False,nonregular_components_classified=False)
    sources=['SYMPLECTIC_TARGET_STRATA.md','tools/verify_symplectic_target_strata.py',
             'tools/verify_symplectic_centralizer_realization.py','tools/verify_type_a_character_degeneration.py']
    report['source_sha256']={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in sources}
    (ROOT/'symplectic-target-strata-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(result='passed',first_block_ideals=len(blocks),even_ranks=len(even),odd_ranks=len(odd),
                         primitive_cases=len(original),trace_pairings=report['total_primitive_trace_pairings'],
                         signed_counts=len(report['signed_root_counts'])),indent=2))


if __name__=='__main__': main()
