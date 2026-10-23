"""Check type C realization using actual contraction-kernel matrices."""
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from sage.all import (QQ, PolynomialRing, matrix, identity_matrix,
                      block_diagonal_matrix, binomial, vector)
from sage.env import SAGE_VERSION
from verify_type_a_character_degeneration import wedge_lie, compound

ROOT=Path(__file__).resolve().parents[1]


def coefficient_model(m):
    R=PolynomialRing(QQ,'z'); z=R.gen()
    index=binomial(4*m-2,2*m-1)-binomial(4*m-2,2*m-3)
    base=sum(z**(2*j) for j in range(m+1))
    variations=[z**k+z**(2*m-k) for k in range(1,m)]+[z**m]
    functions=[(1-z)**j*(1+z)**(2*m-j) for j in range(1,2*m,2)]
    raw=matrix(QQ,[[(f*v)[2*m-1] for v in variations] for f in functions])
    M=(QQ(2)/index)*raw
    # Coerce before division: Sage's binomial may return a Python int.
    assert M.det()==raw.det()*(QQ(2)/index)**m
    assert all((base*f)[2*m-1]==0 for f in functions)
    assert M.det()!=0
    return R,index,base,variations,M


def spectrum(m,p):
    R,index,base,variations,M=coefficient_model(m)
    z=R.gen(); solution=M.solve_right(vector(QQ,p))
    change=sum(a*v for a,v in zip(solution,variations))
    c=QQ(1)
    while True:
        Q=base+c*change
        if Q(1) and Q(-1) and Q.gcd(Q.derivative()).degree()==0: break
        c+=1
    d=2*m
    # Q(z)=det(I+z B), so char_B(x)=x^d Q(-1/x).
    F=sum((-1)**h*Q[h]*z**(d-h) for h in range(d+1))
    B=matrix(QQ,d)
    for j in range(d-1): B[j+1,j]=1
    for j in range(d): B[j,d-1]=-F[j]
    assert B.charpoly().list()==F.list()
    return c,Q,B


def invariant_form(B):
    d=B.nrows(); skew=[]
    for j in range(d):
        for k in range(j+1,d):
            S=matrix(QQ,d); S[j,k]=1; S[k,j]=-1; skew.append(S)
    equations=matrix(QQ,[(B.transpose()*S*B-S).list() for S in skew]).transpose()
    forms=[sum((a*S for a,S in zip(v,skew)),matrix(QQ,d))
           for v in equations.right_kernel().basis()]
    candidates=list(forms)+[sum((QQ(k)**j*S for j,S in enumerate(forms)),matrix(QQ,d))
                            for k in range(1,21)]
    for J in candidates:
        if J.det():
            assert J.transpose()==-J and B.transpose()*J*B==J
            return J
    raise AssertionError('No nondegenerate invariant form among the exact candidates')


def contraction(J,p):
    n=J.nrows(); source=list(combinations(range(n),p))
    target=list(combinations(range(n),p-2)); rows={v:i for i,v in enumerate(target)}
    C=matrix(QQ,len(target),len(source),sparse=True)
    for col,subset in enumerate(source):
        for a in range(p):
            for b in range(a+1,p):
                remainder=tuple(subset[k] for k in range(p) if k not in (a,b))
                C[rows[remainder],col]+=(-1)**(a+b-1)*J[subset[a],subset[b]]
    assert C.rank()==len(target)
    return C


def symplectic_basis(J):
    n=J.nrows(); inverse=J.inverse(); result=[]
    for a in range(n):
        for b in range(a,n):
            S=matrix(QQ,n); S[a,b]=1; S[b,a]=1
            X=inverse*S
            assert X.transpose()*J+J*X==0
            result.append(X)
    return result


def main():
    transforms=[]
    for m in range(1,11):
        R,index,base,variations,M=coefficient_model(m)
        transforms.append(dict(m=m,index=int(index),determinant=str(M.det())))
    cases=[]; pairings=0
    for m in (1,2):
        fixtures=[('power_%s'%r,[QQ(int(j==r)) for j in range(1,2*m,2)])
                  for r in range(1,2*m,2)]
        fixtures += [('mixed',[QQ(j+1) for j in range(m)]),('zero',[QQ(0)]*m)]
        for label,p in fixtures:
            c,Q,B=spectrum(m,p)
            J0=matrix(QQ,2*m)
            for j in range(m): J0[j,m+j]=1; J0[m+j,j]=-1
            J=block_diagonal_matrix(J0,invariant_form(B))
            t=block_diagonal_matrix(identity_matrix(QQ,2*m),B)
            Y=matrix(QQ,2*m)
            for j in range(m-1): Y[j,j+1]=1; Y[m+j+1,m+j]=-1
            Y[m-1,2*m-1]=1
            assert Y.transpose()*J0+J0*Y==0 and Y**(2*m)==0 and Y**(2*m-1)!=0
            u=(identity_matrix(QQ,2*m)+Y)*(identity_matrix(QQ,2*m)-Y).inverse()
            g=t*block_diagonal_matrix(u,identity_matrix(QQ,2*m))
            assert t.transpose()*J*t==J and g.transpose()*J*g==J
            assert t.det()==g.det()==1
            C=contraction(J,2*m)
            V=C.right_kernel().basis_matrix().transpose()
            rows=list(V.transpose().pivots())
            left=V.matrix_from_rows(rows).inverse()
            def restrict(A):
                AV=A*V
                result=left*AV.matrix_from_rows(rows)
                assert V*result==AV
                return result
            group_difference=restrict(compound(g,2*m))-restrict(compound(t,2*m))
            index=binomial(4*m-2,2*m-1)-binomial(4*m-2,2*m-3)
            H=matrix(QQ,4*m); H[0,0]=1; H[m,m]=-1
            dH=restrict(wedge_lie(H,2*m))
            assert (dH*dH).trace()==index*(H*H).trace()
            target=c*sum((p[j]*Y**(2*j+1) for j in range(m)),matrix(QQ,2*m))
            expected=block_diagonal_matrix(target,matrix(QQ,2*m))
            basis=symplectic_basis(J)
            for X in basis:
                dX=restrict(wedge_lie(X,2*m))
                assert (group_difference*dX).trace()==index*(expected*X).trace(),(m,label)
                pairings+=1
            commutators=matrix(QQ,[(t*X-X*t).list() for X in basis]).transpose()
            centralizer_dim=len(basis)-commutators.rank()
            assert centralizer_dim==m*(2*m+1)+m
            cases.append(dict(m=m,label=label,primitive_dimension=V.ncols(),
                              scalar=str(c),spectral_polynomial=str(Q),
                              ambient_tangent_pairings=len(basis),centralizer_dimension=centralizer_dim))
            print('Actual primitive matrices passed:',cases[-1],flush=True)

    paths=['SYMPLECTIC_CENTRALIZER_REALIZATION.md',
           'tools/verify_symplectic_centralizer_realization.py',
           'tools/verify_type_a_character_degeneration.py']
    report=dict(result='passed',date='2026-10-23',sage_version=SAGE_VERSION,
                coefficient_isomorphisms=transforms,actual_primitive_cases=cases,
                actual_primitive_tangent_pairings=pairings,
                contraction_kernels_and_invariance_checked=True,
                trace_indices_and_exact_centralizers_checked=True,
                used_runtime_character_or_borel_api=False,
                full_fiber_component_theorem=False,
                scope='Finite exact primitive-kernel trace corroboration of the separate written all-rank proof; shared Sage arithmetic and signed-wedge helpers; no external audit or novelty certificate.',
                source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'symplectic-centralizer-realization-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('source_sha256','actual_primitive_cases')},indent=2))


if __name__=='__main__':main()
