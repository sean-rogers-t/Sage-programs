"""Exact corroboration of sharp symplectic harmonic localization.

Run with Sage Python. Direct matrix differentiation, independent exterior
trace projection, genuine low-degree harmonic kernels, and conjugacy-fiber
pairs corroborate the uniform proof. No floating point or runtime API edits.
"""
from hashlib import sha256
from pathlib import Path
import json

from sage.all import QQ, IntegerVectors, PolynomialRing, matrix, prod, vector
from sage.env import SAGE_VERSION
from verify_mixed_regularity import (sp_basis, boundary, centralizer_dimension,
    exterior_matrix, exterior_derivative)

ROOT=Path(__file__).resolve().parents[1]


def gram(basis):
    return matrix(QQ,[[(H*K).trace() for K in basis] for H in basis])


def coordinate(X,basis,Gi):
    return Gi*vector(QQ,[(X*H).trace() for H in basis])


def point(n,kind,m=1):
    if kind=='identity':
        return matrix.identity(QQ,2*n)
    if kind=='unipotent':
        return boundary(n,eigenvalue=1)
    if kind=='generic':
        diagonal=matrix.diagonal([QQ(j+2) for j in range(n)]+
                                 [1/QQ(j+2) for j in range(n)])
        S=matrix.identity(QQ,2*n)
        for i in range(n):
            for j in range(n):
                S[i,n+j]=i+j+1
        return S*diagonal*S.inverse()
    assert kind=='critical'
    g=matrix(QQ,2*n,2*n)
    pair_sums=[QQ(2)-QQ(4*i)/(m*n*(n+1)) for i in range(1,n+1)]
    assert len(set(pair_sums))==n and all(0<x<2 for x in pair_sums)
    for j,x in enumerate(pair_sums):
        g[j,n+j]=-1
        g[n+j,j]=1
        g[n+j,n+j]=x
    return g


def theta(g,m):
    n=g.nrows()//2
    I=matrix.identity(QQ,2*n)
    return QQ(2)/4**(n*m)*(g+I).det()**(m-1)*(g-I)*(g+I).adjugate()


def direct_differentials(g,m,basis,Gi):
    """Differentiate D^m (g-I)(g+I)^-1 directly on gH, no determinant lemma."""
    n=g.nrows()//2
    I=matrix.identity(QQ,2*n)
    Mi=(g+I).inverse()
    C=(g-I)*Mi
    D=(g+I).det()
    columns=[]
    cayley=[]
    for H in basis:
        change=g*H
        dC=(change-C*change)*Mi
        relative_dD=(Mi*change).trace()
        Y=dC+m*relative_dD*C
        columns.append(coordinate(Y,basis,Gi))
        cayley.append(coordinate(dC,basis,Gi))
    J=matrix(QQ,columns).transpose()
    JC=matrix(QQ,cayley).transpose()
    scalar=QQ(2)*D**m/4**(n*m)
    return scalar*J,JC,C


def differential_checks():
    rows=[]
    coordinate_rows=[]
    cache={}
    for n in range(2,6):
        basis=sp_basis(n)
        Gi=gram(basis).inverse()
        N=len(basis)
        for m in range(1,4):
            for kind in ['identity','generic','unipotent','critical']:
                g=point(n,kind,m)
                I=matrix.identity(QQ,2*n)
                D=(g+I).det()
                assert D!=0 and g.det()==1
                J,JC,C=direct_differentials(g,m,basis,Gi)
                delta=(I-C).det()
                assert D*delta==4**n
                predicted=(D/4**n)**((n*m-1)*(2*n+1))*(1-n*m+QQ(m)/2*g.trace())
                assert J.det()==predicted
                key=n,tuple(g.list())
                if key not in cache:
                    assert JC.det()==QQ(1)/2**N*delta**(2*n+1)
                    cache[key]=True
                    coordinate_rows.append(dict(rank=n,kind=kind,tensor_power=m))
                if kind=='identity':
                    assert J==matrix.identity(QQ,N)
                if kind=='critical':
                    assert predicted==0 and J.rank()==N-1
                    assert centralizer_dimension(g,basis)==n
                    assert centralizer_dimension(theta(g,m),basis)==n
                    radial=JC.inverse()*coordinate(C,basis,Gi)
                    assert radial!=0 and J*radial==0
                rows.append(dict(rank=n,tensor_power=m,kind=kind,
                                 dimension=N,Jacobian=str(predicted),det_I_plus_g=str(D),
                                 differential_corank=1 if kind=='critical' else N-J.rank()))
        print('Full matrix differentials passed at rank',n,flush=True)
    return rows,coordinate_rows


def exterior_data(g,basis,derivatives):
    n=g.nrows()//2
    A=[exterior_matrix(g,k) for k in range(2*n+1)]
    chi=sum(M.trace() for M in A)
    first=vector(QQ,[sum((M*D).trace() for M,D in zip(A,ds)) for ds in derivatives])
    second=matrix(QQ,[[
        sum((M*DH*DK).trace() for M,DH,DK in zip(A,ds,es))
        for es in derivatives] for ds in derivatives])
    return chi,first,second


def exterior_checks():
    n=2
    basis=sp_basis(n)
    G=gram(basis)
    Gi=G.inverse()
    derivatives=[[exterior_derivative(H,k) for k in range(2*n+1)] for H in basis]
    actual_Gram=matrix(QQ,[[sum((A*B).trace() for A,B in zip(ds,es))
                           for es in derivatives] for ds in derivatives])
    assert actual_Gram==4**(n-1)*G
    cache={}
    rows=[]
    for m in range(1,4):
        for kind in ['identity','generic','unipotent','critical']:
            g=point(n,kind,m)
            key=tuple(g.list())
            if key not in cache:
                cache[key]=exterior_data(g,basis,derivatives)
            chi,v,T=cache[key]
            I_m=m*16**(m-1)*4
            T_m=m*chi**(m-1)*T
            if m>1:
                T_m+=m*(m-1)*chi**(m-2)*v.column()*v.row()
            J=Gi*T_m.transpose()/I_m
            J_direct,_,_=direct_differentials(g,m,basis,Gi)
            assert J==J_direct
            predicted=((chi/16)**((2*m-1)*5)*(1-2*m+QQ(m)/2*g.trace()))
            assert J.det()==predicted
            image=Gi*(m*chi**(m-1)*v)/I_m
            assert image==coordinate(theta(g,m),basis,Gi)
            rows.append(dict(tensor_power=m,kind=kind,representation_dimension=16**m,
                             actual_trace_projection_matches_full_differential=True))
    print('Actual exterior-representation differential checks passed',flush=True)
    return rows,basis,Gi,derivatives,cache


def harmonic_checks(basis,Gi,derivatives,cache):
    N=len(basis)
    R=PolynomialRing(QQ,N,'u')
    x=R.gens()
    G=Gi.inverse()
    q=sum(G[i,j]*x[i]*x[j] for i in range(N) for j in range(N))
    exponents={d:[tuple(v) for v in IntegerVectors(d,N)] for d in range(4)}
    monomials={d:[prod(v**a for v,a in zip(x,e)) for e in exponents[d]] for d in range(4)}
    harmonics={0:monomials[0],1:monomials[1]}
    for d in [2,3]:
        positions={v:i for i,v in enumerate(exponents[d-2])}
        L=matrix(QQ,len(positions),len(exponents[d]))
        for col,e in enumerate(exponents[d]):
            for i in range(N):
                for j in range(N):
                    if not Gi[i,j] or not e[i] or e[j]<(2 if i==j else 1):
                        continue
                    v=list(e)
                    coefficient=e[i]*(e[j]-(1 if i==j else 0))*Gi[i,j]
                    v[i]-=1
                    v[j]-=1
                    L[positions[tuple(v)],col]+=coefficient
        kernel=L.right_kernel()
        assert kernel.dimension()==len(monomials[d])-len(monomials[d-2])
        harmonics[d]=[sum(c*p for c,p in zip(v,monomials[d])) for v in kernel.basis()]
        separated=harmonics[d]+[q*h for h in harmonics[d-2]]
        M=matrix(QQ,[[p.monomial_coefficient(v) for p in separated] for v in monomials[d]])
        assert M.det()!=0
    assert [len(harmonics[d]) for d in range(4)]==[1,10,54,210]
    evaluations=0
    for m in range(1,4):
        for kind in ['generic','critical']:
            g=point(2,kind,m)
            I=matrix.identity(QQ,4)
            C=(g-I)*(g+I).inverse()
            delta=(I-C).det()
            key=tuple(g.list())
            if key not in cache:
                cache[key]=exterior_data(g,basis,derivatives)
            chi,v,_=cache[key]
            image=Gi*(m*chi**(m-1)*v)/(m*16**(m-1)*4)
            c=coordinate(C,basis,Gi)
            for d in range(4):
                for h in harmonics[d]:
                    assert h(*image)==QQ(2)**d*delta**(-m*d)*h(*c)
                    evaluations+=1
    return dict(group='Sp4',degrees=[0,1,2,3],dimensions=[1,10,54,210],
                genuine_invariant_differential_kernel_bases=True,
                separation_matrix_checks=2,actual_exterior_trace_pullback_evaluations=evaluations,
                scope='Finite degrees below the next basic invariant differential order four; full theorem uses known Kostant separation.')


def necessity_pairs():
    rows=[]
    for n in range(2,7):
        g=boundary(n)
        k=matrix.identity(QQ,2*n)
        k[1,1]=k[n+1,n+1]=0
        k[1,n+1]=1
        k[n+1,1]=-1
        g1=k*g*k.inverse()
        assert g1!=g and g1.charpoly()==g.charpoly()
        assert g[1,1]==2 and g1[1,1]==QQ(1)/2
        for m in range(1,4):
            assert theta(g,m)==theta(g1,m)
            rows.append(dict(rank=n,tensor_power=m,distinguishing_entry=[2,2],
                             first_value='2',second_value='1/2',equal_Springer_images=True))
    return rows


def torus_boundary_audit(basis4,Gi4,derivatives):
    """Check the regularity distinction, including actual traces on D=0.

    Unlike inverse-Cayley calculations, the exterior matrices are defined
    at the boundary. Degenerate strata and Sp2 test the rank hypotheses.
    """
    torus_rows=[]
    boundary_rows=[]
    trace_rows=[]
    rank_one=[]
    for n in range(2,6):
        basis=sp_basis(n)
        for pattern in ['positive','negative','reciprocal','mixed_sign']:
            z=[QQ(j+2) for j in range(n)]
            if pattern=='negative': z=[-v for v in z]
            if pattern=='reciprocal': z=[1/v for v in z]
            if pattern=='mixed_sign': z=[(-1)**j*v for j,v in enumerate(z)]
            g=matrix.diagonal(z+[1/v for v in z])
            source_kernel=matrix(QQ,[(H*g-g*H).list() for H in basis]).left_kernel()
            assert source_kernel.dimension()==n
            assert (g+matrix.identity(QQ,2*n)).det()!=0
            for m in range(1,4):
                image=theta(g,m)
                target_kernel=matrix(QQ,[(H*image-image*H).list()
                                        for H in basis]).left_kernel()
                assert source_kernel==target_kernel
                torus_rows.append(dict(rank=n,tensor_power=m,pattern=pattern,
                                       centralizer_subspaces_equal=True))
        for kind in ['single_Jordan_minus_one','semisimple_minus_plane','two_Jordan_minus_one']:
            g=boundary(n)
            if kind=='semisimple_minus_plane': g[0,n]=0
            if kind=='two_Jordan_minus_one':
                g[1,1]=g[n+1,n+1]=-1
                g[1,n+1]=1
            I=matrix.identity(QQ,2*n)
            assert (g+I).det()==0
            omega=matrix(QQ,2*n,2*n)
            for j in range(n):
                omega[j,n+j]=1
                omega[n+j,j]=-1
            assert g.transpose()*omega*g==omega
            source_dimension=centralizer_dimension(g,basis)
            assert (source_dimension==n)==(kind=='single_Jordan_minus_one')
            if n==2: chi,v,_=exterior_data(g,basis4,derivatives)
            for m in range(1,4):
                image=theta(g,m)
                expected_rank=int(kind=='single_Jordan_minus_one' and m==1)
                assert image.rank()==expected_rank
                image_dimension=centralizer_dimension(image,basis)
                assert image_dimension==(2*n*n-n if expected_rank else n*(2*n+1))
                assert image_dimension>n
                boundary_rows.append(dict(rank=n,tensor_power=m,kind=kind,
                    source_centralizer_dimension=source_dimension,image_rank=expected_rank,
                    image_centralizer_dimension=image_dimension))
                if n==2:
                    index=m*16**(m-1)*4
                    projected=Gi4*(m*chi**(m-1)*v)/index
                    assert projected==coordinate(image,basis4,Gi4)
                    trace_rows.append(dict(kind=kind,tensor_power=m,
                        actual_exterior_trace_at_boundary=True))
    basis=sp_basis(1)
    g=boundary(1)
    assert centralizer_dimension(g,basis)==1
    for m in range(1,4):
        image=theta(g,m)
        image_dimension=centralizer_dimension(image,basis)
        assert image_dimension==(1 if m==1 else 3)
        rank_one.append(dict(tensor_power=m,image_centralizer_dimension=image_dimension,
                             globally_positive_control=(m==1)))
    return dict(regular_torus_subspace_checks=torus_rows,boundary_strata_checks=boundary_rows,
                actual_boundary_trace_checks=trace_rows,rank_one_boundary_controls=rank_one)


def main():
    differentials,coordinates=differential_checks()
    exterior,basis,Gi,derivatives,cache=exterior_checks()
    harmonics=harmonic_checks(basis,Gi,derivatives,cache)
    pairs=necessity_pairs()
    audit=torus_boundary_audit(basis,Gi,derivatives)
    paths=['SYMPLECTIC_HARMONIC_LOCALIZATION.md','tools/verify_symplectic_harmonic_localization.py',
           'tools/verify_mixed_regularity.py','MIXED_REGULARITY.md']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        arbitrary_rank_proof='SYMPLECTIC_HARMONIC_LOCALIZATION.md, Sections 2-5 and 7',
        group='Sp_(2n)(C), n>=2',coefficient_field='C; exact rational finite checks',
        harmonic_localization='C[G][D^-1] = C[G]^G[D^-1] tensor theta_m^*Harm(g), module isomorphism',
        minimal_invariant_principal_localization='For nonzero invariant f, iff D divides f',
        full_Cayley_determinant='(D/4^n)^((nm-1)(2n+1)) (1-nm+(m/2)tr(g))',
        full_matrix_differential_checks=len(differentials),
        independent_Cayley_coordinate_Jacobian_checks=len(coordinates),
        regular_critical_corank_one_checks=sum(r['kind']=='critical' for r in differentials),
        actual_exterior_tensor_differential_checks=len(exterior),
        low_degree_harmonic_checks=harmonics,
        explicit_necessity_conjugacy_pairs=len(pairs),
        torus_boundary_audit=audit,
        matrix_differentials=differentials,Cayley_coordinate_checks=coordinates,
        exterior_tensor_checks=exterior,conjugacy_pairs=pairs,
        known_input_corollary=True,novelty_claim=False,runtime_API_changed=False,
        scope='Exact finite corroboration; uniform decomposition follows from known Kostant separation and Cayley coordinates, not finite exhaustion or an external proof review.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'symplectic-harmonic-localization-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    compact={k:v for k,v in report.items() if k not in
        ('source_sha256','matrix_differentials','Cayley_coordinate_checks','exterior_tensor_checks',
         'conjugacy_pairs','torus_boundary_audit')}
    compact['torus_boundary_audit_counts']={k:len(v) for k,v in audit.items()}
    print(json.dumps(compact,indent=2))


if __name__=='__main__':
    main()
