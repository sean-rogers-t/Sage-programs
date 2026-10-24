"""Exact audits of mixed-character B/C regularity; run with Sage Python.

Independent constructions: Freudenthal versus Chebyshev sign orbits;
exterior subsets versus irreducible decompositions and binomial tensor
weights; actual exterior matrices versus the adjugate Springer formula.
Type B witnesses use a finite algebra, not floating point root finding.
"""
from collections import defaultdict
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import json

from sage.all import QQ, PolynomialRing, QuadraticField, WeylCharacterRing, matrix, prod
from sage.env import SAGE_VERSION

ROOT = Path(__file__).resolve().parents[1]


def weights(group,n,labels,rings):
    key=group,n
    if key not in rings:
        rings[key]=WeylCharacterRing(group+str(n),style='coroots')
    result={}
    scale=2 if group=='B' else 1
    for w,m in rings[key](tuple(labels)).weight_multiplicities().items():
        coords=tuple(scale*v for v in w.to_vector())
        assert all(v==int(v) for v in coords)
        result[tuple(int(v) for v in coords)]=int(m)
    return result


def add(*summands):
    result=defaultdict(int)
    for coefficient,w in summands:
        for v,m in w.items():
            result[v]+=coefficient*m
    return dict(result)


def convolution(a,b):
    result=defaultdict(int)
    for v,m in a.items():
        for w,k in b.items():
            result[tuple(x+y for x,y in zip(v,w))]+=m*k
    return dict(result)


def polynomial_character(w,R):
    variables=R.gens()
    maximum=max(max(v) for v in w)
    T=[[R(1),x] for x in variables]
    for j,x in enumerate(variables):
        if maximum>=2:
            T[j].append(x*x-2)
        while len(T[j])<=maximum:
            T[j].append(x*T[j][-1]-T[j][-2])
    F=R(0)
    for v,m in w.items():
        if all(k>=0 for k in v):
            for signs in product([-1,1],repeat=len(v)):
                assert w[tuple(s*k for s,k in zip(signs,v))]==m
            F+=m*prod(T[j][k] for j,k in enumerate(v))
    return F


def trace(w,point,scale=1):
    n=len(point)
    gram=matrix(QQ,n,n,lambda i,j:sum(QQ(m)*v[i]*v[j]/scale**2 for v,m in w.items()))
    S=gram[0,0]
    assert S>0 and gram==S*matrix.identity(QQ,n)
    evaluated=[(v,m,prod(a**k for a,k in zip(point,v))) for v,m in w.items()]
    return (sum(m*ev for v,m,ev in evaluated),S,
            [sum(QQ(m)*v[j]*ev/scale for v,m,ev in evaluated)/S for j in range(n)])


def sp_basis(n):
    result=[]
    for i in range(n):
        for j in range(n):
            H=matrix(QQ,2*n,2*n)
            H[i,j]=1
            H[n+j,n+i]=-1
            result.append(H)
    for i in range(n):
        for j in range(i,n):
            for side in [0,1]:
                H=matrix(QQ,2*n,2*n)
                if side==0:
                    H[i,n+j]=H[j,n+i]=1
                else:
                    H[n+i,j]=H[n+j,i]=1
                result.append(H)
    return result


def boundary(n,eigenvalue=-1):
    g=matrix.diagonal([QQ(1)]*(2*n))
    g[0,0]=g[n,n]=eigenvalue
    g[0,n]=1
    for j in range(1,n):
        g[j,j]=j+1
        g[n+j,n+j]=QQ(1)/(j+1)
    J=matrix(QQ,2*n,2*n)
    for j in range(n):
        J[j,n+j]=1
        J[n+j,j]=-1
    assert g.transpose()*J*g==J and g.det()==1
    return g


def centralizer_dimension(g,basis):
    return len(basis)-matrix(QQ,[(H*g-g*H).list() for H in basis]).rank()


def exterior_matrix(g,k):
    subsets=list(combinations(range(g.nrows()),k))
    return matrix(QQ,[[g.matrix_from_rows_and_columns(list(I),list(J)).det()
                       for J in subsets] for I in subsets])


def exterior_derivative(H,k):
    subsets=list(combinations(range(H.nrows()),k))
    indices={I:j for j,I in enumerate(subsets)}
    result=matrix(QQ,len(subsets),len(subsets))
    for col,I in enumerate(subsets):
        for pos,j in enumerate(I):
            for q in range(H.nrows()):
                if not H[q,j] or q in I[:pos]+I[pos+1:]:
                    continue
                J=I[:pos]+(q,)+I[pos+1:]
                sign=(-1)**sum(J[s]>J[t] for s in range(k) for t in range(s+1,k))
                result[indices[tuple(sorted(J))],col]+=sign*H[q,j]
    return result


def exterior_weights(n):
    basis=[tuple(sign if i==j else 0 for i in range(n))
           for j in range(n) for sign in [-1,1]]
    result=defaultdict(int)
    for k in range(2*n+1):
        for subset in combinations(basis,k):
            result[tuple(sum(v[j] for v in subset) for j in range(n))]+=1
    return dict(result)


def sp_checks(rings):
    K=QuadraticField(-1,'ii')
    ii=K.gen()
    rows=[]
    for n in range(2,5):
        R=PolynomialRing(QQ,n,'x')
        x=R.gens()
        P=prod(v+2 for v in x)
        base=exterior_weights(n)
        # Independently decompose each exterior degree into primitive kernels.
        pieces=[]
        for k in range(2*n+1):
            for r in range(min(k,2*n-k),-1,-2):
                a=tuple(1 if j==r-1 else 0 for j in range(n))
                pieces.append((1,weights('C',n,a,rings)))
        assert base==add(*pieces)
        assert polynomial_character(base,R)==P
        w={(0,)*n:1}
        for m in range(1,4):
            w=convolution(w,base)
            # Binomial model follows directly from (z+2+z^-1)^m.
            binomial={v:prod(comb(2*m,m+k) for k in v)
                      for v in product(range(-m,m+1),repeat=n)}
            assert w==binomial
            F=polynomial_character(w,R)
            assert F==P**m
            for v in x:
                assert F.derivative(v)==m*P**m/(v+2)
            points=[[K(j+2) for j in range(n)]]
            central=[K(v) for j in range(n//2) for v in [j+2,-j-2]]
            if n%2:
                central.append(ii)
            points.append(central)
            points.append([K(-j-2) for j in range(n)])
            for z in points:
                assert all(v*v!=1 for v in z)
                assert all(z[i]!=z[j] and z[i]*z[j]!=1 for i in range(n) for j in range(i+1,n))
                chi,S,image=trace(w,z)
                det=prod(v+2+1/v for v in z)
                assert chi==det**m and S==2*m*4**(n*m-1)
                expected=[QQ(2)/4**(n*m)*det**m*(v-1)/(v+1) for v in z]
                assert image==expected
                assert all(y!=0 for y in image)
                assert all(image[i]!=image[j] and image[i]!=-image[j]
                           for i in range(n) for j in range(i+1,n))
            rows.append(dict(rank=n,tensor_power=m,dimension=sum(w.values()),
                             trace_Gram_scalar=str(S),torus_points=len(points)))
    boundary_rows=[]
    matrix_rows=[]
    for n in range(2,7):
        g=boundary(n)
        basis=sp_basis(n)
        I=matrix.identity(QQ,2*n)
        adj=(g+I).adjugate()
        assert (g+I).det()==0
        theta=QQ(2)/4**n*(g-I)*adj
        D=prod((1+QQ(j))*(1+1/QQ(j)) for j in range(2,n+1))
        expected=matrix(QQ,2*n,2*n)
        expected[0,n]=4*D/4**n
        assert theta==expected and theta.rank()==1 and theta**2==0
        din=centralizer_dimension(g,basis)
        dout=centralizer_dimension(theta,basis)
        assert din==n and dout==2*n*n-n and dout>n
        boundary_rows.append(dict(rank=n,input_centralizer=din,output_centralizer_m1=dout,
                                  output_centralizer_m_ge2=len(basis)))
        # A regular nonsemisimple input inside the good open set stays regular.
        positive=boundary(n,eigenvalue=1)
        assert (positive+I).det()!=0
        assert centralizer_dimension(positive,basis)==n
        open_dims=[]
        for m in range(1,4):
            out=QQ(2)/4**(n*m)*(positive+I).det()**m*(positive-I)*(positive+I).inverse()
            dim=centralizer_dimension(out,basis)
            assert dim==n
            open_dims.append(dim)
        boundary_rows[-1]['good_open_output_centralizers_m_1_2_3']=open_dims
        if n<=3:
            exterior=[exterior_matrix(g,k) for k in range(2*n+1)]
            assert sum(A.trace() for A in exterior)==0
            for H in basis:
                derivatives=[exterior_derivative(H,k) for k in range(2*n+1)]
                rhs=sum((A*B).trace() for A,B in zip(exterior,derivatives))
                assert rhs==4**(n-1)*(theta*H).trace()
                # Every tensor power m>=2 has zero character derivative here.
                assert 2*sum(A.trace() for A in exterior)*rhs==0
            H=basis[0]
            S=sum((exterior_derivative(H,k)**2).trace() for k in range(2*n+1))
            assert S==2*4**(n-1)
            matrix_rows.append(dict(rank=n,Lie_basis_checks=len(basis),
                                    exterior_dimension=4**n,trace_Gram_scalar=str(S)))
    return rows,boundary_rows,matrix_rows


def algebra_add(a,b):
    result=dict(a)
    for mask,c in b.items():
        result[mask]=result.get(mask,0)+c
        if result[mask]==0:
            del result[mask]
    return result


def algebra_scale(c,a):
    return {mask:c*v for mask,v in a.items() if c*v!=0}


def algebra_multiply(a,b,discriminants):
    result={}
    for v,c in a.items():
        for w,d in b.items():
            common=v&w
            coefficient=c*d*prod(discriminants[j] for j in range(len(discriminants)) if common>>j&1)
            mask=v^w
            result[mask]=result.get(mask,0)+coefficient
            if result[mask]==0:
                del result[mask]
    return result


def algebra_power(a,k,discriminants):
    result={0:discriminants[0].parent()(1)}
    for _ in range(k):
        result=algebra_multiply(result,a,discriminants)
    return result


def exact_B_witness(w,F):
    R=F.parent()
    b=R.gens()
    n=len(b)
    Q=F.derivative(b[0])
    # Every possible boundary linear factor is excluded independently.
    boundary_substitutions=0
    for j in range(n):
        for value in [-2,2]:
            assert Q.subs({b[j]:value})!=0
            boundary_substitutions+=1
    for j in range(n):
        for k in range(j+1,n):
            for sign in [-1,1]:
                assert Q.subs({b[j]:sign*b[k]})!=0
                boundary_substitutions+=1
    assert Q.total_degree()>0
    for j,k in combinations(range(n),2):
        values=[-v if t in [j,k] else v for t,v in enumerate(b)]
        assert F(*values)==F
        assert Q(*values)==(-Q if 0 in [j,k] else Q)
    variable=next(j for j,v in enumerate(b) if Q.degree(v)>0)
    U=PolynomialRing(QQ,'X')
    X=U.gen()
    fixed=[j for j in range(n) if j!=variable]
    trials=0
    for values in product(range(3,11),repeat=n-1):
        if len(set(values))<n-1:
            continue
        trials+=1
        substitution=[None]*n
        for j,v in zip(fixed,values):
            substitution[j]=U(v)
        substitution[variable]=X
        q=U(R.hom(substitution,U)(Q))
        if q.degree()<1 or any(q(v)==0 for v in [-2,2]+list(values)+[-v for v in values]):
            continue
        break
    else:
        raise AssertionError('Finite witness grid exhausted')
    q=q.monic()
    A=U.quotient(q,'xx')
    bb=[A(v) if j!=variable else A.gen() for j,v in enumerate(substitution)]
    discriminants=[v*v-4 for v in bb]
    # q is coprime to every boundary factor; all roots give regular inputs.
    assert all(q(v)!=0 for v in [-2,2]+list(values)+[-v for v in values])
    ap=[{0:v/2,1<<j:A(1)/2} for j,v in enumerate(bb)]
    am=[{0:v/2,1<<j:-A(1)/2} for j,v in enumerate(bb)]
    powers={}
    for j in range(n):
        for k in {v[j] for v in w}:
            powers[j,k]=algebra_power(ap[j] if k>=0 else am[j],abs(k),discriminants)
    numerator={}
    for v,m in w.items():
        if not v[0]:
            continue
        ev={0:A(1)}
        for j,k in enumerate(v):
            ev=algebra_multiply(ev,powers[j,k],discriminants)
        numerator=algebra_add(numerator,algebra_scale(QQ(m)*v[0]/2,ev))
    assert not numerator
    return dict(variable=variable+1,fixed_b=[None if j==variable else str(substitution[j]) for j in range(n)],
                root_polynomial=str(q),root_polynomial_degree=int(q.degree()),grid_trials=trials,
                finite_algebra_dimension=int(q.degree())*2**n,
                direct_weight_trace_zero=True,boundary_factor_exclusions=boundary_substitutions)


def B_checks(rings):
    rows=[]
    for n in range(3,6):
        R=PolynomialRing(QQ,n,'b')
        b=R.gens()
        vec=(1,)+(0,)*(n-1)
        spin=(0,)*(n-1)+(1,)
        adj=(0,1)+(0,)*(n-2)
        sym=(2,)+(0,)*(n-1)
        high=(1,)+(0,)*(n-2)+(1,)
        even_spin=(0,)*(n-1)+(2,)
        zero=(0,)*n
        cases=[[(1,vec),(1,spin)],[(1,vec),(2,spin),(3,zero)],
               [(1,adj),(1,spin)],[(1,sym),(1,spin)],
               [(1,vec),(1,high)],[(1,even_spin),(1,spin)]]
        for case in cases:
            w=add(*[(m,weights('B',n,a,rings)) for m,a in case])
            F=polynomial_character(w,R)
            for i in range(n):
                assert F.derivative(b[i])(*([2]*n))>0
            point=[QQ(j+2) for j in range(n)]
            chi,S,image=trace(w,point,scale=2)
            values=[v+1/v for v in point]
            assert chi==F(*values)
            assert image==[(point[j]-1/point[j])*F.derivative(b[j])(*values)/(2*S) for j in range(n)]
            row=exact_B_witness(w,F)
            row.update(rank=n,summands=[dict(multiplicity=m,labels=list(a)) for m,a in case],
                       dimension=sum(w.values()))
            rows.append(row)
        print('Type B mixed exact witnesses passed at rank',n,flush=True)
    # B2 exception in independent B weights and coordinates.
    n=2
    w=add((1,weights('B',n,(1,0),rings)),(2,weights('B',n,(0,1),rings)),
          (3,weights('B',n,(0,0),rings)))
    R=PolynomialRing(QQ,2,'b')
    b=R.gens()
    assert polynomial_character(w,R)==(b[0]+b[1])**2
    K=QuadraticField(-1,'ii')
    point=[K.gen(),K(2)]
    chi,S,image=trace(w,point,scale=2)
    assert all(y!=0 for y in image) and image[0]!=image[1] and image[0]!=-image[1]
    return rows,dict(character='(b1+b2)^2',dimension=sum(w.values()),
                     image_at_a_i_2=[str(y) for y in image])


def main():
    rings={}
    sp,bd,mat=sp_checks(rings)
    print('Symplectic character, torus, boundary and actual exterior-matrix checks passed',flush=True)
    b,exception=B_checks(rings)
    paths=['MIXED_REGULARITY.md','tools/verify_mixed_regularity.py']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        arbitrary_rank_proof='MIXED_REGULARITY.md, Sections 2-6',
        Sp_torus_character_classification='A det(1+g)^m+B, m>=1; actual representation characters',
        Sp_global_regularity_preserving_representations='none at rank >=2',
        Spin_B_torus_regularity_preserving_representations='none at rank >=3; rank two classified via Sp4',
        Spin_B_global_regularity_preserving_representations='none at rank >=2',
        exterior_character_decomposition_ranks=[2,3,4],
        Sp_tensor_character_cases=len(sp),Sp_regular_torus_trace_checks=sum(r['torus_points'] for r in sp),
        independent_tensor_binomial_models=len(sp),Sp_nonsemisimple_boundary_ranks=len(bd),
        Sp_regular_nonsemisimple_good_open_checks=3*len(bd),
        full_exterior_matrix_boundary_ranks=len(mat),
        full_exterior_matrix_Lie_basis_checks=sum(r['Lie_basis_checks'] for r in mat),
        mixed_B_direct_finite_algebra_witnesses=len(b),
        mixed_B_boundary_linear_factor_exclusions=sum(r['boundary_factor_exclusions'] for r in b),
        B2_torus_positive_control=exception,
        symplectic_cases=sp,boundary_centralizers=bd,exterior_matrix_checks=mat,mixed_B_witnesses=b,
        runtime_API_changed=False,novelty_claim=False,
        scope='Exact finite independent corroboration of a written all-representation proof; no external review or priority claim.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'mixed-regularity-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('source_sha256','symplectic_cases','boundary_centralizers','exterior_matrix_checks','mixed_B_witnesses')},indent=2))


if __name__=='__main__':
    main()
