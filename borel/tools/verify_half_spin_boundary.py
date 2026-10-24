"""Independent Jordan-matrix and chiral-trace checks for the half-spin zero fiber."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

from sage.all import (QQ, Partitions, matrix, block_diagonal_matrix, prod)
from sage.env import SAGE_VERSION
from verify_half_spin_cayley import (QI, spin_weights, moment_without_b,
                                    chiral_data, trace_product)

ROOT=Path(__file__).resolve().parents[1]


def orthogonal_partitions(n):
    return [tuple(p) for p in Partitions(n)
            if all(d%2 or m%2==0 for d,m in Counter(p).items())]


def jordan(d):
    return matrix(QI,d,d,lambda i,j: int(j==i+1))


def orthogonal_jordan(p):
    """Construct N and its symmetric form, without an orbit-table API."""
    ns=[]; forms=[]
    for d,m in sorted(Counter(p).items()):
        A=matrix(QI,d,d,lambda i,j: (-1)**i if i+j==d-1 else 0)
        if d%2:
            for _ in range(m): ns.append(jordan(d)); forms.append(A)
        else:
            assert m%2==0
            for _ in range(m//2):
                ns.append(block_diagonal_matrix(jordan(d),jordan(d)))
                B=matrix(QI,2*d,2*d)
                B[:d,d:]=A; B[d:,:d]=A.transpose(); forms.append(B)
    N=block_diagonal_matrix(ns); B=block_diagonal_matrix(forms)
    assert B.transpose()==B and B.det()!=0 and N.transpose()*B+B*N==0
    return N,B


def to_standard(B):
    """An explicit Q(i) congruence for the signed matching form B."""
    n=B.nrows(); R=matrix(QI,n,n); used=set(); ii=QI.gen()
    for i in range(n):
        if i in used: continue
        j=next(j for j in range(n) if B[i,j])
        a=B[i,j]
        if i==j:
            assert a in (1,-1); R[i,i]=1 if a==1 else ii
        else:
            R[i,i]=1; R[j,i]=1/(2*a)
            R[i,j]=ii; R[j,j]=-ii/(2*a)
        used.update([i,j])
    assert R.transpose()*B*R==matrix.identity(QI,n)
    return R


@lru_cache(None)
def standard_nilpotent(p):
    N,B=orthogonal_jordan(p); R=to_standard(B); X=R.inverse()*N*R
    assert X.transpose()==-X
    assert all(X.nrows()-(X**j).rank()==sum(min(j,d) for d in p)
               for j in range(1,max(p)+1))
    return X


def nil_exp(N):
    result=matrix.identity(N.base_ring(),N.nrows()); term=result
    for j in range(1,N.nrows()+1):
        term=term*N/j; result+=term
        if not term: return result
    raise AssertionError('Expected a nilpotent exponential')


def so_centralizer(p):
    columns=[sum(d>=j for d in p) for j in range(1,max(p)+1)]
    return (sum(x*x for x in columns)-sum(d%2 for d in p))//2


def wedge_centralizer(p):
    # sl2 weights on wedge^2 of the natural module: weights 0 and 1
    # count all highest-weight summands, hence the nilpotent centralizer.
    weights=[d-1-2*j for d in p for j in range(d)]
    return sum(weights[i]+weights[j] in (0,1)
               for i,j in combinations(range(len(weights)),2))


def direct_centralizer(X):
    n=X.nrows(); pairs=list(combinations(range(n),2)); columns=[]
    for i,j in pairs:
        E=matrix(QI,n,n); E[i,j]=1; E[j,i]=-1
        C=X*E-E*X
        columns.append([C[a,b] for a,b in pairs])
    return len(pairs)-matrix(QI,columns).rank()


def partition_checks():
    rows=[]; matrix_rows=[]
    for k in range(2,16):
        admissible=[]
        for p in orthogonal_partitions(2*k):
            z=so_centralizer(p); assert z==wedge_centralizer(p)
            if len(p)>=4:
                assert z>=k+4; admissible.append(z)
            if k<=5:
                actual=direct_centralizer(standard_nilpotent(p))
                assert actual==z
                matrix_rows.append(dict(partition=p,centralizer_dimension=actual))
            rows.append(dict(rank=k,partition=p,centralizer_dimension=z,
                             boundary_allowed=len(p)>=4))
        assert min(admissible)==k+4
        assert so_centralizer((2*k-3,1,1,1))==k+4
    gl_checks=0
    for m in range(1,16):
        for p in Partitions(m):
            columns=list(p.conjugate()); assert sum(x*x for x in columns)>=m
            gl_checks+=1
    return rows,matrix_rows,gl_checks


def semisimple_checks():
    rows=[]
    for r in range(3,10):
        for k in range(1,r+1):
            for l in range(r-k+1):
                zs=[QQ(-1)]*k+[QQ(1)]*l+[QQ(j+2)**2 for j in range(r-k-l)]
                expected=(k==r or (k>=2 and l>=2))
                for sigma in (1,-1):
                    actual=moment_without_b(spin_weights(r,sigma),zs)==0
                    assert actual==expected,(r,k,l,sigma)
                    rows.append(dict(rank=r,minus_pairs=k,plus_pairs=l,
                                     chirality=sigma,zero=bool(actual)))
    return rows


def full_boundary_traces(pm,pp,rest=(),conjugate=False):
    """Actual Spin matrices, not evaluation of the proposed boundary formula.

    rest=(m1,m2,...) specifies regular unipotents in GL_mj blocks at
    distinct reciprocal eigenvalue pairs; their auxiliary spin parameters
    are 2,3,... . Empty plus/minus blocks are allowed as negative controls.
    """
    k=sum(pm)//2; l=sum(pp)//2; r=k+l+sum(rest); n=2*r
    N=matrix(QI,n,n); s=matrix.identity(QI,n); rest_planes=[]
    if pm:
        N[:2*k,:2*k]=standard_nilpotent(tuple(pm))
        s[:2*k,:2*k]=-matrix.identity(QI,2*k)
    if pp: N[2*k:2*(k+l),2*k:2*(k+l)]=standard_nilpotent(tuple(pp))
    offset=2*(k+l)
    for j,m in enumerate(rest):
        J=jordan(m); NN=block_diagonal_matrix(J,-J.transpose())
        B=matrix(QI,2*m,2*m)
        B[:m,m:]=matrix.identity(QI,m); B[m:,:m]=matrix.identity(QI,m)
        R=to_standard(B); aa=QQ(j+2); z=aa*aa
        ss=matrix.diagonal(QI,[z]*m+[1/z]*m)
        NN=R.inverse()*NN*R; ss=R.inverse()*ss*R
        N[offset:offset+2*m,offset:offset+2*m]=NN
        s[offset:offset+2*m,offset:offset+2*m]=ss
        for t in range(m): rest_planes.append((offset+t,offset+m+t,aa))
        offset+=2*m
    assert N.transpose()==-N and s*N==N*s
    gam,_,pairs,_,_=chiral_data(r,1); K=2**r; I=matrix.identity(QI,K)
    Bs=[gam[i]*gam[j]/2 for i,j in pairs]
    Bspin=sum((N[i,j]*B for (i,j),B in zip(pairs,Bs)),matrix(QI,K,K))
    S=I
    for j in range(k): S=S*gam[2*j]*gam[2*j+1]
    for i,j,aa in rest_planes:
        alpha=(aa+1/aa)/2; beta=s[i,j]/(2*alpha)
        assert alpha*alpha+beta*beta==1
        S=S*(alpha*I+beta*gam[i]*gam[j])
    g=S*nil_exp(Bspin); a=s*nil_exp(N)
    if conjugate:
        R=QQ(3)/5*I+QQ(4)/5*gam[0]*gam[-1]
        g=R*g*R.inverse()
        aa=matrix.identity(QI,n); aa[0,0]=aa[n-1,n-1]=-QQ(7)/25
        aa[0,n-1]=QQ(24)/25; aa[n-1,0]=-QQ(24)/25
        a=aa*a*aa.transpose()
    # Check the cover itself, not just traces at a matrix assumed to lift a.
    assert all(g*gam[j]==sum((a[i,j]*gam[i] for i in range(n)),matrix(QI,K,K))*g
               for j in range(n))
    assert a.transpose()*a==matrix.identity(QI,n) and a.det()==1
    km=n-(a+matrix.identity(QI,n)).rank()
    kp=n-(a-matrix.identity(QI,n)).rank()
    central=a in (matrix.identity(QI,n),-matrix.identity(QI,n))
    expected=central or (km>=4 and kp>=4)
    rows=[]
    for sigma in (1,-1):
        _,indices,_,Ds,_=chiral_data(r,sigma)
        block=g.matrix_from_rows_and_columns(indices,indices)
        coords=[-QQ(4)/(2**(r-1))*trace_product(block,D) for D in Ds]
        actual=all(x==0 for x in coords)
        assert actual==expected,(pm,pp,rest,sigma,km,kp)
        rows.append(dict(rank=r,minus_partition=pm,plus_partition=pp,
                         regular_GL_blocks=rest,chirality=sigma,
                         minus_eigenvectors=km,plus_eigenvectors=kp,
                         image_zero=actual,conjugated=conjugate))
    return rows


def chiral_checks():
    rows=[]
    for k,l,rest in [(2,2,()),(2,2,(1,)),(2,3,()),(3,2,())]:
        for pm,pp in product(orthogonal_partitions(2*k),orthogonal_partitions(2*l)):
            rows.extend(full_boundary_traces(pm,pp,rest))
        print('Finished chiral boundary block:',k,l,rest,flush=True)
    extra=[((1,)*4,(1,)*4,(2,),False),
           ((1,)*4,(3,1,1,1),(1,),True),
           ((1,)*4,(3,3),(1,),False),
           ((1,)*4,(2,2,2,2),(),False),
           ((1,)*4,(1,)*6,(),True),
           ((1,)*10,(),(),False),
           ((3,1,1,1,1,1,1,1),(),(),False),
           ((1,)*6,(1,)*2,(1,),False)]
    for pm,pp,rest,conj in extra:
        rows.extend(full_boundary_traces(pm,pp,rest,conj))
    return rows


def dimension_checks():
    rows=[]
    for r in range(4,21):
        d=r*(2*r-1); bounds=[]
        for k in range(2,r-1):
            for l in range(2,r-k+1):
                m=r-k-l
                bound=d-(k+4)-(l+4)-m+m
                bounds.append((bound,k,l,m))
        assert max(v[0] for v in bounds)==d-12
        assert all(k==l==2 for v,k,l,m in bounds if v==d-12)
        rows.append(dict(rank=r,zero_fiber_dimension=d-12,
                         maximizing_minus_pairs=2,maximizing_plus_pairs=2))
    spin10=[]
    for k,l,m in [(2,2,1),(2,3,0),(3,2,0)]:
        for pm,pp in product(orthogonal_partitions(2*k),orthogonal_partitions(2*l)):
            if len(pm)<4 or len(pp)<4: continue
            cent=so_centralizer(pm)+so_centralizer(pp)+m
            spin10.append(dict(minus_partition=pm,plus_partition=pp,
                               free_reciprocal_pairs=m,centralizer_dimension=cent,
                               orbit_dimension=45-cent,stratum_dimension=45-cent+m))
    assert sorted(row['stratum_dimension'] for row in spin10)==[24,24,30,30,32,32,33]
    return rows,spin10


def main():
    semisimple=semisimple_checks()
    partitions,centralizers,gl=partition_checks()
    chiral=chiral_checks(); dimensions,spin10=dimension_checks()
    paths=['HALF_SPIN_BOUNDARY.md','tools/verify_half_spin_boundary.py',
           'HALF_SPIN_CAYLEY.md','tools/verify_half_spin_cayley.py',
           'tools/verify_spin_cayley_geometry.py']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        all_rank_underlying_zero_fiber_written_proof=True,
        all_rank_dimension_written_proof=True,Spin10_zero_fiber_dimension=33,
        boundary_scheme_reduced_proved=False,boundary_irreducible_components_classified=False,
        generic_degree_determined=False,novelty_claim=False,externally_reviewed=False,
        runtime_API_changed=False,
        semisimple_weight_checks=len(semisimple),partition_sl2_checks=len(partitions),
        direct_matrix_centralizer_checks=len(centralizers),GL_partition_bounds=gl,
        full_chiral_boundary_checks=len(chiral),finite_rank_dimension_checks=len(dimensions),
        Spin10_Jordan_strata=spin10,semisimple_checks=semisimple,
        partition_checks=partitions,centralizer_checks=centralizers,
        chiral_checks=chiral,dimension_checks=dimensions,
        scope='Exact independent finite corroboration of separate all-rank proofs; the classification is set-theoretic and gives dimension, not boundary scheme structure or generic degree.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'half-spin-boundary-verification.json').write_text(json.dumps(report,indent=2,default=int)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
         ('source_sha256','semisimple_checks','partition_checks','centralizer_checks',
          'chiral_checks','dimension_checks')},indent=2,default=int))


if __name__=='__main__': main()
