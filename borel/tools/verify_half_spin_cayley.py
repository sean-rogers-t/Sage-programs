"""Individual half-spin: weights, chiral matrices and exact zero-fiber ranks."""
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

from sage.all import (QQ, QuadraticField, CyclotomicField, WeylCharacterRing,
                      PolynomialRing, matrix, vector, prod, lcm)
from sage.env import SAGE_VERSION
from verify_spin_cayley_geometry import gamma_matrices
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]
QI=QuadraticField(-1,'ii')


def spin_weights(r,sigma):
    return [e for e in product([1,-1],repeat=r) if prod(e)==sigma]


def moment_without_b(es,zs):
    A=zs[0].parent(); r=len(zs)
    return vector(A,[sum(QQ(e[i])/2*prod(z**((ei-1)//2) for z,ei in zip(zs,e))
                        for e in es) for i in range(r)])


def weight_checks():
    rows=[]
    for r in range(3,9):
        W=WeylCharacterRing('D'+str(r),style='coroots')
        for sigma in [1,-1]:
            labels=[0]*r; labels[-1 if sigma==1 else -2]=1
            ch=W(tuple(labels))
            actual={tuple(int(2*x) for x in w.to_vector()):int(m)
                    for w,m in ch.weight_multiplicities().items()}
            es=spin_weights(r,sigma)
            assert actual=={e:1 for e in es}
            G=matrix(QQ,r,r,lambda i,j:sum(QQ(e[i]*e[j])/4 for e in es))
            assert G==2**(r-3)*matrix.identity(QQ,r)
            aa=[QQ(i+2) for i in range(r)]; zs=[a*a for a in aa]; b=prod(aa)
            xs=[(z-1)/(z+1) for z in zs]; c=prod(a+1/a for a in aa)/2**r
            image=G.inverse()*(b*moment_without_b(es,zs))
            expected=vector(QQ,[2*c*(x+sigma*prod(xs[j] for j in range(r) if j!=i))
                                for i,x in enumerate(xs)])
            assert image==expected
            assert matrix(QQ,[xs,list(image)]).rank()==2
            rows.append(dict(rank=r,chirality=sigma,dimension=2**(r-1),
                             highest_weight_comparison=True,trace_formula=True,not_radial=True))
    return rows


def trace_product(A,B):
    return sum(v*B[j,i] for (i,j),v in A.dict().items())


def pfaffian(A):
    if A.nrows()==0: return A.base_ring().one()
    return sum((-1)**(j+1)*A[0,j]*pfaffian(A.matrix_from_rows_and_columns(
        [k for k in range(A.nrows()) if k not in [0,j]],
        [k for k in range(A.nrows()) if k not in [0,j]]))
        for j in range(1,A.nrows()) if A[0,j])


def pfaffian_cofactor(X):
    N=X.nrows(); r=N//2; out=matrix(QI,N,N)
    for i,j in combinations(range(N),2):
        keep=[k for k in range(N) if k not in [i,j]]
        out[i,j]=QI.gen()**r*(-1)**(i+j)*pfaffian(X.matrix_from_rows_and_columns(keep,keep))
        out[j,i]=-out[i,j]
    assert X*out==QI.gen()**r*pfaffian(X)*matrix.identity(QI,N)
    return out


@lru_cache(None)
def chiral_data(r,sigma):
    N=2*r; gam=gamma_matrices(N,QI)
    indices=[j for j in range(2**r) if (-1)**int(j).bit_count()==sigma]
    pairs=list(combinations(range(N),2))
    Bs=[(gam[i]*gam[j]/2).matrix_from_rows_and_columns(indices,indices) for i,j in pairs]
    dim=2**(r-1)
    G=matrix(QI,[[trace_product(B,C) for C in Bs] for B in Bs])
    assert G==(-QQ(dim)/4)*matrix.identity(QI,len(pairs))
    es=[tuple(1 if not (j>>(r-1-i))&1 else -1 for i in range(r)) for j in indices]
    return gam,indices,pairs,Bs,es


def Clifford_formula_checks():
    rows=[]
    for r in [3,4,5]:
        N=2*r; K=2**r; gam=gamma_matrices(N,QI); I=matrix.identity(QI,K)
        g0=I
        for j in range(r):
            a=QQ(j+2); alpha=(a+1/a)/2; beta=(a-1/a)/(2*QI.gen())
            g0=g0*(alpha*I+beta*gam[2*j]*gam[2*j+1])
        R=QQ(3)/5*I+QQ(4)/5*gam[0]*gam[2]
        B=(gam[0]+QI.gen()*gam[1])*gam[2]/2
        assert B*B==0
        for name,g in [('torus',g0),('conjugated_torus',R*g0*R.inverse()),('unipotent',I+B)]:
            gi=g.inverse()
            transformed=[g*T*gi for T in gam]
            a=matrix(QI,N,N,lambda i,j:trace_product(transformed[j],gam[i])/K)
            assert a.transpose()*a==matrix.identity(QI,N) and a.det()==1
            X=(a-matrix.identity(QI,N))*(a+matrix.identity(QI,N)).inverse()
            c=g.trace()/K; assert c*c*(matrix.identity(QI,N)-X).det()==1
            J=pfaffian_cofactor(X)
            for sigma in [1,-1]:
                _,indices,pairs,Bs,_=chiral_data(r,sigma)
                block=g.matrix_from_rows_and_columns(indices,indices); dim=2**(r-1)
                coords=[-QQ(4)/dim*trace_product(block,D) for D in Bs]
                theta=matrix(QI,N,N)
                for val,(i,j) in zip(coords,pairs): theta[i,j]=val; theta[j,i]=-val
                assert theta==2*c*(X+sigma*J),(r,sigma,name)
                rows.append(dict(rank=r,chirality=sigma,point=name,full_chiral_trace_verified=True))
    return rows


def zero_parameters(r,sigma,j):
    M=int(lcm(4,2*(r-2))); A=CyclotomicField(M); w=A.gen()
    t=w**(M//(2*(r-2))*j); u=t*t; eta=-sigma*(-1)**j
    xs=[t]*(r-1)+[eta*t]
    assert u!=1 and u**(r-2)==1 and prod(xs)==-sigma*u
    zs=[(1+x)/(1-x) for x in xs]
    return A,xs,zs,u


def exact_zero_checks():
    rows=[]
    for r in range(4,10):
        for sigma in [1,-1]:
            for j in range(1,r-2):
                A,xs,zs,u=zero_parameters(r,sigma,j)
                es=spin_weights(r,sigma)
                assert moment_without_b(es,zs)==0
                chi_over_b=sum(prod(z**((e-1)//2) for z,e in zip(zs,eps)) for eps in es)
                assert chi_over_b!=0
                c_over_b=prod(1+z for z in zs)/(2**r*prod(zs))
                assert chi_over_b==2**(r-1)*c_over_b*(1-u)
                assert c_over_b**2*prod(zs)*(1-u)**r==1
                rows.append(dict(rank=r,chirality=sigma,root_index=j,
                                 direct_weight_image_zero=True,half_character_nonzero=True,
                                 two_lifts=True,orbit_dimension=r*(r-1)))
    return rows


def exact_zero_differentials():
    rows=[]
    for r in [4,5,6]:
        for sigma in [1,-1]:
            _,_,_,Bs,es=chiral_data(r,sigma); dim=2**(r-1)
            for j in range(1,r-2):
                A,xs,zs,u=zero_parameters(r,sigma,j)
                ii=A.gen()**(int(A.zeta_order())//4)
                def embed(x):
                    co=QI(x).list(); return A(co[0])+(A(co[1])*ii if len(co)>1 else 0)
                Ds=[{ij:embed(val) for ij,val in B.dict().items()} for B in Bs]
                vals=[prod(z**((e-1)//2) for z,e in zip(zs,eps)) for eps in es]
                def entry(out,inp):
                    return -QQ(4)/dim*sum(vals[a]*v*Ds[out].get((b,a),0)
                                         for (a,b),v in Ds[inp].items())
                J=matrix(A,len(Bs),len(Bs),entry)
                assert J.rank()==r*r,(r,sigma,j,J.rank())
                rows.append(dict(rank=r,chirality=sigma,root_index=j,
                                 actual_half_spin_trace_differential_rank=int(J.rank()),
                                 orbit_tangent_dimension=r*(r-1)))
                print('Half-spin zero differential:',r,sigma,j,r*r,flush=True)
    return rows


def boundary_checks():
    rows=[]
    for r in range(4,9):
        zs=[QQ(-1),QQ(-1),QQ(1),QQ(1)]+[QQ(j+2) for j in range(r-4)]
        for sigma in [1,-1]:
            es=spin_weights(r,sigma)
            assert moment_without_b(es,zs)==0
            assert sum(prod(z**((e-1)//2) for z,e in zip(zs,eps)) for eps in es)==0
            rows.append(dict(rank=r,chirality=sigma,direct_boundary_weight_image_zero=True,
                             family_dimension=r*(2*r-1)-12))
    return rows


def kernel_and_low_rank_checks():
    kernels=[]
    for r in range(3,11):
        for sigma in [1,-1]:
            es=spin_weights(r,sigma); b=QI.one() if r%2==0 else QI.gen()
            actions=[QI(1),QI(-1),sigma*b,-sigma*b]
            assert actions.count(1)==(2 if r%2==0 else 1)
            for bb in [b,-b]:
                assert all(bb*prod(QI(-1)**((e-1)//2) for e in eps)==sigma*bb for eps in es)
            kernels.append(dict(rank=r,chirality=sigma,kernel_order=actions.count(1),
                                standard_cover_kernel_acts=-1))
    # D4 half-spin weights form four opposite orthonormal pairs. Their
    # character sublattice has index two in the actual Spin weight lattice.
    triality=[]
    P=matrix(QQ,[[1,0,0,0],[0,1,0,0],[0,0,1,0],[QQ(1)/2]*4])
    for sigma in [1,-1]:
        es=spin_weights(4,sigma); reps=[e for e in es if e[0]==1]
        M=matrix(QQ,[[QQ(x)/2 for x in e] for e in reps])
        assert M*M.transpose()==matrix.identity(QQ,4)
        coords=M*P.inverse(); assert all(x in QQ and x.denominator()==1 for x in coords.list())
        assert abs(coords.det())==2
        triality.append(dict(chirality=sigma,torus_cover_index=2,vector_degree=16,Spin8_degree=32))
    fibers=[]
    for sigma in [1,-1]:
        r=3; es=spin_weights(r,sigma)
        R=PolynomialRing(QQ,['z1','z2','z3','b','v'],order='degrevlex')
        zs=R.gens()[:3]; b,v=R.gens()[3:]; den=prod(zs)
        eqs=[sum(QQ(e[i])/2*b*prod(z**((ei+1)//2) for z,ei in zip(zs,e)) for e in es)
              -(i+2)*den for i in range(3)]
        eqs += [b*b-prod(zs),v*b*prod(zs)-1]
        gb=R.ideal(eqs).groebner_basis()
        length=staircase_length([tuple(f.lm().exponents()[0]) for f in gb])
        assert length==4
        fibers.append(dict(chirality=sigma,Spin6_original_torus_fiber_length=length))
    return kernels,triality,fibers


def main():
    weights=weight_checks(); clifford=Clifford_formula_checks()
    zeros=exact_zero_checks(); differentials=exact_zero_differentials()
    boundary=boundary_checks(); kernels,triality,fibers=kernel_and_low_rank_checks()
    paths=['HALF_SPIN_CAYLEY.md','tools/verify_half_spin_cayley.py',
           'tools/verify_spin_cayley_geometry.py','tools/monomial_staircase.py']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        uniform_Pfaffian_formula=True,complete_chart_zero_fiber_written_proof=True,
        chart_zero_fiber_smooth_written_proof=True,global_boundary_dimension_lower_bound=True,
        highest_weight_moment_checks=len(weights),full_chiral_Clifford_checks=len(clifford),
        exact_zero_weight_checks=len(zeros),full_zero_differential_rank_checks=len(differentials),
        boundary_family_checks=len(boundary),central_kernel_checks=len(kernels),
        Spin10=dict(chart_noncentral_classes_per_half_spin=4,chart_class_dimension=20,
                    differential_rank=25,boundary_zero_family_dimension=33,
                    generic_degree_determined=False,full_zero_fiber_classified=False),
        weight_checks=weights,Clifford_checks=clifford,zero_checks=zeros,
        zero_differentials=differentials,boundary_checks=boundary,kernel_checks=kernels,
        triality_controls=triality,Spin6_fibers=fibers,
        novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
        scope='Finite independent checks corroborate separate uniform proofs. No half-spin generic degree beyond known low-rank controls or maximal finite target locus is asserted.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'half-spin-cayley-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('source_sha256','weight_checks','Clifford_checks','zero_checks','zero_differentials',
         'boundary_checks','kernel_checks')},indent=2))


if __name__=='__main__': main()
