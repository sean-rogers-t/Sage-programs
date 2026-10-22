"""Full trace Jacobians and finite-model audits for polynomial Spin maps."""
from hashlib import sha256
from pathlib import Path
from itertools import combinations
import json

from sage.all import QQ, PolynomialRing, QuadraticField, matrix, vector, prod, block_diagonal_matrix
from sage.env import SAGE_VERSION
from verify_spin_polynomial_degree import u, examples, mixed_weights, gram, parity, polynomial

ROOT=Path(__file__).resolve().parents[1]


def expected_jacobian(N,F,c,trace_a):
    r=N//2; K=2**r; d=N*(N-1)//2
    k=F.derivative()(K*c)/F.derivative()(K)
    kp=K*F.derivative(2)(K*c)/F.derivative()(K)
    if N==3: return k*k*(c*k+(c*c-1)*kp)
    return c**((N-1)*(N-4)//2)*k**(d-1)*(k+(k+c*kp)*(trace_a-N)/4)


def weight_jacobians():
    rows=[]
    Fs=[u,u**2,u**3,u+u**2,(u+1)**3,(u**2-1)*(u+1)]
    for N in range(3,11):
        r=N//2; K=2**r
        aa=[QQ(i+2) for i in range(r)]; zs=[a*a for a in aa]; b=prod(aa)
        c=prod(1+z for z in zs)/(K*b)
        trace_a=N%2+sum(z+1/z for z in zs)
        for F in Fs:
            W=mixed_weights(r,F); G=gram(r,W)
            v=vector(QQ,[0]*r); Hess=matrix(QQ,r,r)
            for k,mult in W.items():
                p=parity(k); mu=vector(QQ,[QQ(ki)/2 for ki in k])
                term=mult*b**p*prod(z**((ki-p)//2) for z,ki in zip(zs,k))
                v+=term*mu; Hess+=term*mu.column()*mu.row()
            yy=G.inverse()*v; J=(G.inverse()*Hess).det()
            for i,j in combinations(range(r),2):
                for sign in [-1,1]:
                    root_z=zs[i]*zs[j]**sign
                    root_y=yy[i]+sign*yy[j]
                    J*=root_y**2*root_z/(root_z-1)**2
            if N%2:
                for z,y in zip(zs,yy): J*=y*y*z/(z-1)**2
            assert J==expected_jacobian(N,F,c,trace_a),(N,F)
            rows.append(dict(N=N,polynomial=str(F),dimension=2**r,
                             original_weight_Hessian_and_all_root_directions=True))
    return rows


def gamma_matrices(N,A):
    r=N//2; ii=A.gen(); I=matrix.identity(A,2)
    X=matrix(A,[[0,1],[1,0]]); Y=matrix(A,[[0,-ii],[ii,0]])
    Z=matrix(A,[[1,0],[0,-1]])
    def tensor(factors):
        result=matrix.identity(A,1)
        for f in factors: result=result.tensor_product(f)
        return result
    gammas=[]
    for j in range(r):
        for T in [X,Y]: gammas.append(tensor([Z]*j+[T]+[I]*(r-j-1)))
    if N%2: gammas.append(tensor([Z]*r))
    one=matrix.identity(A,2**r)
    assert all(g*g==one for g in gammas)
    assert all(gammas[i]*gammas[j]+gammas[j]*gammas[i]==0
               for i,j in combinations(range(N),2))
    return gammas


def clifford_boundary_checks():
    A=QuadraticField(-1,'ii'); rows=[]
    for N in range(3,9):
        r=N//2; K=2**r; gam=gamma_matrices(N,A); one=matrix.identity(A,K)
        pairs=list(combinations(range(N),2))
        Bs=[gam[i]*gam[j]/2 for i,j in pairs]
        G=matrix(A,[[(B*C).trace() for C in Bs] for B in Bs])
        assert G==(-QQ(K)/4)*matrix.identity(A,len(Bs))
        for zeros in range(1,min(r,2)+1):
            g=one
            for j in range(r):
                biv=gam[2*j]*gam[2*j+1]
                g=g*(biv if j<zeros else QQ(3)/5*one+QQ(4)/5*biv)
            rotation=QQ(3)/5*one+QQ(4)/5*gam[0]*gam[2]
            g=rotation*g*rotation.inverse()
            c=g.trace()/K; assert c==0
            v=vector(A,[(g*B).trace() for B in Bs]); base=G.inverse()*v
            theta=matrix(A,N,N)
            for val,(i,j) in zip(base,pairs): theta[i,j]=val; theta[j,i]=-val
            assert theta.rank()==(2 if zeros==1 else 0)
            JS=matrix(A,len(Bs),len(Bs),lambda j,i:-QQ(4)/K*(g*Bs[i]*Bs[j]).trace())
            invg=g.inverse(); trace_a=sum((g*T*invg*T).trace()/K for T in gam)
            for F in [u,u**2,u+u**2]:
                kval=F.derivative()(K*c)/F.derivative()(K)
                JF=kval*JS+F.derivative(2)(K*c)/F.derivative()(K)*base.column()*v.row()
                actual=JF.det()
                assert actual==expected_jacobian(N,F,c,trace_a),(N,zeros,F)
                rows.append(dict(N=N,zero_rotation_blocks=zeros,polynomial=str(F),
                                 boundary_image_rank=int((kval*theta).rank()),
                                 full_Clifford_trace_Jacobian=str(actual)))
    return rows


def differential_fiber_checks():
    rows=[]
    for r in range(1,5):
        for label,F in examples():
            ys=[QQ(i+2) for i in range(r)]; P,h=polynomial(r,F,ys)
            A=P.parent().quotient(P,'cb'); c=A.gen(); hc=h(c); hp=h.derivative()(c)
            xs=[y/(2*hc) for y in ys]
            R=sum(x*x/(1-x*x) for x in xs); L=hc+c*hp*R
            assert P.derivative()(c)==2**(2*r+1)*hc**(2*r-1)/c**3*L
            rows.append(dict(rank=r,character=label,ramification_identity_verified=True))
    return rows


def nonfinite_boundary_targets():
    A=QuadraticField(-1,'ii'); ii=A.gen(); rows=[]
    nil=matrix(A,[[0,1,ii],[-1,0,0],[-ii,0,0]])
    assert nil.rank()==2 and nil**3==0
    for N in range(5,13):
        r=N//2
        blocks=[nil] if N%2 else [matrix(A,2,2)]
        blocks += [matrix(A,[[0,j+2],[-j-2,0]]) for j in range(r-1)]
        Y=block_diagonal_matrix(blocks)
        assert Y.nrows()==N and Y+Y.transpose()==0 and Y.rank()>=4
        char=Y.charpoly(); assert char[N%2]==0
        rows.append(dict(N=N,matrix_rank=int(Y.rank()),last_even_characteristic_coefficient_zero=True,
                         odd_zero_block_nonsemisimple=bool(N%2)))
    return rows


def regular_critical_invariants():
    rows=[]
    for r in range(2,9):
        for m in range(1,7):
            ss=[QQ(10+i) for i in range(r-1)]
            ss.append(2*r-QQ(4)/m-sum(ss))
            assert len(set(ss))==r and all(s not in [-2,2] for s in ss)
            c_squared=prod(s+2 for s in ss)/4**r
            ts=[(s-2)/4 for s in ss]
            assert c_squared!=0 and 1+m*sum(ts)==0
            assert all((s-2)/(s+2)!=0 for s in ss)
            rows.append(dict(rank=r,tensor_power=m,pair_traces=[str(s) for s in ss],
                             regular_semisimple=True,critical=True,in_finite_target_locus=True))
    return rows


def discriminant_control():
    A=PolynomialRing(QQ,['p','v']); p,v=A.gens()
    R=PolynomialRing(A,'c'); c=R.gen()
    P=c**4-(p+v+2)*c*c/2+(p-v)**2/16
    disc=P.discriminant()
    assert disc==(p-v)**2*(p+1)**2*(v+1)**2
    # At p=v=0 four distinct product-fiber points remain unramified,
    # even though two of them have the same primitive-element value c=0.
    points=[(x,y) for x in [-1,1] for y in [-1,1]]
    assert len(points)==4 and all(4*x*y!=0 for x,y in points)
    assert len({QQ(x+y)/2 for x,y in points})==3
    return dict(monic_scalar_discriminant=str(disc.factor()),
                genuine_branch_factors=['p+1','v+1'],primitive_element_artifact='p-v',
                unramified_product_fiber_points=4,distinct_c_values=3)


def main():
    weight=weight_jacobians(); clifford=clifford_boundary_checks()
    fiber=differential_fiber_checks(); boundary=nonfinite_boundary_targets()
    discrim=discriminant_control(); regular=regular_critical_invariants()
    paths=['SPIN_CAYLEY_GEOMETRY.md','tools/verify_spin_cayley_geometry.py',
           'SPIN_POLYNOMIAL_CAYLEY_DEGREES.md','tools/verify_spin_polynomial_degree.py',
           'SPIN_CAYLEY_DEGREES.md','tools/verify_spin_cayley_degree.py']
    report=dict(result='passed',date='2026-10-22',sage_version=SAGE_VERSION,
        maximal_finite_locus_uniform_written_proof=True,full_global_Jacobian_written_proof=True,
        weight_root_Jacobian_checks=len(weight),Clifford_boundary_Jacobian_checks=len(clifford),
        finite_algebra_ramification_checks=len(fiber),nonfinite_boundary_target_checks=len(boundary),
        regular_semisimple_critical_invariant_checks=len(regular),
        weight_root_checks=weight,Clifford_checks=clifford,ramification_checks=fiber,
        boundary_targets=boundary,Spin4_discriminant_control=discrim,
        regular_critical_invariants=regular,
        novelty_claim=False,externally_reviewed=False,runtime_API_changed=False,
        scope='Independent finite trace/weight corroboration; the uniform finiteness theorem uses a separate boundary-rank and integral-extension proof. No complete exceptional-fiber classification.',
        source_sha256={p:sha256((ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in paths})
    (ROOT/'spin-cayley-geometry-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
        ('source_sha256','weight_root_checks','Clifford_checks','ramification_checks','boundary_targets',
         'regular_critical_invariants')},indent=2))


if __name__=='__main__': main()
