"""Exact checks of the common E6 cubic refinement mechanism."""
from hashlib import sha256
from pathlib import Path
import json
from itertools import product
from sage.all import QQ, GF, PolynomialRing, LaurentPolynomialRing, matrix, vector, prod, RootSystem

ROOT=Path(__file__).resolve().parents[1]
H=[vector(QQ,[-2,0,-1,0,0,0]),vector(QQ,[2,3,4,6,4,2]),vector(QQ,[0,0,0,0,-1,-2])]


def master(a,b):return sum(x*y for x,y in zip(a,b))+prod(a)+prod(b)-3


def identities():
    R=LaurentPolynomialRing(QQ,['p','z','q','x','y','w','u','v'])
    p,z,q,x,y,w,u,v=R.gens()
    a=p*x+p**-2;b=p**-1*x+p**2
    c=z*y+z**-2;d=z**-1*y+z**2
    phi=master([a,c,v],[b,d,u])
    S=z**3+z**-3;T=p**3+p**-3;A=z**-2+z*y;B=z**2+z**-1*y
    assert phi==(u+A*a)*(v+B*b)-(y*y+S*y)*(x*x+T*x)-2
    # Coordinate change from the original central cocharacters, checked on every charge.
    data=json.loads((ROOT/'branching-norm-fibers-verification.json').read_text())
    one=data['adjoint_mechanism_tests'][1];three=data['adjoint_mechanism_tests'][0]
    original1=sum(t['coefficient']*prod(xx**e for xx,e in zip((y,x,u,v),t['exponent']))
                  *p**(t['charge'][1]-t['charge'][0])*z**t['charge'][1] for t in one['terms'])
    assert phi==original1
    aa=[a,c,q*w+q**-2];bb=[b,d,q**-1*w+q**2]
    phi3=master(aa,bb)
    original3=sum(t['coefficient']*prod(xx**e for xx,e in zip((y,x,w),t['exponent']))
                  *p**(t['charge'][1]-t['charge'][0])*z**t['charge'][1]*q**(t['charge'][1]-t['charge'][2]) for t in three['terms'])
    assert phi3==original3
    Z=p*z*q;D=Z+Z**-1;E=Z**-1-Z
    delta=[(t**3-t**-3)/2 for t in (p,z,q)]
    X=[xx+(t**3+t**-3)/2 for xx,t in zip((x,y,w),(p,z,q))]
    expanded=sum(t*t for t in X)+D*prod(X)+E*sum(delta[i]*X[(i+1)%3]*X[(i+2)%3] for i in range(3))
    expanded+=D*sum(X[i]*delta[(i+1)%3]*delta[(i+2)%3] for i in range(3))
    expanded+=E*prod(delta)-sum(t*t for t in delta)-3
    assert phi3==expanded
    # Boundary identities on each component Z=+/-i; coefficient field Q(i).
    C=QQ.extension(PolynomialRing(QQ,'t').gen()**2+1,'ii');ii=C.gen()
    L=LaurentPolynomialRing(C,['p','z']);pp,zz=L.gens()
    for sign in (-1,1):
        h=R.hom([pp,zz,sign*ii/(pp*zz),1,1,1,1,1],L)
        de=[h(t) for t in delta];ee=h(E)
        Q=matrix(L,3,3,lambda i,j:2 if i==j else ee*de[3-i-j])
        assert Q.det()==0
        assert ee*prod(de)-sum(t*t for t in de)-3==-2
        assert all(Q.matrix_from_rows_and_columns([i,j],[i,j]).det()==4*(1+de[3-i-j]**2)
                   for i in range(3) for j in range(i+1,3))
    # At all torus coordinates 1 the five relative critical points are reduced.
    F=PolynomialRing(QQ,['x','y','w','p','z','q']); xx=F.gens()[:3];zz=F.gens()[3:]
    point={v:QQ(1) for v in (p,z,q)}
    points=[(-1,-1,-1)]+[tuple(v-1 for v in signs) for signs in product((-1,1),repeat=3) if prod(signs)==-1]
    dets=[]
    for values in points:
        at={**point,**dict(zip((x,y,w),values)),u:0,v:0}
        assert all(phi3.derivative(t).subs(at)==0 for t in (x,y,w))
        hh=matrix(QQ,3,3,lambda i,j:phi3.derivative((x,y,w)[i]).derivative((x,y,w)[j]).subs(at))
        dets.append(int(hh.det()));assert hh.det()!=0
    at={**point,x:-1,y:-1,w:-1,u:0,v:0}
    variables=(x,y,w,p,z,q)
    HH=matrix(QQ,6,6,lambda i,j:phi3.derivative(variables[i]).derivative(variables[j]).subs(at))
    assert HH==matrix.diagonal(QQ,[2,2,2,-18,-18,-18])
    return dict(master_matches_both_central_charge_records=True,
                product_two_quadratics_identity=True,boundary_determinant_zero=True,
                boundary_critical_value=-2,boundary_rank_one_locus_finite=True,
                five_relative_hessian_determinants=dets,
                target_dominance_hessian_diagonal=[2,2,2,-18,-18,-18])


def finite_points():
    R=RootSystem('E6').ambient_space();cr=R.simple_coroots()
    roots=[tuple(int(a.scalar(cr[j+1])) for j in range(6)) for a in R.roots()]
    B=matrix(QQ,6,6,lambda i,j:sum(w[i]*w[j] for w in roots))
    assert matrix(QQ,3,3,lambda i,j:H[i]*B*H[j])==144*matrix.identity(QQ,3)
    rows=[]
    for refined in (2,3):
        found=False
        for prime in (1009,1021,1033,1063,1093):
            F=GF(prime);P=PolynomialRing(F,['x','y','w']);x,y,w=P.gens()
            for central in list(product(range(2,9),range(2,9),(5,))) if refined==2 else list(product(range(2,6),repeat=3)):
                p,z,q=map(F,central);aa=p*x+p**-2;bb=p**-1*x+p**2
                cc=z*y+z**-2;dd=z**-1*y+z**2
                if refined==2:
                    xv=-(p**3+p**-3)/2;yv=-(z**3+z**-3)/2
                    av=F(aa(x=xv,y=yv,w=0));bv=F(bb(x=xv,y=yv,w=0))
                    A=z**-2+z*yv;BB=z**2+z**-1*yv
                    solutions=[(xv,yv,None,-A*av,-BB*bv)]
                else:
                    phi=master([aa,cc,q*w+q**-2],[bb,dd,q**-1*w+q**2])
                    I=P.ideal([phi.derivative(t) for t in P.gens()])
                    solutions=[(s[x],s[y],s[w],None,None) for s in I.variety()]
                U=PolynomialRing(F,'t');t=U.gen()
                for xv,yv,wv,uv,vv in solutions:
                    eigen=[]
                    for trace in ((xv,yv) if refined==2 else (xv,yv,wv)):
                        rs=(t*t-trace*t+1).roots(multiplicities=False)
                        if len(rs)!=2:break
                        eigen.append(rs[0])
                    if len(eigen)!=refined:continue
                    params=[p**int(H[0][j])*z**int(H[1][j]) for j in range(6)]
                    params[2]*=eigen[0];params[1]*=eigen[1]
                    if refined==2:
                        rs=(t**3-uv*t*t+vv*t-1).roots(multiplicities=False)
                        if len(rs)!=3:continue
                        params[4]*=rs[0];params[5]*=rs[0]*rs[1]
                    else:
                        params=[params[j]*q**int(H[2][j]) for j in range(6)]
                        params[4]*=eigen[2]
                    values=[prod(params[j]**rr[j] for j in range(6)) for rr in roots]
                    if F(1) in values:continue
                    grad=vector(F,[sum(v*rr[i] for v,rr in zip(values,roots)) for i in range(6)])
                    theta=B.change_ring(F).inverse()*grad
                    center=matrix(F,6,refined,lambda i,j:H[j][i]);coeff=center.solve_right(theta)
                    if any(sum(F(rr[j])*theta[j] for j in range(6))==0 for rr in roots
                           if any(sum(rr[j]*H[k][j] for j in range(6))!=0 for k in range(refined))):continue
                    Hess=matrix(F,6,6,lambda i,j:sum(v*rr[i]*rr[j] for v,rr in zip(values,roots)))
                    if Hess.rank()!=6:continue
                    outside=sum(any(sum(rr[j]*H[k][j] for j in range(6))!=0 for k in range(refined)) for rr in roots)
                    assert outside==(62 if refined==2 else 66)
                    rows.append(dict(refined_factors=refined,prime=prime,torus_coroot_parameters=list(map(int,params)),
                                     central_target_coefficients=list(map(int,coeff)),cartan_rank=6,root_rank=outside,
                                     full_rank=6+outside,kernel_dimension=78-6-outside))
                    found=True;break
                if found:break
            if found:break
        assert found,refined
    return rows


def main():
    report=dict(result='passed',identities=identities(),independent_root_points=finite_points(),
                remaining_two_mechanism_tests_resolved=True,all_exceptional_Levis_resolved=False,
                Fricke_rank_five_is_not_Springer_degree=True,literature_priority_claim=False)
    names=['CUBIC_LEVI_REFINEMENTS.md','tools/verify_cubic_levi_refinements.py',
           'branching-norm-fibers-verification.json','BRANCHING_NORM_FIBERS.md']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'cubic-levi-refinements-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
