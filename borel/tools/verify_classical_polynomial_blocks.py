"""Independent ambient-root and original-ideal check; run with sage -python.

Does not import the native audit or its reduced A,B,Q equations.
"""
from pathlib import Path
from hashlib import sha256
from itertools import product
from collections import Counter
import json
from sage.all import QQ,CyclotomicField,PolynomialRing,matrix,vector

ROOT=Path(__file__).resolve().parents[1]

def root_weights(kind,n):
    roots=[]
    for i in range(n):
        for j in range(i+1,n):
            for a,b in product([-1,1],repeat=2):
                r=[0]*n;r[i]=a;r[j]=b;roots.append(tuple(r))
    if kind in 'BC':
        for i in range(n):
            for sign in [-1,1]:
                r=[0]*n;r[i]=sign*(1 if kind=='B' else 2);roots.append(tuple(r))
    weights=Counter(roots);weights[tuple([0]*n)]=n
    return weights

def evaluate(weights,z,K):
    n=len(z);chi=K(0);grad=vector(K,[0]*n);M=matrix(QQ,n,n)
    for r,m in weights.items():
        mon=K(m)
        for a,zi in zip(r,z): mon*=zi**a
        chi+=mon
        grad+=mon*vector(K,r)
        for j in range(n):
            for k in range(n): M[j,k]+=m*r[j]*r[k]
    return chi,M.change_ring(K).solve_right(grad)

def full_ideal_checks():
    """Two-sided ideals, including T=0; no reduced A,B,Q input."""
    checks=[]
    for size in (2,3,4):
        R=PolynomialRing(QQ,names=['d','di','T','W','V']+
                         ['e'+str(j) for j in range(1,size)],order='degrevlex')
        d,di,T,W,V=R.gens()[:5]; es=list(R.gens()[5:])
        e=[R(1)]+es+[d]
        U=e[1]+e[size-1]*di
        U2=e[1]**2-2*e[2]+(e[size-1]*di)**2-2*e[size-2]*di
        for eps in (-1,1):
            phi=((U+W)**2+eps*(U2+V))/2
            full=[phi.derivative(x).subs({W:T-U}) for x in es]
            J=R.ideal([d*di-1]+full)
            tau=-eps*T
            if size==2:
                expected=[(1+di**2)*es[0]-tau*(1+di)]
            elif size==3:
                expected=[es[0]-tau-di,es[1]*di-tau-d]
            else:
                expected=[d+1,es[0]-tau,es[2]*di-tau]
            assert J==R.ideal([d*di-1]+expected),(size,eps)
            if size==2:
                assert J+R.ideal(d*d+1)==R.ideal(d*di-1,d*d+1,T)
            checks.append({'size':size,'epsilon':eps,'two_sided_ideal_equality':True})

    # Full residual ideals: intersections, not just pointwise solutions.
    R=PolynomialRing(QQ,names=('L','T','u','v','w'),order='degrevlex')
    L,T,u,v,w=R.gens()
    rows=[
        ('B1',u*u-1,u*u-1,(u,),[
            R.ideal(u,T-L+1),R.ideal(L+1,T-u*u+2)]),
        ('B2',v,u*u-v-1,(u,v),[R.ideal(u,L-1,T-v-1)]),
        ('C1',u,u*u-1,(u,),[R.ideal(2*u+L,2*T-L)]),
        ('D2',u*v,u*u+v*v-2,(u,v),[
            R.ideal(u,v,T-L),R.ideal(L-2,u+v,T-2+v*v),
            R.ideal(L+2,u-v,T+2-v*v)]),
        ('D3',v,u*w-1,(u,v,w),[R.ideal(u,w,L,T-v)])]
    residual=[]
    for label,trace,chi,coords,branches in rows:
        phi=L*trace+chi
        generators=[T-L-trace]+[phi.derivative(x) for x in coords]
        J=R.ideal(generators)
        intersection=branches[0]
        for branch in branches[1:]:intersection=intersection.intersection(branch)
        assert J==intersection,label
        residual.append({'type':label,'full_ideal_equals_branch_intersection':True})
    return checks,residual

def main():
    ideal_checks,residual_checks=full_ideal_checks()
    K=CyclotomicField(12);w=K.gen();ii=w**3;z3=w**4;z6=w**2
    cases=[('B3 singleton residual B2','B',[z6,-1,2],0,[1,0,0]),
           ('D4 singleton residual D3','D',[ii,1,-1,2],0,[1,0,0,0]),
           ('B2 singleton residual B1','B',[z3,2],1,[1,0]),
           ('C2 two-block constant branch','C',[2,z3/2],1,[1,1]),
           ('C5 four-block residual C1','C',[2,-2,ii/2,-ii/2,ii],-1,[1,1,1,1,0]),
           ('D6 four-block residual D2','D',[2,-2,ii/2,-ii/2,1,-1],-2,[1,1,1,1,0,0])]
    rows=[]
    for label,kind,vals,want,pattern in cases:
        z=list(map(K,vals));n=len(z);weights=root_weights(kind,n)
        for r in weights:
            if not any(r):continue
            mon=K(1)
            for a,zi in zip(r,z):mon*=zi**a
            assert mon!=1,(label,r)
        chi,theta=evaluate(weights,z,K)
        assert chi==want,(label,chi)
        assert theta[0]!=0
        assert theta==theta[0]*vector(K,pattern),(label,theta)
        # Direct tensor-character weight convolution, not the chain formula.
        tensor_checked=n<=4
        if tensor_checked:
            tensor=Counter()
            for r,m in weights.items():
                for t,a in weights.items(): tensor[tuple(x+y for x,y in zip(r,t))]+=m*a
            chi2,theta2=evaluate(tensor,z,K)
            dim=sum(weights.values())
            assert chi2==chi*chi
            assert theta2==chi/dim*theta
        rows.append({'case':label,'rank':n,'adjoint_dimension':sum(weights.values()),'chi':str(chi),'theta':list(map(str,theta)),
                     'all_root_regularities_checked':True,'direct_tensor_convolution':tensor_checked})

    # Retain nilpotents at the D5 GL2 SO6 boundary using the full potential.
    R=PolynomialRing(QQ,names=('d','di','x','u','y','v'),order='degrevlex')
    d,di,x,u,y,v=R.gens()
    phi=d+di+x*x*di+(1+di)*x*y+u*v-1
    J=R.ideal([d*di-1]+[phi.derivative(a) for a in (x,u,y,v)])
    local=J.saturation(R.ideal(y))[0]
    assert (d+1)**2 in local and d+1 not in local
    mu=d*phi.derivative(d)-di*phi.derivative(di)
    assert mu*mu in local and mu not in local
    assert phi+3 in local
    # Thus the support maps to zero although the actual scheme has nilpotents.
    report={'result':'passed','date':'2026-10-27','method':'Full ambient root characters and exact cyclotomic arithmetic; independent Singular ideal membership',
            'root_checks':rows,'D5_boundary':{'nilpotent_retained':True,'target_nilpotent':True,'constant_value':-3},
            'full_block_ideal_checks':ideal_checks,'full_residual_ideal_checks':residual_checks,
            'scope':'Finite independent checks of decisive new branches and boundary; general classification requires the written proof',
            'source_sha256':{name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in ['CLASSICAL_POLYNOMIAL_ADJOINT_BLOCKS.md','tools/verify_classical_polynomial_blocks.py']}}
    (ROOT/'classical-polynomial-blocks-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'result':'passed','root_cases':len(rows),'tensor_convolutions':sum(r['direct_tensor_convolution'] for r in rows),'nonreduced_D5_boundary':True,'block_ideal_equalities':len(ideal_checks),'residual_ideal_decompositions':len(residual_checks)}))

if __name__=='__main__':main()
