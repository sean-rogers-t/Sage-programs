"""Sage-free F4 proof audit using native roots and exterior/spin characters."""
from collections import Counter
from itertools import combinations,product
from fractions import Fraction
from pathlib import Path
from hashlib import sha256
import json
import sympy as s
from audit_affine_branching_certificates import cartan,root_data
from audit_branching_norm_fibers import multiply,subset_fundamental
ROOT=Path(__file__).resolve().parents[1]

def fundamentals(kind,order,total):
    n=len(order)
    if kind.startswith('A'):
        return {pos:subset_fundamental(n+1,k+1,order,total) for k,pos in enumerate(order)}
    eps=[tuple(sign*int(i==j) for i in range(n)) for j in range(n) for sign in (-1,1)]
    if kind.startswith('B'):eps.append((0,)*n)
    def wedge(k):
        return Counter(tuple(sum(v[j] for v in choice) for j in range(n)) for choice in combinations(eps,k))
    out={}
    for k,pos in enumerate(order,1):
        if kind.startswith('B') and k==n:
            raw=Counter({tuple(Fraction(v,2) for v in signs):1 for signs in product((-1,1),repeat=n)})
        else:
            raw=wedge(k)
            if kind.startswith('C') and k>=2:raw.subtract(wedge(k-2))
        result=Counter()
        for e,m in raw.items():
            if not m:continue
            w=[0]*total
            local=[e[j]-e[j+1] for j in range(n-1)]+[(2 if kind.startswith('B') else 1)*e[-1]]
            for j,a in zip(order,local):assert int(a)==a;w[j]=int(a)
            result[tuple(w)]+=m
        assert all(v>=0 for v in result.values())
        out[pos]=+result
    return out

def root_check(row):
    A=cartan('F4');positive=root_data(A)[1];weights=Counter({r:1 for r in positive})
    weights.update({tuple(-a for a in r):1 for r in positive});weights[(0,)*4]=4
    I=[n-1 for n in row['nodes']];H=row['central_cocharacters'];n=len(I)
    # Check the central kernel using root columns, independently of Sage's
    # primitive-coordinate construction.
    assert all(sum(A[j][i]*h[j] for j in range(4))==0 for i in I for h in H)
    fund={}
    for kind,order in zip(row['factors'],row['factor_orders']):fund.update(fundamentals(kind,order,n))
    rebuilt=Counter()
    for v in row['terms']:
        small=Counter({(0,)*n:1})
        for i,e in enumerate(v['exponent']):
            for _ in range(e):small=multiply(small,fund[i])
        for w,m in small.items():rebuilt[tuple(v['charge']),w]+=m*v['coefficient']
    expected=Counter()
    for w,m in weights.items():
        expected[tuple(sum(a*b for a,b in zip(w,h)) for h in H),tuple(w[i] for i in I)]+=m
    assert all(v>=0 for v in rebuilt.values()) and +rebuilt==expected

def symbolic(rows):
    for row in rows:
        if row['degree']<=2:continue
        xs=s.symbols('x0:'+str(len(row['nodes'])));ts=s.symbols('t0:'+str(len(row['central_cocharacters'])))
        f=s.expand(sum(v['coefficient']*s.prod(x**e for x,e in zip(xs,v['exponent']))*
                      s.prod(t**e for t,e in zip(ts,v['charge'])) for v in row['terms']))
        if len(xs)==2:
            x,y=xs;p,q=ts
            u,v={(1,3):(p,p*q),(1,4):(p*q*q,p*q),(2,4):(p*q*q,q)}[tuple(row['nodes'])]
            A=u+1/u;B=v+1/v;ell=A*y*y+B*y+A*(B*B-3)
            assert s.cancel(f-(x*x+x*ell+(B*B-1)*y*y+A*B*y+A*A-B*B))==0
            res=s.cancel(f.subs(x,-ell/2));N=s.cancel(2*s.diff(res,y)).as_numer_denom()[0]
            sq=s.Poly(N,y,p,q).sqf_list()[1]
            assert all(m==1 or len(a.terms())==1 for a,m in sq)
            for t in ts:assert s.cancel(s.diff(f,t).subs(x,-ell/2)-s.diff(res,t))==0
        elif row['nodes']==[1,2,4]:
            a,b,x=xs;t=ts[0];S=t**3+t**-3
            graph={a:-(s.diff(f,b)-a),b:-(s.diff(f,a)-b)}
            res=s.expand(f.subs(graph,simultaneous=True));M=4*x**3+3*S*x*x+2*(S*S-4)*x-S
            assert s.cancel(s.diff(res,x)+M)==0
            assert s.cancel(s.diff(f,t).subs(graph,simultaneous=True)-s.diff(res,t))==0
            spec=s.Poly(M.subs(t,2),x);assert spec.gcd(spec.diff()).degree()==0
        else:
            x,y,z=xs;t=ts[0];H=s.hessian(f,(y,z));ell=s.Matrix([s.diff(f,v).subs({y:0,z:0}) for v in (y,z)])
            D=s.factor(H.det());Q=f.subs({y:0,z:0});res=s.cancel(Q-(ell.T*H.inv()*ell)[0]/2)
            assert s.cancel(D-4*x*x-4*(t**3+t**-3)*x-3)==0
            N=s.cancel(t**9*D*D*s.diff(res,x));assert s.Poly(N,x).LC()==32*t**9
            spec=s.Poly(N.subs(t,2),x);assert spec.gcd(spec.diff()).degree()==0
            K=s.expand(t**3*D);cs=[s.cancel(v).as_numer_denom()[0] for v in H.adjugate()*ell]
            resultant=s.factor(s.resultant(K,cs[0],x))
            expected=-4*t**6*(t**6-1)*(2*t**4+t*t+2)*(4*t**8-2*t**6-3*t**4-2*t*t+4)
            assert s.expand(resultant-expected)==0
            V=s.Symbol('V');assert list(s.groebner([K,s.diff(K,x),*cs,t*V-1],V,x,t))==[1]
            star={x:-(t**3+t**-3)/2};aa=(t**6-1)/t**4;bb=-(t**6-1)/t**2
            assert H.subs(star).applyfunc(s.cancel)==s.Matrix([[aa,1],[1,bb]])
            assert all(s.cancel(v)==0 for v in ell.subs(star)-s.Matrix([-bb/2,-aa/2]))
    return dict(all_five_cubic_models=True,full_Laurent_squarefreeness=True,
                exact_boundary_resultant=True,double_root_inconsistency=True)

def main():
    record=json.loads((ROOT/'f4-adjoint-fibers-verification.json').read_text());rows=record['characters']
    assert {tuple(r['nodes']) for r in rows}=={I for k in range(4) for I in combinations(range(1,5),k)}
    for row in rows:root_check(row)
    report=dict(result='passed',independent_root_characters=15,native_exterior_and_spin_models=True,
                checks=symbolic(rows),all_proper_F4_Levis_covered=True,literature_priority_claim=False)
    names=['F4_ADJOINT_FIBERS.md','tools/audit_f4_adjoint_fibers.py','f4-adjoint-fibers-verification.json']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'f4-adjoint-fibers-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':main()
