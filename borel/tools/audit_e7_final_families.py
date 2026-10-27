"""Independent native-root/SymPy audit; does not import the Sage verifier."""
from pathlib import Path
from hashlib import sha256
import gzip,json
import sympy as s
from audit_branching_norm_fibers import reconstruct
ROOT=Path(__file__).resolve().parents[1]

def run():
    rows=json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    row=next(v['character'] for v in rows if v['character']['type']=='E7' and v['character']['nodes']==[1,2,3,5,6])
    reconstruct(row)
    a,x,b,y,z,p,q,r,t=s.symbols('a x b y z p q r t')
    original=s.expand(sum(s.Rational(v['coefficient'])*s.prod(w**e for w,e in zip((a,x,b,y,z),v['exponent']))*
                     s.prod(w**e for w,e in zip((p,q),v['charge'])) for v in row['terms']))
    aa=-(s.diff(original,b)-a);bb=-(s.diff(original,a)-b)
    graph={a:aa,b:bb};res=s.expand(original.subs(graph,simultaneous=True))
    H=s.Matrix([[-2*(x+r)/(r*t),-(x+r)*(x+1/r)],[-(x+r)*(x+1/r),-2*r*t*(x+1/r)]])
    ell=s.Matrix([-t*(x*x+r*x+r*r-1)-(r*x+1)/(r*r*t*t),-r*t*t*(x+r)-(x*x+x/r+1/r**2-1)/t])
    v=s.Matrix([y,z]);f=s.expand(-2+(v.T*H*v)[0]/2+(ell.T*v)[0])
    assert s.cancel(res-f.subs({r:p**3,t:q*q}))==0
    for u in (x,y,z,p,q):
        assert s.cancel(s.diff(original,u).subs(graph,simultaneous=True)-s.diff(res,u))==0
    D=s.factor(H.det());P=(x+r)*(x+1/r);k=r*x*x+(r*r+1)*x-3*r
    assert s.cancel(D-P*(4-P))==0
    R=s.cancel(-2-(ell.T*H.adjugate()*ell)[0]/(2*D))
    J=s.Poly(s.cancel(r**6*t**3*D**2*s.diff(R,x)),x,r,t)
    assert s.Poly(J.as_expr(),x).LC()==2*r**6*t**3
    spec=s.Poly(J.as_expr().subs({r:2,t:1}),x)
    assert spec.degree()==9 and spec.gcd(spec.diff()).degree()==0
    # Independent full multivariate squarefree check, in addition to the
    # specialization argument used in the written proof.
    assert all(m==1 for _,m in J.sqf_list()[1])
    for val in (-r,-1/r):
        h=H.subs(x,val).applyfunc(s.factor);l=ell.subs(x,val).applyfunc(s.factor)
        i=0 if val==-r else 1
        assert h[i,0]==h[i,1]==0
        target=-(r*r-1)*(r*r*t**3-1)/(r*r*t**(2-i))
        assert s.cancel(l[i]-target)==0
        assert s.factor(s.diff(D,x).subs(x,val))!=0
    consistency=H.adjugate()*ell
    expected=[2*r*(r*r+1)*(r*x-1)*(r*r*t**3+1),2*(x-r)*(r*r+1)*(r*r*t**3+1)]
    for j in range(2):
        polynomial=s.cancel(r**3*t**(j+1)*consistency[j])
        assert s.rem(s.expand(polynomial-expected[j]),k,x)==0
    assert s.cancel(k.subs(x,-(r*r+1)/(2*r))+(r**4+14*r*r+1)/(4*r))==0
    assert s.gcd(r**4+14*r*r+1,(r*r-1)*(r*r+1))==1
    for eps in (-1,1):
        sub={r:eps,x:-eps}
        assert H.subs(sub)==s.zeros(2) and ell.subs(sub)==s.zeros(2,1)
        assert s.expand(f.subs(sub))==-2 and s.expand((t*s.diff(f,t)).subs(sub))==0
    pt={x:-4,y:-3,z:-3,r:1,t:1}
    assert all(s.diff(f,u).subs(pt)==0 for u in (x,y,z))
    diff=lambda f,u:u*s.diff(f,u) if u in (r,t) else s.diff(f,u)
    M=s.Matrix([[diff(diff(f,u),v).subs(pt) for v in (x,y,z,r,t)] for u in (x,y,z,r,t)])
    S=M[3:,3:]-M[3:,:3]*M[:3,:3].inv()*M[:3,3:]
    assert M[:3,:3].det()==-810 and S.det()==s.Rational(11664,5)
    return dict(root_reconstruction=True,full_reduction_and_chain_rule=True,
                numerator_squarefree_over_Q=True,degree_preserving_squarefree_specialization=True,
                all_three_boundary_divisors=True,rank_zero_target_line=True,
                dominance_schur_determinant=str(S.det()))

def main():
    report=dict(result='passed',checks=run(),scope='Exact characteristic-zero identities and independent root reconstruction; geometric argument is in the proof note, no priority claim.')
    names=['E7_FINAL_FAMILIES.md','tools/audit_e7_final_families.py']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'e7-final-families-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':main()
