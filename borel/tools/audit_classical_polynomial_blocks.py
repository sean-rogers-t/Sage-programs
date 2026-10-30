"""Native characteristic-zero checks; no Sage or previous audit imports."""
from pathlib import Path
from hashlib import sha256
import json
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
def zero(x):
    assert S.cancel(x)==0, S.factor(x)

def companion(coeff):
    n=len(coeff)-1
    M=S.zeros(n)
    for j in range(n-1): M[j+1,j]=1
    for i in range(n): M[i,n-1]=-coeff[n-i]
    return M

def main():
    d,T,W,V,q=S.symbols('d T W V q',nonzero=True)
    block_checks=[]
    for size in range(1,5):
        es=S.symbols('e1:'+str(size))
        coeff=[1]+[(-1)**j*es[j-1] for j in range(1,size)]+[(-1)**size*d]
        M=companion(coeff); Mi=M.inv()
        U=S.trace(M)+S.trace(Mi)
        U2=S.trace(M*M)+S.trace(Mi*Mi)
        ss=d+1/d
        for eps in [-1,1]:
            tau=-eps*T
            phi=((U+W)**2+eps*(U2+V))/2
            derivatives=[S.diff(phi,e).subs(W,T-U) for e in es]
            if size==1:
                subs={}; expectedU=ss; expectedU2=ss**2-2
            elif size==2:
                subs={es[0]:tau*(1+1/d)/(1+1/d**2)}
                F=(d+1)**2/(d*d+1)
                expectedU=tau*F;expectedU2=T*T*F-2*ss
            elif size==3:
                subs={es[0]:tau+1/d,es[1]:d*(tau+d)}
                expectedU=2*tau+ss;expectedU2=2*T*T-(ss**2-2)
            else:
                subs={d:-1,es[0]:tau,es[2]:-tau}
                expectedU=2*tau;expectedU2=2*T*T
            for v in derivatives: zero(v.subs(subs,simultaneous=True))
            zero(U.subs(subs,simultaneous=True)-expectedU)
            zero(U2.subs(subs,simultaneous=True)-expectedU2)
            block_checks.append({'size':size,'epsilon':eps,'matrix_inverse_and_critical_equations':True})

    # Symbolic all-q values, computed independently from trace and square trace.
    triples={('B',0):[(0,1,1)],('B',1):[(0,-1,3),(1,1,-1)],
             ('B',2):[(1,-1,1)],('C',0):[(0,0,0)],('C',1):[(1,0,-2)],
             ('D',0):[(0,0,0)],('D',2):[(0,0,4),(1,2,0),(1,-2,0)],('D',3):[(1,0,2)]}
    expected={('B',0):[-q/(2*q-1)],('B',1):[-(3*q-1)/(2*q-1),(2*q-1)/(4*q)],
              ('B',2):[-(2*q+1)/(4*q)],('C',0):[0],('C',1):[-1],
              ('D',0):[0],('D',2):[-2,-1/q,-1/q],('D',3):[-1]}
    constants=[]
    for (kind,m),branches in triples.items():
        eps=1 if kind=='C' else -1
        for index,(c,b,d0) in enumerate(branches):
            t=S.Rational(b)/(1+eps*(2*q+c))
            square_trace=(2*q+c)*t*t+d0
            value=S.factor((t*t+eps*square_trace)/2)
            zero(value-expected[kind,m][index])
            constants.append({'type':kind,'m':m,'triple':[c,b,d0],'T':str(S.factor(t)),'value':str(value)})

    # General reduction with independent trace contributions and free block counts.
    k,p,c,b,d0,eps=S.symbols('k p c b d0 eps')
    s1,s3,s2=S.symbols('s1 s3 s2',nonzero=True)
    for ep in [-1,1]:
        F=1+2/s2
        A=1+ep*(c+2*q+2+F)
        B=b+s1+s3
        Q=d0+(s1*s1-2)-(s3*s3-2)-2*s2
        trace=-ep*(c+2*q+2+F)*T+B
        square=(c+2*q+2+F)*T*T+Q
        zero((T-trace)-(A*T-B))
        zero((T*T+ep*square)/2-(A*T*T+ep*Q)/2)
        zero(((T*T+ep*square)/2).subs(T,B/A)-(B*B/(2*A)+ep*Q/2))

    # All-q constant singular branches with one two-block.
    singular=[]
    for kind,m,c,d0,want in [('C',0,0,0,1/(q+1)),('C',1,1,-2,-(2*q+1)/(2*q+3)),
                           ('D',0,0,0,-1/q),('D',2,0,4,-2-1/q),('D',3,1,2,-1-2/(2*q+1))]:
        ep=1 if kind=='C' else -1
        sval=-S.Rational(2)*ep/(1+ep*(c+2*q+1))
        value=ep*S.Rational(d0,2)-ep*sval
        zero(value-want)
        singular.append({'type':kind,'m':m,'s':str(S.factor(sval)),'value':str(S.factor(value))})

    # Original D5 GL2 x Spin6 potential, with Spin6=SL4 characters u,y,v.
    x,u,y,v=S.symbols('x u y v')
    phi=d+1/d+x*x/d+(1+1/d)*x*y+u*v-1
    eliminate={u:0,v:0,x:-(d+1)*y/2}
    zero(S.diff(phi,x).subs(eliminate))
    zero(S.diff(phi,y).subs(eliminate)+(d+1)**2*y/(2*d))
    reduced=S.factor(phi.subs(eliminate))
    zero(reduced-(d+1/d-1-(d+1)**2*y*y/(4*d)))
    zero(reduced.subs(d,-1)+3)
    zero((d*S.diff(phi,d)).subs({d:-1,x:0,u:0,v:0}))

    report={'result':'passed','date':'2026-10-27','arithmetic':'SymPy exact characteristic zero',
            'block_companion_checks':block_checks,'four_block_constants':constants,
            'one_two_block_constants':singular,'D5_double_zero_boundary':True,
            'scope':'Symbolic identities in arbitrary q plus original-potential boundary; not formal verification of geometric proof',
            'source_sha256':{name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in ['CLASSICAL_POLYNOMIAL_ADJOINT_BLOCKS.md','tools/audit_classical_polynomial_blocks.py']}}
    (ROOT/'classical-polynomial-blocks-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'result':'passed','block_checks':len(block_checks),'symbolic_four_block_values':len(constants),'singular_values':len(singular),'D5_boundary':True}))

if __name__=='__main__': main()
