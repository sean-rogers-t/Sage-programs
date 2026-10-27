"""Independent SymPy identities for the all-rank proof and exact SL3 control."""
from pathlib import Path
from hashlib import sha256
from math import comb
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

ROOT=Path(__file__).resolve().parents[1]
L,a,b=s.symbols('L a b')

def main():
    f=L**3-a*L**2+b*L-1
    rem=lambda x:s.rem(s.expand(x),f,L)
    u=3*L;v=3*L**2-2*a*L;x=a*a*L;y=a*a*L*L
    assert rem(u*v-a*v+3*b*u-9)==2*x
    assert a*a*u==3*x
    assert rem(v*v+3*b*v-3*u+3*a+a*b*u)==y
    gs=[s.Integer(1),u,v,x,y]
    # Solve every pair product, with integer polynomial coefficients.
    closure=0
    for i,g in enumerate(gs):
        for h in gs[i:]:
            z=rem(g*h); c2=s.expand(z).coeff(L,2);c1=s.expand(z).coeff(L,1)
            # c2 lies in (3,a^2), and c1 correction then does too.
            B=s.Poly(c2,a,b); vcoef=s.Integer(0);ycoef=s.Integer(0)
            for (ea,eb),coef in B.terms():
                if coef%3==0:vcoef+=coef//3*a**ea*b**eb
                else:
                    assert ea>=2
                    ycoef+=coef*a**(ea-2)*b**eb
            c1=s.expand(c1+2*a*vcoef)
            ucoef=s.Integer(0);xcoef=s.Integer(0)
            for (ea,eb),coef in s.Poly(c1,a,b).terms():
                if coef%3==0:ucoef+=coef//3*a**ea*b**eb
                else:
                    assert ea>=2
                    xcoef+=coef*a**(ea-2)*b**eb
            const=z.subs(L,0)
            assert s.expand(z-(const+ucoef*u+vcoef*v+xcoef*x+ycoef*y))==0
            closure+=1
    rel=s.Matrix([[0,0],[-a*a,0],[0,-a*a],[3,-2*a],[0,3]])
    inclusion=s.Matrix([[s.expand(g).coeff(L,i) for g in gs] for i in range(3)])
    assert inclusion*rel==s.zeros(3,2)
    diagonal=smith_normal_form(rel.subs(a,3),domain=s.ZZ)
    assert [abs(diagonal[i,i]) for i in range(2)]==[3,3]

    # Universal modular identities; these checks corroborate the residue-class proof.
    modular=[]
    for n in range(3,129):
        c=s.symbols('c1:5');c1,c2,c3,c4=c
        F=lambda j:sum((-1)**k*comb(n-k,j-k)*(1 if k==0 else c[k-1])*L**(j-k) for k in range(j+1))
        for p in s.factorint(n):
            if p!=2 or n%4==0:
                expr=F(2)-c2-c1*L
            else:
                delta=comb(n,4)%2
                expr=F(4)+delta*(F(2)**2+c1*c1*F(2))-(c3+delta*c1**3)*L-(c4+delta*(c2*c2+c1*c1*c2))
            assert s.Poly(s.expand(expr),L,*c,modulus=p).is_zero,(n,p)
            modular.append([n,int(p)])
    # Positive Sp control at augmentation from the original b=L-L^-1 inclusion.
    sp=[]
    h=s.symbols('h')
    for n in range(1,9):
        bclass=s.series((1+h)-1/(1+h),h,0,2*n).removeO()
        coefficients=[s.expand(bclass).coeff(h,i) for i in range(2*n)]
        cols=[[s.Integer(1)]+[s.Integer(0)]*(2*n-1)]
        for j in range(1,2*n):
            cols.append([sum(cols[-1][k]*coefficients[i-k] for k in range(i+1)) for i in range(2*n)])
        M=s.Matrix(2*n,2*n,lambda i,j:cols[j][i])
        assert M.det()==2**(n*(2*n-1))
        sp.append(n)
    report=dict(result='passed',SL3_pair_products=closure,SL3_kernel_factors=[3,3],
                modular_identity_checks=len(modular),rank_range=[3,128],
                symplectic_augmentation_ranks=sp,
                scope='Independent exact identities; all-rank geometric argument is written, not machine-proved.',
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in ['TYPE_A_K_SPECIALIZATION.md','tools/audit_type_a_k_specialization.py']})
    (ROOT/'type-a-k-specialization-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':main()
