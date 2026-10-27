"""Independent exact identities and finite abelian kernels; no Sage imports."""
from pathlib import Path
from hashlib import sha256
from math import comb,prod
import json
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

ROOT=Path(__file__).resolve().parents[1]

def smith(A):
    D=smith_normal_form(A,domain=s.ZZ)
    return [abs(int(D[i,i])) for i in range(min(D.rows,D.cols))]

def main():
    a,v,L=s.symbols('a v L');rows=[]
    for m in [3,4,5,7,8,9,16]:
        p,exponent=next(iter(s.factorint(m).items()))
        vp=lambda j: next(k for k in range(exponent+1) if j % (p**(k+1)))
        r=m-1;epsilon=(-1)**m;F=v**m-a*v**(m-1)-epsilon*a
        rem=lambda x:s.rem(s.expand(x),F,v)
        cs=[s.Integer(comb(m,j))+(comb(m-1,j-1)*a if 1<=j<m else 0) for j in range(m+1)]
        original=sum((-1)**j*cs[j]*L**(m-j) for j in range(m+1))
        assert s.expand(original.subs(L,1+v)-F)==0
        h=lambda j:s.Rational(m,j)*v**j-a*v**(j-1)
        for j in range(1,m):
            fj=sum((-1)**k*comb(m-k,j-k)*cs[k]*L**(j-k) for k in range(j+1))
            unit=comb(m-1,j-1)
            assert unit%p!=0
            assert s.expand(fj.subs(L,1+v)-unit*h(j))==0
        assert rem(h(m))==epsilon*a
        for k in range(1,m-1):
            assert rem(h(1)*h(k+1)-s.Rational(m*(k+2),k+1)*h(k+2)-a*a*v**k)==0
        coefficient=s.Rational(m-2,2*(m-1))
        assert coefficient.p%p and coefficient.q%p
        assert rem(h(2)*h(m-1)-epsilon*s.Rational(m,2*(m-1))*a*h(1)-coefficient*a*a*v**(m-1))==0
        products=0
        for i in range(1,m):
            for j in range(i,m):
                coefficient=s.Rational(m*(i+j),i*j)
                assert coefficient.q%p
                assert s.expand(h(i)*h(j)-coefficient*h(i+j)-a*a*v**(i+j-2))==0
                products+=1
        for j in range(1,m-1):
            assert s.rem(rem(h(m+j))-epsilon*s.Rational(j,m+j)*a*h(j),a*a,a)==0
        # Rescale the cyclic generators by p-local units; only valuations remain.
        w=[exponent-vp(j) for j in range(1,m)]
        orders=[p**w[0]]+[p**(w[j-1]+w[j]) for j in range(1,r)]+[p**w[-1]]
        D=s.diag(*orders);A=s.zeros(r+1)
        for j in range(r-1):A[j+1,j]=p**w[j+1]
        A[r,r-1]=1
        K=s.zeros(r+1,r)
        for j in range(r-1):K[j+1,j]=p**w[j+1]
        K[r,r-1]=1
        for X in (A*D,A*K,K*s.diag(*[p**x for x in w])):
            assert all(int(X[i,j])%orders[i]==0 for i in range(r+1) for j in range(X.cols))
        assert prod(smith(D.row_join(A)))==p**sum(w)
        assert prod(smith(D.row_join(K)))==p**sum(w)
        assert prod(orders)==p**(2*sum(w))
        assert sum(w)==exponent*m-(m-1)//(p-1)
        rows.append({'m':m,'p':p,'exponent':exponent,'product_identities':products,
                     'kernel_order':str(p**sum(w)), 'kernel_factors':sorted(p**x for x in w),
                     'full_image_membership_identities':True})

    # General-module obstruction: the same ordinary exponent p can give p^r.
    controls=[]
    for p in [2,3,5]:
        for r in range(2,9):
            M=p*s.eye(r)
            for i in range(1,r):M[i,i-1]=a
            last=s.zeros(r,1);last[0]=a
            B=M.row_join(last)
            relation=s.Matrix([(-1)**j*a**j*p**(r-j) for j in range(1,r+1)]+[p**r])
            assert B*relation==s.zeros(r,1)
            assert smith(relation.subs(a,0))==[p**r]
            controls.append({'p':p,'rank':r,'ordinary_exponent':p,'defect_exponent':p**r})
    report=dict(result='passed',prime_power_family_checks=rows,general_module_controls=controls,
                scope='Independent exact identities and finite abelian kernels; uniform completeness and quotient transfer are proved in the note.',
                source_sha256={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in ['TYPE_A_K_PRIMARY_LAYERS.md','tools/audit_type_a_k_primary_layers.py']})
    (ROOT/'type-a-k-primary-layers-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'result':'passed','prime_power_families':len(rows),'product_identities':sum(x['product_identities'] for x in rows),'module_controls':len(controls)}))

if __name__=='__main__':main()
