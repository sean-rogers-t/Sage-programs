"""Native lattice, formal trace, and quadratic-coordinate audit."""
from hashlib import sha256
from math import comb
from pathlib import Path
import json
import sympy as s
from audit_type_a_k_first_correction import native_image

ROOT=Path(__file__).resolve().parents[1]


def quadratic(p):
    y,z,v,w=s.symbols('y z v w');d=p+3;r=(p+1)//2
    def even(f):
        f=s.Poly(s.expand((f+f.xreplace({v:-v}))/2),y,z,v)
        return s.Poly(sum(c*y**a*z**b*w**(e//2) for (a,b,e),c in f.terms()),y,z,w,modulus=p)
    def reduce(f):
        f=s.Poly(f,y,z,w,modulus=p)
        for i in range(1,d):
            c=f.coeff_monomial(y**(d-i)*z**i)
            f-=c*even((z+v)**i*(y-v)**(d-i))
        return s.Poly(sum(c*y**a*z**b*w**e for (a,b,e),c in f.terms() if a and b),y,z,w,modulus=p)
    ap=(y+v)**p
    xA=even(ap*(y-v)**2*(z-y))
    xC=even(((y+v)**p+(y-v)**p)/2*(z-y)**3
            +((y+v)*(y-v)**p+(y-v)*(y+v)**p)*(z-y)**2
            +ap*(y-v)**2*(z-y))
    xW=even((y+v)*(y-v)**2*(z-y)**p)
    actual=[reduce(f) for f in [xA,xC,xW]]
    targets=[-4*y*z*w**r,-4*y**p*z*w,-4*y*z**p*w]
    assert all(a==s.Poly(b,y,z,w,modulus=p) for a,b in zip(actual,targets))
    return dict(p=p,reductions=[str(a.as_expr()) for a in actual],distinct_unit_basis_monomials=True)


def trace(p):
    M,L=s.symbols('M L');d=p+3
    a=[M]+list(s.symbols('a1:'+str(d+1)))
    Q=[M]+list(s.symbols('Q1:'+str(d+1)))
    q=[M]+[M*L**j+sum(comb(j,i)*a[i]*L**(j-i) for i in range(1,j+1)) for j in range(1,d+1)]
    formal=sum(comb(3,i)*comb(p,j)*a[i]*a[j]*Q[d-i-j]
               for i in range(4) for j in range(p+1))-M**2*Q[d]
    for t in range(1,d+1):
        formal-=a[t]*sum(comb(3,i)*comb(p,t-i)*Q[3-i]*Q[p-t+i]
                         for i in range(4) if 0<=t-i<=p)
    assert s.expand(formal.subs(dict(zip(Q[1:],q[1:])),simultaneous=True)
                    -M*(q[3]*q[p]-M*q[d]))==0
    remaining=3*(a[1]*Q[2]+a[2]*Q[1])*(a[p]-Q[p])
    remainder=s.Poly(s.expand(formal-remaining).subs(M,0),*(a[1:]+Q[1:]))
    assert all(int(c)%p==0 for c in remainder.coeffs())
    return dict(p=p,variable_mass_exact_trace=True,formal_reduction_mod_p_and_mass=True)


def main():
    cases=[]
    for m,p in [(10,5),(25,5),(14,7)]:
        row=native_image(m,p,offset=3)
        assert row['primary']==([p if m!=25 else 25]*(4 if p==5 else 7))
        cases.append(row);print(json.dumps(row),flush=True)
    roots=[]
    for p in [5,7]:
        row=native_image(p,p,offset=3)
        assert row['primary']==[p]*(2 if p==5 else 5)
        roots.append(row);print(json.dumps(dict(root_boundary=row)),flush=True)
    quadratic_checks=[quadratic(p) for p in [5,7,11]]
    traces=[trace(p) for p in [5,7,11]]
    M,L,a1,a2,a3=s.symbols('M L a1 a2 a3')
    q1=M*L+a1;q2=M*L**2+2*a1*L+a2;q3=M*L**3+3*a1*L**2+3*a2*L+a3
    assert s.expand(q1*q2-M*q3+2*a2*q1-3*a1*a2+M*a3-2*a1**2*L)==0
    names=['TYPE_A_K_NEXT_ELEMENTARY.md','TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md',
           'TYPE_A_K_FIRST_CORRECTION.md','tools/audit_type_a_k_next_elementary.py',
           'tools/audit_type_a_k_first_correction.py']
    report=dict(result='passed',native_elementary_cases=cases,native_root_cases=roots,
                quadratic_reductions=quadratic_checks,
                trace_checks=traces,sage_imported=False,root_boundary_theorem=True,
                variable_mass_literal_root_membership_identity=True,
                all_rank_claim_from_finite_checks=False,
                source_sha256={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names})
    (ROOT/'type-a-k-next-elementary-audit.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
