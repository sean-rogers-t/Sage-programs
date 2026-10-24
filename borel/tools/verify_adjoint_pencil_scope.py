"""Diagram obstruction, exact multiplier transfer and the E7 central target fiber."""
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
from sage.all import (QQ, GF, CartanType, RootSystem, PolynomialRing, matrix,
                      vector, prod)
from verify_branching_norm_fibers import charge_polynomial

ROOT=Path(__file__).resolve().parents[1]


def diagram_candidates(kind):
    A=CartanType(kind).cartan_matrix();n=A.nrows();accepted=[]
    for mask in range(1,2**n-1):
        I={i for i in range(n) if mask>>i&1};left=set(I);parts=[]
        while left:
            first=min(left);left.remove(first);part={first};todo=[first]
            while todo:
                i=todo.pop();new={j for j in left if A[i,j]}
                part.update(new);left-=new;todo.extend(new)
            parts.append(part)
        singles=[p for p in parts if len(p)==1]
        if len(singles)!=1 or len(parts)<2:continue
        if any(sum(A[i,j]!=0 for j in p if j!=i)>2 for p in parts for i in p):continue
        leaves=[p for p in parts if len(p)>1]
        if any(sum(any(A[j,i] for i in p) for p in leaves) not in (0,2)
               for j in range(n) if j not in I):continue
        accepted.append(sorted(i+1 for i in I))
    assert accepted==[[i for i in range(1,n+1) if i!=4]]
    return dict(type=kind,necessary_incidence_candidates=accepted)


def e7_identities(row):
    F=PolynomialRing(QQ,'z').fraction_field();z=F.gen()
    P=PolynomialRing(F,['a','b','x','u','v','w']);a,b,x,u,v,w=P.gens()
    # u is the A3 dual defining character (node 7), w is defining (node 5).
    local=(a,x,b,w,v,u)
    phi=sum(t['coefficient']*z**t['charge'][0]*prod(s**e for s,e in zip(local,t['exponent'])) for t in row['terms'])
    formula=a*b+u*w+x*x-2+x*(z*a*u+z**-1*b*w)+v*(z*z*b+z**-2*a)+x*(z**3*w+z**-3*u)+z**4*a+z**-4*b
    assert phi==formula
    assert phi.derivative(v)==z**2*b+z**-2*a
    assert phi.derivative(a)==b+x*z*u+z**-2*v+z**4
    reduced=phi.subs({a:-z**4*b,v:0})
    up=u+x*(z**-1*b+z**3);wp=w+x*(-z**5*b+z**-3)
    S=z**6-z**-6;y=z*z*b
    assert reduced==up*wp+(x*x-1)*(y*y+S*y)-2
    # All original derivatives, including the original z derivative.
    values={x:0,b:-S/(2*z*z),a:z*z*S/2,u:0,w:0,v:-(z**6+z**-6)/2}
    assert all(phi.derivative(t).subs(values)==0 for t in P.gens())
    # Differentiate Laurent coefficient functions before substituting the critical point.
    central=sum(coeff.derivative(z)*prod(t**e for t,e in zip(P.gens(),exponent))
                for exponent,coeff in phi.dict().items())*z
    assert central.subs(values)==3*(z**12-z**-12)
    # Universal exact transfer with independent S and D(S), on c!=0.
    R=PolynomialRing(QQ,['x','y','S','D','c','h']);xx,yy,ss,dd,c,h=R.gens()
    f=(xx*xx-1)*(yy*yy+ss*yy)
    I=R.ideal([f.derivative(xx),f.derivative(yy),(xx*xx-1)*yy*dd-c,c*h-1])
    J=R.ideal([xx,2*yy+ss,ss*dd-2*c,c*h-1])
    assert I==J
    return dict(original_character_matches_formula=True,multiplier_constraint_exact=True,
                completed_product_identity=True,original_derivatives_checked=True,
                universal_nonzero_target_ideal_equality=True,scalar_target='3*(z^12-z^-12)',
                actual_quotient='3*w^2-c*w-3',exceptional_nonzero_targets='c^2+36=0')


def root_point(row):
    p=1033;F=GF(p);R=PolynomialRing(F,'t');t=R.gen();H=row['central_cocharacters'][0]
    ambient=RootSystem('E7').ambient_space();cr=ambient.simple_coroots()
    roots=[tuple(int(r.scalar(cr[j+1])) for j in range(7)) for r in ambient.roots()]
    B=matrix(QQ,7,7,lambda i,j:sum(r[i]*r[j] for r in roots))
    assert vector(QQ,H)*B*vector(QQ,H)==432
    roots2=(t*t+1).roots(multiplicities=False);assert len(roots2)==2
    for zi in range(2,300):
        z=F(zi);S=z**6-z**-6;c=3*(z**12-z**-12)
        if not c or c*c+36==0:continue
        a=z*z*S/2;b=-S/(2*z*z);v=-(z**6+z**-6)/2
        roots3=(t**3-a*t*t+b*t-1).roots(multiplicities=False)
        roots4=(t**4+v*t*t+1).roots(multiplicities=False)
        if len(roots3)!=3 or len(roots4)!=4:continue
        params=[z**j for j in H]
        params[0]*=roots3[0];params[2]*=roots3[0]*roots3[1];params[1]*=roots2[0]
        for k in range(3):params[4+k]*=prod(roots4[:k+1])
        vals=[prod(ti**e for ti,e in zip(params,r)) for r in roots]
        if F(1) in vals:continue
        grad=vector(F,[sum(vv*r[i] for vv,r in zip(vals,roots)) for i in range(7)])
        theta=B.change_ring(F).inverse()*grad
        assert theta==c/F(432)*vector(F,H)
        Hess=matrix(F,7,7,lambda i,j:sum(vv*r[i]*r[j] for vv,r in zip(vals,roots)))
        assert Hess.rank()==7
        outside=sum(sum(r[i]*H[i] for i in range(7))!=0 for r in roots)
        assert outside==106
        return dict(prime=p,central_parameter=int(z),target_c=int(c),
                    torus_coroot_parameters=list(map(int,params)),trace_norm=432,
                    no_source_root_value_equals_one=True,Cartan_rank=7,root_rank=106,
                    full_rank=113,kernel_dimension=20)
    raise AssertionError('No split regular torus witness in the bounded corroboration search')


def main():
    rows=[charge_polynomial(k,tuple(i for i in range(int(k[1:])) if i!=3)) for k in ('E6','E7','E8')]
    highest=[]
    for row in rows:
        n=int(row['type'][1:]);R=RootSystem(row['type']).root_lattice();theta=R.highest_root()
        coeff=list(theta.to_vector());A=CartanType(row['type']).cartan_matrix()
        labels=list(A*vector(coeff));node=labels.index(1)+1
        assert sum(labels)==1
        assert (n,node,int(coeff[3])) in ((6,2,3),(7,1,4),(8,8,6))
        highest.append(dict(type=row['type'],simple_coefficients=list(map(int,coeff)),
                            highest_weight_node=node,omitted_node_coefficient=int(coeff[3])))
    report=dict(result='passed',diagram_checks=[diagram_candidates(k) for k in ('E6','E7','E8')],
                charge_records=rows,highest_root_obstructions=highest,
                E7=e7_identities(rows[1]),independent_original_root_point=root_point(rows[1]),
                raw_shape_classification='E6, omit node 4; original fundamental coordinates only',
                arbitrary_affine_changes_classified=False,all_exceptional_Levis_resolved=False,
                generic_degree_of_full_Springer_map_claim=False,literature_priority_claim=False)
    names=['ADJOINT_PENCIL_SCOPE.md','tools/verify_adjoint_pencil_scope.py',
           'tools/verify_branching_norm_fibers.py','GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'adjoint-pencil-scope-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='charge_records'},indent=2))


if __name__=='__main__':main()
