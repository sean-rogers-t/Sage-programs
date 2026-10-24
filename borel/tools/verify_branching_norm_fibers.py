"""Charge-preserving checks for the branching norm and fiber criteria."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from pathlib import Path
import json

from sage.all import (QQ, ZZ, CartanType, RootSystem, matrix, vector, lcm, gcd, prod,
                      PolynomialRing, CyclotomicField)
from verify_affine_branching_classification import levi_data, character_weights
from verify_quadratic_branching_fibers import convolve, fundamental_local

ROOT=Path(__file__).resolve().parents[1]


def charge_polynomial(kind,I,highest=None):
    """Exact character inversion independently in each central charge."""
    n=int(kind[1:]); I=tuple(I); J=tuple(i for i in range(n) if i not in I)
    A=CartanType(kind).cartan_matrix(); Ai=A.inverse()
    central=[]
    for j in J:
        v=Ai.column(j); v=v*lcm(x.denominator() for x in v)
        v=v/gcd(list(v)); central.append(v)
    cover=matrix(QQ,n,n,lambda i,j: int(i==I[j]) if j<len(I) else central[j-len(I)][i])
    assert cover.det()!=0
    if highest is None:
        R=RootSystem(kind).ambient_space(); cr=R.simple_coroots()
        weights=Counter({tuple(int(a.scalar(cr[j+1])) for j in range(n)):1 for a in R.roots()})
        weights[(0,)*n]=n
    else:weights=Counter(character_weights(kind,highest))
    groups=defaultdict(Counter)
    for w,m in weights.items():
        q=tuple(int(vector(QQ,w).dot_product(v)) for v in central)
        groups[q][tuple(w[i] for i in I)]+=m
    components,height=levi_data(kind,I); terms=[]
    for q,original in sorted(groups.items()):
        residual=original.copy(); poly={}
        while residual:
            dominant=[w for w in residual if all(v>=0 for v in w)]
            assert dominant
            nu=max(dominant,key=lambda w:(sum(h*v for h,v in zip(height,w)),w))
            coeff=residual[nu]; poly[nu]=poly.get(nu,0)+coeff
            small={(0,)*len(I):1}
            for i,e in enumerate(nu):
                for _ in range(e):small=convolve(small,fundamental_local(kind,I,i))
            for w,m in small.items():
                residual[w]-=coeff*m
                if not residual[w]:del residual[w]
        rebuilt=Counter()
        for exponent,coeff in poly.items():
            if not coeff:continue
            small={(0,)*len(I):1}
            for i,e in enumerate(exponent):
                for _ in range(e):small=convolve(small,fundamental_local(kind,I,i))
            for w,m in small.items():rebuilt[w]+=coeff*m
            terms.append(dict(charge=list(q),exponent=list(exponent),coefficient=int(coeff)))
        assert dict(rebuilt)==dict(original)
    return dict(type=kind,nodes=[i+1 for i in I],factors=[f for f,o in components],
                factor_orders=[list(o) for f,o in components],
                central_cocharacters=[list(map(int,v)) for v in central],
                cover_degree=abs(int(cover.det())),
                highest=None if highest is None else list(highest),
                dimension=sum(weights.values()),terms=terms,
                max_degree=max(sum(t['exponent']) for t in terms),
                reconstructed_every_central_charge=True)


def main():
    cases=[('E6',(1,2,4)),('E6',(1,2,4,5)),('E6',(0,1,2,4,5))]
    rows=[]
    for kind,I in cases:
        row=charge_polynomial(kind,I);rows.append(row)
        print(row['type'],row['nodes'],row['factors'],'degree',row['max_degree'],flush=True)
    for row in rows:
        I=[i-1 for i in row['nodes']]; A=CartanType(row['type']).cartan_matrix()
        Q=A.matrix_from_rows_and_columns(I,I).inverse()
        H=matrix(QQ,row['central_cocharacters']).transpose(); C=(H.transpose()*A*H).inverse()
        for term in row['terms']:
            e=vector(QQ,term['exponent']); q=vector(QQ,term['charge'])
            assert e*Q*e+q*C*q<=2
            if any(q):
                assert e*Q*e<2
                assert all(sum(e[j] for j in order)<=1 for order in row['factor_orders'])
                assert sum(e)<=3
        row['all_charge_norm_budgets_verified']=True
    # A non-affine, genuinely quadratic restriction of the E6 minimal module.
    minimal=charge_polynomial('E6',(1,2,4),(1,0,0,0,0,0))
    assert minimal['max_degree']==2 and minimal['dimension']==27
    report=dict(result='passed',adjoint_mechanism_tests=rows,
                E6_minimal_quadratic_control=minimal,
                exact_suspension_checks=suspension_checks(rows[-1]),
                E6_direct_root_differential=e6_root_point(rows[-1]),
                all_exceptional_adjoint_Levis_resolved=False,
                literature_priority_claim=False,public_API_changed=False)
    names=['BRANCHING_NORM_FIBERS.md','tools/verify_branching_norm_fibers.py',
           'tools/verify_affine_branching_classification.py','tools/verify_quadratic_branching_fibers.py',
           'GENERAL_CAYLEY_FIBER_CRITICAL.md','QUADRATIC_BRANCHING_FIBERS.md']
    report['source_sha256']={name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in names}
    (ROOT/'branching-norm-fibers-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('adjoint_mechanism_tests','E6_minimal_quadratic_control')},indent=2))


def suspension_checks(record):
    F=PolynomialRing(QQ,'c').fraction_field(); c=F.gen()
    R=PolynomialRing(F,['a','b','y','u','v','z','w'])
    a,b,y,u,v,z,w=R.gens(); S=z**3+w**3
    chi=y*y+S*y-2+a*b+u*v+(w*w+z*y)*a*v+(z*z+w*y)*b*u
    coordinates=(a,y,b,u,v)
    reconstructed=R(0)
    for term in record['terms']:
        charge=term['charge'][0]
        reconstructed += term['coefficient']*prod(x**e for x,e in zip(coordinates,term['exponent']))*(z**charge if charge>=0 else w**(-charge))
    assert reconstructed==chi
    equations=[chi.derivative(x) for x in (a,b,y,u,v)]
    equations += [z*chi.derivative(z)-w*chi.derivative(w)-c,z*w-1]
    I=R.ideal(equations)
    J=R.ideal([a,b,u,v,2*y+S,3*(z**6-w**6)+2*c,z*w-1])
    assert I==J and I.dimension()==0
    jac=matrix(R,[[f.derivative(x) for x in R.gens()] for f in equations])
    assert I+R.ideal(jac.det())==R.ideal(1)
    # Full coordinate change, including the base-dependent terms.
    identity=chi-((a+(z*z+w*y)*u)*(b+(w*w+z*y)*v)+(y*y+S*y)*(1-u*v)-2)
    assert identity in R.ideal(z*w-1)
    # The general suspension with TWO hyperbolic pairs retains a free variable.
    T=PolynomialRing(F,['y','u1','v1','u2','v2','free','z','w'])
    y,u1,v1,u2,v2,free,z,w=T.gens(); S=z+w
    phi=(y*y+S*y)*(1-u1*v1-u2*v2)
    eq=[phi.derivative(x) for x in (y,u1,v1,u2,v2,free)]
    eq += [z*phi.derivative(z)-w*phi.derivative(w)-c,z*w-1]
    expected=T.ideal([2*y+S,u1,v1,u2,v2,(z*z-w*w)+2*c,z*w-1])
    assert T.ideal(eq)==expected and expected.dimension()==1
    # Uniform localization at c!=0, including the nonreduced exceptional targets.
    U=PolynomialRing(QQ,['y','u','v','z','w','c','d'])
    y,u,v,z,w,c,d=U.gens(); S=z**3+w**3
    phi=(y*y+S*y)*(1-u*v)
    equations=[phi.derivative(x) for x in (y,u,v)]
    equations += [z*phi.derivative(z)-w*phi.derivative(w)-c,z*w-1,c*d-1]
    expected=U.ideal([2*y+S,u,v,3*(z**6-w**6)+2*c,z*w-1,c*d-1])
    assert U.ideal(equations)==expected
    return dict(E6_original_gradient_ideal_equals_scalar_model=True,
                E6_formula_matches_root_reconstructed_character=True,
                suspension_uniform_on_nonzero_target=True,
                E6_generic_quotient_finite_smooth=True,
                E6_target_equation='3*(z^6-z^(-6))+2*c=0',
                two_pair_suspension_free_coordinate_preserved=True)


def e6_root_point(record):
    # Direct adjoint-weight trace at z=2, SL3 traces zero and SL2 trace -65/16.
    F=CyclotomicField(3); omega=F.gen(); P=PolynomialRing(F,'a'); a=P.gen()
    L=F.extension(a*a+QQ(65)/16*a+1,'lam'); lam=L.gen(); om=L(omega)
    H=vector(QQ,record['central_cocharacters'][0]); params=[L(2)**int(h) for h in H]
    # A2 chains (1,3) and (5,6), A1 node 2.
    for order in ((0,2),(4,5)):
        params[order[0]]*=1; params[order[1]]*=om
    params[1]*=lam
    R=RootSystem('E6').ambient_space(); cr=R.simple_coroots()
    roots=[tuple(int(alpha.scalar(cr[j+1])) for j in range(6)) for alpha in R.roots()]
    values=[prod(params[j]**w[j] for j in range(6)) for w in roots]
    assert all(v!=1 for v in values)
    B=matrix(QQ,6,6,lambda i,j:sum(w[i]*w[j] for w in roots))
    grad=vector(L,[sum(v*w[i] for v,w in zip(values,roots)) for i in range(6)])
    c=-QQ(12285)/128; norm=H*B*H
    assert norm==144 and B.change_ring(L).inverse()*grad==vector(L,H)*c/norm
    Hess=matrix(L,6,6,lambda i,j:sum(v*w[i]*w[j] for v,w in zip(values,roots)))
    torus_rank=Hess.rank(); outside=sum(vector(QQ,w)*H!=0 for w in roots)
    assert torus_rank==6 and outside==58
    return dict(central_trace_norm=int(norm),target_parameter=str(c),
                normalized_target_multiple=str(c/norm),regular_semisimple_source=True,
                cartan_rank=int(torus_rank),root_direction_rank=int(outside),
                full_rank=int(torus_rank+outside),kernel_dimension=78-int(torus_rank+outside))


if __name__=='__main__':main()
