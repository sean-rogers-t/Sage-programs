"""Exact multiplier/E8 certificate with original-root differential corroboration."""
from hashlib import sha256
from pathlib import Path
import json
from sage.all import (QQ, GF, PolynomialRing, LaurentPolynomialRing, RootSystem,
                      matrix, vector, prod)
from verify_branching_norm_fibers import charge_polynomial
from linear_character_elimination import eliminate

ROOT=Path(__file__).resolve().parents[1]


def dz(phi,z):
    return phi.parent()({e:z*v.derivative(z) for e,v in phi.dict().items()})


def helper_checks():
    R=LaurentPolynomialRing(QQ,'z');z=R.gen()
    P=PolynomialRing(R.fraction_field(),['a','b','u','v','y']);a,b,u,v,y=P.gens()
    f=u*(z*a+b-y*y-z)+v*(b/z-y/z)+a*a*b+y*a+z*b
    cert=eliminate(f,(u,v),(a,b),R);graph=cert['graph'];red=cert['reduced']
    assert f.derivative(y).subs(graph)==red.derivative(y)
    assert dz(f,z).subs(graph)==dz(red,z)
    rejected=[]
    for bad in (u*(z-1)*a+a*a,u*u*a+a):
        try:eliminate(bad,(u,),(a,),R)
        except ValueError as error:rejected.append(str(error))
        else:raise AssertionError('Unsupported constraint accepted')
    return dict(coupled_unit_matrix=True,remaining_affine_chain_rule=True,
                central_matrix_derivative_retained=True,negative_cases=rejected)


def symbolic(row):
    R=LaurentPolynomialRing(QQ,'z');z=R.gen()
    P=PolynomialRing(R.fraction_field(),['a','b','x','p','q','r','s']);a,b,x,p,q,r,s=P.gens()
    local=(a,x,b,p,q,r,s)
    phi=sum(t['coefficient']*z**t['charge'][0]*prod(v**e for v,e in zip(local,t['exponent'])) for t in row['terms'])
    formula=a*b+x*x+p*s-2+z*a*x*s+z**2*b*r+z**3*x*q+z**4*a*p+z**5*x*b+z**6*s+x*b*p/z+a*q/z**2+x*r/z**3+b*s/z**4+a*x/z**5+p/z**6
    assert phi==formula
    cert=eliminate(phi,(q,r),(a,b),R);graph=cert['graph']
    assert cert['determinant']==1
    U=z**6*(x*x-1)+z**-9*x;V=z**-6*(x*x-1)+z**9*x
    assert cert['reduced']==p*s-p*V-s*U-2
    pair={p:U,s:V};graph={v:f.subs(pair) for v,f in graph.items()};graph.update(pair)
    psi=-(x*x-1)*(x*x+(z**15+z**-15)*x)-3
    assert phi.subs(graph)==psi
    assert all(phi.derivative(v).subs(graph)==0 for v in (a,b,p,q,r,s))
    assert phi.derivative(x).subs(graph)==psi.derivative(x)
    assert dz(phi,z).subs(graph)==dz(psi,z)
    assert dz(psi,z)==15*x*(1-x*x)*(z**15-z**-15)
    weights=(-10,-20,-15,-24,-18,-12,-6)
    for exponents,coefficient in phi.dict().items():
        for power in R(coefficient).dict():
            assert (sum(i*j for i,j in zip(exponents,weights))+int(power))%30==0
    # Each recovered variable has the required central-kernel weight.
    for var,f in graph.items():
        target=weights[P.gens().index(var)]
        for exponents,coefficient in f.dict().items():
            for power in R(coefficient).dict():
                assert (sum(i*j for i,j in zip(exponents,weights))+int(power)-target)%30==0
    return dict(character_matches=True,unit_constraint_determinant=True,
                all_original_affine_and_central_derivatives=True,
                residual='-(x^2-1)*(x^2+(z^15+z^-15)*x)-3',
                central_action_and_graph_verified=True,no_source_boundary_removed=True)


def quotient_checks():
    R=PolynomialRing(QQ,['Y','k']);Y,k=R.gens()
    P=900*Y*(1-Y)**2*(4*Y-1)*(Y**2-3*Y+1)-k*(1-3*Y)**2
    Delta=4096*k**5+4485075*k**4+1691730000*k**3+292161500000*k*k+23287500000000*k+683437500000000
    K=PolynomialRing(QQ,'k');kk=K.gen();T=PolynomialRing(K,'Y');yy=T.gen()
    f=T([K(P.coefficient({Y:i})) for i in range(7)])
    dd=K(Delta)
    assert f.discriminant()==-30611001600000000*kk*dd
    assert dd.gcd(dd.derivative())==1 and dd(0)!=0
    assert dd.is_irreducible()
    L=K.quotient(dd,'j');FL=PolynomialRing(L,'Y')
    special=FL(f);g=special.gcd(special.derivative())
    assert g.degree()==1 and special.gcd(special.derivative(2))==1
    assert [P.subs({Y:v}) for v in (0,1,QQ(1)/3)]==[-k,-4*k,QQ(400)/81]
    # Rational inverse without taking radicals: Y, (1-Y), (1-3Y) are units.
    C0=PolynomialRing(QQ,'c');c=C0.gen();K0=C0.fraction_field()
    T0=PolynomialRing(K0,'Y');y=T0.gen()
    ff=900*y*(1-y)**2*(4*y-1)*(y*y-3*y+1)-c*c*(1-3*y)**2
    Q=T0.quotient(ff,'y');y=Q.gen()
    # Invariants M=x*z^15 and N=x*z^-15. MN=Y, M-N=c/(15B).
    M=y*(2*y-1)/(1-3*y)+c/(30*(1-y))
    N=y*(2*y-1)/(1-3*y)-c/(30*(1-y))
    assert M*N==y
    assert (1-3*y)*(M+N)==2*y*(2*y-1)
    assert 15*(1-y)*(M-N)==c
    return dict(degree=6,discriminant=str(dd),nonzero_exceptional_targets=10,
                every_exceptional_fiber_one_double_four_simple=True,
                generic_regular_orbits=6,regular_orbit_dimension=28,
                quotient_inverse_checked=True,excluded_coordinate_values_not_roots=True)


def root_point(row):
    H=row['central_cocharacters'][0]
    ambient=RootSystem('E8').ambient_space();cr=ambient.simple_coroots()
    roots=[tuple(int(r.scalar(cr[i+1])) for i in range(8)) for r in ambient.roots()]
    B=matrix(QQ,8,8,lambda i,j:sum(r[i]*r[j] for r in roots))
    assert vector(QQ,H)*B*vector(QQ,H)==1800
    for prime in (1033,3001,32003):
        F=GF(prime);R=PolynomialRing(F,'t');t=R.gen()
        for zi in range(2,min(prime,4000)):
            z=F(zi);S=z**15+z**-15
            for x in (4*t**3+3*S*t*t-2*t-S).roots(multiplicities=False):
                c=15*x*(1-x*x)*(z**15-z**-15)
                if not c:continue
                a=-z**5*x;b=-z**-5*x
                p=z**6*(x*x-1)+z**-9*x;s=z**-6*(x*x-1)+z**9*x
                q=-z*z*(b+z*x*s+z**4*p+z**-5*x)
                r=-z**-2*(a+z**-1*x*p+z**-4*s+z**5*x)
                r2=(t*t-x*t+1).roots(multiplicities=False)
                if len(r2)!=2:continue
                r3=(t**3-a*t*t+b*t-1).roots(multiplicities=False)
                if len(r3)!=3:continue
                r5=(t**5-p*t**4+q*t**3-r*t*t+s*t-1).roots(multiplicities=False)
                if len(r5)!=5:continue
                params=[z**j for j in H]
                params[0]*=r3[0];params[2]*=r3[0]*r3[1];params[1]*=r2[0]
                for i in range(4):params[4+i]*=prod(r5[:i+1])
                vals=[prod(ti**e for ti,e in zip(params,rr)) for rr in roots]
                if F(1) in vals:continue
                grad=vector(F,[sum(v*rr[i] for v,rr in zip(vals,roots)) for i in range(8)])
                assert B.change_ring(F).inverse()*grad==c/F(1800)*vector(F,H)
                Hess=matrix(F,8,8,lambda i,j:sum(v*rr[i]*rr[j] for v,rr in zip(vals,roots)))
                if Hess.rank()!=8:continue
                outside=sum(sum(rr[i]*H[i] for i in range(8))!=0 for rr in roots)
                assert outside==212
                return dict(prime=prime,central_parameter=int(z),x=int(x),target_c=int(c),
                            torus_coroot_parameters=list(map(int,params)),trace_norm=1800,
                            no_source_root_value_equals_one=True,Cartan_rank=8,
                            root_rank=212,full_rank=220,kernel_dimension=28)
    raise AssertionError('No split witness found; this finite search is corroboration only')


def main():
    row=charge_polynomial('E8',(0,1,2,4,5,6,7))
    report=dict(result='passed',helper=helper_checks(),charge_record=row,
                symbolic=symbolic(row),quotient=quotient_checks(),root_point=root_point(row),
                all_exceptional_Levis_resolved=False,full_map_generic_degree_claim=False,
                literature_priority_claim=False)
    names=['MULTIPLIER_BRANCHING_FIBERS.md','tools/linear_character_elimination.py',
           'tools/verify_multiplier_branching_fibers.py','tools/verify_branching_norm_fibers.py',
           'GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report['source_sha256']={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'multiplier-branching-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='charge_record'},indent=2))


if __name__=='__main__':main()
