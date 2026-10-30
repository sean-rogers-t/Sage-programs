"""Exact root and boundary checks for the bounded adjoint nilpotent search.

The prime-ideal argument is written in UNIFORM_ADJOINT_PROOF_TARGET.md.
Checks the entire relative ideal, not just its source-generic localization.
"""
from hashlib import sha256
from itertools import product
from pathlib import Path
import gzip
import json
import sympy as s
from audit_unit_constraint_classification import audit_case

ROOT=Path(__file__).resolve().parents[1]


def independent_screen_check():
    """Independent SymPy modular replay of the Sage specialization certificate."""
    from itertools import permutations
    data=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    for path in (ROOT/'tools/search_adjoint_nilpotents.py',ROOT/'unit-constraint-cases.json.gz'):
        assert data['source_sha256'][path.name]==sha256(path.read_bytes()).hexdigest()
    assert data['generic_relative_algebra_finite_etale_by_spreading']
    assert data['vector_space_dimension']==38
    variables=s.symbols('x u v y')
    names=dict(zip(('x','u','v','y'),variables))
    parse=lambda txt:s.sympify(txt.replace('^','**'),locals=names)
    phi=parse(data['specialized_potential'])
    t,r=s.symbols('t r')
    generic=s.sympify(data['potential'].replace('^','**'),locals=dict(names,t=t,r=r))
    assert s.expand(generic.subs({t:2,r:3})-phi)==0
    rows=json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    row=next(v for v in rows if v['character']['type']=='E8'
             and v['character']['nodes']==[1,2,3,5,6,8])
    rebuilt=0
    for term in row['residual']['terms']:
        ex=term['exponent']; az,bw=term['charge']; bw+=-ex[1]+ex[2]
        assert az%3==bw%3==0
        rebuilt+=s.Rational(term['coefficient'])*t**((az-2*bw)//3)*r**(bw//3)*s.prod(
            vv**ee for vv,ee in zip(variables,ex))
    assert s.expand(generic-rebuilt)==0
    prime=65521
    def mod(f):
        return sum((int(c.p)*pow(int(c.q),-1,prime)%prime)*s.prod(v**e for v,e in
                   zip(variables,ex)) for ex,c in s.Poly(f,*variables,domain=s.QQ).terms())
    partials=[mod(s.diff(phi,v)) for v in variables]
    basis=s.groebner(partials,*variables,modulus=prime,order='grevlex')
    stored=[s.Poly(mod(parse(f)),*variables,modulus=prime) for f in
            data['specialized_groebner_basis']]
    assert set(basis.polys)==set(stored)
    # Independently recompute the Hessian determinant, without trusting the
    # stored Sage remainder, then check its ideal is the unit ideal.
    H=[[mod(s.diff(phi,a,b)) for b in variables] for a in variables]
    determinant=0
    for perm in permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        determinant+=sign*s.prod(H[i][perm[i]] for i in range(4))
    remainder=basis.reduce(s.expand(determinant))[1]
    assert s.Poly(remainder-mod(parse(data['specialized_hessian_remainder'])),
                  *variables,modulus=prime).is_zero
    test=s.groebner(partials+[remainder],*variables,modulus=prime,order='grevlex')
    assert len(test.polys)==1 and test.polys[0].as_expr()==1
    leading=[p.LM(order='grevlex').exponents for p in basis.polys]
    bounds=[min(e[i] for e in leading if e[i] and sum(e)==e[i]) for i in range(4)]
    standard=[e for e in product(*(range(n) for n in bounds))
              if not any(all(a<=b for a,b in zip(lm,e)) for lm in leading)]
    assert len(standard)==38
    return dict(independent_engine='SymPy',prime=prime,
                specialized_basis_recomputed=True,standard_monomials=38,
                hessian_remainder_recomputed=True,jacobian_ideal_unit=True,
                central_boundary_tested=False)


def main():
    rows=json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    selected=[]
    for nodes in ([1,2,3,4,6,8],[1,2,3,5,6,8]):
        row=next(v for v in rows if v['character']['type']=='E8'
                 and v['character']['nodes']==nodes)
        selected.append(row)
        audit_case(row)  # independent full root reconstruction and all derivatives

    a,b,x,y,z,w,r=s.symbols('a b x y z w r')
    H=(-w**2-w**-2-3-x*x*((w*r+1/(w*r))*y+y*y)
       -x*(w*w*r+1/(w*w*r)+(r+1/r)*y*y+2*(w+1/w)*y)
       -(r/w+w/r)*y)
    phi=a*(x+w*y+1/r)+b*(x+y/w+r)+H
    # Residual variables in the archive are a_old,b_old,x,y.
    # a=a_old*z, b=b_old/z, r=z^5*w^2.
    original=sum(s.Rational(t['coefficient'])*s.prod(v**e for v,e in
                 zip((a,b,x,y),t['exponent']))*z**t['charge'][0]
                 *w**t['charge'][1] for t in selected[0]['residual']['terms'])
    pullback=lambda f:s.cancel(f.subs({a:a*z,b:b/z,r:z**5*w*w},simultaneous=True))
    assert s.cancel(original-pullback(phi))==0
    # Chain rule before setting ANY affine derivative to zero.
    assert s.cancel(z*s.diff(original,z)-pullback(
        5*r*s.diff(phi,r)+a*s.diff(phi,a)-b*s.diff(phi,b)))==0
    assert s.cancel(w*s.diff(original,w)-pullback(
        w*s.diff(phi,w)+2*r*s.diff(phi,r)))==0
    xx=-w*y-1/r
    aa=-b-s.diff(H,x).subs(x,xx)
    graph={a:aa,x:xx}
    F=(1/w-w)*y+r-1/r
    h=s.expand(H.subs(x,xx))
    reduced=b*F+h
    assert s.cancel(phi.subs(graph,simultaneous=True)-reduced)==0
    for v in (a,x):
        assert s.cancel(s.diff(phi,v).subs(graph,simultaneous=True))==0
    for v in (b,y,r,w):
        assert s.cancel(s.diff(phi,v).subs(graph,simultaneous=True)
                        -s.diff(reduced,v))==0
    A=1/w-w
    q=s.diff(reduced,y)
    assert s.cancel(q-b*A-s.diff(h,y))==0
    assert s.Poly(h,y).LC()==-w*w
    # F is primitive linear in y over the central Laurent UFD.
    assert s.gcd(1-w*w,r*r-1)==1
    # The open A!=0 is exactly the central open torus graph.
    yy=-(r-1/r)/A
    bb=-s.diff(h,y).subs(y,yy)/A
    assert s.cancel(F.subs(y,yy))==0
    assert s.cancel(q.subs({y:yy,b:bb},simultaneous=True))==0
    boundary=[]
    for epsilon,eta in product((-1,1),repeat=2):
        spec={w:epsilon,r:eta}
        p=s.factor(s.diff(h,y).subs(spec))
        expected=-2*(2*y+epsilon*eta)*(y*y+epsilon*eta*y-1)
        assert s.expand(p-expected)==0
        assert s.degree(p,y)==3
        assert s.degree(s.gcd(p,s.diff(p,y)),y)==0
        u0,v0,gcd=s.gcdex(p,s.diff(p,y),y)
        assert s.expand(u0*p+v0*s.diff(p,y))==1
        jac_minor=s.diff(F,r)*s.diff(q,y)-s.diff(F,y)*s.diff(q,r)
        assert s.cancel(jac_minor.subs(spec)-2*s.diff(p,y))==0
        boundary.append(dict(w=epsilon,r=eta,polynomial=str(p),squarefree=True,
                             bezout=[str(u0),str(v0)],
                             jacobian_minor='2*derivative(boundary_polynomial,y)'))
    # Independent central Jacobian for the finite etale coordinate change.
    assert s.Matrix([[5,2],[0,1]]).det()==5
    # Source-generic screen change in the OTHER selected model.
    assert s.Matrix([[3,0],[6,3]]).det()==9
    for term in selected[1]['residual']['terms']:
        az,bw=term['charge']; ex=term['exponent']
        assert az%3==0 and (bw-ex[1]+ex[2])%3==0
    report=dict(result='passed',search_scope='two selected E8 adjoint models',
        full_characters_independently_reconstructed=True,
        every_original_derivative_replayed=True,
        exact_unit_elimination=True,
        first_model=dict(nodes=selected[0]['character']['nodes'],
            reduced_potential='b*((w^-1-w)*y+r-r^-1)+h(y,w,r)',
            leading_y_coefficient_of_h='-w^2',
            relative_ideal=['(w^-1-w)*y+r-r^-1','b*(w^-1-w)+h_y'],
            boundary=boundary,open_dimension=2,boundary_dimension=1,
            complete_intersection_prime_by_written_argument=True,
            entire_relative_critical_surface_smooth=True,
            all_central_boundaries_retained=True,
            generic_target_reduced_or_empty=True),
        counterexample_found=False,
        uniform_conjecture_proved=False,
        second_model_source_screen_alone_resolves_boundary=False,
        independent_source_screen=independent_screen_check())
    paths=[Path(__file__), ROOT/'unit-constraint-cases.json.gz',
           ROOT/'UNIFORM_ADJOINT_PROOF_TARGET.md',
           ROOT/'tools/search_adjoint_nilpotents.py',
           ROOT/'adjoint-nilpotent-source-screen.json']
    report['sha256']={str(p.relative_to(ROOT)).replace('\\','/'):
                       sha256(p.read_bytes()).hexdigest() for p in paths}
    (ROOT/'adjoint-nilpotent-search-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
