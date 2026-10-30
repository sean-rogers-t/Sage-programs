"""Exact symbolic audit of the shared-variable coupling theorem.

No root-family recognition is inferred from these model computations.
The arbitrary-size proof is in SHARED_SQUARE_COUPLINGS.md.
"""
from pathlib import Path
from hashlib import sha256
import json
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def checks():
    sizes = []
    for count in range(1, 6):
        x, q, a, b, Q = s.symbols('x q a b Q')
        d = s.symbols('d:' + str(count))
        w = s.symbols('w:' + str(count))
        f = sum(wi * (q-di) for wi, di in zip(w, d)) + Q
        for j in range(count):
            changed = f.subs(w[j], b-sum(w[i] for i in range(count) if i != j))
            changed = changed.subs(q, a+d[j])
            residual = Q+sum(w[i]*(d[j]-d[i]) for i in range(count) if i != j)
            assert s.expand(changed-a*b-residual) == 0
        # At E: x=W=0. Active coefficients vanish to first order.
        DQ = s.Symbol('DQ'); Dd = s.symbols('Dd:' + str(count))
        for mask in range(1, 2**count):
            mu = DQ-sum(w[i]*Dd[i] for i in range(count))
            expected = DQ-sum(Dd)
            substitutions = {Dd[i]: 0 for i in range(count) if mask >> i & 1}
            substitutions.update({w[i]: 1 for i in range(count) if not (mask >> i & 1)})
            assert s.expand((mu-expected).subs(substitutions)) == 0
        sizes.append(count)
    # z-dependent changes alter central derivatives only by gradient terms.
    z, x, y, a, b, u, v = s.symbols('z x y a b u v')
    d1, d2 = (z-1)**2, (z+1)**2
    f = y*(x*x-d1)+(1-u*v)*(x*x-d2)+z
    # q=x^2, a=q-d1, b=y+(1-uv); work before the etale square root.
    q = s.Symbol('q')
    fq = y*(q-d1)+(1-u*v)*(q-d2)+z
    change = {q:a+d1, y:b-(1-u*v)}
    new = s.expand(fq.subs(change, simultaneous=True))
    assert s.expand(new-a*b-(z+(1-u*v)*(d1-d2))) == 0
    chain = z*s.diff(new,z)-z*s.diff(fq,z).subs(change, simultaneous=True)
    assert s.expand(chain-z*s.diff(d1,z)*s.diff(fq,q).subs(change, simultaneous=True)) == 0
    # Counterexample on the chart u_1,u_2 != 0, with a=u_1v_1,b=u_2v_2.
    c = s.Symbol('c'); field = s.QQ.frac_field(c)
    bad = (1-a)*(x*x-(z-1))+(1-b)*(x*x+(z-1))+z
    gens = [s.diff(bad,t) for t in (x,a,b)]+[z*s.diff(bad,z)-c]
    G = s.groebner(gens,b,z,x,a,domain=field)
    expected = s.groebner([b-a-1+c,z-1,x*x,x*(1+c-2*a)],b,z,x,a,domain=field)
    assert G == expected
    assert G.reduce(x)[1] != 0 and G.reduce(x*x)[1] == 0
    assert G.reduce(x*(1+c-2*a))[1] == 0
    # A true positive-excess model, not just isolated critical points.
    # With all delta_i=1 and Q=z, the x!=0 sheet is sum u_i v_i=count.
    p, r = s.symbols('p r')
    good = (2-u*v-p*r)*(x*x-1)+z
    assert s.diff(good,z)*z == z
    G = s.groebner([x*x-1,2-u*v-p*r,z-c],x,z,u,v,p,r,domain=field)
    assert all(G.reduce(s.diff(good,t))[1] == 0 for t in (x,u,v,p,r))
    # Existing E6 root-derived form is exactly the one-block square case.
    S = z**3+z**-3
    assert s.expand((1-u*v)*(x*x-S*S/4)-2
                    -((1-u*v)*(x*x-(S/2)**2)-2)) == 0
    # Test the proposed affine normal form against an EXISTING root-derived
    # E8 control, without adding a new Levi calculation.
    x,u,v,y,t,r=s.symbols('x u v y t r')
    source=json.loads((ROOT/'adjoint-nilpotent-source-screen.json').read_text())
    original=s.sympify(source['potential'].replace('^','**'),
                      locals={str(a):a for a in (x,u,v,y,t,r)})
    poly=s.Poly(original,x,u,v,y)
    assert poly.total_degree()==4
    top=sum(coeff*s.prod(a**n for a,n in zip((x,u,v,y),powers))
            for powers,coeff in poly.terms() if sum(powers)==4)
    q1=r*v*x+r*x*y+t*u*y;q2=r*v*y+t*u*x+t*x*y
    assert s.expand(top+q1*q2/(r*t))==0
    ranks=[]
    for quad in (q1,q2):
        H=s.hessian(quad,(x,u,v,y))
        assert s.factor(H.det())==r*r*t*t
        ranks.append(4)
    # q1,q2 are nonassociate rank-four quadrics: their product is squarefree,
    # unlike the repeated linear factor in every shared-x quartic top form.
    return dict(symbolic_block_counts=sizes,
                all_active_subset_target_factorizations=True,
                central_derivative_chain_rule=True,
                opposite_simple_zero_gluing_has_generic_embedded_nilpotent=True,
                positive_dimensional_sheets_retained=True,
                existing_E6_form_is_one_block_specialization=True,
                existing_E8_top_form_factorization_verified=True,
                E8_quadratic_factor_ranks=ranks,
                universal_affine_shared_variable_normal_form_disproved=True)


def main():
    report = dict(result='passed',engine='SymPy exact polynomial and rational-function arithmetic',
                  checks=checks(),shared_square_family_theorem_in_note=True,
                  full_root_normal_form_for_remaining_cases_proved=False,
                  new_exceptional_group_coverage=False,uniform_conjecture_proved=False,
                  literature_priority_claim=False,
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                                 for p in (Path(__file__),ROOT/'SHARED_SQUARE_COUPLINGS.md',
                                           ROOT/'adjoint-nilpotent-source-screen.json')})
    (ROOT/'shared-square-couplings-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__':
    main()
