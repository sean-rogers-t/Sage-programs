"""Bounded search, NOT a uniform theorem or an exhaustive exceptional census.

The source-generic screen intentionally tests only a dense central-torus open.
A negative screen cannot discard boundary components of the actual generic
TARGET fiber. No source-generic result is reported as target smoothness.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import gzip
import json
import time
from sage.all import QQ, PolynomialRing, prod

ROOT = Path(__file__).resolve().parents[1]


def selected_record():
    rows = json.loads(gzip.decompress((ROOT/'unit-constraint-cases.json.gz').read_bytes()))
    return next(row for row in rows if row['character']['type'] == 'E8'
                and row['character']['nodes'] == [1,2,3,5,6,8])


def screen(specialized=False):
    row = selected_record()
    T = PolynomialRing(QQ, names='t,r')
    K = T.fraction_field()
    t,r = map(K,T.gens())
    if specialized:
        K=QQ
        t,r=QQ(2),QQ(3)
    P = PolynomialRing(K,names='x,u,v,y')
    x,u,v,y = P.gens()
    # Original residual variables (x,u_old,v_old,y), central z,w.
    # u_old=u/w, v_old=v*w, t=z^3, r=z^6*w^3.
    # This is a finite etale central change, and unit scaling of affine vars.
    terms = []
    for term in row['residual']['terms']:
        ex = term['exponent']
        a,b = term['charge']
        b += -ex[1]+ex[2]
        assert a % 3 == b % 3 == 0
        terms.append(QQ(term['coefficient'])*t**((a-2*b)//3)*r**(b//3)
                     *prod(vv**ee for vv,ee in zip(P.gens(),ex)))
    phi = sum(terms)
    J = P.ideal([phi.derivative(vv) for vv in P.gens()])
    start = time.monotonic()
    print('source-specialized' if specialized else 'source-generic','ideal started',flush=True)
    gb = J.groebner_basis()
    print('Groebner basis',len(gb),'dimension',J.dimension(),
          'elapsed',time.monotonic()-start,flush=True)
    report = dict(status='source_specialized_screen' if specialized else 'source_generic_screen_only', type='E8',
                  nodes=row['character']['nodes'], potential=str(phi),
                  dimension=int(J.dimension()),
                  vector_space_dimension=int(J.vector_space_dimension()),
                  radical_test_completed=False,
                  target_generic_claim=False,
                  source_boundary_not_tested=True,
                  elapsed_seconds=time.monotonic()-start,
                  source_sha256={p.name:sha256(p.read_bytes()).hexdigest() for p in
                      (Path(__file__),ROOT/'unit-constraint-cases.json.gz')})
    path=ROOT/('adjoint-nilpotent-specialized-screen.json' if specialized else 'adjoint-nilpotent-source-screen.json')
    path.write_text(json.dumps(report,indent=2)+'\n')
    if specialized:
        radical = J.radical()
        outside = [f for f in radical.gens() if J.reduce(f) != 0]
        report.update(radical_test_completed=True,
                      radical_generators_outside_ideal=[str(f) for f in outside],
                      specialized_relative_ideal_radical=not outside)
    else:
        # Spread the monic generic basis to t=2,r=3. Equal specialized reduced
        # bases provide a free rank-38 algebra on a parameter open. A unit
        # Jacobian at this point makes its norm nonzero over QQ(t,r).
        from sage.all import matrix
        Q=PolynomialRing(QQ,names='x,u,v,y')
        def specialize(f):
            values={}
            for exponent,c in f.dict().items():
                den=c.denominator()(2,3)
                assert den!=0
                values[exponent]=c.numerator()(2,3)/den
            return Q(values)
        specialized_phi=specialize(phi)
        I=Q.ideal([specialized_phi.derivative(vv) for vv in Q.gens()])
        bg=I.groebner_basis()
        assert all(f.lc()==1 for f in gb)
        assert [specialize(f) for f in gb]==list(bg)
        assert [f.lm().exponents() for f in gb]==[f.lm().exponents() for f in bg]
        H=matrix(Q,[[specialized_phi.derivative(a).derivative(b)
                     for b in Q.gens()] for a in Q.gens()])
        det=I.reduce(H.det())
        unit=(I+Q.ideal(det)).is_one()
        assert unit
        report.update(jacobian_test_completed=True,specialized_jacobian_ideal_unit=True,
                      monic_generic_basis_specializes_exactly=True,
                      every_basis_coefficient_regular_at_specialization=True,
                      generic_relative_algebra_finite_etale_by_spreading=True,
                      specialization=dict(t=2,r=3),
                      specialized_potential=str(specialized_phi),
                      specialized_hessian_remainder=str(det),
                      specialized_groebner_basis=[str(f) for f in bg])
    report['elapsed_seconds']=time.monotonic()-start
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in
          ('potential','specialized_potential','specialized_hessian_remainder',
           'specialized_groebner_basis')},indent=2),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source-generic-screen',action='store_true')
    parser.add_argument('--specialized-screen',action='store_true')
    args=parser.parse_args()
    if args.source_generic_screen:
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
        screen()
    if args.specialized_screen:
        screen(specialized=True)
