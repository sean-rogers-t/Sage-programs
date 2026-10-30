"""Exact controls for the fiber review; no Sage or cached polynomial basis.

These controls test scheme/target distinctions, not the general theorem or
the open simple-adjoint conjecture. The E7 membership flags must stay false.
"""
from hashlib import sha256
import gzip
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def same_ideal(left, right, variables, field):
    a = sp.groebner(left, *variables, domain=field, order='lex')
    b = sp.groebner(right, *variables, domain=field, order='lex')
    assert a == b, (a, b)
    return a


def main():
    x, u, v, z, c = sp.symbols('x u v z c')
    field = sp.QQ.frac_field(c)
    D = lambda f: sp.expand(z * sp.diff(f, z))
    controls = {}

    # Nonreduced incidence of large dimension, but a single bad target.
    phi = (z-1)**2*x+z
    same_ideal([sp.diff(phi, x), D(phi)-c], [1], (x, z), field)
    assert D(phi).subs(z, 1) == 1
    controls['nonreduced_source_nondominant_bad_image'] = True

    # The rank-drop locus gives the ENTIRE nonempty generic fiber.
    phi = (z-1)*x*x/2
    equations = [sp.diff(phi, x), D(phi)-c]
    same_ideal(equations, [z-1, x*x-2*c], (z, x), field)
    assert sp.discriminant(x*x-2*c, x) == 8*c
    assert sp.Matrix([[sp.diff(f, t) for t in (x, z)]
                      for f in equations]).det().subs(z, 1) == -x*x
    controls['rank_drop_branch_retained_and_generically_reduced'] = True

    # Generic smoothness cannot be invoked on an arbitrary source scheme.
    phi = x**3/3+z
    basis = same_ideal([sp.diff(phi, x), D(phi)-c],
                       [x*x, z-c], (x, z), field)
    assert basis.reduce(x)[1] == x and basis.reduce(x*x)[1] == 0
    controls['persistent_double_structure_detected'] = True

    # Identical support with a good/bad collision coefficient.
    good = x*x+x+z+(x+1)*u*v
    bad = x*x+2*x+z+(x+1)*u*v
    assert set(sp.Poly(good, x,u,v,z).monoms()) == set(sp.Poly(bad, x,u,v,z).monoms())
    basis = same_ideal([sp.diff(bad, t) for t in (x,u,v)]+[D(bad)-c],
                       [2*x+2+u*v, u*u*v, u*v*v, z-c], (x,u,v,z), field)
    assert basis.reduce(u*v)[1] == u*v
    assert basis.reduce((u*v)**2)[1] == 0
    controls['same_support_persistent_collision_nilpotent'] = True

    # A smooth excess direction must not be declared singular by oversized minors.
    phi = x*x/2+z  # u is an unused affine coordinate
    relative = [sp.diff(phi,t) for t in (x,u)]
    jac = sp.Matrix([[sp.diff(f,t) for t in (x,u,z)] for f in relative])
    assert jac.rank() == 1
    assert all(jac[:, list(cols)].det() == 0 for cols in ((0,1),(0,2),(1,2)))
    same_ideal(relative+[D(phi)-c], [x,z-c], (x,u,z), field)
    controls['smooth_excess_not_a_maximal_minor_singularity'] = True

    # Smooth source does not mean every target fiber is smooth.
    phi = x**3/3+(z-1)*x
    assert sp.diff(sp.diff(phi,x),z) == 1
    same_ideal([sp.diff(phi,x),D(phi)-c],
               [z+x*x-1,x**3-x+c], (z,x), field)
    assert sp.discriminant(x**3-x+c,x) == 4-27*c*c
    controls['smooth_source_has_special_branch_targets'] = True

    # A base-dependent variable change must retain the target chain rule.
    a, b = sp.symbols('a b')
    phi = a*(z*b-x)+b*b+z
    graph = {b:x/z,a:-2*x/z**2}
    reduced = x*x/z**2+z
    assert sp.expand(phi.subs(graph, simultaneous=True)-reduced) == 0
    assert sp.expand(D(phi).subs(graph, simultaneous=True)-D(reduced)) == 0
    controls['unit_elimination_preserves_original_central_derivative'] = True

    # Independent defining-SL2 trace normalization, in a matrix chart.
    a,b,d,t = sp.symbols('a b d t')
    g = sp.Matrix([[a,b],[d,(1+b*d)/a]])
    theta = g-sp.trace(g)*sp.eye(2)/2
    basis = [sp.Matrix([[0,1],[0,0]]),sp.Matrix([[0,0],[1,0]]),sp.diag(1,-1)]
    for H in basis:
        velocity = g*H
        derivative = sum(sp.diff(sp.trace(g),q)*dq for q,dq in
                         zip((a,b,d),(velocity[0,0],velocity[0,1],velocity[1,0])))
        assert sp.cancel(derivative-sp.trace(theta*H)) == 0
    torus = sp.diag(z,1/z)
    torus_theta = torus-sp.trace(torus)*sp.eye(2)/2
    assert sp.expand(D(sp.trace(torus))-2*t-sp.trace((torus_theta-sp.diag(t,-t))*sp.diag(1,-1))) == 0
    controls['trace_gradient_and_target_normalization'] = True

    record = json.loads((ROOT/'e7-quartic-boundary-exploration.json').read_text())
    payload = (ROOT/record['archive']).read_bytes()
    assert sha256(payload).hexdigest() == record['archive_sha256']
    candidate = json.loads(gzip.decompress(payload))
    assert not candidate['reverse_ideal_inclusion_independently_certified']
    assert not record['companion_generic_smoothness_promoted']
    pencil = json.loads((ROOT/'e7-symmetric-pencil-verification.json').read_text())
    assert pencil['resolved_proper_standard_Levis'] == 101
    assert sum((10,12,3,1)) == pencil['remaining_cubic_embeddings'] == 26
    stats = [{'degree':max(sum(e) for e,c in p),'terms':len(p)}
             for p in candidate['selected_members']]
    assert len(candidate['original_generators']) == 19 and len(stats) == 6
    names = ['CAYLEY_FIBER_REVIEW.md','tools/audit_cayley_fiber_review.py',
             'GENERAL_CAYLEY_FIBER_CRITICAL.md','AFFINE_BRANCHING_FIBERS.md',
             'QUADRATIC_BRANCHING_FIBERS.md','MULTIPLIER_BRANCHING_FIBERS.md',
             'SPECTRAL_BRANCHING_FIBERS.md','E7_MULTIPLIER_BOUNDARY.md',
             'E7_SYMMETRIC_PENCIL.md','e7-quartic-boundary-exploration.json',
             'e7-symmetric-pencil-verification.json']
    report = dict(result='passed',date='2026-10-24',uses_Sage=False,
                  controls=controls,control_count=len(controls),
                  selected_certificate_polynomials=stats,
                  E7_resolved=101,E7_unresolved=26,
                  selected_memberships_verified=False,new_group_classification=False,
                  all_simple_adjoint_conjecture_proved=False,literature_priority_claim=False,
                  scope='Exact illustrative ideal and trace checks plus archive/status consistency; general proof and literature review are written separately.',
                  source_sha256={n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names})
    (ROOT/'cayley-fiber-review-audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__':
    main()
