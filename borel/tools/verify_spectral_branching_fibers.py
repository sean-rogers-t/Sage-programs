"""Research-only recognition and exact scheme checks for the spectral criterion."""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

from sage.all import QQ, PolynomialRing, matrix, prod, CyclotomicField, vector

ROOT = Path(__file__).resolve().parents[1]


def recognize_shape(record, center):
    """Read exact central-charge support; return only a shape certificate.

    This does NOT certify the spectral or boundary hypotheses of the theorem.
    Indices refer to the retained-node order in the existing charge records.
    """
    terms = {}
    for t in record['terms']:
        key = (tuple(t['charge']), tuple(t['exponent']))
        terms[key] = terms.get(key, 0) + t['coefficient']
    terms = {k: v for k, v in terms.items() if v}
    quadratic = [(charge, coeff) for (charge, e), coeff in terms.items()
                 if e[center] == 2 and sum(e) == 2]
    if len(quadratic) != 1:
        return dict(accepted=False, reason='quadratic center coefficient is not a Laurent unit')
    edges = set()
    for (charge, e), coeff in terms.items():
        leaf = [(i, n) for i, n in enumerate(e) if i != center and n]
        if not leaf:
            if e[center] > 2:
                return dict(accepted=False, reason='central degree exceeds two')
            continue
        if len(leaf) != 2 or any(n != 1 for i, n in leaf) or e[center] > 1:
            return dict(accepted=False, reason='non-pair leaf monomial', exponent=list(e))
        edges.add(tuple(sorted(i for i, n in leaf)))
    graph = {}
    for i, j in edges:
        graph.setdefault(i, set()).add(j)
        graph.setdefault(j, set()).add(i)
    colors = {}; components = []
    for start in sorted(graph):
        if start in colors:
            continue
        colors[start] = 0; stack = [start]; sides = [[], []]
        while stack:
            i = stack.pop(); sides[colors[i]].append(i)
            for j in sorted(graph[i]):
                if j in colors and colors[j] == colors[i]:
                    return dict(accepted=False, reason='odd cycle')
                if j not in colors:
                    colors[j] = 1-colors[i]; stack.append(j)
        components.append(sides)
    for flips in product((0, 1), repeat=len(components)):
        left = sorted(i for pair, f in zip(components, flips) for i in pair[f])
        right = sorted(i for pair, f in zip(components, flips) for i in pair[1-f])
        if len(left) == len(right):
            return dict(accepted=True, left=left, right=right,
                        free=[i for i in range(len(record['nodes'])) if i != center and i not in graph],
                        pair_count=len(left), edges=[list(e) for e in sorted(edges)],
                        quadratic_center_coefficient_is_Laurent_unit=True,
                        full_central_charges_retained=True,
                        spectral_and_boundary_hypotheses_certified=False)
    return dict(accepted=False, reason='no balanced bipartition')


def character_checks():
    records = json.loads((ROOT/'branching-norm-fibers-verification.json').read_text())['adjoint_mechanism_tests']
    results = []
    for rec in records:
        center = rec['nodes'].index(2)
        result = recognize_shape(rec, center)
        results.append(dict(nodes=rec['nodes'], **result))
    assert [r['accepted'] for r in results] == [False, False, True]
    rec = records[-1]; shape = results[-1]; center = rec['nodes'].index(2)
    F = PolynomialRing(QQ, 'z').fraction_field(); z = F.gen()
    P = PolynomialRing(F, ['r'+str(i) for i in range(len(rec['nodes']))]); variables = P.gens()
    phi = sum(t['coefficient'] * z**t['charge'][0] * prod(v**e for v, e in zip(variables, t['exponent']))
              for t in rec['terms'])
    x = variables[center]
    M = matrix(P, [[phi.derivative(variables[i]).derivative(variables[j])
                    for j in shape['right']] for i in shape['left']])
    g = phi.subs({variables[i]: 0 for i in shape['left']+shape['right']})
    assert phi == g + sum(variables[i]*M[a,b]*variables[j]
                         for a,i in enumerate(shape['left']) for b,j in enumerate(shape['right']))
    A = M.subs({x: 0}); B = M.apply_map(lambda f: f.derivative(x))
    assert B.det() != 0
    L = -B.inverse()*A; S = z**3+z**-3
    assert g == x*x+S*x-2 and M.det() == -x*(x+S)
    assert L*(L+S*matrix.identity(P,2)) == 0
    assert L.rank() == 1 and (L+S*matrix.identity(P,2)).rank() == 1
    assert g.subs({x: 0}) == g.subs({x: -S}) == -2
    return dict(shape_tests=results, E6_matrix=[list(map(str,row)) for row in M.rows()],
                determinant=str(M.det()), semisimple_on_S_nonzero=True,
                eigenvalue_multiplicities=[1,1], spectral_target_maps_zero=True)


def exact_ideals():
    field = PolynomialRing(QQ, ['b','lam','c']).fraction_field(); b,lam,c = field.gens()
    R = PolynomialRing(field, ['x','u1','u2','v1','v2','z']); x,u1,u2,v1,v2,z = R.gens()
    pairing = u1*v1+u2*v2
    phi = x*x+b*x+z+(x-lam)*pairing
    I = R.ideal([phi.derivative(t) for t in (x,u1,u2,v1,v2)]+[z-c])
    ordinary = R.ideal([x+b/2,u1,u2,v1,v2,z-c])
    spectral = R.ideal([x-lam,pairing+2*lam+b,z-c])
    assert I == ordinary.intersection(spectral)
    assert ordinary+spectral == R.ideal(1)
    # Finite model m=2 corroborates, but the note proves all block sizes.
    S = PolynomialRing(QQ, ['X','u','v']); X,u,v = S.gens()
    bad = S.ideal([2*X+u*v,X*u,X*v])
    assert bad == S.ideal([2*X+u*v,u*u*v,u*v*v])
    assert u*v not in bad and (u*v)**2 in bad
    # Uniform boundary identity, including nilpotents, in the E6 completion.
    T = PolynomialRing(QQ, ['x','u','v','S','dS']); x,u,v,s,ds = T.gens()
    phi = x*(x+s)*(1-u*v)
    boundary = T.ideal([phi.derivative(t) for t in (x,u,v)]+[s])
    assert x*ds*(1-u*v) in boundary
    # One local active coordinate suffices to get the exact spectral ideal.
    V = PolynomialRing(field, ['x','u','v','h']); x,u,v,h = V.gens()
    original = V.ideal([2*x+b+u*v,(x-lam)*u,(x-lam)*v,u*h-1])
    expected = V.ideal([x-lam,u*v+2*lam+b,u*h-1])
    assert original == expected
    return dict(two_pair_generic_ideal_is_disjoint_branch_intersection=True,
                active_coordinate_localization_exact=True,
                collision_nilpotent_uv_nonzero=True, collision_nilpotent_uv_square_zero=True,
                E6_boundary_target_zero_in_quotient_ring=True,
                arbitrary_m_claim_comes_from_written_proof=True)


def main():
    report = dict(result='passed', character=character_checks(), ideals=exact_ideals(),
                  reductive_realization=reductive_realization(),
                  all_exceptional_classification=False, new_exceptional_case_claim=False,
                  simple_group_counterexample_claim=False,
                  public_API_changed=False, literature_priority_claim=False)
    names = ['SPECTRAL_BRANCHING_FIBERS.md','tools/verify_spectral_branching_fibers.py',
             'branching-norm-fibers-verification.json','GENERAL_CAYLEY_FIBER_CRITICAL.md']
    report['source_sha256'] = {n:sha256((ROOT/n).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for n in names}
    (ROOT/'spectral-branching-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


def reductive_realization():
    v = [(1,0,0,0),(-1,0,0,0)]
    w = [(0,1,0,0),(0,-1,1,0),(0,0,-1,0)]
    dual = [tuple(-a for a in ww) for ww in w]
    tensor = lambda aa,bb:[tuple(a+b for a,b in zip(t,r)) for t in aa for r in bb]
    ww = tensor(w,dual)
    weights = tensor(v,v)+v+v+[(0,0,0,1)]+ww+tensor(v,ww)
    assert len(weights) == 36
    B = matrix(QQ,4,4,lambda i,j:sum(t[i]*t[j] for t in weights))
    assert B == matrix(QQ,[[30,0,0,0],[0,36,-18,0],[0,-18,36,0],[0,0,0,1]])
    K = CyclotomicField(3); omega = K.gen()
    F = PolynomialRing(K,'c').fraction_field(); c = F.gen()
    params = [F(omega),F(1),F(omega),c]
    vals = [prod(t**e for t,e in zip(params,wgt)) for wgt in weights]
    grad = vector(F,[sum(vv*ww[i] for vv,ww in zip(vals,weights)) for i in range(4)])
    assert B.change_ring(F).inverse()*grad == vector(F,[0,0,0,c])
    Hess = matrix(F,4,4,lambda i,j:sum(vv*ww[i]*ww[j] for vv,ww in zip(vals,weights)))
    assert Hess.rank() == 2
    eigen2 = [omega,omega**-1]; eigen3 = [K(1),omega,omega**-1]
    assert len(set(eigen2)) == 2 and len(set(eigen3)) == 3
    # All root directions have zero derivative: target is central and source
    # regular semisimple, so equivariance divides by nonzero root(g)-1.
    return dict(group='SL2 x SL3 x Gm',representation_dimension=36,
                trace_Gram=[list(map(int,row)) for row in B.rows()],
                regular_semisimple_source=True,central_target_symbolic=True,
                full_differential_rank=2,regular_fiber_local_dimension=9,
                regular_fiber_tangent_dimension=10,group_is_simple=False,
                representation_is_adjoint=False)


if __name__ == '__main__':
    main()
