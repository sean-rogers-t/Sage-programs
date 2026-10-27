"""Fresh rational audit of the final E8 degree proof's local geometry.

No Sage, F4 engine, or earlier audit helpers are imported. The finite-field
ideal memberships remain a separately disclosed dependency, not rechecked
by this script. Geometry-to-group implications are audited in the note.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations, permutations
import gzip
import json
import sympy as sp
from sympy.polys.rings import ring

ROOT = Path(__file__).resolve().parents[1]


def determinant(matrix):
    n = len(matrix)
    value = 0
    for perm in permutations(range(n)):
        sign = (-1) ** sum(perm[i] > perm[j] for i in range(n) for j in range(i + 1, n))
        term = sign
        for i, j in enumerate(perm):
            term *= matrix[i][j]
        value += term
    return value


def evaluate(poly, point):
    return sum(coefficient * sp.prod(v ** e for v, e in zip(point, exponent))
               for exponent, coefficient in poly.items())


def main():
    P, x, y, w, a, A, B, C, h = ring('x,y,w,a,A,B,C,h', sp.QQ)
    variables = P.gens
    p, q, r = x*y+w*a, x*w+y*a, x*a+y*w
    f = (x*x+y*y+w*w+a*a-p*p-q*q-r*r+x*y*w*a
         -A*q*r-B*p*r-C*p*q-B*C*p-A*C*q-A*B*r-A*B*C-4)
    equations = [f.diff(v) for v in (x, y, w)] + [A*A+B*B+C*C-A*B*C-4]
    jacobian = [[g.diff(v) for v in variables[:7]] for g in equations]
    originals = list(equations)
    for cols in combinations(range(7), 4):
        minor = determinant([[row[j] for j in cols] for row in jacobian])
        if minor and minor not in originals:
            originals.append(minor)
    assert len(originals) == 38
    for prime in (32003, 32009):
        data = json.loads(gzip.decompress((ROOT/f'e8-last-{prime}-certificate.json.gz').read_bytes()))
        decoded = [P.from_dict({tuple(e)+(0,): sp.QQ(c) for e, c in row})
                   for row in data['original_affine_generators']]
        assert set(decoded) == set(originals)

    def homogenize(poly):
        degree = max(map(sum, poly))
        return P.from_dict({e[:7]+(degree-sum(e),): c for e, c in poly.items()})

    homogeneous = [homogenize(g) for g in originals]
    orders = [min(sum(e[4:]) for e in g) for g in homogeneous]
    assert orders[:3] == [1, 1, 1] and min(orders[3:]) == 3
    # Extract the linear normal terms from freshly constructed equations.
    def linear_coefficient(poly, normal_index):
        return P.from_dict({e[:4]+(0,)*4: c for e, c in poly.items()
                            if sum(e[4:]) == 1 and e[normal_index] == 1})

    matrix = [[linear_coefficient(g, j) for j in (4, 5, 6)] for g in homogeneous[:3]]
    linear_h = [linear_coefficient(g, 7) for g in homogeneous[:3]]
    D = determinant(matrix)
    assert D == -2*p*q*r*(-a**3+a*(x*x+y*y+w*w)-2*x*y*w)
    adjugate = [[(-1)**(i+j)*determinant(
        [[matrix[k][l] for l in range(3) if l != i] for k in range(3) if k != j])
        for j in range(3)] for i in range(3)]
    ell = [sum(adjugate[i][j]*linear_h[j] for j in range(3)) for i in range(3)]
    c = D*sum(v*v for v in ell)+ell[0]*ell[1]*ell[2]-4*D**3
    witness = (1, 2, 3, 4, 0, 0, 0, 0)
    d_at = evaluate(D, witness)
    c_at = evaluate(c, witness)
    assert d_at == 61600 and c_at != 0
    normal_cubic = sp.cancel(c_at/d_at**3)
    # A nonzero cubic coefficient makes H=h^3 times a unit after formal
    # elimination. All remaining generators have order >=3, so the local
    # algebra is exactly F[[h]]/(h^3), independently of the global bound.

    witnesses = [
        ('X_A', (1, -6, 2, 3, 0, 16, -9, 0)),
        ('X_B', (1, sp.Rational(-2, 3), 2, 3, 16, 0, -5, 0)),
        ('X_C', (1, 2, 3, -6, 16, -9, 0, 0)),
    ]
    checks = []
    for label, point in witnesses:
        assert all(evaluate(g, point) == 0 for g in homogeneous), label
        rows = [[evaluate(g.diff(v), point) for v in variables] for g in homogeneous]
        rank = sp.Matrix(rows).rank()
        assert rank == 4, (label, rank)
        checks.append(dict(component=label, point=list(map(str, point)), full_jacobian_rank=rank))

    result = dict(
        result='passed', original_generators_rebuilt=38,
        independent_original_input_matches_both_certificates=True,
        normal_orders=orders,
        normal_determinant_witness=str(d_at), cubic_numerator_witness=str(c_at),
        eliminated_H_cubic_coefficient_witness=str(normal_cubic),
        generic_linear_component_length=3,
        sextic_smooth_witnesses=checks,
        finite_field_memberships_independently_replayed=False,
        scope='Exact rational local audit; the group transfer is a written geometric argument. No new generic nonemptiness or map degree claim.')
    names = ['tools/audit_e8_degree_transfer.py', 'E8_DEGREE_TRANSFER_AUDIT.md',
             'E8_LAST_ADJOINT_FIBER.md', 'ADJOINT_GENERIC_FIBER_THEOREM.md']
    result['source_sha256'] = {n: sha256((ROOT/n).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                               for n in names}
    result['raw_sha256'] = {f'e8-last-{p}-certificate.json.gz': sha256((ROOT/f'e8-last-{p}-certificate.json.gz').read_bytes()).hexdigest()
                            for p in (32003, 32009)}
    (ROOT/'e8-degree-transfer-audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'raw_sha256')}, indent=2))


if __name__ == '__main__':
    main()
