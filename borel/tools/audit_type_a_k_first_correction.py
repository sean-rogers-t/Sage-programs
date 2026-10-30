"""Native characteristic-zero checks for the first specialization correction.

Independent column lattices and symbolic coefficient identities. No Sage
or first-correction verifier is imported. Finite checks do not replace
the written uniform proof.
"""
from functools import lru_cache
from hashlib import sha256
from math import comb
from pathlib import Path
import json
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from binary_tensor_certificates_native import binary_tensor_certificate_native

ROOT = Path(__file__).resolve().parents[1]


@lru_cache(None)
def monomials(weights, n):
    if not weights:
        return ((),) if n == 0 else ()
    return tuple((a,)+tail for a in range(n//weights[0]+1)
                 for tail in monomials(weights[1:], n-a*weights[0]))


def bounded_primary(A, p, exponent_bound):
    """Exact p-local Smith factors using the established annihilator bound.

    Work modulo p^(bound+1), so every genuine pivot is detected. At each
    stage eliminate unit pivots; if none remain, divide the entire residual
    matrix by p. This avoids integer Smith coefficient explosion.
    """
    modulus=p**(exponent_bound+1)
    rows=[[int(x)%modulus for x in row] for row in A.tolist()]
    rank=A.rows
    depth=0
    factors=[]
    while rows:
        pivot=next(((i,j) for i,row in enumerate(rows) for j,x in enumerate(row) if x%p),None)
        if pivot is None:
            assert depth<exponent_bound, 'Annihilator bound or full-rank hypothesis failed'
            rows=[[x//p for x in row] for row in rows]
            modulus//=p
            depth+=1
            continue
        i,j=pivot
        rows[0],rows[i]=rows[i],rows[0]
        for row in rows:
            row[0],row[j]=row[j],row[0]
        inv=pow(rows[0][0],-1,modulus)
        top=[x*inv%modulus for x in rows[0]]
        rows=[[(x-row[0]*y)%modulus for x,y in zip(row[1:],top[1:])]
              for row in rows[1:]]
        if depth:
            factors.append(p**depth)
        rank-=1
    assert rank==0
    return factors


def native_image(m, p, offset=2):
    d = p+offset
    count = min(d, m-1)
    variables = s.symbols('L '+ ' '.join('c%d' % i for i in range(1, count+1)))
    L = variables[0]
    cs = [s.Integer(1)]+list(variables[1:])
    if d >= m:
        cs.append((-1)**(m+1)*sum((-1)**i*cs[i]*L**(m-i) for i in range(m)))
    fs = [sum((-1)**i*comb(m-i, j-i)*cs[i]*L**(j-i) for i in range(j+1))
          for j in range(1, min(d, m-1)+1)]
    generators = cs[1:]+fs
    degrees = list(range(1, len(cs)))+list(range(1, len(fs)+1))
    weights = tuple([1]+list(range(1, count+1)))
    bases, images = {0: [(0,)*len(variables)]}, {0: s.Matrix([[1]])}
    for n in range(1, d+1):
        basis = monomials(weights, n)
        bases[n] = basis
        index = {e: i for i, e in enumerate(basis)}
        source, coefficient = [], []
        for gi, (g, degree) in enumerate(zip(generators, degrees)):
            if degree > n:
                continue
            terms = s.Poly(g, *variables).terms()
            op = s.zeros(len(basis), len(bases[n-degree]))
            for col, e in enumerate(bases[n-degree]):
                for f, value in terms:
                    op[index[tuple(a+b for a, b in zip(e, f))], col] += value
            block = op*images[n-degree]
            source.extend(block.col(j) for j in range(block.cols))
            if gi < len(cs)-1:
                coefficient.extend(block.col(j) for j in range(block.cols))
        # mL=f_1+c_1 puts m^r times every ambient L^r monomial in J.
        # The diagonal sublattice gives this certified determinant multiple;
        # modular HNF avoids coefficient explosion without a Smith prediction.
        determinant_multiple = m**sum(e[0] for e in basis)
        J = hermite_normal_form(s.Matrix.hstack(*source), D=determinant_multiple)
        images[n] = J
        C = s.Matrix.hstack(*coefficient)
    coordinates = J.inv()*C
    assert all(x.q == 1 for x in coordinates)
    coordinates = s.Matrix([[int(x) for x in row] for row in coordinates.tolist()])
    # The ordinary monomial is last, and column HNF is upper triangular.
    # Hence the first n-1 source columns are an integral ordinary-kernel basis.
    # The previously proved m-annihilator bounds its coefficient-ideal index
    # by m^(n-1). This avoids uncontrolled intermediate HNF coefficients;
    # it does not assume the multiplicities under test.
    if d < m:
        assert all(x == 0 for x in coordinates[-1, :])
        relative = hermite_normal_form(coordinates[:-1, :], D=m**(J.rows-1))
    else:
        # Above the root boundary the ordinary homogeneous image is zero.
        relative = hermite_normal_form(coordinates, D=m**J.rows)
    bound=0
    temp=m
    while temp%p==0:
        bound+=1
        temp//=p
    factors=bounded_primary(relative,p,bound)
    return dict(m=m, p=p, degree=d, ambient_rank=len(basis), primary=factors,
                method='Independent elementary lattices and exact modular Smith through p^(v_p(m)+1); uses the earlier m-annihilator bound')


def symbolic_trace_controls():
    records = []
    M, L = s.symbols('M L')
    for p in [3, 5, 7]:
        d = p+2
        aa = [M]+list(s.symbols('a1:'+str(d+1)))
        q = [M]+[M*L**j+sum(comb(j, i)*aa[i]*L**(j-i) for i in range(1, j+1))
                  for j in range(1, d+1)]
        trace = sum(comb(2, i)*comb(p, j)*aa[i]*aa[j]*q[d-i-j]
                    for i in range(3) for j in range(p+1))-M**2*q[d]
        correction = 0
        for r in range(1, d+1):
            correction += aa[r]*sum(comb(2, i)*comb(p, r-i)*q[2-i]*q[p-r+i]
                                    for i in range(3) if 0 <= r-i <= p)
        assert s.expand(trace-correction-M*(q[2]*q[p]-M*q[d])) == 0
        # In the FORMAL certificate modulo (p,M), only these two terms survive.
        Q = [M]+list(s.symbols('Q1:'+str(d+1)))
        formal = sum(comb(2, i)*comb(p, j)*aa[i]*aa[j]*Q[d-i-j]
                     for i in range(3) for j in range(p+1))-M**2*Q[d]
        for r in range(1, d+1):
            formal -= aa[r]*sum(comb(2, i)*comb(p, r-i)*Q[2-i]*Q[p-r+i]
                                for i in range(3) if 0 <= r-i <= p)
        remainder = s.Poly(s.expand(formal-2*aa[1]*Q[1]*(aa[p]-Q[p])).subs(M, 0),
                          *(aa[1:]+Q[1:]))
        assert all(int(c) % p == 0 for c in remainder.coeffs())
        records.append(dict(p=p, symbolic_mass_trace_identity=True,
                            formal_certificate_reduction_mod_p_and_mass=True))
    return records


def high_degree_certificate():
    # Replay the dimension-independent elimination on its four bottom rows.
    x2, x3, x4, y1, y2, y3 = s.symbols('x2 x3 x4 y1 y2 y3')
    rows = [y3+2*y2, y3+3*y2, y2+2*y1,
            x4+3*x3+3*y2, x4+4*x3+6*y2,
            x3+3*x2+3*y1,
            x4+4*x3+4*x2+2*y2+4*y1]
    variables = (x2, x3, x4, y1, y2, y3)
    A = s.Matrix([[s.expand(row).coeff(v) for v in variables] for row in rows])
    smith = smith_normal_form(A, domain=s.ZZ)
    nonunits = [abs(int(smith[i, i])) for i in range(6) if abs(smith[i, i]) > 1]
    assert A.rank() == 6 and all(x & (x-1) == 0 for x in nonunits)
    return dict(rank=6, nonunit_smith=nonunits,
                subsequent_rows_eliminate_recursively_with_unit_pivots=True)


def two_coordinate_injectivity():
    ux, uy, uz, vx, vy, vz = s.symbols('ux uy uz vx vy vz')
    a, b, c = s.symbols('a b c')
    f = a*ux*vy*vz+b*uy*vx*vz+c*uz*vx*vy
    f = s.Poly(s.expand(f.subs({ux: uy+uz, vx: vy+vz})), uy, uz, vy, vz)
    A = s.Matrix([[term.coeff(t) for t in (a, b, c)] for term in f.coeffs()])
    smith = smith_normal_form(A, domain=s.ZZ)
    assert [abs(int(smith[i, i])) for i in range(3)] == [1, 1, 1]
    return dict(rank=3, smith=[1, 1, 1], all_characteristics=True)


def main():
    cases = []
    for m, p in [(3, 3), (6, 3), (9, 3), (5, 5), (10, 5)]:
        cases.append(native_image(m, p))
        print(json.dumps(cases[-1]), flush=True)
    assert [row['primary'] for row in cases] == [[], [3], [3, 9], [5, 5], [5, 5, 5]]
    tensors = [binary_tensor_certificate_native(k, 2) for k in [2, 3]]
    assert all(t['gluing_obstruction_dimension'] == 0 for t in tensors)
    assert all(not t['restriction_nonunits'] for t in tensors)
    assert tensors[1]['coefficient_kernel_dimension'] == 1
    report = dict(result='passed', native_elementary_cases=cases,
                  tensor_certificates=tensors,
                  formal_trace=symbolic_trace_controls(),
                  high_coefficient_elimination=high_degree_certificate(),
                  overlap_injectivity=two_coordinate_injectivity(),
                  sage_imported=False, external_review=False,
                  all_rank_claim_from_finite_checks=False)
    names = ['TYPE_A_K_FIRST_CORRECTION.md', 'tools/audit_type_a_k_first_correction.py',
             'tools/binary_tensor_certificates_native.py']
    report['source_sha256'] = {name: sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                             for name in names}
    (ROOT/'type-a-k-first-correction-audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
