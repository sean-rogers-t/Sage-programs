"""Actual elementary images and a separate moment/Witt lattice comparison.

Finite corroboration of TYPE_A_K_FIRST_CORRECTION.md, not an all-rank
proof. No public API is added. Run from borel/ with Sage Python.
"""
from hashlib import sha256
from math import comb
from pathlib import Path
import json

from sage.all import ZZ, QQ, PolynomialRing, matrix, vector, binomial, lcm
from verify_type_a_k_stable_multiplicities import exponents, expected, increasing_partitions
from binary_tensor_certificates import binary_tensor_certificate

ROOT = Path(__file__).resolve().parents[1]


def build_image(m, d, model):
    if model == 'moment':
        count = d
    elif model == 'auxiliary':
        count = d
    elif model == 'unrooted':
        count = min(m, d)
    else:
        assert model == 'actual'
        count = min(m-1, d)
    R = PolynomialRing(ZZ, names=['L']+['c%d' % i for i in range(1, count+1)])
    L, *coefficients = R.gens()
    inputs = list(coefficients)
    cs = [R(1)]+inputs
    if model == 'actual' and d >= m:
        cm = (-1)**(m+1)*sum((-1)**i*cs[i]*L**(m-i) for i in range(m))
        cs.append(cm)
        inputs.append(cm)
    if model == 'moment':
        source = [m*L**j+sum(comb(j, i)*cs[i]*L**(j-i) for i in range(1, j+1))
                  for j in range(1, d+1)]
    else:
        top = d if model == 'auxiliary' else min(d, m-(model == 'actual'))
        source = [sum((-1)**i*binomial(m-i, j-i)*cs[i]*L**(j-i)
                      for i in range(j+1)) for j in range(1, top+1)]
    generators = inputs+source
    weights = list(range(1, len(inputs)+1))+list(range(1, len(source)+1))
    bases = {0: [(0,)*(count+1)]}
    images = {0: matrix(ZZ, [[1]])}
    for degree in range(1, d+1):
        basis = exponents([1]+list(range(1, count+1)), degree)
        bases[degree] = basis
        index = {e: i for i, e in enumerate(basis)}
        source_rows, ideal_rows = [], []
        for gi, (g, w) in enumerate(zip(generators, weights)):
            if w > degree:
                continue
            op = matrix(ZZ, len(bases[degree-w]), len(basis), sparse=True)
            for row, e in enumerate(bases[degree-w]):
                for z, value in g.dict().items():
                    op[row, index[tuple(x+y for x, y in zip(e, z))]] += value
            block = images[degree-w]*op
            source_rows.extend(block.rows())
            if gi < len(inputs):
                ideal_rows.extend(block.rows())
        J = matrix(ZZ, source_rows).hermite_form(include_zero_rows=False)
        C = matrix(ZZ, ideal_rows)
        images[degree] = J
    return R, basis, J, C, source


def primary(J, C, p):
    coords = J.transpose().solve_right(C.transpose()).transpose()
    assert all(x.denominator() == 1 for x in coords.list())
    return [int(p**ZZ(x).valuation(p)) for x in matrix(ZZ, coords).elementary_divisors()
            if x and ZZ(x).valuation(p) > 0]


def predicted(m, p, root=False):
    a = expected(p+2)
    if root:
        assert m == p
        return [p]*(a-2)
    s = ZZ(m).valuation(p)
    return ([int(p**(s-1))] if s > 1 else [])+[int(p**s)]*(a-1)


def elementary_controls():
    records = []
    for p, m in [(3, 3), (3, 6), (3, 9), (5, 5), (5, 10), (5, 25),
                 (7, 7), (7, 14), (7, 49), (11, 22)]:
        R, basis, J, C, source = build_image(m, p+2, 'actual')
        factors = primary(J, C, p)
        assert factors == predicted(m, p, root=m == p), (m, p, factors)
        q = [R(m)]
        d = p+2
        for j in range(1, d+1):
            q.append(sum((-1)**(i-1)*source[i-1]*q[j-i]
                         for i in range(1, min(j, len(source)+1)))
                     + ((-1)**(j-1)*j*source[j-1] if j <= len(source) else 0))
        representatives = [q[i]*q[d-i]-m*q[d] for i in range(3, (p+1)//2+1)]
        if m > p:
            g = R((q[p]-q[1]**p)/p)
            representatives += [q[1]*q[p+1]-m*q[d], q[2]*g-((m-m**p)//p)*q[d]]
        for k in range(3, d+1):
            remainder = d-(2**k-1)
            if remainder < 0:
                break
            for part in increasing_partitions(remainder, k):
                representatives.append(R.prod(q[2**i+part[i]] for i in range(k))-m**(k-1)*q[d])
        index = {e: i for i, e in enumerate(basis)}
        extra = []
        for f in representatives+([q[d]] if m > p else []):
            row = [ZZ(0)]*len(basis)
            for e, value in f.dict().items():
                row[index[tuple(e)]] = value
            extra.append(row)
        spanning = matrix(ZZ, C.rows()+extra)
        coords = J.transpose().solve_right(spanning.transpose()).transpose()
        assert all(x.denominator() == 1 for x in coords.list())
        diagonal = matrix(ZZ, coords).elementary_divisors()
        assert len([x for x in diagonal if x]) == J.nrows()
        assert all(ZZ(x).valuation(p) == 0 for x in diagonal if x), (m, p, diagonal)
        row = dict(m=m, p=p, degree=p+2, ambient_rank=len(basis), primary=factors,
                   explicit_representatives_generate=True,
                   representative_count=len(representatives))
        records.append(row)
        print(json.dumps(row), flush=True)
    return records


def root_controls():
    records = []
    for p in [3, 5, 7]:
        stages = {}
        for model in ['auxiliary', 'unrooted', 'actual']:
            R, basis, J, C, source = build_image(p, p+2, model)
            stages[model] = primary(J, C, p)
        assert stages['auxiliary'] == [p]*(expected(p+2)-1)
        assert stages['unrooted'] == stages['actual'] == [p]*(expected(p+2)-2)
        records.append(dict(p=p, stages=stages,
                            loss_occurs_at_coefficient_deletion=True,
                            root_equation_preserves_remaining_defect=True))
    return records


def moment_witt_controls():
    records = []
    for p, m in [(3, 3), (3, 9), (5, 10), (5, 25), (7, 49)]:
        d = p+2
        R, basis, J, C, source = build_image(m, d, 'moment')
        # Remove scalar coefficient coordinates before adjoining the Witt direction.
        keep = [i for i, e in enumerate(basis) if e[0] > 0]
        J = J.matrix_from_columns(keep).hermite_form(include_zero_rows=False)
        C = C.matrix_from_columns(keep).hermite_form(include_zero_rows=False)
        basis = [basis[i] for i in keep]
        index = {e: i for i, e in enumerate(basis)}
        K = R.change_ring(QQ)
        L, *aa = K.gens()
        a = [K(m)]+aa
        q = [K(m)]+[K(f) for f in source]
        alpha = (m-m**p)//p
        H = alpha*L**p+sum((comb(p, i)//p)*(a[i]-m**(p-i)*a[1]**i)*L**(p-i)
                          for i in range(1, p))
        b = (a[p]-a[1]**p)/p
        V = b*(q[2]-a[2])
        v, w1, w2, w3 = m*b*a[1]*L, a[2]*H, a[1]**2*H, a[1]*q[1]*H
        E = q[2]*H-alpha*q[d]
        R1, R2 = q[1]*q[p+1]-m*q[d], q[2]*q[p]-m*q[d]

        def vec(f):
            out = [QQ(0)]*len(basis)
            for e, c in f.dict().items():
                if e[0]:
                    out[index[tuple(e)]] = c
            return vector(QQ, out)

        def order_in(f, lattice):
            coords = lattice.transpose().solve_right(vec(f))
            return int(ZZ(lcm(x.denominator() for x in coords)).valuation(p))

        assert primary(J, C, p) == [int(p**ZZ(m).valuation(p))]*expected(d)
        image_orders = [order_in(f, J) for f in [v, w1, w2, w3, E]]
        assert image_orders == [0, 1, 0, 0, 1]
        for f in [v-QQ(m)/(2*p)*R2,
                  w2-QQ(m)/p*(R1-R2/2), w3+QQ(m)/(2*p)*R2,
                  p*E-R2]:
            assert order_in(f, C) == 0
        # No nonzero combination of the two Bockstein directions is in J.
        for x in range(p):
            for y in range(p):
                if x or y:
                    assert order_in(x*w1+y*E, J) == 1
        new_source = [p*row for row in J.rows()]+[p*vec(f) for f in [V, v, w1, w2, w3, E]]
        new_ideal = [p*row for row in C.rows()]+[p*vec(f) for f in [V, v, w1, w2, w3]]
        JJ = matrix(ZZ, new_source).hermite_form(include_zero_rows=False)
        CC = matrix(ZZ, new_ideal).hermite_form(include_zero_rows=False)
        factors = primary(JJ, CC, p)
        assert factors == predicted(m, p)
        records.append(dict(m=m, p=p, old_moment_factors=primary(J, C, p),
                            new_witt_factors=factors,
                            new_direction_image_order_valuations=image_orders,
                            full_ideal_class_identities=True,
                            two_bockstein_classes_independent=True))
    return records


def main():
    tensors = [binary_tensor_certificate(k, 2) for k in [2, 3]]
    assert all(t['gluing_obstruction_dimension'] == 0 for t in tensors)
    assert all(not t['restriction_nonunits'] for t in tensors)
    assert tensors[0]['compatibility_nonunits'] == []
    assert tensors[1]['compatibility_nonunits'] == [2]
    assert tensors[1]['coefficient_kernel_dimension'] == 1
    actual = elementary_controls()
    roots = root_controls()
    comparison = moment_witt_controls()
    names = ['TYPE_A_K_FIRST_CORRECTION.md', 'tools/verify_type_a_k_first_correction.py',
             'tools/verify_type_a_k_stable_multiplicities.py',
             'tools/binary_tensor_certificates.py']
    report = dict(result='passed', actual_elementary_cases=actual,
                  root_boundary_stages=roots, moment_witt_comparisons=comparison,
                  tensor_certificates=tensors,
                  all_rank_claim_from_finite_checks=False,
                  source_sha256={name: sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                                 for name in names})
    (ROOT/'type-a-k-first-correction-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(result='passed', actual_cases=len(actual),
                         root_stage_cases=len(roots), moment_comparisons=len(comparison))))


if __name__ == '__main__':
    main()
