"""Independent native apolar constraints and symbolic moment identities.

No Sage, jet-lattice verifier, or precomputed source basis is imported.
Finite tests support the written general argument; they are not a formal
verification of that argument or an external mathematical review.
"""
from pathlib import Path
from hashlib import sha256
from itertools import permutations, combinations
from math import factorial, prod, gcd, comb
import json
import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

ROOT = Path(__file__).resolve().parents[1]


def partitions(n, length, ceiling=None):
    if length == 0:
        return [()] if n == 0 else []
    ceiling = n if ceiling is None else min(ceiling, n)
    return [(a,) + tail for a in range(ceiling, -1, -1)
            for tail in partitions(n-a, length-1, a)]


def compositions(n, length):
    if length == 0:
        return [()] if n == 0 else []
    return [(a,) + tail for a in range(n+1)
            for tail in compositions(n-a, length-1)]


def dimension(k, n):
    return sum(len(partitions(j, k)) for j in range(n+1)) if n >= 0 else 0


def rank_mod(rows, p):
    if not rows:
        return 0
    pivots = {}
    for source in rows:
        row = [int(x) % p for x in source]
        while any(row):
            j = next(i for i, x in enumerate(row) if x)
            if j not in pivots:
                inv = pow(row[j], -1, p)
                pivots[j] = [x*inv % p for x in row]
                break
            multiple = row[j]
            row = [(a-multiple*b) % p for a, b in zip(row, pivots[j])]
    return len(pivots)


def constraints(k, d):
    # Symmetric dual basis z0^(d-|lambda|) m_lambda(z1,...,zk).
    basis = [(d-n, part) for n in range(d+1) for part in partitions(n, k)]
    j_rows, c_rows = {}, {}
    for column, (power, part) in enumerate(basis):
        for size in range(k+1):
            values = {}
            for e in set(permutations(part)):
                for alpha in compositions(power, size):
                    coefficient = factorial(power) // prod(factorial(a) for a in alpha)
                    target = tuple(e[i] + (alpha[i] if i < size else 0)
                                   for i in range(k))
                    values[target] = values.get(target, 0) + coefficient
            for e, coefficient in values.items():
                key = (size, e)
                j_rows.setdefault(key, [0]*len(basis))[column] += coefficient
                if size < k and any(e[size:]):
                    c_rows.setdefault(key, [0]*len(basis))[column] += coefficient
    return basis, list(j_rows.values()), list(c_rows.values())


def audit_apolar_spaces():
    records = []
    for k, degrees in [(1, range(1, 9)), (2, range(1, 11)),
                       (3, range(3, 11)), (4, [6, 7, 8])]:
        for d in degrees:
            basis, J, C = constraints(k, d)
            def rational_rank(rows):
                if not rows:
                    return 0
                return len(DomainMatrix.from_Matrix(s.Matrix(rows)).convert_to(s.QQ).rref()[1])
            jrank, crank = rational_rank(J), rational_rank(C)
            # Constraint ranks are source-space ranks by the perfect pairing.
            expected_j = len(basis) - dimension(k, d-2**k)
            if k == 1:
                expected_c = len(basis) - dimension(k, d-1)
                defect = 0
            else:
                expected_c = len(basis) - 1 - dimension(k, d-(2**k-1))
                defect = len(partitions(d-(2**k-1), k)) if d >= 2**k-1 else 0
            assert (jrank, crank) == (expected_j, expected_c), (k, d, jrank, crank)
            assert jrank-crank-1 == defect
            primes = [int(s.nextprime(max(d, k))), int(s.nextprime(s.nextprime(max(d, k))))]
            for p in primes:
                assert (rank_mod(J, p), rank_mod(C, p)) == (jrank, crank)
            records.append(dict(k=k, degree=d, source_rank=jrank, ideal_rank=crank,
                                sector_rank=defect, primes=primes))
    return records


def audit_trace():
    M, Y, L, t = s.symbols('M Y L t')
    aa = s.symbols('a1:9')
    def q(j, variable):
        return M if j == 0 else M*variable**j + sum(s.binomial(j, i)*aa[i-1]*variable**(j-i)
                                                  for i in range(1, j+1))
    for j in range(1, 9):
        assert s.expand(q(j, Y+t)-sum(s.binomial(j, r)*q(j-r, Y)*t**r for r in range(j+1))) == 0
    expressions = [q(1, Y)*q(2, Y)-M*q(3, Y),
                   q(2, Y)*q(4, Y)-M*q(6, Y),
                   q(1, Y)*q(2, Y)*q(4, Y)-M**2*q(7, Y)]
    for H in expressions:
        H = s.expand(H)
        assert s.expand(H.subs({a: 0 for a in aa})) == 0
        T = sum(coefficient*q(exponent[0], L) for exponent, coefficient in s.Poly(H, Y).terms())
        rhs = M*H.subs(Y, L) + sum(aa[r-1]*s.diff(H, Y, r).subs(Y, L)/factorial(r)
                                  for r in range(1, s.degree(H, Y)+1))
        assert s.expand(T-rhs) == 0
    return dict(hasse_identities=8, symbolic_mass_kernel_trace_identities=len(expressions))


def audit_first_generator_square():
    m, L, c1, c2 = s.symbols('m L c1 c2')
    f1 = m*L-c1
    f2 = m*(m-1)*L**2/2-(m-1)*c1*L+c2
    q2 = f1**2-2*f2
    a1, a2 = -c1, c1**2-2*c2
    assert s.expand(f1**2-m*q2-(a1**2-m*a2)) == 0
    return dict(symbolic_every_rank_newton_identity=True,
                combined_with_proved_m_annihilator=True,
                first_generator_squared_annihilates_entire_defect=True)


def audit_primitive_classes():
    records = []
    for k in range(2, 5):
        z = s.symbols('z0:'+str(k+1))
        factors = [z[0]-sum(z[i+1] for i in B)
                   for size in range(k) for B in combinations(range(k), size)]
        assert len(factors) == 2**k-1
        # The full-hyperplane restriction can be built independently as
        # the product of NONEMPTY subset sums.
        restricted = s.Poly(1, *z[1:])
        for size in range(1, k+1):
            for B in combinations(range(k), size):
                restricted *= s.Poly(sum(z[i+1] for i in B), *z[1:])
        target = tuple(2**i for i in range(k))
        coefficient = restricted.coeff_monomial(target)
        assert coefficient == 1
        assert all(all(e[i] >= 1 for i in range(k)) for e, value in restricted.terms())
        pairing = prod(factorial(i) for i in target)
        # For k=2,3 explicitly check the substitution as well.
        if k <= 3:
            actual = s.Poly(s.prod(factors).subs(z[0], sum(z[1:])), *z[1:])
            assert actual == restricted
        records.append(dict(k=k, threshold=2**k-1,
                            binary_exponents=list(target), coefficient=int(coefficient),
                            differential_pairing=pairing))
    return records


def audit_triangular_basis():
    records = []
    for k, degrees in [(2, range(3, 11)), (3, range(7, 11)), (4, [15, 16, 17])]:
        z = s.symbols('z1:'+str(k+1))
        product_poly = s.Poly(1, *z)
        for size in range(1, k+1):
            for B in combinations(range(k), size):
                product_poly *= s.Poly(sum(z[i] for i in B), *z)
        for d in degrees:
            # Sort by the claimed leading monomial order, largest variable first.
            tuples = [tuple(reversed(part)) for part in partitions(d-(2**k-1), k)]
            tuples.sort(key=lambda part: tuple(reversed(part)))
            rows = []
            for mu in tuples:
                monomial_symmetric = s.Poly(sum(s.prod(z[i]**e[i] for i in range(k))
                                               for e in set(permutations(mu))), *z)
                poly = product_poly*monomial_symmetric
                rows.append([int(poly.coeff_monomial(tuple(2**i+lam[i] for i in range(k))))
                             for lam in tuples])
            assert all(rows[i][i] == 1 for i in range(len(rows)))
            assert all(rows[i][j] == 0 for i in range(len(rows)) for j in range(i+1, len(rows)))
            assert s.Matrix(rows).det() == 1
            records.append(dict(k=k, degree=d, rank=len(rows),
                                normalized_pairing_determinant=1))
    return records


def audit_cubic_generation():
    m, L, c1, c2, c3 = s.symbols('m L c1 c2 c3')
    f1 = m*L-c1
    f2 = m*(m-1)*L**2/2-(m-1)*c1*L+c2
    f3 = m*(m-1)*(m-2)*L**3/6-(m-1)*(m-2)*c1*L**2/2+(m-2)*c2*L-c3
    R = (m-2)*f1*f2-3*m*f3
    scalar = -(m-2)*c1*c2+3*m*c3
    assert s.expand((f1-c1)*R-2*(m-2)*(c1**2*f2-c2*f1**2)-scalar*(f1-c1)) == 0
    assert s.expand((m-1)*f1**2-2*m*f2-((m-1)*c1**2-2*m*c2)) == 0
    residue_families = []
    for r in range(12):
        u = 2 if r % 2 == 0 else 1
        g = gcd(r-2, 6)
        h = gcd(r, 3)
        numerator, denominator = 3*u, g*h
        assert numerator % denominator == 0
        multiple = numerator // denominator
        assert multiple == (3 if r % 3 == 1 else 1)
        assert gcd(multiple, r if r else 12) == 1
        residue_families.append(dict(residue_mod_12=r, generator_multiple=multiple))
    return dict(symbolic_source_identities=2, residue_families=residue_families,
                order_input='IMG-39 degree-three functional and maximal-minor proof',
                new_claim='R_m generates the complete cubic defect; f_1 kills it')


def weighted_exponents(weights, degree):
    if not weights:
        return [()] if degree == 0 else []
    return [(a,) + tail for a in range(degree//weights[0]+1)
            for tail in weighted_exponents(weights[1:], degree-a*weights[0])]


def first_prime_piece(m, p, elementary, degree=None, correction_class=False):
    top = p if degree is None else degree
    assert top <= m+1 or not elementary
    # This small native implementation is separate from the Sage verifier.
    # Only the last p-primary piece is returned; all preceding image
    # lattices are built exactly to obtain the full coefficient ideal.
    count = m-1 if elementary and m <= top else top
    variables = s.symbols('L ' + ' '.join('v%d' % i for i in range(1, count+1)))
    L = variables[0]
    coefficients = [s.Integer(1)] + list(variables[1:])
    if elementary and m <= top:
        coefficients.append((-1)**(m+1)*sum((-1)**i*coefficients[i]*L**(m-i)
                                           for i in range(m)))
    if elementary:
        source = [sum((-1)**i*comb(m-i, j-i)*coefficients[i]*L**(j-i)
                      for i in range(j+1)) for j in range(1, min(m, top+1))]
    else:
        source = [m*L**j+sum(comb(j, i)*coefficients[i]*L**(j-i)
                             for i in range(1, j+1)) for j in range(1, top+1)]
    generators = coefficients[1:] + source
    degrees = list(range(1, len(coefficients))) + list(range(1, len(source)+1))
    weights = [1] + list(range(1, count+1))
    bases = {0: [(0,)*len(variables)]}
    images = {0: s.Matrix([[1]])}
    answer = None
    for d in range(1, top+1):
        basis = weighted_exponents(weights, d)
        bases[d] = basis
        index = {e: j for j, e in enumerate(basis)}
        columns, ideal_columns = [], []
        for gi, (g, w) in enumerate(zip(generators, degrees)):
            if w > d:
                continue
            op = s.zeros(len(basis), len(bases[d-w]))
            for j, e in enumerate(bases[d-w]):
                for z, value in s.Poly(g, *variables).terms():
                    op[index[tuple(a+b for a, b in zip(e, z))], j] += value
            block = op*images[d-w]
            new = [block[:, j] for j in range(block.cols)]
            columns.extend(new)
            if gi < len(coefficients)-1:
                ideal_columns.extend(new)
        H = hermite_normal_form(s.Matrix.hstack(*columns))
        images[d] = H
        if d == top:
            coordinates = H.inv()*s.Matrix.hstack(*ideal_columns)
            assert all(x.q == 1 for x in coordinates)
            coordinates = s.Matrix([[int(x) for x in row] for row in coordinates.tolist()])
            diagonal = smith_normal_form(hermite_normal_form(coordinates), domain=s.ZZ)
            answer = []
            for i in range(min(diagonal.shape)):
                value = abs(int(diagonal[i, i]))
                if value == 0:
                    continue
                power = 1
                while value % p == 0:
                    power *= p
                    value //= p
                if power > 1:
                    answer.append(power)
            if correction_class:
                assert elementary and top==p+2
                q = [s.Integer(m)]
                for j in range(1,top+1):
                    q.append(s.expand(sum((-1)**(i-1)*source[i-1]*q[j-i]
                                         for i in range(1,min(j,len(source)+1)))
                                      +((-1)**(j-1)*j*source[j-1] if j<=len(source) else 0)))
                polynomial = s.Poly(q[1]*q[p+1]-m*q[p+2],*variables)
                vector = s.Matrix([polynomial.coeff_monomial(prod(v**a for v,a in zip(variables,e)))
                                   for e in basis])
                ideal_basis = hermite_normal_form(s.Matrix.hstack(*ideal_columns))
                solution, parameters = ideal_basis.gauss_jordan_solve(vector)
                assert parameters.rows==0
                order = s.ilcm(*[int(x.q) for x in solution])
                assert (m//p)%order==0
                answer = dict(p_primary_factors=answer,correction_class_order=int(order),
                              annihilator=m//p)
    return answer


def audit_first_prime_transition():
    cases = []
    for m, p in [(2, 2), (3, 3), (6, 3), (9, 3), (5, 5), (25, 5)]:
        elementary = first_prime_piece(m, p, True)
        moment = first_prime_piece(m, p, False)
        multiplicity = sum(len(partitions(p-(2**k-1), k)) for k in range(2, p+1)
                           if p >= 2**k-1)
        power, rest = 1, m
        while rest % p == 0:
            power *= p
            rest //= p
        assert elementary == moment == [power]*multiplicity, (m, p, elementary, moment)
        cases.append(dict(m=m, p=p, degree=p, elementary=elementary, moment=moment,
                          root_equation_in_elementary_model=m==p))
    N, L, b = s.symbols('N L b')
    symbolic = []
    for p in [2, 3, 5, 7]:
        M = p*N
        aa = s.symbols('a1:'+str(p+1))
        q1 = M*L+aa[0]
        qp = M*L**p+sum(comb(p, i)*aa[i-1]*L**(p-i) for i in range(1, p+1))
        g = s.expand((qp.subs(aa[-1], aa[0]**p+p*b)-q1**p)/p)
        expected_g = b+(M-M**p)*L**p/p+sum(s.Rational(comb(p, i), p)
                        *(aa[i-1]-M**(p-i)*aa[0]**i)*L**(p-i) for i in range(1, p))
        assert s.expand(g-expected_g) == 0
        assert all(coefficient.q == 1 for _, coefficient in s.Poly(g, N, L, b, *aa[:-1]).terms())
        lifted = s.expand((M**(p-1)*(qp-aa[-1])-(q1-aa[0])**p)/p)
        rhs = sum(s.Rational(comb(p, i), p)*M**(i-1)*aa[i-1]*(q1-aa[0])**(p-i)
                  for i in range(1, p))
        assert s.expand(lifted-rhs) == 0
        symbolic.append(dict(p=p, all_m_divisible_by_p_symbolic=True,
                             integral_first_coordinate_and_frobenius_lift=True))
    return dict(cases=cases, symbolic_families=symbolic,
                multiplicity_formula_at_d_equals_p_claimed=True)


def audit_prime_boundary_source_ranks():
    """Direct moment words, with coefficient degree; no differential pairing.

    Set M=1 only after homogenizing each word to coefficient degree k.
    Ambient monomials of coefficient cost <=k retain that information.
    """
    records = []
    for p, d in [(p, p) for p in [2, 3, 5, 7]] + [(p, p+1) for p in [3, 5, 7]] + [(3,5)]:
        L = s.Symbol('L')
        aa = s.symbols('a1:'+str(d+1))
        variables = (L,) + aa
        q = [L**j+sum(comb(j, i)*aa[i-1]*L**(j-i)
                      for i in range(1, j+1)) for j in range(1, d+1)]
        generators = list(aa) + q
        words = []
        for e in weighted_exponents(list(range(1, d+1))*2, d):
            polynomial = prod(g**a for g, a in zip(generators, e))
            words.append((sum(e), any(e[:d]),
                          dict(s.Poly(polynomial, *variables).terms())))
        ambient = weighted_exponents([1]+list(range(1, d+1)), d)
        for k in range(1, d):
            basis = [e for e in ambient if sum(e[1:]) <= k]
            J, C = [], []
            for cost, has_input, polynomial in words:
                if cost > k:
                    continue
                row = [polynomial.get(e, 0) for e in basis]
                J.append(row)
                if has_input:
                    C.append(row)
            rational = tuple(len(DomainMatrix.from_Matrix(s.Matrix(rows))
                                 .convert_to(s.QQ).rref()[1]) for rows in (J, C))
            modular = rank_mod(J, p), rank_mod(C, p)
            predicted_j = len(basis)-dimension(k, d-2**k)
            predicted_c = (len(basis)-dimension(k, d-1) if k == 1 else
                           len(basis)-1-dimension(k, d-(2**k-1)))
            assert rational == modular == (predicted_j, predicted_c), (p, k, rational, modular)
            records.append(dict(p=p, degree=d, coefficient_degree=k,
                                ambient_rank=len(basis), rational_ranks=rational,
                                modular_ranks=modular, both_images_p_saturated=True))
    return records


def audit_first_witt_quadratic_relations():
    N, L, b = s.symbols('N L b')
    records = []
    for p in [3, 5, 7, 11]:
        n, r, m = p+1, (p+1)//2, p*N
        aa = [m] + list(s.symbols('a1:'+str(n+1)))
        aa[p] = aa[1]**p+p*b
        q = [m]+[m*L**j+sum(comb(j, i)*aa[i]*L**(j-i)
                            for i in range(1, j+1)) for j in range(1, n+1)]
        g = s.expand((q[p]-q[1]**p)/p)
        coefficients = {i: s.Rational((-1)**i*comb(n, i), p*(2 if i == r else 1))
                        for i in range(2, r+1)}
        assert all(c.q == 1 for c in coefficients.values())
        assert sum(coefficients.values()) == 1
        geometric = sum(q[1]**(2*(r-1-h))*aa[1]**(2*h) for h in range(r))
        identity = (N*(q[n]-aa[n])-n*N*(q[2]-aa[2])*geometric
                    -n*(q[1]*g-aa[1]*b)
                    +sum(c*(q[i]*q[n-i]-aa[i]*aa[n-i]) for i, c in coefficients.items()))
        assert s.expand(identity) == 0
        records.append(dict(p=p, symbolic_all_m_divisible_by_p=True,
                            integral_coefficients={str(i): int(c) for i,c in coefficients.items()},
                            full_polynomial_identity_verified=True))
    return records


def audit_first_mixed_comparison():
    M, m, L = s.symbols('M m L')
    symbolic = []
    for p in [3, 5, 7, 11]:
        aa = [M]+list(s.symbols('a1:'+str(p+1)))
        qp = M*L**p+sum(comb(p,i)*aa[i]*L**(p-i) for i in range(1,p+1))
        q1 = M*L+aa[1]
        W = aa[1]*sum(s.Rational(comb(p,i),p)*(aa[i]-M**(p-i)*aa[1]**i)*L**(p-i)
                     for i in range(1,p))
        source_word = aa[1]*(qp-q1**p-aa[p]+aa[1]**p)
        assert s.expand(source_word-(M-M**p)*aa[1]*L**p-p*W) == 0
        divided_difference = sum(M**(p-1-j)*m**j for j in range(p))
        assert s.expand(m-m**p-M+M**p+(M-m)*(1-divided_difference)) == 0
        c = [s.Integer(1)]+list(s.symbols('c1:'+str(p+2)))
        fs = [sum((-1)**i*s.binomial(p-i,j-i)*c[i]*L**(j-i) for i in range(j+1))
              for j in range(1,p+2)]
        assert s.expand(fs[-1]-(-1)**(p+1)*c[-1]) == 0
        assert s.expand(fs[-2]-sum((-1)**i*c[i]*L**(p-i) for i in range(p+1))) == 0
        symbolic.append(dict(p=p, bockstein_source_identity=True,
                             bockstein_divided_difference=True,
                             root_boundary_redundant_scalar_coordinate=True))
    actual = []
    for m,p in [(3,3),(6,3),(9,3)]:
        elementary = first_prime_piece(m,p,True,p+1)
        moment = first_prime_piece(m,p,False,p+1)
        power, rest = 1,m
        while rest % p == 0:
            power *= p
            rest //= p
        assert elementary == moment == [power]
        actual.append(dict(m=m,p=p,degree=p+1,elementary=elementary,moment=moment,
                           root_equation_in_elementary_model=m==p))
    reduced = [dict(m=m,p=3,degree=5,**first_prime_piece(m,3,True,5,True)) for m in [6,9]]
    return dict(symbolic_families=symbolic, independent_actual_images=actual,
                first_reduced_order_class_checks=reduced)


def audit_tensor_gluing():
    """Build tensor substitutions as integer coefficient dictionaries.

    No polynomial-substitution builder or Sage certificate is imported.
    Indices encode u_i*v_j*v_l, with j<=l.
    """
    perms = list(permutations((1,2,3)))
    subsets = [tuple(B) for n in range(4) for B in combinations((1,2,3),n)]
    local = [(i,(j,l)) for i in (1,2,3) for j in (1,2,3) for l in range(j,4)]
    ambient = [(i,(j,l)) for i in range(4) for j in range(4) for l in range(j,4)]
    def permute(mon, pm):
        f = lambda x: 0 if x==0 else pm[x-1]
        return f(mon[0]),tuple(sorted(f(j) for j in mon[1]))
    def orbits(items, group):
        seen, result = set(), []
        for mon in items:
            if mon in seen:
                continue
            orbit = {permute(mon,pm) for pm in group}
            seen.update(orbit)
            result.append(orbit)
        return result
    def transform(polynomial, substitutions):
        result = {}
        for (i,(j,l)), coefficient in polynomial.items():
            for a,ca in substitutions.get(i,[(i,1)]):
                for b,cb in substitutions.get(j,[(j,1)]):
                    for c,cc in substitutions.get(l,[(l,1)]):
                        mon = a,tuple(sorted((b,c)))
                        result[mon] = result.get(mon,0)+coefficient*ca*cb*cc
        return {mon:c for mon,c in result.items() if c}
    families, representatives = [], []
    for size in range(4):
        B0 = tuple(range(1,size+1))
        group = [pm for pm in perms if set(pm[:size])==set(B0)]
        for orbit in orbits(local,group):
            family = {B:{} for B in subsets}
            for B in subsets:
                if len(B)==size:
                    pm = next(pm for pm in perms if set(pm[:size])==set(B))
                    family[B] = {permute(mon,pm):1 for mon in orbit}
            families.append(family)
            representatives.append((B0,next(mon for mon in local if mon in orbit)))
    constraints = {}
    for B,C in combinations(subsets,2):
        coefficients = {i:int(i in B)-int(i in C) for i in (1,2,3)}
        pivot = next(i for i,c in coefficients.items() if c)
        substitution = {pivot:[(i,-c//coefficients[pivot]) for i,c in coefficients.items()
                               if i!=pivot and c]}
        for column,family in enumerate(families):
            difference = family[B].copy()
            for mon,c in family[C].items():
                difference[mon] = difference.get(mon,0)-c
            for mon,c in transform(difference,substitution).items():
                constraints.setdefault((B,C,mon),[0]*28)[column] += c
    global_polynomials = [{mon:1 for mon in orbit} for orbit in orbits(ambient,perms)]
    restricted = {(B,i):transform(poly,{0:[(j,1) for j in B]})
                  for B in subsets for i,poly in enumerate(global_polynomials)}
    restriction = s.Matrix([[restricted[B,i].get(mon,0) for i in range(11)]
                            for B,mon in representatives])
    compatibility = s.Matrix(list(constraints.values()))
    ideal_constraints = {}
    for B in subsets:
        if len(B)==3:
            continue
        for i in range(11):
            for mon,c in restricted[B,i].items():
                if mon[0] not in B or any(j not in B for j in mon[1]):
                    ideal_constraints.setdefault((B,mon),[0]*11)[i] = c
    ideal = s.Matrix(list(ideal_constraints.values()))
    assert compatibility*restriction == s.zeros(compatibility.rows,restriction.cols)
    def nonzero_smith(mat):
        diagonal = smith_normal_form(mat,domain=s.ZZ)
        return [abs(int(diagonal[i,i])) for i in range(min(diagonal.shape)) if diagonal[i,i]]
    cd,ad,idg = map(nonzero_smith,(compatibility,restriction,ideal))
    assert cd == [1]*16+[2] and ad == [1]*11 and idg == [1]*10
    assert ideal.nullspace()==[s.Matrix([1]+[0]*10)]
    return dict(family_dimension=28,ambient_invariant_dimension=11,
                compatibility_shape=list(compatibility.shape),compatibility_rank=17,
                compatibility_nonunit_smith=[2],ambient_restriction_nonunit_smith=[],
                coefficient_condition_rank=10,coefficient_condition_nonunit_smith=[],
                coefficient_kernel='span(u0*v0^2)',
                method='Native coefficient-dictionary tensor substitutions and SymPy integer Smith forms.')


def main():
    spaces = audit_apolar_spaces()
    trace = audit_trace()
    primitive = audit_primitive_classes()
    triangular = audit_triangular_basis()
    cubic = audit_cubic_generation()
    first_prime = audit_first_prime_transition()
    square = audit_first_generator_square()
    boundary_ranks = audit_prime_boundary_source_ranks()
    witt_quadratics = audit_first_witt_quadratic_relations()
    mixed_comparison = audit_first_mixed_comparison()
    tensor_gluing = audit_tensor_gluing()
    sequence = [sum(len(partitions(d-(2**k-1), k)) for k in range(2, d+1)
                    if d >= 2**k-1) for d in range(21)]
    assert sequence[:11] == [0, 0, 0, 1, 1, 2, 2, 4, 4, 6, 7]
    names = ['TYPE_A_K_STABLE_MULTIPLICITIES.md',
             'tools/audit_type_a_k_stable_multiplicities.py']
    report = dict(result='passed', apolar_spaces=spaces, moment_trace=trace,
        primitive_classes=primitive, complete_cubic_piece=cubic,
        explicit_partition_basis_pairings=triangular,
        first_prime_degree_comparison=first_prime,
        square_zero_first_generator=square,
        prime_boundary_direct_moment_ranks=boundary_ranks,
        first_witt_quadratic_relations=witt_quadratics,
        first_mixed_degree_comparison=mixed_comparison,
        tensor_gluing_certificate=tensor_gluing,
        multiplicity_sequence_through_20=sequence,
        method='Native Python/SymPy restriction constraints, separate rational and modular row elimination, symbolic mass identities, subset-product coefficients.',
        sage_or_elementary_jet_builder_imported=False,
        all_rank_claim_from_finite_checks=False,
        uniform_cubic_annihilator_uses_exact_parameter_independent_tensor_certificate=True,
        full_degree_p_plus_two_multiplicity_formula=False, external_review=False,
        source_sha256={name: sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                       for name in names})
    (ROOT/'type-a-k-stable-multiplicities-audit.json').write_text(
        json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(result='passed', apolar_spaces=len(spaces),
                         moment_trace=trace, primitive_classes=primitive,
                         multiplicity_sequence=sequence), indent=2))


if __name__ == '__main__':
    main()
