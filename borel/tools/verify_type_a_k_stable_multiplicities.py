"""Finite elementary-jet controls for the written stable-range theorem.

This deliberately builds the actual homogeneous image over ZZ, without
the moment-coordinate or apolar calculation used in the proof. No public
API is changed. The finite checks are not the arbitrary-rank proof.
"""
from pathlib import Path
from hashlib import sha256
from math import comb
from itertools import permutations, combinations
import json
from sage.all import ZZ, PolynomialRing, matrix

ROOT = Path(__file__).resolve().parents[1]


def exponents(weights, degree):
    if not weights:
        return [()] if degree == 0 else []
    return [(a,) + tail for a in range(degree // weights[0] + 1)
            for tail in exponents(weights[1:], degree - a * weights[0])]


def partitions_bounded(k, n):
    if n < 0:
        return 0
    counts = [1] + [0] * n
    for part in range(1, k + 1):
        for j in range(part, n + 1):
            counts[j] += counts[j - part]
    return counts[n]


def expected(d):
    return sum(partitions_bounded(k, d - (2**k - 1))
               for k in range(2, d + 1))


def increasing_partitions(total, length, lower=0):
    if length == 0:
        return [()] if total == 0 else []
    return [(a,) + tail for a in range(lower, total//length + 1)
            for tail in increasing_partitions(total-a, length-1, a)]


def in_theorem_range(m, p, d):
    return d <= p or (p > 2 and d == p+1)


def run_case(m, p, top, boundary_expected=None):
    assert p in ZZ(m).prime_divisors() and top <= m+1
    assert in_theorem_range(m, p, top) or boundary_expected is not None
    # At top=m, eliminate c_m by F(L)=0, retaining the actual image.
    count = min(top, m-1)
    ring = PolynomialRing(ZZ, names=['L'] + ['c%d' % i for i in range(1, count + 1)])
    variables = ring.gens()
    L = variables[0]
    cs = [ring(1)] + list(variables[1:])
    if top >= m:
        cs.append((-1)**(m+1)*sum((-1)**i*cs[i]*L**(m-i) for i in range(m)))
    fs = [sum((-1)**i * comb(m-i, j-i) * cs[i] * L**(j-i)
              for i in range(j+1)) for j in range(1, min(top,m)+1)]
    q = [ring(m)]
    for j in range(1, top+1):
        q.append(sum((-1)**(i-1)*fs[i-1]*q[j-i] for i in range(1, min(j,len(fs)+1)))
                 + ((-1)**(j-1)*j*fs[j-1] if j<=len(fs) else 0))
    generators = cs[1:] + fs
    weights = list(range(1, len(cs))) + list(range(1,len(fs)+1))
    ambient_weights = [1] + list(range(1, count+1))
    bases = {0: [(0,) * len(variables)]}
    images = {0: matrix(ZZ, [[1]])}
    records = []
    counter = None
    correction = None
    for d in range(1, top+1):
        basis = exponents(ambient_weights, d)
        bases[d] = basis
        index = {e: i for i, e in enumerate(basis)}
        source_rows, ideal_rows = [], []
        for gi, (g, w) in enumerate(zip(generators, weights)):
            if w > d:
                continue
            op = matrix(ZZ, len(bases[d-w]), len(basis), sparse=True)
            for row, e in enumerate(bases[d-w]):
                for z, coefficient in g.dict().items():
                    op[row, index[tuple(x+y for x, y in zip(e, z))]] += coefficient
            block = images[d-w] * op
            source_rows.extend(block.rows())
            if gi < len(cs)-1:
                ideal_rows.extend(block.rows())
        image = matrix(ZZ, source_rows).hermite_form(include_zero_rows=False)
        images[d] = image
        ideal = matrix(ZZ, ideal_rows)
        coordinates = image.transpose().solve_right(ideal.transpose()).transpose()
        assert all(x.denominator() == 1 for x in coordinates.list())
        factors = [abs(int(x)) for x in matrix(ZZ, coordinates).elementary_divisors()
                   if abs(x) > 1]
        primary = [int(p**ZZ(x).valuation(p)) for x in factors if x % p == 0]
        predicted = ([int(p**ZZ(m).valuation(p))] * expected(d)
                     if in_theorem_range(m, p, d) else boundary_expected[d])
        assert primary == predicted, (m, p, d, primary, predicted)
        basis_check = None
        if in_theorem_range(m, p, d):
            representatives = []
            for k in range(2, d+1):
                remaining = d-(2**k-1)
                if remaining < 0:
                    break
                for part in increasing_partitions(remaining, k):
                    if d == p+1 and k == 2 and part == (0, p-2):
                        g = ring((q[p]-q[1]**p)/p)
                        representatives.append(q[1]*g-((m-m**p)//p)*q[d])
                    else:
                        representatives.append(ring.prod(q[2**i+part[i]] for i in range(k))
                                               - m**(k-1)*q[d])
            assert len(representatives) == expected(d)
            extra = []
            ordinary_lift = ring(0) if d>=m else (fs[d-1] if d==p else q[d])
            for polynomial in representatives + [ordinary_lift]:
                vector = [ZZ(0)]*len(basis)
                for e, coefficient in polynomial.dict().items():
                    vector[index[tuple(e)]] = coefficient
                extra.append(vector)
            span = matrix(ZZ, ideal.rows()+extra).hermite_form(include_zero_rows=False)
            assert span.nrows() == image.nrows()
            index_valuation = abs(span.det()).valuation(p)-abs(image.det()).valuation(p)
            assert index_valuation == 0
            basis_check = dict(integral_representatives=len(representatives),
                               source_generation_index_valuation=int(index_valuation))
        records.append(dict(degree=d, ambient_rank=len(basis),
                            full_integer_factors=factors, p_primary=primary,
                            partition_multiplicity=expected(d), in_theorem_range=in_theorem_range(m,p,d),
                            explicit_partition_basis=basis_check))
        print(json.dumps(dict(m=m, p=p, degree=d, primary=primary)), flush=True)
        if d == 7 and p > 7:
            # Newton moments from the ELEMENTARY source, independently
            # of the moment-polynomial construction in the proof.
            target = q[1] * (q[2]*q[4] - m*q[6])
            row = [ZZ(0)] * len(basis)
            for e, coefficient in target.dict().items():
                row[index[tuple(e)]] = coefficient
            h = ideal.hermite_form(include_zero_rows=False)
            solution = h.transpose().solve_right(matrix(ZZ, [row]).transpose())
            order = ZZ(1)
            for x in solution.list():
                order = order.lcm(x.denominator())
            assert order.valuation(p) == ZZ(m).valuation(p)
            counter = dict(degree=7, full_order=int(order),
                           p_primary_order=int(p**order.valuation(p)))
        if p>2 and d==p+2:
            target = q[1]*q[p+1]-m*q[p+2]
            row = [ZZ(0)]*len(basis)
            for e,coefficient in target.dict().items():
                row[index[tuple(e)]] = coefficient
            h = ideal.hermite_form(include_zero_rows=False)
            solution = h.transpose().solve_right(matrix(ZZ,[row]).transpose())
            order = ZZ(1)
            for x in solution.list():
                order = order.lcm(x.denominator())
            assert (m//p)%order==0
            correction = dict(degree=d,full_order=int(order),proved_annihilator=m//p,
                              full_smith_formula_in_this_degree_claimed=False)
    if boundary_expected is not None:
        assert records[-1]['p_primary'] != [int(p**ZZ(m).valuation(p))] * expected(top)
    return dict(m=m, p=p, v_p_m=int(ZZ(m).valuation(p)), degrees=records,
                root_boundary_included=top>=m,
                boundary_control=boundary_expected is not None,
                first_reduced_order_class=correction,
                nonannihilation_by_first_generator=counter)


def tensor_gluing_certificate():
    k=3;r=2
    R=PolynomialRing(ZZ,names=['u'+str(i) for i in range(k+1)]+['v'+str(i) for i in range(k+1)])
    u=R.gens()[:k+1];v=R.gens()[k+1:]
    perms=list(permutations(range(1,k+1)))
    subsets=[tuple(B) for n in range(k+1) for B in combinations(range(1,k+1),n)]
    mons=[(i,(j,l)) for i in range(1,k+1) for j in range(1,k+1) for l in range(j,k+1)]
    def poly(mon):
     i,jj=mon
     return u[i]*R.prod(v[j] for j in jj)
    def permmon(mon,pm):return (pm[mon[0]-1],tuple(sorted(pm[j-1] for j in mon[1])))
    def orbits(items,ps):
     seen=set();out=[]
     for mon in items:
      if mon in seen:continue
      orb=set(permmon(mon,pm) for pm in ps);seen.update(orb);out.append(orb)
     return out
    families=[]
    for size in range(k+1):
     B0=tuple(range(1,size+1))
     stab=[pm for pm in perms if set(pm[:size])==set(B0)]
     for orb in orbits(mons,stab):
      family={B:R(0) for B in subsets}
      for B in subsets:
       if len(B)!=size:continue
       pm=next(pm for pm in perms if set(pm[:size])==set(B))
       family[B]=sum(poly(permmon(mon,pm)) for mon in orb)
      families.append(family)
    constraints={}
    for B,C in combinations(subsets,2):
     coeff={i:int(i in B)-int(i in C) for i in range(1,k+1)}
     pivot=next(i for i in coeff if coeff[i])
     sub={u[pivot]:-sum(coeff[i]*u[i] for i in coeff if i!=pivot)/coeff[pivot],v[pivot]:-sum(coeff[i]*v[i] for i in coeff if i!=pivot)/coeff[pivot]}
     for col,family in enumerate(families):
      dif=R((family[B]-family[C]).subs(sub))
      for exp,c in dif.dict().items():constraints.setdefault((B,C,tuple(exp)),[0]*len(families))[col]+=c
    C=matrix(ZZ,list(constraints.values()))
    # Ambient invariant orbit sums, including fixed index 0.
    ambientmons=[(i,(j,l)) for i in range(k+1) for j in range(k+1) for l in range(j,k+1)]
    def apm(mon,pm):
     f=lambda x:0 if x==0 else pm[x-1]
     return(f(mon[0]),tuple(sorted(f(j) for j in mon[1])))
    seen=set();ambient=[]
    for mon in ambientmons:
     if mon in seen:continue
     orb=set(apm(mon,pm) for pm in perms);seen.update(orb);ambient.append(sum(poly(t) for t in orb))
    # Expand all family vectors in the common direct sum of local monomials.
    rowkeys=[(B,mon) for B in subsets for mon in mons]
    expof={mon:next(iter(poly(mon).dict())) for mon in mons}
    F=matrix(ZZ,[[fam[B].monomial_coefficient(poly(mon)) for fam in families] for B,mon in rowkeys])
    Araw=matrix(ZZ,[[R(a.subs({u[0]:sum(u[i] for i in B),v[0]:sum(v[i] for i in B)})).monomial_coefficient(poly(mon)) for a in ambient] for B,mon in rowkeys])
    A=matrix(ZZ,F.solve_right(Araw));assert C*A==0
    K=C.right_kernel_matrix();coords=K.transpose().solve_right(A)


    ideal_rows={}
    for B in subsets:
     if len(B)==k:continue
     for j,a in enumerate(ambient):
      f=R(a.subs({u[0]:sum(u[i] for i in B),v[0]:sum(v[i] for i in B)}))
      for mon in mons:
       if mon[0] in B and all(i in B for i in mon[1]):continue
       c=f.monomial_coefficient(poly(mon))
       if c:ideal_rows.setdefault((B,mon),[0]*len(ambient))[j]=c
    I=matrix(ZZ,list(ideal_rows.values()))
    assert C.rank()==17 and K.nrows()==11 and A.rank()==11
    assert [int(x) for x in C.elementary_divisors() if abs(x)>1]==[2]
    assert not [x for x in A.elementary_divisors() if abs(x)>1]
    assert I.rank()==10 and not [x for x in I.elementary_divisors() if abs(x)>1]
    assert list(I.right_kernel_matrix()[0])==[1]+[0]*10
    return dict(family_dimension=28,ambient_invariant_dimension=11,
                compatibility_shape=list(map(int,C.dimensions())),compatibility_rank=17,
                compatibility_nonunit_smith=[2],ambient_restriction_nonunit_smith=[],
                coefficient_condition_rank=10,coefficient_condition_nonunit_smith=[],
                coefficient_kernel='span(u0*v0^2)',
                method='Independent Sage polynomial substitution, invariant orbit sums, exact integer Smith forms.')


def cubic_prime_three_certificate():
    # M=1 only encodes homogenization to coefficient degree three:
    # a word of cost c<=3 represents M^(3-c) times that word.
    ring = PolynomialRing(ZZ,names=['L']+['a%d'%i for i in range(1,6)])
    L,*aa = ring.gens()
    q = [L**j+sum(comb(j,i)*aa[i-1]*L**(j-i) for i in range(1,j+1))
         for j in range(1,6)]
    basis = [e for e in exponents([1,1,2,3,4,5],5) if sum(e[1:])<=3]
    J,C = [],[]
    for e in exponents([1,2,3,4,5]*2,5):
        if sum(e)>3:
            continue
        polynomial = ring.prod(g**a for g,a in zip(aa+q,e))
        row = [polynomial.monomial_coefficient(ring.prod(x**a for x,a in zip(ring.gens(),b)))
               for b in basis]
        J.append(row)
        if any(e[:5]):
            C.append(row)
    J,C = matrix(ZZ,J),matrix(ZZ,C)
    assert (J.rank(),C.rank())==(16,15)
    assert [int(x) for x in J.elementary_divisors() if abs(x)>1]==[2]
    assert [int(x) for x in C.elementary_divisors() if abs(x)>1]==[2]
    return dict(prime=3,weight=5,coefficient_degree=3,source_rank=16,ideal_rank=15,
                source_nonunit_smith=[2],ideal_nonunit_smith=[2],three_saturated=True)


def main():
    # Prime-power, mixed, and higher-exponent controls, not a rank census.
    cases = [run_case(25, 5, 7,{7:[5,25,25,25]}), run_case(14, 7, 8),
             run_case(22, 11, 12), run_case(121, 11, 12), run_case(7, 7, 8), run_case(3,3,4),
             run_case(10, 5, 7, {5: [5, 5], 6: [5, 5], 7: [5, 5, 5]})]
    names = ['TYPE_A_K_STABLE_MULTIPLICITIES.md',
             'tools/verify_type_a_k_stable_multiplicities.py']
    report = dict(result='passed', cases=cases,
        tensor_gluing_certificate=tensor_gluing_certificate(),
        cubic_prime_three_certificate=cubic_prime_three_certificate(),
        method='Exhaustive homogeneous elementary-jet lattices over ZZ; full coefficient ideal in each checked degree; integral Smith quotients.',
        moment_or_apolar_builder_imported=False,
        all_rank_claim_from_finite_checks=False,
        uniform_cubic_annihilator_uses_exact_parameter_independent_tensor_certificate=True,
        full_degree_p_plus_two_multiplicity_formula=False,
        unrestricted_degrees_above_p_claimed=False, prime_boundary_claimed=True,
        degree_p_plus_one_claimed_for_odd_p_including_root_boundary=True,
        public_api_changed=False, external_review=False,
        source_sha256={name: sha256((ROOT/name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
                       for name in names})
    (ROOT/'type-a-k-stable-multiplicities-verification.json').write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(result='passed', cases=len(cases),
                         checked_degrees=sum(len(c['degrees']) for c in cases))))


if __name__ == '__main__':
    main()
