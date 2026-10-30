"""Independent Sage checks for critical-value geometry and adjoint witnesses.

Run with sage -python. This implementation does not import the native
SymPy audit, its derivatives, root enumerator, or compact E8 potential.
Finite calculations corroborate the written all-rank/geometric proofs.
"""
from pathlib import Path
from hashlib import sha256
from itertools import product
from collections import Counter
import json
from sage.all import QQ, PolynomialRing, LaurentPolynomialRing, NumberField, matrix, vector, WeylCharacterRing, RootSystem

ROOT = Path(__file__).resolve().parents[1]


def main():
    # Scheme controls use polynomial ideals with an explicit torus inverse.
    R = PolynomialRing(QQ, names=('x', 'y', 'z', 'zi', 'c'), order='lex')
    x, y, z, zi, c = R.gens()
    phi = x*x+y*y+(z*z+zi*zi)*x*y-1
    inverse = z*zi-1
    # Log derivative on the Laurent quotient: D z=z, D zi=-zi.
    D = lambda f: z*f.derivative(z)-zi*f.derivative(zi)
    original = R.ideal([inverse, phi.derivative(x), phi.derivative(y)])
    # On z=1, x=-y, the original critical ideal is contained in this
    # reduced component, but transverse multiplicity is two away from y=0.
    wall = R.ideal([inverse, z-1, x+y])
    assert all(f in wall for f in original.gens())
    assert D(phi) in wall
    # Change x=(u+v)/2,y=(u-v)/2 and saturate at v(z+1).
    S = PolynomialRing(QQ, names=('u', 'v', 'z', 'zi'), order='lex')
    u, v, z, zi = S.gens()
    changed = (u*u*(z+zi)**2-v*v*(z-zi)**2)/4-1
    J = S.ideal([z*zi-1, changed.derivative(u), changed.derivative(v)])
    localized = J.saturation(S.ideal(v*(z+1)))[0]
    expected = S.ideal([z*zi-1, u, (z-1)**2])
    assert localized == expected

    # Independent exact inversion of the nonlinear conormal constraint.
    R = PolynomialRing(QQ, names=('u', 'z', 'w', 'c1', 'c2'), order='lex')
    u, z, w, c1, c2 = R.gens()
    I = R.ideal([z+w-1, u*z-c1, u*w-c2])
    assert u-c1-c2 in I
    assert (c1+c2)*z-c1 in I and (c1+c2)*w-c2 in I
    # A nonclean component can even have a nonzero affine target image.
    phi_bad = u*u*(z+w-1)**2/2+z+w
    nonclean_support = R.ideal([z+w-1])
    assert phi_bad.derivative(u) in nonclean_support
    assert z*phi_bad.derivative(z)+w*phi_bad.derivative(w)-1 in nonclean_support

    # Original composed critical ideal, using a second CAS and genuine
    # ideal membership. The target is nilpotent but is not set to zero.
    C = PolynomialRing(QQ, names=('a','b','u','z','zi','c'), order='degrevlex')
    a, b, u, z, zi, c = C.gens()
    F = z+zi
    q = a*b+F*u+F*F
    potential = q*q
    central = z*potential.derivative(z)-zi*potential.derivative(zi)
    composed_ideal = C.ideal([potential.derivative(v) for v in (a,b,u)]+[central-c,z*zi-1])
    assert c*c in composed_ideal
    assert c not in composed_ideal
    # Actuality is checked through the full representation ring, not only
    # an abstract polynomial with assigned positive coefficients.
    W = WeylCharacterRing(['E',6], style='coroots')
    adj = W((0,1,0,0,0,0))
    tensor_character = (adj+2)**2
    assert int(tensor_character.degree()) == 6400
    assert all(coefficient >= 0 for coefficient in tensor_character.monomial_coefficients().values())

    # Rational Hessian criterion: compute support using reduced fractions
    # and ranks over the rational-function field, without SymPy.
    P = PolynomialRing(QQ, names=('z', 'w', 'v'))
    z, w, v = P.gens()
    F = P.fraction_field()
    functions = [F((z+w)/(z-w)), F((z*z-w*w)/(z-w)),
                 F((z+w+v)/(z*w*v+1)), F(1/(z+w)), F(z/w+v)]
    rational_ranks = []
    for f in functions:
        support = set(f.numerator().dict()) | set(f.denominator().dict())
        first = next(iter(support))
        M = matrix(QQ, [[a-b for a,b in zip(e, first)] for e in support])
        hessian = matrix(F, [[vv*(uu*f.derivative(uu)).derivative(vv)
                              for vv in P.gens()] for uu in P.gens()])
        assert M.rank() == hessian.rank()
        rational_ranks.append(int(M.rank()))

    # E8 full archived character, termwise exact derivatives over Q(i).
    # Ordinary derived variables and logarithmic central variables are
    # implemented by differential operators on a Laurent polynomial ring.
    K = NumberField(PolynomialRing(QQ, 'j').gen()**2+1, 'ii')
    ii = K.gen()
    row = next(r for r in json.loads((ROOT/'e8-final-character-models.json').read_text())
               if r['seed'] == [1,2,3,5,7])
    L = LaurentPolynomialRing(K, names=('x1','x2','x3','x5','x7','s','z0','z1'))
    variables = L.gens()
    term_variables = variables[:5]+(variables[5]+variables[5]**-1,)+variables[6:]
    full = L.zero()
    for term in row['terms']:
        monomial = L(term['coefficient'])
        for variable, power in zip(term_variables, term['exponent']+term['charge']):
            monomial *= variable**power
        full += monomial
    point = tuple(K(a) for a in (-3,0,-3,0,0,ii,1,1))
    def derivative(f, j):
        result = f.derivative(variables[j])
        return result if j<5 else variables[j]*result
    gradient = [derivative(full,j) for j in range(8)]
    assert all(g(*point) == 0 for g in gradient)
    hessian = matrix(K, [[derivative(g,j)(*point) for j in range(8)] for g in gradient])
    derived = hessian[:5,:5]
    schur = hessian[5:,5:]-hessian[5:,:5]*derived.inverse()*hessian[:5,5:]
    assert derived.det() == 216
    assert schur.det() == -69120
    assert hessian.det() == -14929920

    # Independent polynomial differentiation in SIMPLE-ROOT torus
    # coordinates. These differ from the native coroot-coordinate root
    # sums by an invertible Cartan congruence, so ranks must agree.
    involution_checks = []
    for n in (6,7,8):
        root_space = RootSystem(['E',n]).root_lattice()
        positive = list(root_space.positive_roots())
        roots = positive+[-a for a in positive]
        ring = LaurentPolynomialRing(QQ, names=tuple('a'+str(j) for j in range(n)))
        variables = ring.gens()
        character = ring(n)
        for root in roots:
            term = ring.one()
            for variable, exponent in zip(variables, root.to_vector()):
                term *= variable**int(exponent)
            character += term
        gradient = [v*character.derivative(v) for v in variables]
        hessians = [[v*g.derivative(v) for v in variables] for g in gradient]
        ranks = Counter()
        for point in product((QQ(1),QQ(-1)), repeat=n):
            assert all(g(*point) == 0 for g in gradient)
            H = matrix(QQ, [[f(*point) for f in row] for row in hessians])
            ranks[int(H.rank())] += 1
        assert ranks == ({6:28,1:36} if n==6 else {n:2**n})
        involution_checks.append(dict(type=f'E{n}', hessian_rank_counts=dict(ranks),
                                      characters_checked=2**n))

    # Check the actual adjoint character from an independent Weyl-character
    # implementation, then the all-type zero-fiber character identities.
    records = []
    for kind, n in [('A',n) for n in range(1,7)]+[('C',n) for n in range(2,7)]+[('B',3),('B',4),('D',4),('D',5),('G',2),('F',4),('E',6),('E',7),('E',8)]:
        ring = WeylCharacterRing([kind,n], style='coroots')
        space = RootSystem([kind,n]).root_lattice()
        highest = space.highest_root()
        labels = tuple(int(highest.scalar(a)) for a in space.simple_coroots())
        adjoint = ring(labels)
        fundamental = [ring(tuple(int(i==j) for i in range(n))) for j in range(n)]
        if kind == 'A':
            expression = fundamental[0]*fundamental[-1]-1
            label = 'x1*xr-1' if n>1 else 'x1^2-1'
        elif kind == 'C':
            expression = fundamental[0]**2-fundamental[1]-1
            label = 'x1^2-x2-1'
        else:
            assert sum(labels) == 1
            expression = fundamental[labels.index(1)]
            label = 'x'+str(labels.index(1)+1)
        assert expression == adjoint
        records.append(dict(type=f'{kind}{n}', adjoint_labels=list(labels), character=label,
                            dimension=int(adjoint.degree())))

    # Independently check the all-type-A symbolic elimination in a rational
    # function field. The original equations, not their radical, are used.
    R = PolynomialRing(QQ, names=('X','A','B','z','zi','zj','U','V','t','r'))
    F = R.fraction_field()
    X,A,B,z,zi,zj,U,V,t,r = F.gens()
    phi = (z*X+A)*(X/z+B)-1
    xc = -(z*B+A/z)/2
    value = -1-(z*B-A/z)**2/4
    assert phi.subs({X:xc}) == value
    assert phi.derivative(X).subs({X:xc}) == 0
    assert zi*(zi*V+U/zi)-zj*(zj*V+U/zj) == (zi**2-zj**2)*V
    # Under z=u^t, v_j=u^(-2q), put r=u^(2q+t). The character value
    # is this nonconstant rational function for every t>=1.
    assert value.subs({B:t*r/z, A:t*z/r}) == -1-t**2*(r-1/r)**2/4
    assert (-1-t**2*(r-1/r)**2/4).derivative(r) != 0

    # Reciprocal-polynomial Newton sums, independent of the native
    # companion-matrix computation. The quartic's constant term is a unit.
    R = PolynomialRing(QQ, names=('tau','L','q'))
    tau, L, q = R.gens()
    F = R.fraction_field()
    P = PolynomialRing(F, 'u')
    u = P.gen()
    p = u**4-tau*u**3+L*u**2+tau*u-1
    reciprocal = -P(list(reversed(p.list())))
    assert reciprocal == u**4-tau*u**3-L*u**2+tau*u-1
    def first_two_newton(poly):
        e1 = -poly[3]
        e2 = poly[2]
        return e1,e1**2-2*e2
    assert first_two_newton(p) == (tau,tau**2-2*L)
    assert first_two_newton(reciprocal) == (tau,tau**2+2*L)
    inv_u = u**3-tau*u**2+L*u+tau
    assert (u*inv_u-1).mod(p) == 0
    assert (tau*(u-inv_u)-(u**2-inv_u**2)-L).mod(p) == 0
    TB = F(1)/(1-2*q)
    assert (TB**2-(2*q*TB**2+1))/2 == -q/(2*q-1)
    assert (2*q-1)*(-q/(2*q-1))+q == 0

    report = dict(result='passed',
        SL4_actual_saturation='(u,(z-1)^2,z*zi-1)',
        nonclean_affine_target_control='c1+c2=1',
        nonlinear_conormal_inverse=True,
        E6_composed_target_squared_in_original_ideal=True,
        E6_composed_target_not_in_original_ideal=True,
        E6_shifted_adjoint_tensor_square_dimension=6400,
        type_A_composition_completed_square_and_unit_elimination=True,
        type_A_composition_nonconstant_restriction=True,
        BCD_four_block_reciprocal_Newton_check=True,
        BCD_four_block_constant_values=['-q/(2*q-1)','0','0'],
        rational_support_hessian_ranks=rational_ranks,
        E8_field='Q(i)', E8_original_derived_hessian_determinant=216,
        E8_original_central_schur_determinant=-69120,
        E8_original_full_hessian_determinant=-14929920,
        exceptional_involution_character_checks=involution_checks,
        adjoint_zero_character_checks=records,
        scope='Independent Sage exact controls, original-character E8 witness and Weyl-character zero-fiber identities; no mechanical verification of the geometric proof.',
        literature_priority_claim=False)
    names = ['tools/verify_critical_value_symmetry.py','ADJOINT_NONEMPTINESS_CRITICAL_VALUES.md',
             'CHARACTER_COMPOSITION_FIBERS.md','e8-final-character-models.json']
    report['source_sha256'] = {name:sha256((ROOT/name).read_bytes().replace(b'\r\n',b'\n')).hexdigest() for name in names}
    (ROOT/'critical-value-symmetry-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'}, indent=2))


if __name__ == '__main__':
    main()
