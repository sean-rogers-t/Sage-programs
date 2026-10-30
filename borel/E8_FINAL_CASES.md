# The last five E8 adjoint Levi families: exact schemes and boundary certificates

26 October 2026. This note records case proofs, computational certificates,
and the remaining obstruction separately. No claim of literature priority.
The group is the simple, simply connected complex group of type E8; the
representation is its full adjoint representation, with Killing-form
normalization

    B(theta(g),X)=tr(Ad(g) ad(X)).

For a proper Levi K put Z_K^o={h in z(k):Z_G(h)=K}. The question concerns
the actual scheme theta^(-1)(h) intersect G_reg, including its possible
nilpotents and its regular nonsemisimple points.

## 1. Scope and current theorem

**Theorem.** For each of the four types in the proved rows below there is
a dense open U_K in Z_K^o such that every actual regular-source fiber over
U_K is smooth, allowing emptiness. A nonempty such fiber has the indicated
pure dimension. Nonemptiness is not asserted by this note.

| Derived Levi type | Retained-node seed | Standard subsets | Fiber dimension if nonempty | Status |
| --- | --- | --- | --- | --- |
| A1^4 | 1257 | 7 | 8 | Homogeneous Jacobian certificate |
| A2^2+A1 | 12356 | 8 | 14 | Schur polynomial and every rank boundary |
| A2+A1^3 | 12357 | 8 | — | Remaining affine singular-image problem |
| A3+A1^2 | 12457 | 10 | 16 | Multiplier and homogeneous Jacobian certificate |
| A3+A2+A1 | 123567 | 4 | 20 | Binary criterion including its full boundary |

The already certified [Weyl transfers](E8_WEYL_TRANSFERS.md) apply to the
actual schemes and centralizer-exact targets. These four rows add 29
standard subsets to the preceding 218, giving **247/255**, with one
eight-member class remaining. This is not yet the all-E8 theorem.

Published inputs are the adjoint trace construction, fundamental-character
coordinates for the simply connected derived groups, characteristic-zero
generic smoothness, and the usual Jacobian criterion. The existing
[regular-quotient reduction](GENERAL_CAYLEY_FIBER_CRITICAL.md) supplies the
passage from the quotient critical scheme to the regular group fiber.
The explicit models and the boundary deductions below are repository work.
The dimension certificates are finite, reproducible exact computations;
they are not estimates from sampled target fibers.

## 2. Shared geometric and arithmetic arguments

On a finite central cover use derived fundamental characters x and central
torus coordinates z. If Phi is the full restricted character, use precisely

    Y=V(Phi_x),                  mu_j=z_j Phi_zj.

Unit-determinant paired elimination is an equality of quotient rings. It
retains all remaining affine and central derivatives by the chain rule.
The coefficient matrix may depend on the retained affine variables,
provided its determinant is a Laurent unit and it is independent of the
eliminated variables. Its adjugate then gives an inverse over the original
ring. No division by a nonunit is hidden in the eliminations in this note.

**Jacobian dimension lemma.** Let Y be defined by m polynomials in a smooth
ambient space of dimension m+s. Let S be the closed scheme defined by
those polynomials and all m-rowed Jacobian minors. If dim S<s, then every
morphism mu:Y->A^s has smooth, zero-dimensional generic fibers, possibly
empty. Indeed Y outside S is smooth of pure dimension s by the Jacobian
criterion. Remove the closure of mu(S), whose dimension is less than s,
and then apply characteristic-zero generic smoothness to each remaining
component. A nondominant component also has a proper closed target image
closure. There are finitely many components. This argument never replaces
Y by its reduction: the surviving open subscheme itself is smooth. It
does not assert that all relative critical schemes have dimension s; the
known E6 excess example is outside this lemma's sufficient hypotheses.

**Homogeneous specialization certificate.** Start with polynomials over Q,
clear rational denominators, and homogenize every specified generator
using a new variable h. Let I_h be the ideal they generate, without
silently saturating by h. At a denominator-safe prime p, for each degree
d the multiplication matrix defining (I_h)_d is a finite matrix over
Z localized away from p. Its rank over Q is at least its rank modulo p.
Consequently

    Hilb(Q[X,h]/I_h,d) <= Hilb(F_p[X,h]/(I_h mod p),d).

Thus the characteristic-zero homogeneous dimension is at most the
finite-field dimension. A set of homogeneous members of the finite-field
ideal is enough: the ideal of their leading monomials is contained in the
initial ideal of I_h mod p and supplies an upper bound on its dimension.
Membership is checked by replaying the original-ideal calculation, not by
assuming that a candidate list generates the original ideal. A truncated
standard-basis calculation suffices for such membership and leading-term
bounds; it is not asserted to be a complete Groebner basis.

The open set h!=0 is the product of the original affine scheme with G_m.
Hence a homogeneous bound d+1 proves an affine bound d. Extra components
at infinity can weaken this test, but cannot make a successful bound
incorrect. An affine finite-field dimension alone would not justify this
argument: components can disappear on reduction at a prime. This is why
homogenization is essential here.

### Central covers and ramification

The [full root-coordinate archive](e8-final-character-models.json) records
eight integral cocharacter columns for each model. Their nonzero
determinants give finite etale torus covers. The initial columns are the
retained simple coroots. Next come simple coroots of the subsystem
orthogonal to the retained roots. The remaining columns are orthogonal
to that entire semisimple subsystem. The commuting root types in the
five rows are A1^4, A2, A1, A1^2, and A1, respectively.

Replacing the commuting root torus by its fundamental characters is
finite. Its reflection walls are roots of E8 outside K. The full
character is invariant under each of these reflections while the derived
K characters stay fixed. Its normal derivative therefore vanishes on
the wall, so any target there satisfies B(h,alpha^vee)=0. Such a target
is not in Z_K^o. This excludes all reflection walls and their intersections
from the actual target fibers, not only from their reduced supports.
Off those walls the finite map is etale and the critical-scheme arguments
transfer. No discriminant of the *derived K torus* is inverted.

The additional one-dimensional trace substitution used in section 4 has
its own explicit target-hyperplane argument. It is not inferred to be an
E8 root wall.

## 3. A3+A2+A1: the binary case

The exact paired elimination for seed 123567 leaves variables x=x2 and
x3, a commuting A1 character a, and a torus variable z. Put y=z^2 x3 and
t=z^6. This is a unit change of the derived variable together with a
finite etale torus map. The full residual potential is

    f=alpha(x)y^2+beta(x)y+gamma(x),
    alpha=x^2-(t+t^-1)ax+a^2-1,
    beta=ax^3-2ta^2x^2+a^3x+(t-t^-1)x^2-ax
         +(t-t^-1)(a^2-1),
    gamma=-ta^3x^3+tax^3+a^2x^2+ta^3x-(2t+t^-1)ax-3.

The base C[a,t,t^-1] has dimension two. Define

    P=t^5(4alpha^2 gamma_x-2alpha beta beta_x+beta^2 alpha_x).

It has degree seven in x, leading coefficient -4a^2t^5, and content t^2,
a Laurent unit. Also gcd(alpha,beta)=1. At a=3,t=2 it is

    -1152(x^7-85x^6/8+133x^5/8+2355x^4/32
          -151x^3-645x^2/8+292x-130),

which is squarefree and retains degree seven. Primitivity and this
specialization prove squarefreeness over the entire Laurent base, not
only over its function field.

For the exceptional binary boundary, alpha_x=0 gives
x=(t+t^-1)a/2. The resultant in a of the numerators of alpha and beta at
this value, up to a nonzero rational factor, is

    t^4(t-1)^2(t+1)^2(t^2+1)^4(t^4+1)^2.

At t=+/-1 the numerator of alpha is a nonzero constant. At every other
possible nonzero t it is a nonzero quadratic in a. Therefore the common
zero set of alpha,beta,alpha_x is finite, and so is the smaller exceptional
set that also imposes beta_x=gamma_x=0. The full binary boundary has
dimension at most one. The [binary theorem](BINARY_BRANCHING_FIBERS.md)
now proves the asserted generic smoothness, with every boundary retained.

## 4. A3+A1^2: a multiplier without a denominator boundary

After eliminating the two endpoint characters of A3, let x,y be the A1
characters, v its middle character, a,b the two commuting A1 characters,
and put c=z^2+z^-2. The exact residual is

    f=-x^2y^2-cxy^2a-y^2a^2-cx^2yb-3xyab-cya^2b
      -x^2b^2-cxab^2-a^2b^2+v(ax+by+c)
      +x^2+y^2+cxa+a^2+cyb+b^2-4.

Its v equation is ax+by+c=0. To study the *total* relative critical
scheme, solve this equation for c globally and substitute it into f_x
and f_y. This is a polynomial ring isomorphism, even when a=b=0.
Denote the two resulting polynomials in x,y,v,a,b by g1,g2.

Form their 2-by-5 full Jacobian and its ten two-rowed minors. Homogenize
these twelve generators in six variables. The resulting homogeneous
ideal has dimension three modulo 32003. Independent SymPy computation
at 32009 gives the same dimension (297 basis polynomials). The preceding
homogeneous specialization lemma gives a characteristic-zero affine
singular-locus bound of two. Since the central target has dimension
three, the Jacobian dimension lemma applies.

The trace substitution c=z^2+z^-2 ramifies at z^4=1. On the actual
relative critical scheme, including the eliminated-variable graph,

    z Phi_z=(2z^2-2z^-2) f_c.

Thus this entire ramification locus maps into the proper target
hyperplane with this logarithmic coordinate zero. Excluding that target
hyperplane is legitimate for a generic theorem. The two commuting A1
trace walls are handled by the root-wall argument of section 2.

## 5. A2^2+A1: all Schur boundaries

Write the remaining derived characters as x,y,z, and the normalized
central coordinates as a,b,t. In terms of the root archive, these are

    x=x2, x5=z0*y, x6=z0^-1*z,
    a0=z0*a, a1=z0^-1*b, t=z0^3.

The central variable normalization is an isomorphism on affine character
coordinates over the torus, followed by a finite etale isogeny. It must
not be confused with a statement that the target coordinates are simply
rescaled: the actual target morphism is retained throughout.

The full residual has the form

    f=Q+(y,z)H(y,z)^t/2+ell^t(y,z),
    H11=-2b(x+t), H22=-2a(x+t^-1),
    H12=H21=-x^2-ab-(t+t^-1)x,
    ell1=-ax^2-b^2x-atx-at^2-b^2/t+a,
    ell2=-a^2x-bx^2-a^2t-bx/t+b-b/t^2,
    Q=(1-ab)(x^2+(t+t^-1)x)-3.

Set D=det H and N=QD-ell^t adj(H)ell/2. On D!=0 exact elimination
gives f=N/D and the single equation

    J=t^6(N_x D-N D_x)=0.

J has degree nine in x with Laurent-unit leading coefficient 2t^6.
At a=2,b=3,t=2 its degree remains nine and gcd(J,J_x)=1. Hence J is
primitive and squarefree over the full base. The surviving hypersurface
is reduced of pure dimension three; its nonsmooth locus has dimension
at most two.

It remains to cover all of D=0. A nonzero symmetric rank-one 2-by-2
matrix over C has a nonzero diagonal entry. On the chart Hjj!=0,
consistency of H(y,z)^t+ell=0 is

    D=0, C_j=Hjj ell_(1-j)-Hj,(1-j) ell_j=0.

The base coordinates are x,a,b,t with t!=0. Exact rational Groebner
calculations, independently recomputed in SymPy, give

| Chart | Dimension of V(D,C_j), Hjj*t invertible | Dimension after also imposing D_x=0 |
| --- | --- | --- |
| H11!=0 | 2 | 1 |
| H22!=0 | 2 | 1 |

The certificates include an inverse variable u with u*t*num(Hjj)=1;
they do not count the discarded denominator boundary as covered by
localization. Each chart leaves one free kernel coordinate. In the
remaining equation f_x=0 its quadratic coefficient on D=0 is a unit
multiple of D_x. Off D_x=0 this cuts the kernel coordinate finitely.
Over D_x=0 the base has dimension at most one, so even an identically
vanishing final equation gives total dimension at most two.

Finally H=0 and ell=0 have precisely the following support, with all
intersections retained:

* a=b=0 and x=0 or x=-t-t^-1;
* ab=1, t=+/-1, x=-t.

This follows from an exact equality of the consistency ideal with the
ten-generator ideal recorded in the verifier, not from a generic-rank
argument. In the first branch

    f_x=(2x+t+t^-1)(1-yz).

Except at the two values t=+/-i it imposes a genuine equation in y,z;
at those two values the base is zero-dimensional and both free variables
still give total dimension two. In the second branch the coefficients
of y^2,z^2 in f_x are -b,-a, both nonzero. Again total dimension is two.
Thus the entire Schur boundary has dimension at most two. Its target
image cannot dominate the three-dimensional center. This completes all
rank strata, including rank zero and their intersections.

## 6. A1^4: a global Jacobian certificate

The orthogonal A1^8 root cover has determinant 16. Let the four derived
characters be x1,x2,x5,x7 and the four central A1 characters a0,a1,a2,a3.
The full character is the sum of their eight squares, minus eight, plus
the following fourteen tetrads:

    x1*x2*x5*a0 + x1*x2*x7*a1 + x5*x7*a0*a1
    + x1*x5*x7*a2 + x2*x7*a0*a2 + x2*x5*a1*a2
    + x1*a0*a1*a2 + x2*x5*x7*a3 + x1*x7*a0*a3
    + x1*x5*a1*a3 + x2*a0*a1*a3 + x1*x2*a2*a3
    + x5*a0*a2*a3 + x7*a1*a2*a3.

Expanding the eight SL2 characters reconstructs all 240 root weights
and eight zero weights, independently in Sage and native Cartan code.
Take the four derived derivatives and the seventy maximal minors of
their full 4-by-8 Jacobian. Homogenization gives 74 generators in nine
variables. At primes 32003 and 32009 a degree-12 truncated standard-basis
calculation gives 41 selected original-ideal members whose leading
monomial ideal has dimension four. This is checked both by CAS dimension
and by enumerating all 512 subsets of the nine variables. The two
compressed certificates contain the selected polynomials, not merely a
claimed dimension. Replaying the original-ideal computation checks
their membership.

The homogeneous specialization lemma therefore bounds the actual
characteristic-zero affine Jacobian-deficiency locus by three. The
central base has dimension four. Apply the Jacobian dimension lemma.
All four central trace walls a_i=+/-2 map into their respective normal
target hyperplanes. Their intersections are included in this statement.
This completes the class without any Schur denominator or derived
root-discriminant deletion.

## 7. The remaining A2+A1^3 model and exact gap

The full character and its unit A2 elimination have been reconstructed
and verified; that is not yet a smoothness proof. Let x,y,w be the three
remaining derived A1 characters, a the commuting A1 character, and

    U=z0^3/z1, V=z0^6/z1,
    A=U+U^-1, B=V+V^-1, C=UV+(UV)^-1.

Then h=A^2+B^2+C^2-ABC-4=0. The residual potential is

    f=x^2+y^2+w^2+a^2
      -x^2y^2-x^2w^2-y^2w^2-x^2a^2-y^2a^2-w^2a^2
      -C(x^2yw+xy^2a+xw^2a+ywa^2)
      -B(xy^2w+x^2ya+yw^2a+xwa^2)
      -A(xyw^2+x^2wa+y^2wa+xya^2)-5xywa
      -BC(xy+wa)-AC(xw+ya)-AB(yw+xa)-ABC-4.

The exact critical ideal is (f_x,f_y,f_w,h). The outstanding sufficient
certificate is a dimension bound of two for the full Jacobian-deficiency
scheme, or a direct proof that its actual central target image is
nondense. The finite torus-to-(A,B,C) map also requires its four ramification
points to be handled in the target. Naive homogenization has extra
components at infinity and has not supplied the required bound. Neither
finite-field affine evidence nor the already solved a=0 specialization
is treated as a proof for arbitrary a.

The [reproducible affine experiment](tools/explore_e8_last_affine.py)
computes dimension one over F_32003 (2,384 basis polynomials). This is
supporting evidence only. The tested powers one through four of
(a^2-4)(2C-AB), (a^2-4)(A^2-4)(B^2-4)(C^2-4), and
(2C-AB)(2B-AC)(2A-BC) do not belong to that finite-field ideal.
Thus these tests do not furnish the desired simple ramification
annihilator; they do not rule out a different power or a different
annihilator. The complete basis and rational partial calculations are
preserved in a local research cache. The missing certificate is a
characteristic-zero original-ideal dimension bound or a justified
annihilator of the actual singular target image, including embedded
structure. The experiment is summarized in
[e8-last-case-exploration.json](e8-last-case-exploration.json).

## 8. Audit scope and controls

[Sage verification](tools/verify_e8_final_cases.py) checks the full root
characters, exact unit eliminations, original derivative chain rules,
binary identities and every Schur boundary. The [native audit](tools/audit_e8_final_cases.py)
uses independently generated Cartan roots and SymPy. The
[A1^4 certificate replay](tools/verify_e8_a1_four_certificate.py) validates
the two homogeneous bounds separately. None is a novelty certificate.

The Sp4 control remains consistent because all conclusions are generic,
not statements about every centralizer-exact target. The E6 excess
control remains admissible: zero-dimensional transverse fibers are a
sufficient mechanism here, not a new uniform hypothesis. The altered
character control is also retained: every calculation starts from the
full adjoint root character, including linear and constant coefficients,
and uses the original scheme equations rather than a radical ideal.
