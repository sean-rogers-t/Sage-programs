# Exceptional adjoint fibers: structural exclusions and an E6 excess theorem

Date: 2026-10-23. Written deductions from standard root, character and
adjoint-quotient theory, with independent exact checks. No priority claim.
This note advances the exceptional extension by locating a genuine
obstruction to extending the classical dimension formula unchanged.

## 1. Conclusions and conventions

G is simple complex simply connected, rho=Ad, and theta is normalized by

    B_Ad(theta(g),X)=tr(Ad(g) ad(X)),
    B_Ad(X,Y)=tr(ad(X) ad(Y)).                                  (1)

We study F_h^reg=theta^(-1)(h) intersect G_reg for a semisimple h whose
centralizer is a specified proper Levi K. This is the regular-source
part of the fiber as a scheme, not the full fiber or its reduction.

The new conclusions are:

1. A root-length obstruction excludes regular fiber points for every
   target in many Levi strata, uniformly across root systems. In a
   simply laced G, any D_m (m>=4) or E-type simple factor in [K,K]
   is enough for the exclusion, even with other Levi factors present.
2. A linear-coupling transfer formula identifies entire quotient-critical
   schemes, including their unused coordinates and nonreduced structure.
3. For the E6 contact Levi with derived group SL6, every nonzero central
   target has a nonempty smooth regular fiber of pure dimension 32.
   Here dim K-rank G=30, so the quotient contributes two extra dimensions.
   Its actual quotient-critical scheme is A^2, after the central cover
   is correctly divided out.

Thus the [classical adjoint dimension theorem](CLASSICAL_ADJOINT_FIBERS.md)
does not extend to all simple groups with the same dimension conclusion.
This is an adjoint-representation obstruction, distinct from the earlier
[spin-representation excess](SPIN_FIBER_EXCESS.md). It does not refute
generic smoothness.

## 2. A root-length obstruction for Levi factors

Let S be a simple factor of [K,K]. Give root spaces their inherited
W-invariant Euclidean metric. Suppose:

- the highest root theta_S is a fundamental weight of S, so the adjoint
  character of S is one fundamental-character coordinate y;
- ||theta_S||^2 is the largest root length squared in G.

**Theorem 1.** Under these hypotheses, F_h^reg is empty for every h
with Z_G(h)=K.

Proof. For a root alpha of G not belonging to S, its projection to
the span of S has norm strictly smaller than ||theta_S||. Indeed,
projection decreases norm strictly unless alpha lies in that span, and
a standard Levi root subsystem contains every ambient root in its span.
All roots of G have norm at most ||theta_S||. Thus every highest weight
nu of an external S-constituent has ||nu||<||theta_S||.

Its character polynomial cannot involve y. To see this without assuming
that the constituent is fundamental, use triangular inversion between
irreducible characters and products of fundamental characters. Every
monomial arising has dominant highest weight mu<=nu. For dominant
weights the norm is nondecreasing in dominance order:

    ||nu||^2-||mu||^2=(nu-mu,nu+mu)>=0.

Any monomial involving y has highest weight theta_S+eta, with eta
dominant, and hence norm at least ||theta_S||. This is impossible.
The argument works coefficientwise with the other Levi characters and
central Laurent characters retained. No specialization or cancellation
of central coefficients is used.

The restriction of Ad_G contains Ad_S once, with trivial central
character. Other derived factors and the central trivial summands do
not involve y. Consequently

    chi_Ad|K = y + R(other fundamental characters,central variables).

Its y derivative is 1. The quotient critical ideal is therefore the
unit ideal. By the [regular critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md),
there are no regular source points over this target stratum.

**Consequences.** In simply laced groups every root has the same norm;
the adjoints of D_m (m>=4) and E6,E7,E8 are fundamental. Thus the first
conclusion of section 1 holds for any such Levi factor. The same theorem
applies to ambient-long B_m factors with m>=3, including B3 in F4.
It is a sufficient obstruction, not a classification of the remaining
A-type factors. It asserts regular-part emptiness, not whole-fiber emptiness.

## 3. Exact transfer through a linearly coupled character

Let v=(v_1,...,v_d), y be affine coordinates and z in C*. Suppose

    chi(v,y,z)=P(v)+a(z)y+q(z).                                (2)

Write D=z partial_z. Assume the zero scheme of a in C* is finite and
reduced. Empty zero scheme is allowed. For c in C the critical equations
of W_c=chi-c log z are

    P_i(v)=0,  a(z)=0,  (Da)(z)y+(Dq)(z)=c.                  (3)

**Theorem 2 (linear-coupling transfer).** For every c, the actual
critical scheme is the disjoint union, over a(zeta)=0, of

    Crit(P) with z=zeta and y=(c-Dq(zeta))/Da(zeta).           (4)

This is a scheme isomorphism, with no radical taken. The Hessian adds
a nonsingular (y,z) block to the Hessian of P at each such point.
In particular it preserves excess dimensions and any transverse
nonreduced structure; it does not automatically make Crit(P) smooth.

Proof. Modulo a(z), Da is invertible because all zeros are simple.
The final equation in (3) eliminates y, and the Chinese remainder
theorem separates the distinct zeta. The mixed Hessian entry in
coordinates (y,log z) is Da(zeta), so its two by two determinant is
-(Da(zeta))^2. The v equations and Hessian block are independent.

In the contact application, a=z+z^(-1), q=z^2+z^(-2) plus an optional
constant. Its zeros are zeta=+/-i, Da=2zeta, Dq=0 there. Thus y=c/(2zeta).
Finite central-cover descent must still be performed before counting
components or identifying the quotient on the actual K.

## 4. The exceptional contact grading as one structural family

Let theta_G be the highest long root, H=theta_G^vee, and
K=Z_G(H). The contact grading is

    g=g_(-2)+g_(-1)+g_0+g_1+g_2,
    g_0=s+C H,  dim g_(+/-2)=1,  g_(+/-1)=V.                 (5)

The last equality means isomorphic S-modules; V is self-dual. Use
the central cocharacter z -> theta_G^vee(z) on the cover S times C*.
The adjoint character is

    chi_Ad=chi_(Ad S)+1+z^2+z^(-2)+(z+z^(-1))chi_V.           (6)

These gradings and modules are published representation theory; see
Cahen--Schwachhofer, section 2 and Table 1 of
[Special Symplectic Connections](https://arxiv.org/pdf/math/0402221).
The following table states that input, not a newly discovered branching.

| G | Omitted Bourbaki node | S | V | dim V | B_Ad(H,H) |
| --- | --- | --- | --- | --- | --- |
| G2 | 2 | SL2 on the short root | Sym^3 C^2 | 4 | 16 |
| F4 | 1 | Sp6 | primitive exterior cube, V(omega_3) | 14 | 36 |
| E6 | 2 | SL6 | exterior cube, V(omega_3) | 20 | 48 |
| E7 | 1 | Spin12 | one half-spin module | 32 | 72 |
| E8 | 8 | E7 | V(omega_7) | 56 | 120 |

The normalization in the last column follows directly from the grading:
tr(ad(H)^2)=2 dim V+8. For h=tH set c=B_Ad(h,H); in E6 this is c=48t.
Every t!=0 has Z_G(h)=K. The symbols H and c are not normalized by an
unspecified proportionality constant.

For F4,E6,E7,E8 the module V is one fundamental coordinate y, absent
from the adjoint character polynomial of S. Theorem 2 applies:

- F4: chi_(Ad Sp6)=x_1^2-x_2-1. After removing y=x_3, the residual
  polynomial has derivative -1 in x_2, so the regular fiber is empty.
- E7: chi_(Ad Spin12) is its adjoint fundamental coordinate, distinct
  from either half-spin coordinate. Its derivative is 1, so the regular
  fiber is empty. Theorem 1 independently gives this exclusion.
- E8: chi_(Ad E7) is x_1, distinct from y=x_7. Again the regular fiber
  is empty, also independently covered by Theorem 1.
- E6: chi_(Ad SL6)=x_1 x_5-1 and y=x_3. The residual critical scheme
  is x_1=x_5=0 with x_2,x_4 free. It is A^2, not a finite set.

These assertions hold for every nonzero central target in this contact
Levi stratum. They do not classify other exceptional Levi strata.

## 5. E6: exact quotient scheme, central descent and dimension

Let x_j=chi_(exterior^j C^6), j=1,...,5. Formula (6) reads

    chi_Ad=x_1 x_5+z^2+z^(-2)+(z+z^(-1))x_3.                (7)

The critical ideal, in the Laurent ring with parameter c, is exactly

    (x_1,x_5,z^2+1,2z x_3-c).                               (8)

Thus on the central cover it is two copies of A^2_(x_2,x_4). Its
Hessian has rank four in the six quotient coordinates: the x_1,x_5
pair and x_3,log z pair are nonsingular; x_2,x_4 are free.

The central cover is

    SL6 times C* -> K,   kernel=<(-I_6,-1)> of order two.     (9)

This can be checked in the simply connected cocharacter lattice.
In the E6 simple-coroot basis,

    H=(1,2,2,3,2,1).

Replacing the omitted second coroot by H has determinant two. At z=-1,
the remaining exponents on the A5 chain (1,3,4,5,6) are (1,0,1,0,1)
modulo two, which is -I_6 in SL6. This proves (9), including the kernel
element rather than just its order.

The deck transformation sends z to -z and x_j to (-1)^j x_j. It
exchanges the two copies in (8), fixes x_2,x_4, and preserves x_3=c/(2z).
Consequently the critical scheme on the actual K//K is one A^2.
There are not two distinct components arising from z=+i and z=-i.

For h=tH, t!=0, every point of this quotient scheme has a K-regular
lift. The critical equations imply theta(g)=h at such a lift.
Equivariance forces Z_G(g) subset Z_G(h)=K, so K-regularity becomes
G-regularity. Conversely every regular point of this Springer fiber
arises in this way. The quotient map is smooth there, of relative
dimension dim K-rank G=36-6=30. Therefore F_h^reg is nonempty smooth
of pure dimension 30+2=32. In particular d theta has rank 78-32=46
at every point of F_h^reg.

At h=0, K is not the centralizer of the target, so this argument does
not identify the G-regular part. E6 adjoint is fundamental, and its zero
fiber has no G-regular points by the existing zero-fiber theorem.
The K-relative critical scheme at c=0 is not a counterexample to that
theorem: its K-regular lifts are not G-regular.

## 6. The G2 member and the limit of the transfer hypothesis

Here the relative module has character x^3-2x, not an independent
fundamental coordinate. Formula (6) is

    chi=x^2+z^2+z^(-2)+(z+z^(-1))(x^3-2x).                  (10)

Set u=z+z^(-1), v=z-z^(-1). The relative critical curve X and target
map mu are

    2x+u(3x^2-2)=0,
    mu=v(2u+x^3-2x).                                       (11)

The denominator 3x^2-2 never vanishes on X, so u=-2x/(3x^2-2).
The curve is smooth: away from z=+/-1 its z derivative is nonzero;
at z=+/-1 the two x-roots are distinct (the corresponding quadratic
discriminant is 100), so its x derivative is nonzero. Also

    mu^2=-4 x^6(3x^2-8)^2(9x^2-4)(x^2-1)/(3x^2-2)^4.       (12)

This is nonconstant. The double cover in x has four distinct branch
points x=+/-1,+/-2/3, so X is irreducible. It follows that mu is
dominant and its generic fibers are nonempty finite and smooth.
Central-cover descent and the regular quotient presentation give
nonempty smooth generic regular Springer fibers of dimension
dim K-rank G=4-2=2. No generic point count is asserted here.

This completes the generic regular-fiber outcomes for the five
exceptional contact Levi strata, using the single grading formula
and precisely identifying the member to which linear transfer does
not apply. It is not an all-Levi exceptional classification.

## 7. Independent checks, sources and next gate

The [Sage verifier](tools/verify_exceptional_adjoint_contact.py) compares
each grade with native Weyl characters, checks the E6 ideal and deck
action, and proves the G2 generic curve assertions by exact ideals.
An independent root-basis computation evaluates all 72 E6 roots at
an explicit regular torus point: it checks the normalized target,
Cartan differential rank four, 42 nonzero root-direction derivatives,
and hence full differential rank 46 and kernel dimension 32.
See the [verification report](exceptional-adjoint-contact-verification.json).

The separate [Sage-free auditor](tools/audit_exceptional_adjoint_contact.py)
reconstructs ambient roots from explicit Cartan matrices and reflections,
and reconstructs V from exterior, primitive exterior, half-spin or
minuscule weight models. It checks all five graded decompositions,
trace normalizations and cover degrees, rather than trusting a branching
table or interpolating a rank pattern. Its
[audit report](exceptional-adjoint-contact-audit.json) states the shared
conventions and source bindings.

Known inputs: the published contact gradings and modules above;
fundamental-character polynomial rings, triangular character inversion
and regular adjoint quotients; and the
[Kostant--Michor trace-projection construction](https://www.mat.univie.ac.at/~michor/kostant.pdf).
The root-length exclusion and exact critical-scheme transfer are written
deductions here, as are the E6 fiber consequence and G2 curve argument.
They use elementary elimination once the standard input is supplied.
No claim is made that these deductions are absent from the literature.
Targeted contact-grading and generalized-Cayley E6-fiber searches on
2026-10-23 located the published grading source, but did not establish
priority for these fiber statements. A complete novelty audit remains separate.

The next structural gate is to classify the residual critical schemes
for the remaining A-type Levi factors in E6/E7/E8. Theorem 1 already
disposes of factors with a maximal-length fundamental adjoint; Theorem 2
shows why unused character coordinates must be retained. A proposed
uniform theorem must allow excess dimension, and must distinguish that
from nonreduced structure. The next proof should identify when the
restricted adjoint character has a smooth residual critical scheme,
with its dimension computed from that scheme, rather than assume
isolated critical points or start another unrestricted table.
