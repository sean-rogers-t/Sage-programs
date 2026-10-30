# A single irreducible family realizes all type A translated nilpotent maps

Date: 2026-10-23. Written characteristic-zero theorem with independent finite
matrix corroboration. Literature priority is not established.

## 1. Precise theorem

Fix m>=2, G=SL_(2m)(C), and the irreducible fundamental representation

    rho = exterior^m C^(2m),     highest weight omega_m,
    D = binomial(2m-2,m-1),      B_rho(X,Y)=D tr(XY).

The normalized Springer map theta is defined by

    B_rho(theta(g),Y)=tr(rho(g) d rho(Y)).

It has identity differential at the identity. We use the simply connected
group and its full weight lattice. The representation is infinitesimally
faithful but has finite central kernel mu_m. If a faithful group
representation is required, it descends to SL_(2m)/mu_m; the Lie algebra
and the formulas below are unchanged. No integral lattice comparison,
parabolic cohomology, or Borel characteristic homomorphism is asserted.

**Theorem.** Given any polynomial

    P(x)=p_1 x+...+p_(m-1) x^(m-1), including P=0,

there are a noncentral semisimple t in G and c in C*, such that

    H=Z_G(t)^0 is isomorphic to S(GL_m x (C*)^m),
    [H,H]=SL_m,

and, under the block inclusion iota: sl_m -> sl_(2m),

    theta(t diag(I_m+N,I_m))-theta(t) = c iota(P(N))        (1)

for every nilpotent m-by-m matrix N. Thus the entire translated
centralizer map from the unipotent variety of H to its nilpotent cone
is exactly cP in these coordinates. The equality also restricts to
the upper triangular unipotent subgroup of H.

The representation and group are fixed for each m; only the semisimple
point changes with P. The scalar c may be chosen outside an explicitly
defined finite set. It cannot simply be omitted from this theorem.

Every SL_m-equivariant algebraic map from its nilpotent cone to itself
is a polynomial P of this form. Consequently the theorem realizes all
such maps, up to a nonzero target scalar, inside one irreducible
Springer-map family. In particular all orders 1,...,m-1 and the zero
map occur. There is no rank-independent bound on nonzero orders in
irreducible ambient representations at noncentral points.

This is stronger than the [previous exterior-mixture result](TYPE_A_CHARACTER_DEGENERATION.md):
there we changed a reducible representation and used a central point.
Here we keep one irreducible representation and change a noncentral point.

## 2. General exterior-power block formula

For clarity, first let G=SL_(m+q), rho=exterior^p C^(m+q),
1<=p<=m+q-1, and

    t=diag(alpha I_m,beta_1,...,beta_q),
    alpha^m product(beta_a)=1,
    Q_beta(z)=product_a(1+beta_a z),
    D_p=binomial(m+q-2,p-1).

Assume the beta_a differ from alpha; they need not be pairwise distinct
for this formula. For nilpotent N the translated map is supported entirely
in the first m-by-m block and equals sum_(j=1)^(m-1) b_j N^j, where

    b_j = (-1)^(j-1) alpha^j / D_p *
          [z^(p-j)] (1+alpha z)^(m-j-1) Q_beta(z).          (2)

A negative coefficient index means zero. The normalization D_p is the
quadratic trace index in the ambient group, not the index of a single
summand of its restriction to the centralizer.

To prove (2), use the established identity

    sum_p z^p chi_(exterior^p)(g)=det(I+zg).

Its trace gradient is zg adj(I+zg), followed by traceless projection.
On the first block its coefficient of N^j, j>=1, is

    (-1)^(j-1) (alpha z)^j (1+alpha z)^(m-j-1) Q_beta(z).

This follows by expanding the inverse in
det(I+zg) zg(I+zg)^(-1), an identity in rational functions which is
polynomial after using the adjugate. The other blocks depend on N only
through det(I+alpha z(I+N))=(1+alpha z)^m, so their translated differences
are zero. The first block's scalar term is also unchanged, and the
positive powers of N have trace zero. Subtract theta(t), extract z^p,
and divide by D_p. This proves an equality in the full ambient Lie algebra,
not merely its projection to sl_m.

## 3. Invertible coefficient transform

Now set q=p=m and choose alpha=i, where i^2=-1. Write

    beta_a=i gamma_a,
    Q(z)=product_a(1+gamma_a z)=sum_(h=0)^m e_h z^h,
    e_0=1, e_m=(-1)^m.

The last condition ensures det(t)=i^(2m)e_m=1. Substituting into (2),
with h=m-j, gives

    b_(m-h) = (-1)^(m-h-1) i^m / D * S_h,
    S_h = sum_(k=1)^h binomial(h-1,k-1) e_k,
    1<=h<=m-1.                                            (3)

There is no e_0 term, because the coefficient sought has degree h while
(1+z)^(h-1) has smaller degree. Neither does e_m occur. The Pascal
matrix in (3) is triangular with diagonal one and inverse

    e_h = sum_(k=1)^h (-1)^(h-k) binomial(h-1,k-1) S_k.     (4)

Thus e_1,...,e_(m-1) can prescribe all b_1,...,b_(m-1) freely, before
checking the semisimple spectrum. For a prescribed P, solve (3)-(4)
with b_j=p_j to obtain coefficients E_h. For b_j=c p_j use e_h=c E_h.
Equivalently,

    E_h = D i^(-m) (-1)^(m+h-1) *
          sum_(k=1)^h binomial(h-1,k-1) p_(m-k).            (5)

This formula is valid for P=0 as well.

## 4. The spectrum exists, with exactly the required centralizer

Consider

    Q_c(z)=1+(-z)^m+c sum_(h=1)^(m-1) E_h z^h.             (6)

Its constant coefficient is one and its leading coefficient is the
nonzero fixed number (-1)^m. Therefore it factors over C as
product(1+gamma_a z), with all gamma_a nonzero. At c=0, Q_0=1+(-z)^m
is squarefree, and Q_0(-1)=2. Hence

    disc_z(Q_c) Q_c(-1)

is a nonzero polynomial in c. Choose c nonzero outside its finite zero
set. Then the gamma_a are distinct and none equals one. Consequently
the diagonal matrix

    t=diag(i,...,i, i gamma_1,...,i gamma_m)

has exactly one repeated eigenvalue, i of multiplicity m, and m other
distinct eigenvalues. It belongs to SL_(2m), is semisimple and noncentral,
and has precisely the stated centralizer. Combining (3)-(6) proves (1).

This argument is existential over C, but effective: equations (5)-(6)
construct the spectrum polynomial, and only a discriminant and one
evaluation must be avoided. For algebraic P, c and the eigenvalues may
be chosen algebraic. No positivity of artificial multiplicities is needed:
the ambient representation is the actual irreducible exterior power;
the variable complex numbers are eigenvalues of t.

For P=0, take any c!=0; Q_c=1+(-z)^m is already squarefree and avoids
gamma=1. The translated map vanishes on the whole unipotent variety of H.
In particular theta^(-1)(theta(t)) contains the closed subvariety
t U(H), of dimension m^2-m. Here U(H) means the full unipotent variety,
not the unipotent radical of a chosen Borel. No assertion that it is a
whole fiber component, or the entire fiber, is made.

## 5. Why this covers all equivariant nilpotent maps

Let F: N(sl_m) -> N(sl_m) be SL_m-equivariant and let J be a regular
nilpotent Jordan block. The subgroup exp(sJ) fixes J, so F(J) commutes
with J. The commutant of J consists of polynomials in J. Tracelessness
removes the constant term, giving F(J)=sum_(j=1)^(m-1) p_j J^j.
Equivariance gives the same identity on the regular nilpotent orbit,
which is dense in the nilpotent cone. Regularity extends it everywhere.
This familiar dense-orbit argument is included to specify the scope;
we do not claim the abstract polynomial classification as new.

If P is nonzero with first coefficient p_r!=0, write P(x)=p_r x^r v(x),
v(0)=1. The polynomial coordinate change Y=N v(N)^(1/r), truncated
modulo N^m, has a polynomial inverse. Thus P(N)=p_r Y^r. This recovers
the order classification and identifies its fiber schemes with matrix-power
fiber schemes, after scaling the target. Neither this change nor theorem (1)
classifies the fibers of theta on all of G.

## 6. Published inputs, deduction, checks and limitations

[Kostant--Michor, 4.12](https://arxiv.org/pdf/math/0109066) already define
the equivariant translated centralizer map and identify it as
theta(tu)-theta(t). Their 4.13 gives an isomorphism under an invertible
full-differential hypothesis. We study degenerations outside that
hypothesis; this is not a correction to their theorem and is not described
as an explicitly posed open problem in their paper.

The determinant character identity is standard; see
[van Ritbergen--Schellekens--Vermaseren, (49), Section 5.3](https://arxiv.org/pdf/hep-ph/9802376).
Both primary sources were inspected directly, including KM's assumption
of infinitesimal faithfulness. The nilpotent polynomial framework is also
discussed in [Sobaje, Section 3.1](https://arxiv.org/pdf/1708.04153) and
the [earlier assessment](CAYLEY_FIBER_ASSESSMENT.md).

The construction (3)-(6) and the resulting fixed-irreducible-family
realization theorem are deductions proved here. A targeted literature
search did not locate this formulation, but that does not establish
novelty; the proof uses elementary identities. In particular no claim
to have resolved an established open problem is made.

The [independent verifier](tools/verify_irreducible_centralizer_realization.py)
checks the Pascal inverse symbolically, constructs semisimple companion
blocks over Q(i), and computes the actual middle-exterior group matrices
and Lie actions. It checks all ambient tangent directions, not just a
single moment or the chosen centralizer projection. Tests cover pure
orders, a mixed polynomial, and the zero map for m=2,3,4, including one
fully symbolic upper-nilpotent block for m=3. Companion blocks are
diagonalizable because their polynomials are checked squarefree; conjugating
them to a torus gives the theorem's diagonal setting. The verifier does
not call the runtime character or root-Borel APIs. It shares Sage arithmetic
and signed-wedge helpers with the earlier audit; this is an independent
calculation route, not an external review or a separate-CAS proof.

The result settles representation-realization for this growing ambient
family, including irreducibility. It does not classify all semisimple
strata for a fixed arbitrary representation, exceptional groups, the
full-group generic degree, or the complete fiber over theta(t).

## 7. Next step and its purpose

Stop asking whether higher orders can occur: this theorem answers that
without a bound across ranks. The next useful question is what the
realization says about the full generalized Cayley map. For the zero
model, determine whether t U(H) is a fiber component and compute the
transverse differential and local fiber equations. Success would turn
the uniform construction into a theorem about non-finite fibers of theta,
rather than another existence example inside a centralizer. The distinction
between source critical points and target branch/non-finite loci must be
retained throughout.
