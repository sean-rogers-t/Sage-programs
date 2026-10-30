# Ordinary Springer images: an adapted basis and a divided-power subring

Research checkpoint, 29 October 2026. This concerns ordinary integral
cohomology, not equivariant cohomology or K-theory. The basis construction
below is a constructive consequence of the existing repository prefix
theorem. It is not a new proof of that theorem or an external review.
The divided-power argument uses classical symmetric-function identities;
no literature-priority claim is made.

## 1. The baseline and the exact source

The ledger already records complete normalized single-tail B and D Smith
formulas: [IMG-03 and IMG-05](RESEARCH_LEDGER.md), including even-tail
gluing. There is no remaining even-tail multiplicity gap to close.
Type C supplies the primitive reference image. The open two-step ambient
B/D problem is different. These statements concern selected normalized
generators, not all integral Levi characters for every representation.

Fix a,k>=1, n=a+k. In type B use X_B=OG(a,2n+1); in type D use
X_D=OG(a,2n). For D, k=1 denotes the parabolic omitting both fork nodes.
The selected square inputs are q_i=e_i(x_(a+1)^2,...,x_n^2).
They are the normalized defining/vector Springer inputs, with the Borel
convention beta(omega_i)=sigma_(s_i). Let I_B and I_D^+ be their images,
and S_B and S_D^+ their saturations in integral cohomology. All comparisons
are inside the selected rational image, not all of H*(X;Z).

Use head coordinates y_i=x_i^2. The common source modulo its kernel is

    R = Z[e_1(y),...,e_a(y)]/(h_(k+1)(y),...,h_(k+a)(y)).

Its Schur basis is indexed by lambda contained in k^a. The head/tail
change sends s_lambda(tail^2) to (-1)^|lambda| s_(lambda')(head^2).
Weight m in R is polynomial degree 2m and cohomological degree 4m.
Write phi_B and phi_D for the induced injections of R into cohomology.

The essential earlier inputs are the
[primitive prefix theorem](TYPE_B_GENERAL_SMITH.md), its
[actual-head extension](TYPE_D_HEAD_EULER_B_COMPARISON.md#2-primitivity-for-the-actual-head),
and [B/D coefficient-row equality](TYPE_D_SQUARE_THEOREM.md).

## 2. One explicit basis for all prefixes

Put r=floor(n/2), c=ceil(n/2), epsilon=n mod 2. Pad every partition
gamma in k^a to length a and define

    beta_i = gamma_i+a-i,       1<=i<=a,
    d(gamma) = #{i: beta_i<r},
    e(gamma) = #{i: beta_i>=c+d(gamma)}.

Set t=e(gamma), b=a-t. Choose the least integer j satisfying

    max(0,b-r) <= j <= min(b,c,floor((b+epsilon)/2)),
    (gamma_(t+1),...,gamma_a)
        <= ((c-j)^j,(r-b+j)^(b-j))

coordinatewise. Such a j exists by the earlier corner criterion.
Split gamma into its first t parts mu, next j parts alpha, and remaining
b-j parts beta, and set

    v_gamma = h_mu s_alpha s_beta in R.

Treat this as a coefficient ROW in the Schur-coordinate lattice of R_m,
not as the source class whose Borel image will be divided.

**Simultaneous basis theorem.** The rows v_gamma form an integral basis
of R_m. Moreover the rows with e(gamma)<=t form an integral basis of
the coefficient prefix P_t, for every t at once.

Proof. Order the complete-product coordinates h_gamma by increasing
sum gamma_i^2. Jacobi--Trudi gives v_gamma=h_gamma plus terms of strictly
higher energy; terms with an index exceeding k vanish. Thus the FULL
matrix in complete coordinates is unit triangular. The change between
the rectangular complete and Schur bases is unimodular. Each v_gamma
belongs to P_(e(gamma)) by the earlier unbalanced-incidence argument.
For any t the selected rows are part of an integral basis, hence span
a primitive sublattice. Their number equals rank P_t by the earlier
corner/count theorem. Inclusion and equal rank force equality. This
uses the previously proved integral prefix theorem; it does not replace it.

Let V_m have rows v_gamma in Schur coordinates. Define f_gamma by the
gamma-th COLUMN of V_m^(-1), interpreted in the source Schur basis.
Equivalently, the coefficient pairing satisfies

    <v_delta,f_gamma> = delta_(delta,gamma).

This is the ordinary coordinate pairing (the Schur basis is orthonormal),
not the Poincare pairing between complementary cohomological degrees.
The inverse is integral; it can be computed by finite unit-triangular
back substitution and the integral complete/Schur change of basis.

## 3. Explicit integral classes and multiplication

**Theorem.** In every weight m the classes

    b_gamma^B = phi_B(f_gamma)/2^e(gamma),
    b_gamma^D = phi_D(f_gamma)/2^e(gamma)

are integral bases of S_B,m and S_D,m^+, respectively, and

    I_B,m = direct_sum_gamma 2^e(gamma) Z b_gamma^B,
    I_D,m^+ = direct_sum_gamma 2^e(gamma) Z b_gamma^D.

Proof. The earlier row formula and the simultaneous basis give

    row(phi_B) = sum_t 2^t P_t
               = direct_sum_gamma 2^e(gamma) Z v_gamma.

A rational source vector z has integral Borel image exactly when every
row pairs integrally with z. The dual of this row lattice therefore has
basis f_gamma/2^e(gamma). This proves both integrality and that the
displayed classes exhaust the saturation. The B/D row equality gives
the identical statement for D. No guessed division of an individual
Schur generator is used.

This also gives a direct membership test: phi_B(z) is integral exactly
when 2^e(gamma) (V_m z)_gamma is an integer for every gamma. Once expressed
in the b basis, membership in the original image means that its gamma
coordinate is divisible by 2^e(gamma).

Write the product in the integral source basis as

    f_alpha f_beta = sum_gamma C_(alpha,beta)^gamma f_gamma.

The constants C are obtained from ordinary rectangular Schur multiplication
and the displayed unimodular matrices. Then

    b_alpha b_beta = sum_gamma
       2^(e(gamma)-e(alpha)-e(beta)) C_(alpha,beta)^gamma b_gamma.

Every displayed coefficient is integral: the saturation of a subring
in a torsion-free ring is again a subring. In particular

    2^max(0,e(alpha)+e(beta)-e(gamma)) divides C_(alpha,beta)^gamma.

This is an exact constructive multiplication rule, not a positive
Littlewood--Richardson rule or a claim that the b_gamma are individual
Schubert classes. Their ambient Schubert coordinates are obtained by
applying the existing characteristic map to f_gamma and dividing as above.

**Stronger B/D comparison.** Restriction along X_D -> X_B induces a
graded ring isomorphism S_B -> S_D^+, carrying I_B onto I_D^+ and every
b_gamma^B onto b_gamma^D. Indeed it carries the same rational source
polynomial to itself, and equality of row lattices identifies exactly
which rational polynomials have integral images. This is a formal
consequence of the earlier row equality, stronger than merely equal
Smith lists. It makes no claim about the additional D Euler sector.

### The canonical multiplicative Smith filtration

Define D_t={f in R:phi_B(f) is divisible by 2^t in H*(X_B;Z)}, and let
F_t be its image in R/2R. The basis theorem gives

    D_t = direct_sum_gamma 2^max(t-e(gamma),0) Z f_gamma,
    F_t = span_F2 {f_gamma mod 2 : e(gamma)>=t}.

This definition is independent of the chosen basis. It satisfies
F_i F_j contained in F_(i+j), F_0=R/2R, and F_(ell+1)=0, where
ell=min(a,k). Thus the Smith multiplicities are the graded dimensions
of a multiplicative filtration on the ordinary Grassmannian ring
over F2. This follows directly by multiplying integral divisible
images; the structure constants above compute the filtered products.
The same filtration is obtained from D by row equality.

The first step F_1 is the kernel of Frobenius f -> f^2 on R/2R.
This is the existing zero-negative-entry argument in the
[unit-factor proof](TYPE_B_UNIT_SMITH.md#2-the-zero-negative-entry-comparison-matrix):
the surviving rows are ordinary Grassmannian restriction of f(x^2),
which equals f(x)^2 modulo two. This identifies kernels, not just ranks.
Higher F_t are additional divisibility data, not iterates of Frobenius
and not asserted to be powers of F_1. In particular an element of F_1
has square zero in R/2R, while its integral divided powers can occupy
higher Smith levels.

## 4. A uniform divided-power subring

For this section put the smaller block at the head: a=ell<=k=A.
For the original opposite block order use the established integral
block transfer; individual original Schubert rows are not preserved.
Thus c=ceil((A+ell)/2)<=A and 2c>=n=A+ell.

Let nu have at most ell parts, each in [c,A], and let m_nu(y) be the
monomial symmetric function. Write l(nu) for its number of parts.

**Corner theorem.** For these partitions,

    e(nu)=l(nu),             f_nu=m_nu(y) in R.

Consequently

    z_nu = phi_B(m_nu)/2^l(nu)

is a member of the integral saturation basis. The same holds in the
D square sector. Their span, including z_empty=1, is a primitive graded
subring, with multiplication

    z_lambda z_mu = [product_(j=c)^A
          binomial(m_j(lambda)+m_j(mu),m_j(lambda))] z_(lambda union mu)

when l(lambda)+l(mu)<=ell, and zero otherwise. Here m_j is multiplicity
of the PART j, not the monomial symmetric function.
Equivalently this is the divided-power algebra on generators of weights
c,c+1,...,A, truncated by total divided-power length greater than ell.
Its intersection with the original normalized image is exactly

    direct_sum_nu 2^l(nu) Z z_nu.

### Proof of the basis identification

The small beta numbers of nu are precisely the ell-l(nu) padded-zero
ones; the remaining beta numbers exceed the threshold. Hence e(nu)
=l(nu), and its chosen row is h_nu. To compute its dual, recall the
classical Hall duality <h_lambda,m_nu>=delta_(lambda,nu).

The Schur expansion of m_nu is supported on partitions dominated by nu,
so its largest part is at most nu_1<=A. The full Schur expansion of
v_delta=h_mu s_alpha s_beta has length at most ell (Pieri/LR). Thus
discarding Schur terms outside the rectangle in either factor does
not change their Hall pairing. This support check is essential: the
Hall pairing does not descend to an arbitrary quotient without it.

In Jacobi--Trudi for alpha, each nonzero complete index is at most
alpha_1+l(alpha)-1<=c-1. For beta it is at most r-1<c. Therefore a term
h_nu, all of whose positive indices are at least c, can occur in
h_mu s_alpha s_beta only if alpha and beta are empty. In that case
delta=mu, and its coefficient is one precisely for delta=nu. Hall
duality proves that m_nu is the required f_nu. This proof also shows
why arbitrary Schur polynomials divided according to e(nu) would fail.

### Proof of the product and truncation

Pass to the complete ordinary flag splitting ring on n variables y_i.
The Grassmannian pullback is an integral injection. In the complete
coinvariant ring, product_i(T-y_i)=T^n, hence every y_i^n=0. A collision
of two exponents from [c,A] has sum at least 2c>=n and vanishes. The
only surviving products of monomial symmetric functions assign their
parts to disjoint head variables. Counting repeated equal parts gives
the displayed product of binomial coefficients. More than ell parts
cannot be assigned disjointly, proving truncation. The powers of two
in the B/D normalization multiply by total length, so the same product
law holds for the z classes. Primitivity and the intersection statement
follow because these are a subset of the simultaneous integral basis.

### The first class and every first layer

For 0<=e<=ell set

    z_e=z_((c^e))=phi_B(m_((c^e)))/2^e.

Then z_1=phi_B(p_c(y))/2, z_1^e=e! z_e, and
z_i z_j=binomial(i+j,i) z_(i+j), with z_(ell+1)=0.
Each z_e is primitive and has exact order 2^e modulo the normalized
image. It realizes the first occurrence (weight ec) of the exponent-e
Smith factor. For e=1 this is the first discrepancy between image and
saturation: there is none below c, and at c the quotient is Z/2.

These are explicit expressions in the original Springer generators.
Newton identities express p_c in the elementary q_i, and the identity
of total squared Chern classes gives p_c(head)=-p_c(tail). Thus, in
the original smaller-head convention, the first class is precisely
-beta(p_c(tail))/2. Higher repeated-part classes are
beta(p_c(head)^e)/(2^e e!); the proof above justifies every division
inside integral cohomology, rather than inferring it from rational
Chern-character notation. The total source relation also forces
p_(2c)(head)=0, so the factorial identities already hold integrally
in the rectangular source after expressing them as monomial functions.

This exhibits an integral reason why multiplying the first new class
need not produce the higher primitive classes: the factorial is real,
not a change in normalization. It does NOT imply indecomposability in
the entire saturation ring; products involving other classes must be
checked separately. Nor does this divided-power corner exhaust every
Smith layer. At length ell it includes all largest-exponent classes,
by the earlier highest-factor partition criterion.

### A concrete modular consequence

The divided-power corner is primitive in S_B, and S_B is primitive in
H*(X_B;Z), so its reduction injects into H*(X_B;F2). Write u_(j,s) for
the reduction of z_((j^(2^s))), where c<=j<=A and 2^s<=ell. Its ring is

    F2[u_(j,s)]/(u_(j,s)^2, all squarefree monomials
                with sum of their 2^s values greater than ell).

Indeed the binomial product rule and binary carry criterion give square
zero for each generator. The binary expansion of every part multiplicity
expresses every basis class as a product with odd coefficient, hence
coefficient one in F2. Those products are distinct basis classes until
the length bound is exceeded. This is the standard mod-two structure
of a truncated divided-power algebra, applied to the explicit integral
classes above; it is not a claim of a new abstract divided-power identity.

For the literal defining/vector polynomial-character form, the original
image in weight m is multiplied by 2^(2m). Its corner intersection is
therefore direct_sum_nu 2^(2|nu|+l(nu)) Z z_nu. The saturation and its
ring are unchanged, but the normalized Smith exponents are not the
literal-character exponents.

## 5. What is established and what is not

Published inputs: rectangular Grassmannian presentations, Schur/complete
and monomial Hall duality, Jacobi--Trudi, Pieri/LR, and complete flag
splitting. The earlier repository prefix and B/D row theorems remain
substantial dependencies, not newly verified published results.

The new deductions here are a single simultaneous explicit basis, its
dual integral source/ambient basis, the multiplication/divisibility rule,
and the divided-power corner with its exact image intersection. The
B/D saturated ring comparison is explicitly a corollary of existing work.

The [full D Euler-glued basis](TYPE_D_EXPLICIT_INTEGRAL_BASIS.md) is now
constructed by binary row reduction on the unit sector, with every
higher-level basis vector retained. The [geometric denominator formula](CLASSICAL_SATURATION_GENERATORS.md)
also identifies the saturation with sum_e 2^(-e) ker(R->H*(W_(e-1);Z)).

Unresolved: a geometric cycle interpretation or positivity for the full
b basis; a compact
generator-and-relation presentation for the entire saturation; the
ambient B/D two-step Smith problem; and literature priority. No new
API, PDF, arbitrary-representation theorem, or equivariant claim is made.

[Sage verification](type-b-adapted-basis-verification.json) checks every
prefix and multiplication in eight rectangles, 96 recorded actual B
matrix embeddings, and ten divided-power-corner cases.
[The native audit](type-b-adapted-basis-audit.json) separately uses
permutation Jacobi--Trudi, Pieri chains, integer matrix inversion, and
Hermite forms. It replays the 96 embeddings and checks five product
cases and eight corner cases. It shares the actual B matrix input;
it is not an independent reconstruction of the BGG engine. The proofs
above, not those finite counts, supply the arbitrary-rank statements.

A further [native Chevalley audit](type-b-basis-chevalley-audit.json)
DOES recompute actual embeddings independently of BGG: signed permutations,
positive coroots (including the short-coroot coefficient 2), and divisor
multiplication reproduce all coefficient-row multisets in 26 recorded
B components. The same independent method checks the D components below.

### Literature comparison performed in this session

[Morrison--Sottile, *Two Murnaghan--Nakayama rules in Schubert calculus*,
Section 2](https://franksottile.github.io/research/pdf/Murnaghan.pdf)
records Jacobi--Trudi, the power-sum hook expansion, and the classical
Schubert quotient. Their results put power-sum calculations in established
Schubert calculus; the divided-power multiplication here is an elementary
consequence of that quotient, not a new Murnaghan--Nakayama rule.

The previously identified
[Chari--Schneider--Shereen--Wand paper](https://arxiv.org/abs/1310.5191)
studies current-algebra modules and Demazure flags. The earlier proof note
already matches its character multiplicities to our counting polynomial.
That match is not a construction of our integral ambient classes.

The familiar prime-power arithmetic of binomial coefficients, if used to
analyze divided-power generators, is classical; see
[McTague, *On the Greatest Common Divisor of Binomial Coefficients*](https://arxiv.org/abs/1510.06696).
This focused source check is not an exhaustive novelty audit of the
integral orthogonal subring. In particular, no absence-of-prior-work
claim is made for the basis or the corner embedding.
