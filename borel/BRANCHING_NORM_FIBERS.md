# Structural Springer-fiber criteria from weight norms and character couplings

Date: 2026-10-23. Characteristic-zero deductions with exact independent
corroboration. Published ingredients are identified below; priority is not
established. This implements the theorem-first continuation of GEO-50/51.

Subsequent progress: the [cubic refinement theorem](CUBIC_LEVI_REFINEMENTS.md)
now resolves the A1^3 and A1^2+A2 mechanism tests identified below.
It does not complete the all-exceptional classification.

## 1. Scope, dependency audit, and the question

G is simple complex simply connected, initially simply laced. V is an
actual nontrivial rational representation, with character chi and trace
form B_V. Normalize its Springer map by

    B_V(theta_V(g),X)=tr(rho(g) d rho(X)).                         (1)

For a proper Levi K and h in z(k) with Z_G(h)=K, we study the actual
scheme F_h^reg=theta_V^(-1)(h) intersect G_reg. These are not Springer
resolution fibers, cohomology images, or whole group fibers. Normalize
root lengths squared to two when G is simply laced.

The dependency audit revisited the following arguments rather than
assuming a new generic-smoothness theorem:

| Input | Audit finding and limitation |
| --- | --- |
| [Critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md) | Equivariance gives scheme containment in K; the trace differential identifies the fiber ideal in a Maurer-Cartan frame. Smooth regular quotients give the formal product. Nonregular points do not permit cancellation of the quotient differential. |
| [Quadratic theorem](QUADRATIC_BRANCHING_FIBERS.md) | The dimension bound on the image of the nonsmooth incidence locus is essential. It includes rank-drop loci and nonreduced incidence schemes. Generic smoothness is applied only after excluding their nondominant images. |
| [Square substitution](CLASSICAL_ADJOINT_FIBERS.md#2-a-criterion-for-square-substitutions) | It cancels a square derivative only on its invertible locus; the omitted ramification locus needs a separate nondominance proof. It is not an arbitrary cubic theorem. |
| [Contact transfer and root-length obstruction](EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md) | Simple coefficient zeros are needed for scheme elimination. The E6 central cover identifies the two apparent sheets. Root projection strictness uses the Levi hypothesis. |

No correction to these particular dependencies was found. This is an
internal dependency audit, not an independent human review of every older
proof. Known controls (rank-drop contributions, persistent cubic double
structure, E6 excess, and quotient ramification) remain in the bound
verification reports.

The all-simple-adjoint generic-smoothness question remains open here.
The goal below is a verifiable branching criterion and a precise reduction
of what is left, not an assertion that all residual schemes are smooth.

## 2. A charge-sensitive character degree bound

Work on the finite central cover S times Z -> K, S=[K,K], with
fundamental characters x_i of S and torus coordinates z. The inherited
orthogonal splitting of the real weight space is E_S plus E_Z. Write a
central weight gamma as a vector in E_Z, not merely its coordinate list.
Let R^2 be the largest squared highest-weight norm among the irreducible
summands of V. For V(lambda), R^2=||lambda||^2.

**Theorem A (monomial norm budget).** If z^gamma x^a occurs with nonzero
coefficient in chi|K, and mu_a=sum_i a_i varpi_i, then

    ||mu_a||^2+||gamma||^2 <= R^2.                              (2)

All norms are inherited from G. For d>=1 put

    M_d(S)=min_(i_1,...,i_d) ||varpi_(i_1)+...+varpi_(i_d)||^2,   (3)

allowing repetitions. If R^2-||gamma||^2<M_d(S), the entire gamma
coefficient has character degree at most d-1. In particular the strict
inequality for d=3 in every central charge proves joint quadraticity.
For a torus Levi there are no x variables and no norm test is necessary.

Proof. Every weight of V(lambda) lies in conv(W lambda), so its norm is
at most ||lambda||. A Levi constituent with central weight gamma and
highest S-weight nu comes from an ambient weight nu+gamma. Thus
||nu||^2+||gamma||^2<=R^2. Triangular character inversion expands its
irreducible S-character into products of fundamental characters with
dominant highest weights mu<=nu. Dominance between dominant weights
does not increase the norm, since (nu-mu,nu+mu)>=0. This proves (2).
Fundamental-weight Gram entries are nonnegative. A monomial of degree
at least d contains a degree-d submonomial of no larger norm, proving
the assertion from (3). Work separately in each central charge; no
specialization z=1 or cancellation of distinct charges is used.

This extends the earlier [affine norm-gap criterion](AFFINE_BRANCHING_CLASSIFICATION.md)
to actual character-polynomial degrees. It is an elementary deduction
from standard highest-weight theory, not a new weight-polytope theorem.

For simply laced root systems each fundamental weight has squared norm
at least 1/2. In A_m the least norm is m/(m+1); in D_m it is at least
one; in E6/E7/E8 the least norms are 4/3, 3/2, 2. Orthogonal factors
have zero cross pairing and within-factor pairings are positive. Hence

    M_3(S)>=3/2.                                               (4)

**Corollary A1 (beyond adjoint representations).** If G is simply laced
and R^2<3/2, every Levi restriction is jointly quadratic. For each K,
the regular Springer fibers over a dense open set of central targets
are smooth, possibly empty. On a contributing quadratic incidence
component their dimension is

    dim K-rank G + ell-k-q,                                   (5)

where ell=rank S, k is its generic quadratic-matrix rank, and q is the
codimension of its image in the central torus. Nonemptiness is a
separate dominance question; different components may have different
excess dimensions.

In particular this applies to BOTH 27-dimensional E6 representations
V(omega_1), V(omega_6), because their squared norms are 4/3. It is an
all-Levi written argument, not an inference from 63 separate calculations.
It also applies whenever the same numerical hypothesis holds in A or D,
for example their defining representations. No conclusion for the E7
56-dimensional representation follows from this strict bound: its norm
is exactly 3/2. Failure of the test is not a singularity counterexample.

## 3. The adjoint cubic obstruction has a restricted shape

For the adjoint representation of simply laced G, any Levi containing
a D_m (m>=4) or E-type simple factor already has empty regular fibers
by the root-length obstruction. Consider the remaining case

    S=product_j SL_(m_j+1).

**Theorem B (adjoint branching shape).** Its restricted character is
the sum of the Levi adjoint characters, a central constant, and terms
z^gamma times products of fundamental characters on DISTINCT factors.
Each external term has at most three factors. If it has three, at
least one is A1. Consequently the full character is jointly quadratic
if there are at most two simple factors, or if no factor is A1.

Proof. A root outside the Levi has nonzero projection gamma to E_Z,
so its S-projection has norm strictly less than two. In A_m, every
dominant weight whose fundamental coefficients sum to at least two
has norm at least two: it suffices to use

    (varpi_i,varpi_j)=min(i,j)(m+1-max(i,j))/(m+1)

on pairs; the minimum pair norm is two. For N=m+1 and i<=j, this
quadratic expression is concave in j, so its minimum lies at j=i or
j=N-1. The first endpoint is 4i(N-i)/N>=2; the second is
2+(i-1)*(1-(i-1)/N)>=2. Thus each external constituent
is trivial or a single fundamental on each A-factor. Its character is
the corresponding product, with positive central multiplicity.
Four nontrivial factors have norm at least 4/2=2 and are impossible.
Three factors all of rank at least two have norm at least 3*(2/3)=2
and are impossible. Finally the intrinsic A1 adjoint character is
x_1^2-1, and for A_m, m>=2, it is x_1 x_m-1. These have degree two.

This proves generic regular-fiber smoothness for all the quadratic
cases just specified, by the existing quadratic theorem. It does not
assert that quadratic fibers have zero excess: E6/A5 still has excess
two. Nonemptiness and the values k,q in (5) still require proof.

In particular the unresolved simply laced ADJOINT smoothness question
is confined to all-A Levis with at least three factors, including A1,
and with an actual cubic term. This is a sufficient reduction, not a
claim that every such Levi produces singular fibers. The shapes of the
three-factor terms, not a rank sweep, motivate the tests in section 5.

## 4. An exact suspension transfer beyond quadratic characters

Let z be coordinates on T=(C*)^s, s>=1, and S(z) a Laurent polynomial.
Allow additional free coordinates w. Consider a character potential of
the following form after an invertible change in affine coordinates
over the Laurent ring:

    chi=sum_i a_i b_i + (y^2+S(z)y)(1-sum_j u_j v_j)+constant.   (6)

The target correction is -sum_l c_l log z_l. Affine coordinate changes
may depend on z; transform the ENTIRE potential before differentiating.

**Theorem C (square-discriminant suspension).** Over c!=0 the entire
critical scheme is isomorphic to

    {z in T: D_l(-S(z)^2/4)=c_l for every l} times A_w,
    a_i=b_i=u_j=v_j=0, y=-S(z)/2,  D_l=z_l partial_z_l.          (7)

This is an equality after scheme elimination, not a radical equality.
It holds on every target open c_l!=0, and hence on their union.
There is no assertion about c=0.

Proof. Put P=y^2+Sy and Q=1-sum u_j v_j. The equations, after the
a,b pairs are eliminated, are

    (2y+S)Q=0, P u_j=P v_j=0, y Q D_l S=c_l.                  (8)

Where c_l is invertible, the final equation makes y and Q invertible
in the quotient ring. Thus 2y+S=0 and P=-y^2 is invertible. All u,v
vanish, Q=1, and the remaining equations are (7). Conversely (7)
satisfies (8); its nonzero target equation already makes S invertible.
This proves the exact scheme isomorphism and retains the free w.

The torus map Psi=D(-S^2/4):T->A^s has smooth source. Characteristic-zero
generic smoothness gives generic smooth fibers, possibly empty; if
Psi is dominant, those fibers are nonempty and zero-dimensional.
Consequently a representation satisfying (6) has generic regular fibers
of dimension dim K-rank G+dim w when nonempty.

A concrete sufficient condition for dominance is: S has nonnegative
real coefficients, is nonzero, and its nonconstant exponent vectors
span R^s. At positive real z, S>0 and its logarithmic Hessian is positive
definite. The Hessian of S^2/4 is
(dS dS^t+S Hess S)/2, also positive definite. Thus det dPsi is not
identically zero. This is a written criterion for nonemptiness, not
an assumption that an unspecified residual scheme is smooth.

The suspension includes degree-four potentials and is not the earlier
quadratic theorem relabeled. It is a sufficient normal form; no claim
is made that all representation characters have it.

## 5. Three mechanism tests and an E6 theorem

Use Bourbaki node numbering. The exact charge-preserving tests are:

| G and retained nodes | Derived Levi | Purpose and outcome |
| --- | --- | --- |
| E6, {2,3,5} | A1+A1+A1 | Genuine three-factor cubic; subsequently resolved by the Fricke boundary argument in the cubic refinement note. |
| E6, {2,3,5,6} | A1+A1+A2 | Mixed cubic; subsequently resolved by the product transfer in the cubic refinement note. |
| E6, {1,2,3,5,6} | A2+A1+A2 | Cubic terms pair into the suspension below, giving an exact theorem. |

These are three discriminating tests, not an exceptional all-Levi table.
Their cubic degrees did not prove failure: the subsequent refinement
note establishes nonempty smooth generic regular fibers in both cases.

For the third Levi, set H=(2,3,4,6,4,2) in the E6 simple-coroot basis,
K=Z_G(H), and h=tH, t!=0. On the cover SL3 times SL2 times SL3 times C*,
let (a,b) and (u,v) be the two pairs of SL3 fundamental characters, y
the SL2 defining character, and z its central cocharacter coordinate.
The restricted adjoint character is exactly

    chi=ab+uv+y^2-2+(z^3+z^(-3))y
          +(z^(-2)+zy)av+(z^2+z^(-1)y)bu.                    (9)

Put S=z^3+z^(-3), A=z^(-2)+zy, B=z^2+z^(-1)y. Since
AB=1+Sy+y^2, the invertible substitution

    a'=a+B u, b'=b+A v

turns (9) into

    a'b'+(y^2+Sy)(1-uv)-2.                                  (10)

Theorem C gives a=b=u=v=0, y=-S/2, and

    3*(z^6-z^(-6))+2c=0,  c=B_Ad(h,H)=144t.                  (11)

The trace norm is 144 from the root sum, so this specifies the actual
Springer normalization. Here dim K=20 and rank G=6.

The cocharacter columns for the five retained simple coroots and H
have determinant six. H is primitive. Therefore the multiplication
cover has kernel of order six, projecting isomorphically onto mu6 in
the central z-coordinate. Its action preserves the critical equations;
the zero SL3 coordinates are fixed and y transforms as required by
y=-S/2. On the actual quotient set w=z^6. Its critical scheme is

    Spec C[w,w^(-1)]/(3*w^2+2*c*w-3).                        (12)

For c!=0 and c^2+9!=0, (12) consists of two reduced points. Each has
one regular K-conjugacy class of lifts. Equivariance gives Z_G(g)
subset K for theta(g)=h, so K-regularity equals G-regularity there.
Thus F_h^reg is the disjoint union of two smooth K-orbits, each of
dimension 20-6=14. In particular its excess dimension is zero.
This orbit count concerns the regular part only, not the whole fiber.
At c^2+9=0 the displayed critical scheme is a double point; no claim
about c=0 follows because the target centralizer changes there.

An independent direct adjoint-root evaluation takes z=2, SL3 eigenvalues
(1,omega,omega^2) on both factors, and SL2 eigenvalues lambda,lambda^(-1)
with lambda^2+(65/16)lambda+1=0. It gives

    c=-12285/128, theta(g)=(-1365/2048)H,
    rank dtheta=6+58=64, dim ker dtheta=14.                  (13)

Every root value at this point differs from one. This direct trace
check bypasses character elimination; it corroborates but does not
replace the all-parameter scheme proof.

## 6. Published inputs, verification, and remaining gate

Published inputs and what they establish:

- [Kostant--Michor](https://www.mat.univie.ac.at/~michor/kostant.pdf),
  sections 2 and 4: the trace construction, equivariance, and its
  compatibility with centralizers and semisimple Jordan parts.
- [Steinberg](https://www.numdam.org/item/PMIHES_1965__25__49_0/),
  sections 6--8, especially Theorem 1.5: fundamental-character
  quotient coordinates and regular-quotient smoothness.
- [Kostant, Root Systems for Levi Factors](https://arxiv.org/pdf/0711.2809),
  Theorems 0.1 and 1.9: irreducibility and multiplicity-one of each
  nonzero central-weight space in the adjoint restriction. This is
  published branching structure, not a new claim here.
- Standard highest-weight convexity and triangular character inversion;
  characteristic-zero generic smoothness and the Jacobian criterion.

The norm budget, its quadratic-fiber consequences, the restricted
cubic shape, and the suspension application are deductions written
here. Targeted searches on 2026-10-23 for generalized Cayley fibers,
adjoint Levi branching, regular quotients, and quadratic critical loci
located these inputs; they did not establish priority for our fiber
statements. We do not claim to solve a named Kostant--Michor open problem.

The [Sage verifier](tools/verify_branching_norm_fibers.py) reconstructs
every central charge in the three adjoint tests and the E6 minimal
control, checks monomial norm budgets, exact ideals, and (13).
The [native auditor](tools/audit_branching_norm_fibers.py) independently
enumerates roots by explicit Cartan reflections and builds A-factor
fundamental characters from subset weights. It also checks the
suspension identities and the actual mu6 quotient algebra. See the
[verification](branching-norm-fibers-verification.json) and
[audit](branching-norm-fibers-audit.json) reports. Finite checks are
corroboration; the norm and transfer arguments are the general proofs.

The next gate is the interaction of THREE or more A-factors with an
A1 factor. The subsequent cubic refinement note closes the two tests
above; seek an embedding-sensitive condition explaining when these
product and Fricke reductions extend to other branching configurations.
The all-adjoint generic-smoothness statement remains a conjectural
target, and excess dimensions remain allowed. Do not replace this
gap with an unrestricted exceptional table. Nonregular source points,
nonsemisimple targets, generic degrees, and integral images remain
separate projects; the public API and PDFs are unchanged.
