# A rank formula for individual half-spin Springer degrees

Updated 2026-10-22. This note gives a written characteristic-zero
deformation proof. It is an internal research deduction, not an externally
reviewed or formally verified theorem, and carries no literature-priority
claim. It extends the [Spin_10 calculation](HALF_SPIN_DEGREE.md) and uses
the product character behind the [zero-fiber work](HALF_SPIN_BOUNDARY.md).

## 1. Precise theorem

Let r>=3, G=Spin_(2r)(C), and S_sigma one of its two irreducible
half-spin representations. Define the Springer/generalized Cayley map
theta_sigma:G -> so_(2r) by trace projection:

    tr(dS_sigma(theta_sigma(g)) dS_sigma(Y))
        = tr(S_sigma(g) dS_sigma(Y))  for all Y in so_(2r).

**Theorem.** Its generic degree is

    d_r = 2^(r-1) (r^2-3r-2) + 4r.

The same degree applies to the restriction T -> t and the induced map
T/W -> t/W. It counts distinct points over a general complex target.
Both chiralities have the same degree, since the diagram automorphism
interchanging the spin nodes intertwines the maps.

| r | Group | d_r |
| --- | --- | --- |
| 3 | Spin_6 | 4 |
| 4 | Spin_8 | 32 |
| 5 | Spin_10 | 148 |
| 6 | Spin_12 | 536 |
| 7 | Spin_14 | 1692 |
| 8 | Spin_16 | 4896 |

The rank-seven value also matches a completed original-fiber calculation;
the rank-eight entry is an evaluation of the written formula.
This is a degree of a specified morphism, not the minimum Cayley degree.
For even r the half-spin representation has kernel of order two; the
induced map on its faithful image group has degree d_r/2. For odd r
the half-spin representation is faithful. The theorem itself always
uses the simply connected group and its actual character lattice.

In particular, the requested Spin_12 correction is

    536 = 1376 - 840 = 1376 - 12(30 + 2*20).

Thirty curve strata in each coordinate facet contribute one point
each under a generic deformation. Twenty triple junctions contribute
two each. Deeper faces contribute no escaping solutions; this last
assertion is proved below, rather than assumed from the facet counts.

## 2. Lattice, normalization, and the total intersection

Use the weight coordinates epsilon_1,...,epsilon_r. The weights are

    mu=e/2,  e_i in {1,-1},  product(e_i)=sigma,

each with multiplicity one. Their second moment is 2^(r-3)I. Thus

    theta_sigma(t) = 2^(3-r) sum_mu mu e^mu(t).

It suffices to count the logarithmic gradient of chi_sigma=sum e^mu;
the target scaling is invertible. The lattice is

    M=P(D_r)=Z^r+Z delta,  delta=(1/2,...,1/2),

with basis epsilon_1,...,epsilon_(r-1),delta and covolume 1/2 in the
epsilon coordinates. A weight's exponent vector in this basis is

    ((e_1-e_r)/2,...,(e_(r-1)-e_r)/2,e_r).

Let P_r be the half-spin weight polytope and X_r its normal projective
toric variety in M, with ample Cartier line bundle L_r. The cube
[-1/2,1/2]^r loses 2^(r-1) corner simplices of volume 1/r!, giving

    integral_X_r c_1(L_r)^r = Vol_M(P_r) = 2 r! - 2^r.       (1)

For r>=5 the nonsimplex facets are the 2r coordinate facets
mu_i=+1/2 or -1/2. The other facets are simplices. Any face not lying
in a simplex facet is a smaller demihypercube obtained by fixing some
coordinates. Faces with at most three free coordinates are simplices.
The ranks three and four have only simplex facets and give degrees
4 and 32 directly by the argument of the preceding note.

For a proper face F, the initial gradient equations are

    sum_(mu in F) mu e^mu = 0.

Since its supporting affine hyperplane misses zero, these equations
are equivalent to the face polynomial and all its tangent derivatives
vanishing. On a simplex face, linear independence of the vertices
rules this out.

## 3. A generic deformation and conservation of number

Choose a Laurent polynomial h with independent general coefficients
on the same half-spin weights, and deform chi to chi+eta h. Keep a
general target y fixed. These are computational deformations of
sections; h need not be a representation character.

For general nonzero eta all face polynomials are nonsingular on their
face tori. The compactified gradient equations therefore have no
boundary solutions. They have the full number (1) of solutions, all
simple for general choices. The standard toric top-intersection
formula supplies this count. Generic smoothness follows also from
the incidence equations: weight coefficients give full-rank gradient
evaluation, and characteristic zero gives generic reducedness.

At eta=0 the torus map is generically finite: its normalized derivative
at the identity is the identity. Choose y so its torus fiber is finite
and reduced. Properness of X_r implies that each solution of the
deformation specializes either to one of these torus points or to a
boundary common zero. Each torus point contributes one branch. It
remains to count, without omission or duplication, the branches
specializing to the boundary. One may work with convergent germs or
over an algebraic closure of C((eta)); the latter makes orders precise.

We choose y with all epsilon coordinates nonzero and y_i != +/- y_j
for i!=j. Additional genericity conditions mentioned below are finite
in number at each fixed rank.

## 4. All face critical strata and their local form

On a face with m free coordinates, its character polynomial, after
a monomial unit, is

    F(z) = product_(i=1)^m (1+z_i)
             + tau product_(i=1)^m (1-z_i),  tau in {1,-1}.       (2)

The simultaneous equations F=dF=0 hold exactly when at least two
z_i equal -1 and at least two equal +1. Indeed, if both products are
nonzero, F=0 makes any partial derivative equal to a nonzero multiple
of 2/(1-z_i^2). If both products vanish, a unique zero factor in either
product gives a nonzero derivative in that coordinate. Having at
least two zero factors in both products makes every derivative zero.

Stratify by exactly k coordinates equal to -1, exactly l equal to +1,
and f=m-k-l remaining coordinates q_j in C* minus {1,-1}; k,l>=2.
There are m!/(k!l!f!) such strata for a fixed (k,l,f).

In coordinates transverse to one such stratum, (2) is a unit times
the sum of two products on disjoint variable sets:

    A(a)+B(b),   A=a_1...a_k,   B=b_1...b_l.                (3)

For example, use a_i=(1+z_i)/(1-z_i) near -1 and
b_i=(1-z_i)/(1+z_i) near +1. Factors from the q_j are units. They can
be absorbed into one a and one b coordinate, with invertible Jacobian,
uniformly in q. This gives the displayed local form. For Spin_12 the
facet has m=5: (k,l,f)=(2,2,1) gives a smooth critical curve with a
nondegenerate transverse quadratic form; (2,3,0) and (3,2,0) give
the germ a_1 a_2+b_1 b_2 b_3, whose critical locus has three branches.

We will use a precise ideal observation. If F=U_A A+U_B B with U_A,U_B
units, then

    (A,B) is contained in m_trans J_trans(F),                (4)

where J_trans uses the a,b derivatives and m_trans=(a,b).
Apply the two radial derivatives sum a_i partial_(a_i) and
sum b_j partial_(b_j). Their coefficients on (A,B) form a matrix
with invertible diagonal (k U_A,l U_B) modulo m_trans. Inverting it
proves (4). Derivatives in q or in any normal parameters lie in
(A,B) as well. This observation handles vanishing orders even when
the critical locus is not smooth.

## 5. The free-coordinate calculation

Let H(q) be a general polynomial with every vertex of [0,1]^f in its
support. For general constants kappa_j the system

    q_j partial_j H - kappa_j H = 0,  j=1,...,f,             (5)

has exactly f! simple solutions in the torus; H is nonzero at each.
We use 0!=1. The assertion also holds for a translated unit cube,
with kappa translated accordingly.

Here is the needed proof, including boundary exclusion. Compactify
in (P^1)^f, whose cube line bundle has top self-intersection f!.
On a proper face of the cube, some exponent is fixed at 0 or 1.
The corresponding equation is a nonzero scalar multiple of H_face,
provided kappa_j is neither 0 nor 1. The other equations would then
make H_face singular on its face torus. A general H has no such
singularity on any face. Thus there are no boundary solutions.
The same nonsingularity excludes H=0 in the full torus. Evaluation
of the coefficient parameters in (5) has rank f because the cube
vertices affinely span; the generic incidence fiber is reduced and
zero-dimensional. The intersection count is f!.

A general H also makes these solutions avoid q_j=+1,-1. This follows
from the same incidence rank argument on each of these divisors.

For the restriction of our perturbation h to a stratum, H really has
independent coefficients on a full f-cube. For each pattern of the
free exponents, the fixed k+l coordinates have one parity of binary
exponents. Those exponent vectors affinely span R^(k+l), since
k+l>=4. Evaluating at the fixed nonzero coordinates therefore allows
the value and all k+l transverse first derivatives to be prescribed
independently. The coefficient sets for distinct free patterns are
disjoint. This proves both generality of H and the first-jet freedom
used in the next step; it does not assume that the gradient equations
have independent arbitrary coefficients.

## 6. Contribution of a coordinate facet

On a coordinate facet the transverse primitive toric parameter t
has order one for the constant target monomial; weights on the
opposite facet have order two. Locally the character section has
the form F(z)+t^2 G(z). With logarithmic tangent derivatives D_j,
the deformed equations can be written

    F + eta H - t^2(G+eta K) - y_b t = 0,
    D_j F + eta D_j H + t^2(D_j G+eta D_j K) - y_j t = 0.     (6)

Here y_b is nonzero; monomial units have been included in F,H,G,K.
These are the actual Spin lattice coordinates, as in the Spin_10
proof, so there is no extra factor of two in this local count.

Use the first equation to eliminate t. Observation (4), applied to
the transverse equations, shows along any branch approaching the
stratum that F and its q derivatives are o(eta+t). Hence t=O(eta),
and the free-coordinate equations to first order in eta become (5)
with kappa_j=y_j/y_b (up to the stated monomial translation). A
solution with H=0 would also have dH=0 and is excluded by Section 5.
Thus at each of the f! free-coordinate solutions, t is a nonzero
constant times eta to leading order.

There are (k-1)(l-1) transverse branches over each such solution.
To justify this count, eliminate t and use (3). At eta=0 the
remaining tangent ideal is the Jacobian ideal of A+B: subtraction
of a multiple of F, and the terms of order F^2 from eliminating t,
do not change this ideal. Weighted Euler expresses F in its
Jacobian ideal with coefficients vanishing at the point, so the
resulting change of generators is an invertible matrix. The same
argument works after the unit and coordinate changes above.

For the deformation, the transverse equations consequently have
leading form

    product_(h!=i) a_h = eta alpha_i,
    product_(h!=j) b_h = eta beta_j,                         (7)

with general nonzero alpha_i,beta_j. More precisely their right
sides are eta times units depending on all the local variables.
First-jet freedom from Section 5 guarantees nonzero constant terms.
The free-coordinate Jacobian in (5) is invertible.

The k equations of the first group of (7) have k-1 solutions. If
P=product a_i, then a_i=P/(eta alpha_i), and

    P^(k-1)=eta^k product alpha_i.

All solutions are simple after scaling a_i by eta^(1/(k-1)). The
second group similarly has l-1 solutions. There are no other
valuation patterns: dividing two equations forces all a_i to have
the same order, hence order 1/(k-1), and likewise for b_j.
After a common finite base change in eta, the rescaled system has
invertible Jacobian. Its leading block for q is (5), followed by
the two product systems. The formal implicit function theorem gives
exactly (k-1)(l-1) branches, with no missing higher-order branches.

The contribution of one stratum is therefore

    f! (k-1)(l-1).                                          (8)

## 7. No solutions escape to faces of codimension at least two

This is essential: summing (8) alone would be incomplete. Consider
a coordinate face fixing c>=2 coordinates, and a critical stratum
of its face polynomial. Work at a point of that stratum. On a finite
local toric cover, let the normal parameters be n_1,...,n_c; then
the squared normal torus coordinates are n_i^2 and the constant
target section is

    Z = (unit) n_1...n_c.

The covering group is the even-sign subgroup of {+1,-1}^c. For c=2
the quotient is the ordinary toric surface singularity XY=Z^2.
Any Puiseux branch on the quotient lifts after a finite base change,
so proving nonexistence on the cover suffices; no multiplicity is
being inferred by counting cover points.

The undeformed character section is

    U_A(n,q,a,b) A(a) + U_B(n,q,a,b) B(b),

with both coefficients units at the point. In the product character
formula their normal factors are products of (1+n_i^2) and
(1-n_i^2). Thus (4) applies. Along a hypothetical deformed solution,
the transverse derivative equations imply

    J_trans(F) = O(eta,Z),
    A, B, F, partial_q F, normal derivatives of F = o(eta,Z). (9)

The first assertion follows by absorbing the F terms using (4);
the strict order in the second follows because a,b tend to zero.
These statements are comparisons of Puiseux orders: the order of
the second line is strictly greater than min(ord eta,ord Z).

Let lambda_i=+1/2 or -1/2 be the fixed weight coordinates. The normal
equations now have leading terms

    lambda_i eta H(q) - y_i Z + o(eta,Z) = 0,
          i among the fixed coordinates.                   (10)

Choose two of them. Our generic target condition makes the matrix
with rows (lambda_i,-y_i) invertible. It follows that ord Z>ord eta
and H(q_0)=0. Indeed, order Z<order eta contradicts (10), and at
equal order its invertible two-by-two leading system forces both
leading coefficients to vanish.

The remaining free tangent equations, using (9), then give
partial_j H(q_0)=0 for every free coordinate. But the restriction
of general h to this stratum is a general full-cube polynomial,
up to a monomial unit, and has no singular zero in its torus. If
there are no free coordinates, its constant value is nonzero.
Both cases are contradictions.

This excludes every deeper coordinate face, including junctions
between closures of facet strata. Simplex faces were already
excluded in Section 2. Points where an extra free coordinate
becomes +/-1 belong to a different (k,l,f) stratum and are counted
there; they are not silently discarded as ends of the curves.

## 8. Summation and the degree

Put m=r-1. The total contribution of one coordinate facet is

    C_m = sum_(k,l>=2, k+l<=m)
              m!/(k! l!) (k-1)(l-1)
        = m! - (m-2) 2^(m-1) - 2.                          (11)

For a direct identity, the generating function of C_m/m! is

    [((x-1)e^x+1)^2]/(1-x)
        = (1-x)e^(2x) - 2e^x + 1/(1-x).

There are 2r coordinate facets. Sections 6--7 account for every
escaping branch, so conservation of number gives

    d_r = 2r! - 2^r - 2r C_(r-1)
        = 2^(r-1)(r^2-3r-2)+4r.

For r=3,4 the same expression agrees with the direct simplex-face
calculation; C_2=C_3=0. For r=6, C_5=70: the thirty (2,2,1)
strata give 30, and the twenty (2,3,0)/(3,2,0) strata give 40.
The sixty codimension-two demihypercube faces have six critical
points each, but Section 7 proves that none receives a deformed
solution for the chosen general target. One must not subtract those
360 points again.

The group, torus and quotient degrees agree by the argument in the
Spin_10 note: a preimage of a regular semisimple y commutes with y
and lies in T, and the two Weyl quotient degrees cancel. This also
agrees with Kostant--Michor's Corollary 3.3.

## 9. Verification and published inputs

Run from borel:

```bash
PYTHONPATH=src sage -python tools/verify_half_spin_rank_degrees.py
```

The [verifier](tools/verify_half_spin_rank_degrees.py) and
[report](half-spin-rank-degrees-verification.json) separate the proof
from the finite evidence. They check the exact Spin_12 polytopes
in P(D6), the representation's index-two weight sublattice, all
facet and codimension-two face counts, the two-product local
identities, transverse first-jet ranks, product-gradient lengths,
free-cube counts, and the combinatorial sum through rank forty.
The generating-function identity is checked symbolically as well.

Original Spin_12 fibers are checked in two equation models: a weight
enumeration in the independent P(D6) basis, and the separate product
character in coordinates z_1,...,z_6,b with b^2=product z_i. In the
latter model,

    F=(product(1+z_i)+sigma product(1-z_i))/2,
    chi=F/b,
    theta_i=(2z_i partial_i F-F)/(16b).

Both ideals are saturated to keep all torus coordinates nonzero.
The native staircase count and Sage normal basis agree. The report
records whether a multiplication characteristic polynomial is
squarefree for each fiber; failure of that test alone would not
prove nonreducedness. These models share Sage/Singular arithmetic
but independently construct the defining equations. Finite-field
counts do not establish the characteristic-zero or all-rank theorem.

A separate [rank-seven replay record](experiments/half-spin-rank7-fiber.json)
records a completed original product-character fiber over F_32003:
chirality +, normalized target (1,2,3,5,7,11,17), length 1692. Its
[reproducible exploratory script](tools/explore_half_spin_rank6_product.py)
accepts the rank as its third argument. The computation took about
438 seconds; its output and source hash are retained, but its full
Groebner basis was not saved and reducedness was not separately
certified. The main verifier checks this recorded source binding;
it does not silently rerun the longer rank-seven experiment.

The trace-projection construction and quotient degree reduction are
published in [Kostant--Michor](https://arxiv.org/pdf/math/0109066).
Their full-spin representation in even dimension is S_+ direct sum
S_-; its degree 2r is a different result. The toric intersection
input is standard (Cox--Little--Schenck, *Toric Varieties*, Theorem
13.4.3, [book](https://doi.org/10.1090/gsm/124)). Formal implicit
functions, conservation of number and generic smoothness are also
standard inputs. Sections 4--8 are the specific deformation argument
recorded here, rather than a theorem quoted from those references.
Targeted searches did not identify this rank formula explicitly;
they do not establish novelty or settle its publication status.

## 10. Next goal and remaining limits

The degree now counts the generic inverse branches. It does not
describe their ramification or determine the maximal open target
over which the original map is finite. The next useful goal is to
derive an effective generic-fiber presentation from the half-spin
Pfaffian formula and recover this degree from it. That would turn
the intersection count into usable inversion data and give equations
for branch and non-finiteness loci. Progress means an exact
presentation with the correct generic rank, a separate check against
the original weights, and a proof of its domain of validity.
The source critical locus, target branch locus, non-finite target
locus, and minimum Cayley degree remain distinct questions.
