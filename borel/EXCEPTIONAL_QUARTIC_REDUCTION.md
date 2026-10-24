# From Levi branching to quartic systems: scope and the E6 completion

Date: 2026-10-23. Written uniform reductions, an exhaustive finite
classification, and an E6 consequence of earlier proofs. Independent
exact checks are provided; literature priority remains unverified.

## 1. Setup and the uniform reduction

Let G be a simple complex simply connected group of simply laced type,
rho=Ad, and K a proper Levi whose derived factors are A_m. Normalize

    B_Ad(theta(g),X)=tr(Ad(g)ad(X)).

We study theta^(-1)(h) intersect G_reg as a scheme, for h in z(k) with
Z_G(h)=K. Work on the actual finite central cover [K,K] times Z -> K.
Use fundamental characters on each derived factor and Laurent central
coordinates z. Retain ALL logarithmic target equations D_j chi=c_j,
with c_j=B_Ad(h,H_j). No parabolic cohomology or whole-fiber assertion
is involved.

The [adjoint norm theorem](BRANCHING_NORM_FIBERS.md#3-the-adjoint-cubic-obstruction-has-a-restricted-shape)
proves that every nonzero central-charge term uses at most one fundamental
character from each factor and at most three factors. The intrinsic
adjoint character of A_m is x_1 x_m-1 for m>=2, and x^2-1 for A1.

**Theorem 1 (uniform factor elimination).** For every A_m factor with
m>=2, put a=x_1, b=x_m and let u=(x_2,...,x_(m-1)). The character has
the form

    chi=ab+aF(y,z)+bG(y,z)+H(u,y,z),                        (1)

where y consists of all other factors' characters. The degrees of F,G
in y are at most two; neither depends on u. H has degree at most three
and is affine linear jointly in u. The entire critical scheme, including
the central target equations, is isomorphic to that of

    Psi=H-FG, a=-G, b=-F.                                (2)

Its degree is at most four and its homogeneous quartic part is
-F_2 G_2. All middle characters u remain linear multipliers with their
original coefficient equations. The reduction is defined over the
entire central torus and is compatible with central descent.

If an A1 factor is chosen instead, the character has the form

    chi=x^2+xF(y,z)+H(y,z).

Eliminating x=-F/2 gives Psi=H-F^2/4, also of degree at most four,
with quartic part -F_2^2/4. This includes the case of all-A1 Levis.

Proof. The norm theorem forbids an external term containing two
characters of the same factor. Thus ab is the only term using both
endpoints, and no coefficient of a or b uses a middle character. This
proves (1). Complete the product:

    chi=(a+G)(b+F)+H-FG.

The affine change is globally invertible; the derivatives of a,b set
the two new coordinates to zero. The remaining derivatives, including
D_j, are exactly those of H-FG by the chain rule. The A1 argument is
the same completion of a square. Center weights of a,b are dual, so
G,F transform as a,b, respectively; the substitution graphs are
equivariant. The A1 center acts by sign and preserves its square
completion. This proves the scheme and descent assertions.

**Corollary 1.1.** A nonempty original-coordinate linear-multiplier
system with Laurent-unit constraint determinant exists if and only if
the all-A Levi has at least one A_m factor with m>=2. The forward
implication excludes all-A1 Levis because each coordinate has its own
intrinsic square; the reverse implication uses one endpoint as multiplier
and the other as paired variable, with coefficient 1. The torus Levi
has no such coordinates. This is an existence criterion, not a claim
that every possible elimination order or coordinate change succeeds.

These statements hold in simply laced types in arbitrary rank. They
are deductions from the branching shape and elementary elimination,
not new smoothness theorems for arbitrary quartic polynomials.

## 2. Complete first-step classification in E6, E7, and E8

Here an admissible simultaneous system means exactly the form in
[Theorem 1 of the multiplier note](MULTIPLIER_BRANCHING_FIBERS.md#1-exact-simultaneous-elimination):
choose disjoint lists tau,a in the ORIGINAL fundamental coordinates,
of the same size k; require chi to be jointly affine in tau, its tau
coefficients to be C(z)a-f, and det C to be a Laurent unit. No coordinate
changes precede this test. All such choices are exhausted, including
matrices with several nonzero determinant terms before cancellation.
Thus the test is stronger than restricting to unique matchings.

**Theorem 2 (finite classification).** The proper standard subsets give:

| Ambient type | All-A Levis | Quadratic character | Cubic character | Non-A Levis |
| --- | --- | --- | --- | --- |
| E6 | 60 | 49 | 11 | 3 |
| E7 | 120 | 88 | 32 | 7 |
| E8 | 240 | 147 | 93 | 15 |

The torus subset is included; the full-group subset is excluded.
For the cubic cases the maximum simultaneous system size is:

| Ambient type | k=0 | k=1 | k=2 | k>=3 |
| --- | --- | --- | --- | --- |
| E6 | 5 | 6 | 0 | 0 |
| E7 | 12 | 19 | 1 | 0 |
| E8 | 28 | 60 | 5 | 0 |

Precisely the all-A1 cubic cases have k=0. Among these finite E-type
cases, k=2 occurs exactly for derived type A1+A2+A3 or A1+A2+A4.
All other cubic cases with an A-factor of rank at least two have k=1.
This last maximum-size assertion is an exhaustive finite result,
not an arbitrary-rank diagram theorem.

For a reproducible residual, choose a maximum-size system, breaking ties
lexicographically by multiplier-node list and then paired-node list.
If k=0, complete the square at the smallest retained node. These choices
give exact quartic systems with the following numbers of remaining
AFFINE variables; all central torus variables remain:

| Ambient type | Two affine variables | Three | Four |
| --- | --- | --- | --- |
| E6 | 10 | 1 | 0 |
| E7 | 23 | 9 | 0 |
| E8 | 53 | 34 | 6 |

In particular every genuinely cubic all-A E6/E7/E8 case reduces globally
to a degree-four potential in at most four affine variables. This does
not say its full critical locus is finite or smooth. The number of
central variables can exceed one; it is not included in the bound four.
The factored quartic description in Theorem 1 concerns its specified
single-factor reduction. We do not assume the quartic part of every
maximum-system residual is a single product without checking it.

The six four-variable cases are E8 retained node sets

    {1,2,3,4,6,8}, {1,2,4,5,6,8},
    {1,2,5,6,7,8}, {2,3,5,6,7,8}       (A4+A1+A1),
    {1,2,3,5,6,8}, {1,2,3,5,7,8}       (A2+A2+A1+A1).

This is a list of maximal residual sizes under the specified first step,
not a list of exactly the unsolved geometries. Further nonlinear changes
or additional eliminations may simplify these and other residuals.

### Proof certificate and independent verification

There are 63+127+255=445 proper standard subsets. Classify their induced
Dynkin components and reconstruct the entire adjoint character in each
central charge. The 25 non-A cases are excluded from the all-A census;
the earlier root-length theorem makes their regular fibers empty.
The remaining characters have degree at most three by the norm theorem.

For each of the 136 cubic characters, enumerate all disjoint equal-size
coordinate sets. Second derivatives test joint affinity in tau and
whether the cross-coefficient matrix lies in the Laurent base ring.
Compute its exact determinant, retaining central exponents; it is a
unit exactly when it has one nonzero Laurent monomial. This finite
enumeration proves both the positive and negative maximum-size claims.
Apply the elimination to the chosen system and verify all original
affine and central derivatives. The residual terms prove the last table.
Each fundamental character is an eigenfunction of the finite central
kernel. The chosen derivative equations are therefore equivariant; their
unique solved coordinates make the graph equivariant as well. Thus the
scheme isomorphisms descend with the induced action retained, even when
that action has not been written as individual generators in the report.

The [Sage verifier](tools/verify_unit_constraint_classification.py) uses
full fundamental-character inversion and sparse permutation determinants.
The [native audit](tools/audit_unit_constraint_classification.py) uses
Cartan reflections, exterior subset models, SymPy matrix determinants,
and separate sparse rational polynomial substitution. It reconstructs
every cubic character, repeats the exhaustive coordinate search, and
checks EVERY original derivative after substitution. It also enumerates
all 445 subsets independently to verify completeness.

The [summary and Weyl words](unit-constraint-classification-verification.json),
[compressed exact cases](unit-constraint-cases.json.gz), and
[independent audit](unit-constraint-classification-audit.json) retain the
actual central cocharacters, cover degrees, constraints, residual
coefficients and scope. There is no specialization z=1 or generic-source
localization. The compressed archive is hash-bound to the summary.

## 3. All proper E6 adjoint Levi strata: a geometric consequence

**Theorem 3.** Let G be complex simply connected E6 and theta=theta_Ad
with the trace normalization above. For EVERY proper Levi subgroup K
there is a dense open subset U_K of z(k), consisting of targets whose
centralizer is K, such that for every h in U_K the actual scheme

    theta^(-1)(h) intersect G_reg

is smooth, possibly empty. For derived types A1^3, A1^2+A2 and A1+A2^2
these fibers are nonempty of pure dimensions 6, 10 and 14 respectively.
For Levis with a D-type factor the regular part is empty for every
target whose centralizer is K.

This is generic WITHIN EACH centralizer stratum. It asserts neither
global regularity preservation nor smoothness at every target. It does
not determine all nonemptiness cases, all orbit counts, whole fibers,
or the zero-target fiber (whose centralizer is G).

Proof. Conjugate K to a standard Levi. For the 49 quadratic cases apply
the [quadratic incidence theorem](QUADRATIC_BRANCHING_FIBERS.md), which
includes rank-drop components and proves generic smoothness or emptiness.
The three non-A cases have D-type factors, so apply the
[root-length obstruction](EXCEPTIONAL_ADJOINT_CONTACT_FIBERS.md#2-a-root-length-obstruction-for-levi-factors).
The eleven cubic cases fall into three Weyl-conjugacy classes:

| Derived type | Seed retained nodes | All standard cubic embeddings |
| --- | --- | --- |
| A1^3 | {2,3,5} | {1,2,5}, {1,2,6}, {1,4,6}, {2,3,5}, {2,3,6} |
| A1^2+A2 | {2,3,5,6} | {1,2,3,5}, {1,2,3,6}, {1,2,4,6}, {1,2,5,6}, {2,3,5,6} |
| A1+A2^2 | {1,2,3,5,6} | {1,2,3,5,6} |

This is verified by explicit Weyl words mapping the ENTIRE seed root
subsystem onto the target subsystem, not by matching Dynkin types alone.
The certificate lists words applied left to right. Each is checked by
the reflection formula r -> r-r_i*column_i(Cartan) in fundamental-weight
coordinates. For example [1,3] transports the seed {2,3,5} to {1,2,5}.

A normalizer representative of each word conjugates the Levi groups,
their centers and the target open sets. Adjoint character, trace form
and Springer map are equivariant under that conjugation, so it identifies
the fiber SCHEMES. The [two cubic refinement proofs](CUBIC_LEVI_REFINEMENTS.md)
give the first two seeds; the [E6 suspension theorem](BRANCHING_NORM_FIBERS.md)
gives the last. Their dimensions and nonemptiness therefore transfer to
every embedding in the table. This covers all 63 proper standard cases,
and conjugation covers arbitrary proper Levis.

This completes the E6 portion of the generic-smoothness question by
organizing existing proofs and adding exact conjugacy coverage. It is
not eleven new local fiber calculations.

## 4. What is known and what to do next

Published inputs include [Kostant's irreducible Levi central-root spaces](https://arxiv.org/pdf/0711.2809),
standard type-A fundamental-character identities and Weyl conjugation,
and the [Kostant--Michor trace-map construction](https://www.mat.univie.ac.at/~michor/kostant.pdf).
Completing products/squares is elementary algebra. The norm, quadratic,
root-length and three seed-fiber results are earlier written repository
deductions; the links above state their hypotheses and proof status.

The additions here are the uniform quartic reduction, finite exhaustive
constraint classification, and the E6 all-proper-Levi consequence. A
targeted literature search on 2026-10-23 did not identify this precise
package, but was not a complete priority investigation. No claim of
originality or external peer review follows.

The E7/E8 all-Levi smoothness question remains open here. The cubic
case counts 32 and 93 include previously resolved cases; they are not
counts of unresolved problems. The finite classification does not cover
arbitrary affine changes or every recursive elimination sequence.
G2 and F4 are not simply laced and are outside this theorem.

**Next step:** prove a generic reducedness/smoothness criterion for the
factorized quartic residuals, retaining their central target derivatives
and coefficient rank-drop loci. Begin by grouping the remaining E7
systems under Weyl conjugacy and exact residual equivalence, then select
one obstruction shared by a family. Success means a theorem resolving
that family or an actual adjoint counterexample; another table of sample
points would not answer the remaining question. The six E8 four-variable
systems give a precise later test of whether such a criterion scales.
