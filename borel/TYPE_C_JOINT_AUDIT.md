# Type C joint images: proof and literature audit

Audit date: 2026-10-19. Proof baseline: `22d9969a55dee0e28e4ecea6cb6650890a78edef`.
Scope: [the normalized type C joint-image theorem](TYPE_C_JOINT_PRIMITIVITY.md),
including its integral kernel and rank corollaries. This is a mathematical
and source audit, not external peer review or a certificate of originality.

**Verdict.** No mathematical gap was found under the stated assumptions
a,b,k>=1, c_j=e_j(M), q_i=e_i(T^2), and integral singular cohomology of the
compact symplectic flag. The joint image is saturated; the sum of the
separate images is primitive; their intersection is the constants.
The kernel and rank formulas survive the audit unchanged.

**Literature classification.** These are corollaries of classical integral
homogeneous-space and Schubert theory in the selected Springer/Borel
coordinates. The exact combined statement was not located in the sources
and searches below. That absence is not evidence sufficient to claim a new
theorem. In particular, the kernel presentation and joint saturation should
not be promoted as new cohomological phenomena.

## 1. Check of the integral argument

Use X, Y, Q, R, L, M_0 and J as in the proof. The delicate steps are the
identification of the actual integral image and the simultaneous basis,
not the rational Hilbert series.

| Step | What must be established | Audit finding |
| --- | --- | --- |
| Quaternionic base | R is the integral block-characteristic ring, free over Z | Borel's hypotheses apply; no division by a Weyl-group order is used. |
| Identification with a complex flag ring | The isomorphism preserves both selected block subrings | It matches characteristic generators after ordering the blocks M,A,T. The relevant complex flag has steps b and n-k. |
| Common base basis | R_M and R_T are coordinate submodules in the same basis | Grassmannian pullback preserves opposite Schubert classes with minimal right-coset labels. The allowed descent sets are {b} and {n-k}. |
| Global fiber lifts | Fiber Schubert representatives exist integrally on Y | The polynomials in the Chern classes of E_M dual exist globally. The absolute Lagrangian Schubert theorem suffices. A flag of isotropic bundles on Q is unnecessary. |
| Joint generation | The c_j and q_i generate all of H*(Y;Z) | Squared middle classes are integral polynomials in c. Whitney inversion eliminates the head classes. The remaining fiber lifts are also polynomials in c. |
| Exact separate image | The middle subalgebra equals, rather than merely contains, the claimed module | The Cartesian square with Z_M identifies the same lifts on both sides. Its image is the direct sum of u_lambda R_M. |
| Primitive sum | The two separate images have a common integral basis | In the basis u_lambda sigma_w, the middle image uses all lambda and middle-only descent labels; the tail uses empty lambda and tail-only labels. Their sum is a coordinate subgroup. |
| Ambient saturation | J is a unit summand in H*(X;Z) | The head Lagrangian fiber has global integral Chern-polynomial lifts including 1. Leray–Hirsch gives a free complement to J. |
| Integral kernel | Elimination introduces no torsion or hidden relations | Integral free base change and triangular elimination establish the displayed ideal. Borel's general theorem independently gives the same pre-elimination ring. |

For clarity, the common-basis argument is the following purely integral
calculation. Set D={b,n-k}, and transfer the opposite Schubert basis
sigma_w of Fl(b,n-k;n;C) to R. Then

$$\mathcal B=\{u_\lambda\sigma_w:\lambda\text{ strict},\lambda_1\le b,
                    \operatorname{Des}_R(w)\subset D\}$$

is a Z-basis of J. The subset for L has Des_R(w) contained in {b}; the
subset for M_0 has lambda empty and Des_R(w) contained in {n-k}.
These are exactly the images, by the compatible-lift argument, not just
subsets of them. The quotient by their sum has the complementary basis.
Since b<n-k, the intersection of these two index sets is just
(lambda empty,w identity).

The quotient by a sum of two primitive subgroups need not be torsion-free
in general: Z(1,0) and Z(1,2) in Z^2 are a counterexample. The proof avoids
that invalid inference by constructing this common basis.

The conventions also check out. The p_i are signed even Chern classes
with squared roots; C(v)C(-v)=sum (-1)^i p_i v^(2i).
Changing dual-bundle conventions multiplies homogeneous representatives by
signs, not by nonunits. Polynomial degree d corresponds to H^(2d).
Leray–Hirsch is used as a module isomorphism, not a tensor-product ring
isomorphism. All maps involving Y and Q are maps of compact homogeneous
spaces; the argument does not assume that Y is an algebraic partial flag.

## 2. What the literature already supplies

| Reference and checked location | Published result relevant here | Status of our application |
| --- | --- | --- |
| [Borel (1953)](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Borel-Sur.pdf), Prop. 30.2, journal p.201 | Integral characteristic presentation for equal-rank homogeneous spaces under torsion-freeness hypotheses | Directly gives the ring of Y; the joint kernel follows by eliminating head classes. |
| Borel, Prop. 31.1 and numbered case (31.2), p.202 | Complex and quaternionic partial-flag presentations | Gives R and its degree-doubled complex analogue. Case (31.2) is not Proposition 31.2. |
| [Pragacz–Ratajski (1997; linked preprint)](https://arxiv.org/pdf/alg-geom/9602019), Thm. 2.1(i), Props. 4.2, 4.3, 4.7 | Integral Q-tilde Schubert representatives, pair-removal identities and symmetric-polynomial basis | Supplies compatible fiber lifts and the free base change. Q-tilde is essential; the rescaled P-tilde polynomials are different. |
| Pragacz–Ratajski, Appendix A, §§10.1–10.5 | Quaternionic Schubert theory and the quaternionic-span Lagrangian fibration | Already contains the geometric mechanism behind the one-block lift used here. |
| [Brion (2004)](https://arxiv.org/pdf/math/0410240), §1.2, Example 1.2.3(6), Definition 1.2.4; §1.3 | Integral Schubert bases and compatibility with parabolic projections | Taking unions and intersections of the two basis subsets gives the primitive-sum and constant-intersection statements. |
| [Tamvakis, A Giambelli formula for classical G/P spaces](https://www.math.umd.edu/~harryt/papers/giamcl.pdf), Prop. 6 and Cor. 2 | Integral symplectic partial-flag presentation and Schubert representatives | Covers the ambient cohomology; not identified here as a source explicitly stating our selected joint-image quotient. |

For Borel's direct application take G=Sp(n), K=Sp(a) x U(b) x Sp(k).
The group K has rank n. Both groups and both full flag spaces G/T,K/T
have torsion-free integral cohomology. The restrictions of universal
characteristic classes are the coefficients of P_A P_M P_T, with P_M
expressed through the middle Chern classes. Thus this is an application
of the integral theorem with verified hypotheses, not the rational Borel
presentation silently used over Z.

The relative Chow-basis paragraph after Pragacz–Ratajski's Theorem 2.1
assumes a flag on the base. Our proof does not invoke that paragraph on Q:
it uses the absolute fiber theorem, global Chern polynomials, and
topological Leray–Hirsch. This distinction has now been made explicit.

Also inspected: [Mare, Equivariant cohomology of quaternionic flag
manifolds](https://arxiv.org/pdf/math/0605539), Theorem 1.1, and
[Panin–Walter, Quaternionic Grassmannians and Borel classes](https://arxiv.org/pdf/1011.0649),
Theorems 11.2–11.3. These corroborate the established quaternionic
framework. The algebraic oriented-cohomology statements in the latter
are not being applied directly to our compact bases without a comparison.

## 3. Connection to the Kumar papers

[Kumar, Part I](https://kumar.math.unc.edu/papers/kumar76.pdf) constructs
the Springer-induced map for general parabolics over C. It provides the
framework for selecting these factor generators; its complex-coefficient
statement alone does not settle their integral lattices.

[Kumar–Rogers, Part II](https://kumar.math.unc.edu/papers/kumar81.pdf),
Proposition 11 and Theorem 16, treats the symplectic factor on isotropic
Grassmannians. The proof of Theorem 16, journal p.353, explicitly chooses
the polynomial Z-form generated by elementary symmetric functions of the
squared normalized Springer coordinates. Our q_i use that normalization.
The paper's stable injectivity result is different from the finite-rank,
two-factor image and primitive-sum assertions here.

[Kumar–Xie, Part III](https://kumar.math.unc.edu/papers/kumar86.pdf) treats
the even orthogonal case and explicit exceptional-group calculations.
It is not the source of a general type C two-step integral-image theorem.

What this note adds to the repository is the explicit identification of
the two-factor image with H*(Y), its simultaneous integral basis, and
formulas usable by the calculations. These are useful extensions of our
computational treatment of the Springer/Borel map. The audit does not
establish a new mathematical priority claim for them.

Suggested wording for Shrawan:

> Applying classical integral flag-cohomology results to the normalized
> Springer generators, we obtain a saturated joint image on symplectic
> two-step flags. A compatible Schubert and Leray–Hirsch basis shows that
> the sum of the separate images is primitive and that their intersection
> consists of constants. We regard this as a corollary of classical theory;
> we have not located this exact formulation or established novelty.

## 4. Hypotheses that cannot simply be dropped

**Zero head.** If a=0 and b=k=1, the space becomes IG(1,4)=CP^3. With
c=c_1, its integral ring is Z[c]/(c^4), and the Whitney relation gives
q_1=-c^2. Consequently L=J is the whole ring and M_0=Z[c^2] is contained
in it. Their intersection has a nonzero degree-two class. Thus the
constant-intersection assertion, and the correction term in the quotient
rank formula, genuinely require the distinct descent positions supplied
by a>0. This does not contradict the theorem's stated range.

**Scaled input.** For a,b,k>0, replacing c_1 by 2c_1 while leaving the
other generators in their degrees gives degree-one image 2Z beta(c_1).
The original beta(c_1) is primitive and nonzero. The new image therefore
has 2-torsion in its saturation quotient. The theorem is about the fixed
normalized polynomial lattice, not every integral character convention.

The proof also does not settle joint images in types B or D, where the
integral Chern-polynomial fiber-basis step requires separate analysis.

## 5. Search scope and verification

The search combined exact Springer terminology with broader descriptions
of the underlying spaces and subrings. Representative queries on the
audit date included:

- `"Springer" "joint image" cohomology`
- `"joint image" "Borel"`
- `"Representation ring of Levi subgroups versus cohomology ring of flag varieties" integral`
- `"quaternionic flag" "Chern" "Lagrangian"`
- `"quaternionic flag" "cohomology" integral`
- `"symplectic" "flag" "saturated" "image" cohomology`
- `"Schubert" "pullbacks" "sum" cohomology`
- `"flag varieties" "primitive sublattice"`
- `"Sp(n)" "Sp(a)" "U(b)" cohomology`

The positive conclusions above come from reading primary sources, not
from search snippets. Searches using “Springer”, “primitive”, or
“saturation” also return unrelated Springer resolutions, primitive Hodge
cohomology and tensor-product saturation; these are not evidence for this
specific lattice statement. This was a targeted web and primary-source
search, not an exhaustive MathSciNet/Zentralblatt citation review.

The [geometric-basis verifier](experiments/type_c_joint_primitivity.py)
is rerun against the clarified proof. Its [report](experiments/type_c_joint_primitivity_verification.json)
records the exact ranks, block sizes, integral lattice checks and source
hashes; its independent Cartesian BGG check covers C3.

An additional [formal-algebra audit](experiments/type_c_joint_algebra_audit.py)
imports no library code. For 1<=b<=5 and weighted degrees 0 through 14 it
checks pair-removal identities and that the proposed p-monomial times
strict-Q-tilde basis has determinant +/-1 against ordinary Chern monomials.
This tests precisely the potential denominator issue in the integral
base-change step. Its [report](experiments/type_c_joint_algebra_audit_verification.json)
records finite coverage and hashes. Neither experiment substitutes for
the all-rank proof.

Changes from the audited baseline are expository and attributional:
explicit global-lift and Cartesian-square explanations, the triangular
elimination calculation, a second source-based check of the ring, and
this literature classification. No theorem, Smith factor, API behavior,
or historical calculation was changed.
