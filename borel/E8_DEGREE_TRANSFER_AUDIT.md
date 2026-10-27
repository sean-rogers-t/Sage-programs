# Independent audit of the final E8 degree argument and group transfer

27 October 2026. GEO-92. This is a fresh internal audit of
[the final E8 proof](E8_LAST_ADJOINT_FIBER.md) and the common scheme reduction
in [the consolidated theorem](ADJOINT_GENERIC_FIBER_THEOREM.md). It is not
external peer review or an independent reimplementation of every finite-field
ideal membership. The earlier review package remains pinned to its historical
proof commit; this note is an addendum.

## 1. Verdict and exact scope

No fatal mathematical gap was found in these two arguments. The conclusion
remains generic smoothness of the ACTUAL regular-source fibers of the adjoint
Springer map over each proper Levi central stratum. It allows empty fibers,
regular nonsemisimple points and, in other cases, smooth positive excess.

The new rational audit strengthens the infinity calculation: the four known
components have generic multiplicities exactly 3,1,1,1. These are multiplicities
of the auxiliary homogenized Jacobian-deficiency scheme. They are not
multiplicities or degrees of a Springer map or of its generic fibers.

There are two qualifications to the audit verdict:

1. The two finite-prime ideal-membership replays use the same pinned F4 engine.
   Their agreement is not independent software verification of membership.
   Explicit multiplier vectors or a second membership implementation would
   reduce this remaining computational dependency. Neither is supplied here.
2. The short returning-to-the-group paragraph in the original E8 note must be
   read with the finite-pullback clarification in the consolidated theorem.
   The central target map generally does NOT descend to the trace surface.
   Taking a target image on that surface without this clarification would be
   an unjustified step. Section 4 below gives the valid argument.

This audit addresses the final degree proof and the shared geometric transfer.
It does not newly audit every earlier exceptional seed proof or establish
generic nonemptiness, a uniform root-theoretic proof, or literature priority.

## 2. Fresh characteristic-zero checks

The new [standalone auditor](tools/audit_e8_degree_transfer.py) imports neither
Sage nor any previous audit helper. Starting from the displayed full residual
potential, it uses SymPy sparse rational arithmetic to differentiate, form all
maximal minors by the determinant permutation formula, and homogenize afresh.
It recovers 38 distinct nonzero original generators and matches both archived
certificate inputs exactly. It checks normal orders from these newly computed
polynomials, rather than from a serialized list alone.

At the generic point of L=V(A,B,C,h), set F=Q(x,y,w,a). The three homogenized
gradient equations have the form

    g = M v + h L0 + E,  v=(A,B,C),  E in (A,B,C,h)^3.

The determinant D of M is the nonzero polynomial in the original proof.
Formal elimination gives

    v = -h ell/D + O(h^3),  ell=adj(M)L0.

Substitute in H=h(A^2+B^2+C^2)-ABC-4h^3. Then

    H = (c/D^3) h^3 + O(h^5),
    c = D(ell_A^2+ell_B^2+ell_C^2)+ell_A ell_B ell_C-4D^3.

The fresh check at (x,y,w,a)=(1,2,3,4) gives

    D=61600,
    c=79325621760000,
    c/D^3=57489/169400 != 0.

Thus c is nonzero in F. All other homogenized generators have normal order
at least three. After eliminating v, H generates (h^3), and every remaining
generator belongs to (h^3). Consequently the completed local algebra is
EXACTLY F[[h]]/(h^3). This proves generic length three without using the
global dimension upper bound first. The previous lower-bound argument was
valid, but this separates the local multiplicity proof more cleanly.

For each sextic component the auditor evaluates all 38 equations and their
full 38-by-8 Jacobian at a rational point, in coordinate order
(x,y,w,a,A,B,C,h):

| Component | Point | Jacobian rank |
| --- | --- | --- |
| X_A | (1,-6,2,3,0,16,-9,0) | 4 |
| X_B | (1,-2/3,2,3,16,0,-5,0) | 4 |
| X_C | (1,2,3,-6,16,-9,0,0) | 4 |

The contained integral sextic has dimension four. Rank four bounds the
ambient scheme's tangent dimension at the point by four; its local dimension
is at least four. Hence the actual homogenized scheme is regular there, and
the sextic component has generic multiplicity one. These witness checks do
not replace the original irreducibility proof: that proof uses the irreducible
quadric p=0, the one-dimensional locus p=q=r=0, and the line bundle on its
complement. Complete-intersection purity excludes an extra component on the
complementary locus, and generic reducedness plus unmixedness gives integrality.

The [machine-readable report](e8-degree-transfer-audit.json) records the exact
checks and their scope. The script does not report the geometric arguments
of the next sections as mechanically proved.

## 3. Audit of the degree budget

The integer homogeneous input is essential. In every ordinary degree, reduction
modulo a prime can only lower the rank of its coefficient matrix. Therefore
the characteristic-zero quotient Hilbert function is bounded ABOVE by its
finite-field counterpart. Leading monomials of actual finite-field members
give another upper bound. They need not form a complete Groebner basis.
The weighted monomial ordering causes no conflict: all polynomials are
homogeneous in the ordinary grading used for this comparison.

Subject to the recorded ideal memberships, the monomial bound is dimension
four and degree 21. The four distinct integral infinity components have
dimensions four, degrees 1,6,6,6 and multiplicities 3,1,1,1. Hence their total
contribution is exactly 21. Positivity of top-dimensional contributions
excludes every additional dimension-four component, in particular one meeting
h!=0. The underlying multiplicities are lengths at generic points, as in
[Stacks, cycles of coherent sheaves](https://stacks.math.columbia.edu/tag/02QV).

The logical order matters: containment of these components gives dim V(J)>=4;
the Hilbert-function upper bound gives dim V(J)<=4. Only after equality of
dimensions is known does comparison of leading Hilbert coefficients give
deg V(J)<=21. An upper bound by a HIGHER-dimensional Hilbert function alone
would not bound degree. The original argument has the required containment,
so this is an ordering clarification, not a missing hypothesis.

Localizing the individually homogenized ideal at h gives the original affine
ideal times the Laurent polynomial ring in h. Consequently

    dim V(J)[h!=0] <= 3,
    dim V(I) <= 2.

The second inequality has a dimension shift of ONE from the localized cone;
the total shift of TWO is from the original upper bound four. No projective
and affine dimensions are being interchanged. Saturating J at h before
counting the four infinity components would destroy this particular proof.

Embedded components cannot cancel degree. A lower-dimensional embedded or
nonreduced locus off infinity is still part of V(I) and is removed by its
bad target image. The argument does not infer reducedness merely from generic
reducedness at the minimal primes. Nor does it infer a characteristic-zero
dimension bound from an affine finite-field computation alone.

## 4. Transfer to the actual group fiber

Write Y_tr for the actual four-equation trace model. The full-rank locus of
its 4-by-7 Jacobian is smooth of dimension three. This follows directly from
the Jacobian criterion, without assuming height or a regular sequence in
advance. V(I) contains every point where that sufficient criterion fails.

First pull Y_tr and V(I) back to the original central torus. This map is
finite: the exponent substitution has determinant three, each one-variable
trace has a monic quadratic equation, and the simultaneous-inversion trace
map is finite. Thus the pulled-back bad locus still has dimension at most
two. Only on this pullback do we use the ORIGINAL logarithmic target map
mu to the three-dimensional central target. Its bad-image closure is proper.
No target map on Y_tr is asserted.

The finite map is etale after removing its stated ramification loci:

- For a=s+s^-1, at s=+/-1 the original target coordinate is
  (s-s^-1)f_a=0. Its image lies in a proper target hyperplane.
- For (U,V) modulo simultaneous inversion, the only fixed points are
  U,V=+/-1. At each, both logarithmic trace derivatives vanish, so both
  corresponding target coordinates are zero. This is a proper target line.

The full potential, including its central summand, and the chain rule after
the exact unit elimination are required here. Invertible linear changes of
central target coordinates preserve properness of these subsets. No derived
root discriminant is removed. The finite group quotient on the complement
of its fixed points is etale, so the pulled-back good incidence is smooth.

Remove the closures of these bad target images and the images of nondominant
components. Apply characteristic-zero generic smoothness to the remaining
smooth incidence. This yields a dense target open with smooth zero-dimensional
transverse fibers, possibly empty. Properness is unnecessary; the applicable
statement is [generic smoothness in the target, Theorem 3.3](https://math.stanford.edu/~vakil/0506-216/216class5152.pdf).

Now let h in z(k) satisfy Z_G(h)=K. The identity

    (Ad(g)-1)h = -(Ad(g)-1)(theta(g)-h)

puts the stabilizer equations inside the fiber ideal. Thus the group fiber
lies in K as a scheme. On K, the nondegenerate restricted Killing form and
the Maurer--Cartan frame identify its defining ideal with the coefficient
ideal of d chi_Ad minus the logarithmic central target form. This gives the
original critical equations, not their radical.

On the finite central cover S x Z(K)^0, the regular adjoint quotient is smooth.
The derived group S is simply connected, and Steinberg's independence of
fundamental-character differentials applies to regular elements, including
nonsemisimple ones; see [Steinberg, Theorems 6.1 and 8.1](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf).
Pulling back the transverse fiber along this smooth quotient preserves
smoothness. At a fiber point, equivariance gives Z_G(g) contained in K,
so Z_G(g)=Z_K(g): the regular loci agree ON THE FIBER. The finite central
cover is etale and smoothness descends. The relative quotient dimension is
dim K-rank K=12 in this case.

Every step either preserves the actual ideal or removes a proper subset of
the TARGET. Intersecting the resulting dense target open with Z_K^o retains
a dense open. These facts justify the passage to actual group fibers, without
discarding nilpotents by definition or silently deleting nonsemisimple points.

## 5. An immediate extension, and limits

**Central-isogeny corollary (formal deduction).** The completed theorem holds
for every connected complex algebraic group with simple Lie algebra, including
the adjoint form and intermediate central quotients, with the same Killing-form
normalization and the corresponding proper Levi central strata.

Indeed, let pi:G_sc -> G be a finite central isogeny and identify Lie algebras
by d pi. The adjoint representations and Killing forms agree, hence

    theta_G composed with pi = d pi composed with theta_G_sc.

The map pi is finite etale over C. Its inverse image of the regular locus
is the regular locus, since the Lie centralizer dimensions agree. The inverse
image of a Levi centralizer is the corresponding Levi centralizer upstairs.
For every h the regular fiber upstairs is therefore the finite etale
surjective pullback of the regular fiber downstairs. Smoothness is equivalent
by etale descent. This uses standard isogeny facts and the completed theorem;
it is not a new exceptional calculation or a claim of novelty.

A product statement follows when EVERY simple-factor Levi is proper, by
taking products of the good target opens and then central-isogeny descent.
A merely proper Levi in a semisimple product may contain an entire simple
factor. Its target in that factor is zero, which the present theorem excludes.
Thus a blanket all-proper-Levi semisimple extension needs a separate zero-fiber
argument. For a reductive group with positive-dimensional center, the Killing
form is degenerate; the adjoint normalization used here does not define the
desired full Lie-algebra-valued map without an additional choice.

## 6. Recommended direction

Stop the Levi completion sequence: the requested smoothness question is
settled by the recorded case argument, with its stated computational
dependencies. If stronger certification is desired, a bounded next task is
to export explicit multiplier identities for the selected finite-field
members, checked by a small independent polynomial verifier. This strengthens
reviewability without claiming a new mathematical theorem.

For new mathematics, return to OPEN-14 in the [research ledger](RESEARCH_LEDGER.md):
equivariant integral images in the classical families. In the symplectic
projective family, the ordinary Smith formula and the balanced formulas for
all Fitting ideals are already proved. The next question is what happens
after the actual equivariant moment substitution and then augmentation:
determine the kernel of base change from the embedded equivariant image to
the ordinary image, and the resulting primary arithmetic, especially at two.
This targets a family theorem connecting the extensions back to the original
Springer/Borel images, rather than another list of Levi computations.

Fitting ideals of a finitely presented cokernel commute with arbitrary base
change ([Stacks, Lemma 15.8.4](https://stacks.math.columbia.edu/tag/07Z6));
that standard fact is not the proposed new theorem. The unresolved
part is the specialized embedded image, its possible Tor contribution, and
the primary decomposition, which need not survive substitution unchanged.
More precisely, put R=R(T), M=K_T^0(G/P), J=the specified embedded image,
and C=M/J. Since M is free over R, the standard Tor sequence already gives

    0 -> Tor_1^R(C,Z) -> J tensor_R Z -> image(J tensor_R Z -> M tensor_R Z) -> 0.

The map R -> Z here is augmentation. The research target is to COMPUTE
this Tor module and the primary arithmetic uniformly for the existing
symplectic family, and identify the right-hand image with the correctly
normalized ordinary image. Writing the exact sequence alone is not progress
on that target. The ordinary Smith theorem should be a consistency check,
not recomputed as a new goal.

If the program stays exclusively with generalized Cayley maps, the more
substantial next question is a structural nonemptiness/excess criterion for
proper Levi targets. Smoothness alone permits emptiness, and the E6 contact
example already rules out a universal zero-excess assertion. Such a project
should begin with a proposed structural criterion and stopping rule, not
another unbounded exceptional census.
