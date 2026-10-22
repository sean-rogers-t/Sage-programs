# Classical equivariant characteristic images

Working API, 22 October 2026. The first backend treats the defining
representation of SL_m on every standard partial flag, in equivariant
ordinary cohomology or equivariant K-theory. The second treats the full
literal defining Sp_(2n)/P1 K-image and its integral cokernel.

These are finite classical groups. The Spin central-parity classification
remains shelved. No affine or general representation backend is asserted.

## 1. Coefficient rings and actual source

Use the same omitted-node convention as BorelMap. The retained nodes
generate the Levi Weyl group. Coordinates x_i are the epsilon weights,
with x_m=-sum_(i<m) x_i, rather than fundamental weights. Thus
omega_i=x_1+...+x_i.

The defining SL_m torus source is Sym_Z(Q_root), where Q_root has basis
x_i-x_m. A polynomial must be Levi invariant, and its actual pullback

\[
 x_i\longmapsto t_i-\frac1m\sum_jt_j,\qquad \prod_jt_j=1,
\]

must have integral Laurent coefficients. The default ZZ context checks
both conditions. Merely having an integral ordinary Borel image does
not pass the literal-source check.

Theory H has coefficient ring Z[u_1,...,u_(m-1)], with u_m=-sum u_i.
Theory K has coefficient ring Z[t_1^+-1,...,t_(m-1)^+-1], with
t_m=(prod t_i)^-1. Rational normalization is available only by explicitly
selecting QQ. There is no automatic denominator clearing.

## 2. Type A examples

```python
from sage.all import ZZ, QQ
from sage_borel import type_a_equivariant

h = type_a_equivariant(4, omitted_nodes=1, theory="H")
E1, E2, E3 = h.source_generators()
assert h.expand(E1).specialize() == {(2, 1, 3, 4): 4}
assert h.is_integral_source(E2)

k = type_a_equivariant(4, omitted_nodes=1, theory="K")
F1, F2, F3 = k.source_generators()
assert k.expand(F1 * F2) == k.expand(F1) * k.expand(F2)
assert k.expand(F1).restrictions() == k.restrictions(F1)
assert not k.is_integral_source(k.weights[0])

# General partial flag of steps 1 and 3, not a maximal parabolic.
p = type_a_equivariant(5, omitted_nodes=(1, 3), theory="K")
assert p.block_sizes == (1, 2, 2) and p.schubert_rank == 30
f = 5*p.weights[0]
assert p.is_integral_source(f)
p.expand(f)

# A normalized rational input must be requested explicitly.
r = type_a_equivariant(3, omitted_nodes=1, theory="K", coefficients=QQ)
r.expand(r.weights[0])
```

An expansion stores coefficients in the opposite Schubert basis.
For H, forgetting equivariance sets u_i=0. For K it sets t_i=1;
the resulting coordinates remain structure-sheaf Schubert coordinates,
not Chern-character coordinates. Products are reconstructed from
pointwise products of restrictions and checked for exact divisibility.

| API | Meaning |
| --- | --- |
| weights, cartan_ring | Rational epsilon-weight input coordinates |
| linear_lattice_matrix | Root basis x_i-x_m in the independent epsilon coordinates |
| source_generators() | Complete integral source generators on G/B and SL_m/P1 |
| is_integral_source(f), pullback(f) | Literal character membership and the actual Laurent pullback |
| restrictions(f), expand(f) | Fixed-point values and exact Schubert coordinates |
| restriction_matrix() | Rows fixed points, columns opposite Schubert classes |
| schubert(w) | Basis class indexed by a minimal permutation |
| expansion.specialize() | Ordinary H or K Schubert coordinates |

All standard parabolics support input membership and evaluation.
Automatic enumeration of a full integral invariant-ring generating
set is supplied only on G/B and SL_m/P1. Elsewhere source_generators()
raises NotImplementedError. A supplied list of invariant inputs is
not advertised as the full Levi source.

Likewise the type A backend does not yet compute a complete module
presentation of the coefficient-span image or its specialization
kernel. The known SL4 cohomological specialization defect remains
in [the pilot](EQUIVARIANT_IMAGE_PILOT.md). Evaluating its generators
is different from classifying that module.

## 3. Why the restriction calculation works

The minimal right cosets are permutations increasing inside each Levi
block. They are generated as ordered partitions of the values, not by
enumerating the full symmetric group. The default guard is 720 fixed
points, configurable with max_fixed_points. This bounds enumeration;
it is not a promise that every dense 720-class symbolic calculation
will be quick.

For a reduced word of w, take the sequence of prefix roots beta_j.
In our signs the additive factor is u_b-u_a and the K factor is
1-t_a/t_b for beta_j=epsilon_a-epsilon_b. This agrees with
beta(omega_i)=epsilon_si and the positive projective line restricting
to u_i or t_i.

Additive restrictions sum products over reduced subwords. K restrictions
sum over subwords with the specified 0-Hecke/Demazure product, weighting
each repeated, non-length-increasing step by a minus sign. If q_j is
the factor at the final step i and v has right descent i, the K recursion is

\[
 R(v,j)=(1-q_j)R(v,j-1)+q_jR(vs_i,j-1).
\]

Otherwise R(v,j)=R(v,j-1). The additive formula replaces the right
side by R(v,j-1)+q_jR(vs_i,j-1) in the descent case.
The resulting square restriction matrix is triangular by Bruhat order.
Forward substitution recovers Schubert coefficients. Computation in a
fraction field is temporary: each coefficient must coerce back into
the specified integral polynomial/Laurent ring.

These are applications of published restriction formulas:
[Billey](https://pmc.ncbi.nlm.nih.gov/articles/PMC19230/),
[Willems](https://www.numdam.org/item/BSMF_2004__132_4_569_0/) and
the [Griffeth-Ram characteristic framework](https://arxiv.org/abs/math/0405333).
The backend is not a new construction of ambient equivariant theories.

## 4. Symplectic projective K-images

```python
from sage_borel import symplectic_projective_k_image

c = symplectic_projective_k_image(3, spectral=True)
s1, s2, s3 = c.parameters
F = (s1+s2)*(s1+s3)*(s2+s3)
assert c.presentation_matrix().nrows() == 3
assert c.cyclic_cokernel_factors() == (2, 2*F, 2*F)
assert c.fitting_generators(1) == (4*F,)
assert c.fitting_generators(2) == (2,)
assert c.ordinary_smith_factors() == (1, 1, 1, 2, 128, 128)

torus = symplectic_projective_k_image(2)
assert torus.ordinary_smith_factors() == (1, 1, 8, 8)

c4 = symplectic_projective_k_image(4, spectral=True)
assert c4.half_presentation_matrix().nrows() == 2
assert c4.block_presentation_matrix().nrows() == 4
assert len(c4.fitting_generators(1, max_minors=1)) == 2
assert len(c4.fitting_generators(2, max_minors=1)) == 3
assert c4.ordinary_smith_factors() == (1, 1, 1, 1, 16, 16, 1024, 1024)

c6 = symplectic_projective_k_image(6, spectral=True)
# Three initial-column minors replace all nine size-two minors.
assert len(c6.half_minor_generators(2, max_minors=3)) == 3
assert len(c6.fitting_generators(2, max_minors=1)) == 6
assert c6._matrix is None
assert c6.half_hankel_presentation_matrix().is_symmetric()

c8 = symplectic_projective_k_image(8, spectral=True)
# Balancing needs four entry generators and ten unordered products.
assert len(c8.fitting_generators(6, max_minors=14)) == 10
```

The default coefficient ring is the actual torus ring Z[t_i^+-1].
Spectral mode uses Z[s_1,...,s_n], with s_i=t_i+t_i^-1 under base
change. The torus ring is finite free over this spectral ring, so
Fitting ideals base change exactly.

The image is R[b], b=L-L^-1, in the projective bundle K-ring.
Its monic relation is image_relation. The separate cokernel presentation
uses cokernel_polynomial and multiplication by twice odd; its matrix
has size n, rather than 2n. Those are module relations for the
class of L in the cokernel, not equations defining the image algebra.

The [even/odd proof](SYMPLECTIC_K_FITTING.md) now gives a further uniform
splitting over R. Write P(s)=E(s^2)+sO(s^2), n=2r+epsilon and
let H multiply O modulo E in even rank, or E modulo O in odd rank.
Both monic quotients have rank r. The entire cokernel is two copies
of coker(2H), plus one R/(2) in odd rank. This is an R-module
decomposition, not a splitting compatible with the source coordinate b.

| Type C method | Meaning |
| --- | --- |
| presentation_matrix() | Original reduced rank-n multiplication presentation |
| half_presentation_matrix() | Single r-by-r block 2H, occurring twice |
| half_hankel_presentation_matrix() | Integrally equivalent r-by-r Hankel block from the coefficient pairing |
| block_presentation_matrix() | Equivalent rank-n block diagonal matrix, with optional scalar 2 |
| half_minor_generators(j) | Exact ideal I_j(2H), using only the first j columns |
| fitting_generators(i) | Exact ideal generators; one balanced product of compressed half-block minor ideals |
| ordinary_smith_factors() | Full integer inclusion factors after augmentation, including n units |

For rank four put F=product_(i<j)(s_i+s_j), J=(e1(s),e3(s)).
The complete list is (16F^2), 8F J, 4J^2, 2J, then R. These are
closed ideal identities over the integral coefficient ring; no minors
are enumerated for these cases. The [Schur-minor proof](SYMPLECTIC_K_MINORS.md)
also gives explicit rank-five and cubic-block rank-six/seven formulas.
For n=5, set J=(e1*e4-e5,e1*e2-e3); the list is (32F^2),16FJ,8J^2,4J,(2),R.
For n=6,7, the companion note displays J and a three-generator ideal K
explicitly, and gives every Fitting ideal in terms of J,K,F.

In every rank, integral Schur cyclicity proves that only the first j
columns of the half block are needed for its size-j minor ideal.
Thus half_minor_generators(j) computes binomial(floor(n/2),j) determinants
instead of its square. The general fallback convolves these compressed
ideals. The subsequent [integral Hankel proof](SYMPLECTIC_K_BALANCING.md)
reduces the full convolution to one product in every rank. If k=n-i,
epsilon=n mod 2 and ell=k-epsilon, the ideal is
2^k I_floor(ell/2)(H) I_ceil(ell/2)(H). A coefficient pairing makes H
integrally equivalent to a Hankel matrix; published Hankel product
theorems and a monic integer lift prove the balancing, including at two.
The API enumerates unordered generator products for squares.
The max_minors guard bounds the required initial-column
minor and generator-product operations; it is conservative
and is not a count of distinct final generators.

These are Fitting ideals over a multivariable ring, not a global
Smith diagonalization. Global cyclic decompositions are proved only
for n<=3; cyclic_cokernel_factors() raises NotImplementedError for n>3.
Explicit closed ideal lists now extend through n=9, with six polynomial
generators for the quartic two-minor ideal and four for its three-minor
ideal, displayed in the Hankel proof. The all-rank generators are also the Schur coordinates of a product
of q at the roots of a universal factor of p, as proved in the companion
note. They need not be a minimal family. Further minimal-generator and
specialized primary arithmetic classifications remain open; the fallback
returns exact, possibly redundant generators. Generic Hankel primary
decompositions are not automatically primary decompositions after moment
substitution. Large spectral-polynomial expansions may still be costly.

ordinary_smith_factors() augments before constructing the half-size
integer matrix. It duplicates that block's factors and adds a single
2 when n is odd, then n units for the full inclusion. It never builds
either symbolic presentation. This matters in rank nine and above.

## 5. Checks and limits

The new tests compare projective restrictions with Newton bundle
polynomials, full-flag K divisors with line-bundle formulas, H
specializations with the existing independent BGG engine, and the
0-Hecke recursion with a direct nonreduced-subword sum. The C
presentation is checked against the original 2n projective inclusion
at regular and singular integer spectral points through rank eight and
at rank nine after ordinary augmentation. The audit verifies explicit
unimodular changes over independent elementary parameters through rank
eight and the complete rank-four minor ideals in independent coefficients.
Half-minor convolution checks against the original inclusion in ranks
five and six, symbolic determinants and torus base change also agree.
The new verifier supplies 94 explicit integral Schur matrix identities,
72 original-inclusion checks of all rank-five/six/seven Fitting formulas,
and integer and finite-field checks including characteristic two.
The Hankel verifier supplies 70 universal quartic integer polynomial
identities, 68 original-inclusion checks of every rank-eight/nine ideal,
and polynomial-fiber checks over characteristics two and three.

The proof note supplies the uniform mathematical argument; finite
matrix checks are corroboration. Full integral type A image modules,
general isotropic Grassmannians, B/D quadric backends and a closed
minimal or arithmetic classification for every larger-block Fitting
ideal remain work. The uniform compressed generator formula is proved.
