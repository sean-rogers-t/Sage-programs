# E6 relative Smith factors and representation comparisons

Version 0.22 adds E6 to `exceptional_root_weight_comparison`, covering all
63 proper standard parabolics and G/G. It packages the existing
[60-case theorem](E6_ROOT_WEIGHT_COMPARISON.md) and
[three-case completion](E6_REMAINING_CASES.md); it adds no image theorem.

## Maps and integral sources

Let G be simply connected E6, X=X(T) its weight lattice, Q its root
lattice, and I the retained Bourbaki Levi nodes. Define

\[
 M_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dQ)^{W_I}),\qquad
 N_I^d=\beta^{P_I}((\operatorname{Sym}_{\mathbb Z}^dX)^{W_I}).
\]

These are full invariant-source images in integral cohomology. The
computer-assisted theorem determines

\[
 N_I^d/M_I^d\cong
 (\mathbb Z/3)^{m_{I,1,d}}\oplus
 (\mathbb Z/9)^{m_{I,2,d}}\oplus
 (\mathbb Z/27)^{m_{I,3,d}}
\]

in every degree. The local proofs at three determine the global relative
group, because Q=X at every other prime. This relative quotient does not
determine the absolute quotient H^(2d)(G/P_I;Z)/N_I^d.

For the normalized trace Springer morphisms, the full literal character
sources are Sym_Z(6Q) for V(omega1), dimension 27, and Sym_Z(24X) for
the adjoint V(omega2), dimension 78. Under the papers' maps
xi_lambda^P=beta^P after inverse Springer pullback, their images are
6^d M_I^d and 24^d N_I^d. The normalization is the
[proved representation comparison](E6_REPRESENTATION_COMPARISON.md).

## Usage and conventions

```python
from sage_borel import exceptional_root_weight_comparison

# Omit node 2; the retained Levi is A5.
c = exceptional_root_weight_comparison("E6", retained_nodes=(1,3,4,5,6))
assert c.dimension == 21 and c.schubert_rank == 72
row = c.degree(3)  # cohomology H^6
assert row.relative_smith.multiplicities == {1: 1, 27: 1}
assert row.nonunit_multiplicities == {27: 1}
assert row.relative_smith.torsion_factors() == (27,)
assert c.smith_multiplicity_polynomials[27][3] == 1
assert row.as_dict()["cartan_type"] == "E6"
```

Retained nodes generate W_I: () means G/B, all six nodes mean G/G.
Polynomial degree d means cohomological degree 2d. Node order does not
matter; duplicates and invalid nodes/degrees are rejected. Above the flag
dimension, the image rank and relative quotient are zero. The returned
objects and coefficient tuples are immutable; dictionary views are fresh.

The common E6/E7 properties are:

| Property | Meaning |
| --- | --- |
| `cohomology_ranks[d]` | Rank of H^(2d), also of M_I^d and N_I^d |
| `smith_multiplicity_polynomials` | Factor-to-coefficient tuple dictionary: factors 3,9,27 for E6 and 2 for E7 |
| `multiplicity_coefficients[d]` | Total number of nonunit factors in degree d; retains its previous E7 meaning |
| `degree(d).nonunit_multiplicities` | Nonzero factor-to-count dictionary in one degree |
| `degree(d).relative_smith` | Compressed Smith inclusion M -> N, with unit slots retained |
| `proof_family`, `proof_references`, `certificate_sha256` | Provenance for the finite classification |

The older `multiplicity_2` degree field remains for compatibility and is
zero for E6; use `nonunit_multiplicities` for group-independent code.
`as_dict()` is JSON-safe. Factor lists and integer indices are expanded
only by explicitly calling budgeted methods, as described in the
[E7/shared API guide](E7_RELATIVE_IMAGE_API.md).

## Scalar and representation comparisons

In a simultaneous relative Smith basis of N, M has diagonal entries
s_i in {1,3,9,27}. For the degree-d images of Sym(aQ) and Sym(bX), use
A_i=a^d s_i and B_i=b^d. Sum and intersection entries are gcd(A_i,B_i)
and lcm(A_i,B_i). This proves the four quotient formulas implemented by
`compare_scaled_sources`; it also handles scales sharing factors of three.
No actual simultaneous Schubert basis is claimed or returned.

```python
representations = c.compare_representations(degree=3)
assert representations.left_scale == 6
assert representations.right_scale == 24
assert representations.sum_mod_left.torsion_factors() == (27,)
assert representations.sum_mod_right.torsion_factors() == (64, 64)

# M has relative entries (1,27); compare M with 3^3 N.
shared = c.degree(3).compare_scaled_sources(root_scale=1, weight_scale=3)
assert shared.right_in_left and not shared.left_in_right
assert shared.sum_mod_right.torsion_factors() == (27,)
```

Named selectors are `minimal` and `adjoint`, in either order, including
comparison of a representation with itself. The four properties
`sum_mod_left`, `sum_mod_right`, `left_mod_intersection` and
`right_mod_intersection` describe degreewise subgroup quotients. Their
sums do not include mixed products from other degrees. Embedded image
membership/preimages still require the [character image API](MINUSCULE_CHARACTER_API.md).

## Verification and exact scope

The [canonical manifest](experiments/e6-all-parabolics-summary.json) binds
the previous proofs and independent audits; the
[packager](tools/package_e6_relative_images.py) copies all 63 records into
the installed wheel. The [new tests](tests/test_e6_relative_images.py)
independently reconstruct cohomology polynomials with Sage Levi root data,
the 60 uniform formulas, the closed A5 formula and the common-basis
extension of the audited A2+A2+A1 base. That base coefficient list remains
a shared, explicitly hashed mathematical input.

Fresh checks compare actual minimal/adjoint E6/P1 image lattices in degrees
0,1,2,3,5,8. A separate fresh A5 degree-three calculation constructs both
full invariant sources as integer reflection kernels and obtains Z/27
after BGG. Transformed integer lattices test all four scalar quotients,
including shared three-primary scales. These validate the implementation;
they do not replace the original source-completeness proofs or establish
literature priority. See the [API audit](relative-image-verification.json)
and [current release binding](release-0.22-verification.json).
