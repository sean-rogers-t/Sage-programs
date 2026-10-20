# Extensions after the orthogonal formulas

The full-paper project is on hold. The priorities below concern useful library
features and mathematical experiments. Proposed extensions are not claims of
new theorems or novelty; each needs its own proof and literature check.

## Implemented: fast B/D multiplicities

[orthogonal_smith](ORTHOGONAL_API.md) exposes the proved closed formulas,
including separate square/Euler sectors and explicit character scaling.
The results are compressed, so a large image rank does not itself require
a large list of Smith factors. Explicit Schubert membership and preimages
continue to use [integral_degree](INTEGRAL_IMAGES.md).

## Implemented: compare and combine integral images

Version 0.9 adds degreewise comparison of selected generator algebras on the
same flag variety: containment, equality, sum, intersection and exact
quotient groups. Free parts and finite torsion are reported separately.
See the [API and mathematical algorithms](INTEGRAL_COMPARISONS.md).

These are integer-lattice operations on the existing exact
matrix reports. They are more informative than comparing Smith lists alone:
two sublattices can have the same Smith factors in different directions.
They expose the interactions needed for the next research question.

## Joint Levi factors and more flag steps

The generic `BorelMap` already accepts arbitrary standard parabolics, and
`GeneratorAlgebra` already analyzes chosen invariant inputs. The current
`two_step_symplectic(n,b,k)` helper builds the middle and tail factors of
`GL(a) x GL(b) x Sp(2k)`, with `a=n-b-k`, and their joint algebra.
The extension is to supply equally explicit classical multi-block input
helpers and systematic integral comparisons, rather than merely allowing
more omitted nodes (which is already supported).

The [first joint-image experiment](TWO_STEP_INTEGRAL.md) now compares the
middle-only, tail-only and joint images in five symplectic cases and three
orthogonal cases, covering 96 degree rows. The tested symplectic joint
images are saturated; the orthogonal ones can have factors 2 and 4. Every
tested separate-image sum is primitive inside the joint image. The
[general type C argument](TYPE_C_JOINT_PRIMITIVITY.md) now proves joint
saturation and primitivity of that sum, using an integral basis compatible
with both separate images. It also gives their constant-only intersection,
the integral kernel and the quotient's graded free ranks.
The [B/D joint-image proof](ORTHOGONAL_JOINT_IMAGES.md) now establishes
primitive separate sums and constant-only intersections in both orthogonal
families, with integral kernels and closed rank formulas. It excludes
odd-prime ambient Smith factors, but the joint images can have factors
2,4,... . Their closed multiplicities and more general multi-block
configurations remain further work. Rational kernel presentations alone
cannot establish integral primitivity; the D proof supplies an explicit
free integral decomposition.

Published ambient formulas exist for arbitrary classical partial flags:
[Tamvakis, Giambelli and degeneracy locus formulas for classical G/P spaces](https://math.umd.edu/~harryt/papers/giamdegcl.pdf).
Their existence does not automatically identify these chosen integral
Springer-image sublattices.

## Exceptional integral images

**Current focus.** The [first complete G2/F4 tables](EXCEPTIONAL_INTEGRAL.md)
now cover every maximal parabolic and every degree for two explicit source
forms, with the `exceptional_algebra` helper and independent matrix checks.
The [E6/P1 and E7/P7 extension](EXCEPTIONAL_MINUSCULE.md) now covers every
degree too, checked by fixed-point interpolation and published integral
presentations with ordinary Chevalley multiplication. Integer generator
changes identify their primitive-image algebras as Z[h,16s] and
Z[h,81s,19683t]. The [E6/P1 literal character form](E6_CHARACTER_FORM.md)
is now identified as Sym_Z(6Q)^(W_L), and its entire image is compared with
the two selected forms. The [E7/P7 character form](E7_CHARACTER_FORM.md)
is also identified, as Sym_Z(12Q)^(W_L) for the E7 root lattice, with exact
comparisons including incomparable images. Version 0.13 packages both
certified character-form calculations in the [public API](MINUSCULE_CHARACTER_API.md).
Further exceptional parabolics still require an independent ordinary-ring check.

Version 0.10 now extracts a [general source criterion](CHARACTER_SOURCES.md)
and implements normalized representation inputs, exact linear lattices,
all-primes sufficient certificates and bounded modular-obstruction searches.
The [G2 application](G2_CHARACTER_FORM.md) proves the complete seven-dimensional
character source, both maximal-Levi source and image presentations, and the
adjoint torus form. The [F4 case](F4_CHARACTER_FORM.md) now proves the
26-dimensional torus source Sym_Z(6 Q_long), computes all four maximal-Levi
images, and packages the degreewise character-source intersection in the API.
The E6/P1 and E7/P7 helpers now compute these images without loading saved tables.
The [E6 adjoint comparison](E6_REPRESENTATION_COMPARISON.md) proves the
full torus source Sym_Z(24P), all-Levi sum/intersection formulas, and
complete degreewise incomparability on E6/P1.

The exceptional catalogue and generic integral matrix method are already
available. The completed tables retain explicit integral input forms,
generator contents, every source-monomial matrix and the primes in the
saturation quotient. Further E cases need equally explicit normalization
and independent verification.

A rational Cartan polynomial can have integral Schubert coefficients, while
another can fail integrality. Scaling to an integer image changes the input
lattice and must be explicit. A general exceptional classification or a
restriction to powers of two is not assumed. See the existing
[catalogue](EXCEPTIONAL.md) and [normalization audit](INTEGRAL_LITERATURE.md).

## Representation and character-lattice changes

Compare explicitly specified normalized, SO-character, and Spin-character
sources, or sources arising from other representations. For the currently
specified SO form, uniform degree scaling is proved and implemented.
The [trace-source theorem](TRACE_SOURCE_THEOREM.md) now proves degree-one
generation for every supported normalized trace representation map.
The [irreducible classification](IRREDUCIBLE_CHARACTER_LATTICES.md) now
determines the full linear lattice for every nontrivial irreducible:
Lambda_lambda=I_lambda*(Q+Z*lambda)^*. The fast API computes it without
representation weights. Integral Levi invariants and cohomology images
still need calculation. Direct sums require their combined weighted
coefficient lattice, with a same-index F4 counterexample to the
irreducible shortcut recorded in the proof. A change
of root/weight coordinates alone does not identify another representation's
character lattice.

## Equivariant calculations

An equivariant backend would retain torus parameters, evaluate fixed-point
restrictions, multiply equivariant Schubert classes, and specialize those
parameters to zero to recover ordinary calculations. Published inputs include
[double theta polynomials](https://math.umd.edu/~harryt/papers/dbltheta.pdf)
and [double eta polynomials](https://math.umd.edu/~harryt/papers/dbleta.pdf).

This is a substantial backend extension, not an unimplemented proof of
those published formulas. For integral image questions, the coefficient
ring becomes a polynomial ring over Z, which is generally not a PID.
Ordinary Smith normal form therefore does not apply unchanged; Fitting
ideals, carefully specified specializations, or finite graded Z-lattices
would be appropriate objects to investigate.

## Later: quantum products

Quantum multiplication adds curve-degree parameters and Gromov-Witten
corrections to the ordinary cup product. A focused implementation could
start with a Grassmannian family and published
[quantum Pieri rules](https://arxiv.org/abs/0809.4966).
The ordinary Springer-image formulas do not automatically define a quantum
Springer map or prove corresponding integral image results.

## Suggested order

1. Completed: F4's 26-dimensional character source, independent lattice audit,
   and complete comparisons with the existing images. See
   [the F4 proof and API](F4_CHARACTER_FORM.md).
2. Completed: the [F4 minimal/adjoint comparison](F4_REPRESENTATION_COMPARISON.md),
   including all-degree 3^d scaling and the representation selector.
   Also completed: E6/P1 and E7/P7 public character helpers, and the
   [E6 minimal/adjoint comparison](E6_REPRESENTATION_COMPARISON.md).
   Also completed: the [E7 minimal/adjoint comparison](E7_REPRESENTATION_COMPARISON.md),
   with its own all-primes certificate, integral weight coordinates and
   complete P7 image audit. The E6/E7 comparisons illustrate a general
   coprime-index lemma for certified sources.
   Also completed: the [general trace-source theorem](TRACE_SOURCE_THEOREM.md)
   proves degree-one generation for every supported representation, and a
   separate Python audit reconstructs all 59 surveyed representation cases.
   Next: calculate new linear lattices and their integral Levi/Borel images,
   particularly representation comparisons where scalar and lattice-index
   coprimality fails, and extend coverage to further exceptional parabolics.
3. Equivariant computations after choosing the coefficient-ring API.

Further classical multi-block work, packaging the proved joint formulas and
closed ambient B/D joint multiplicities are deferred. Equivariant and quantum
work require separate mathematical specifications and checks.
