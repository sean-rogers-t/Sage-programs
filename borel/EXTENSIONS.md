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

The exceptional catalogue and generic integral matrix method are already
available. Extend their coverage with explicit integral input forms and
independent small-degree tables, starting with G2 and F4 before larger E
cases. Record which generators are integral, which denominators remain,
and which primes occur in the resulting saturation quotient.

A rational Cartan polynomial can have integral Schubert coefficients, while
another can fail integrality. Scaling to an integer image changes the input
lattice and must be explicit. A general exceptional classification or a
restriction to powers of two is not assumed. See the existing
[catalogue](EXCEPTIONAL.md) and [normalization audit](INTEGRAL_LITERATURE.md).

## Representation and character-lattice changes

Compare explicitly specified normalized, SO-character, and Spin-character
sources, or sources arising from other representations. For the currently
specified SO form, uniform degree scaling is proved and implemented.
Additional weights or representations require a separately justified source
and integrality analysis. A change of root/weight coordinates alone does
not certify a statement about a different representation-ring lattice.

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

1. Package the proved B/C/D joint rank and primitive quotient formulas;
   then investigate closed ambient B/D joint Smith multiplicities.
2. Add reusable classical multi-block input helpers.
3. Exceptional integral tables with explicit normalizations.
4. Equivariant computations after choosing the coefficient-ring API.

The first two steps most directly build on the existing code and the
single-factor B/D deductions. Equivariant and quantum work should be
separate projects with their own mathematical specifications and checks.
