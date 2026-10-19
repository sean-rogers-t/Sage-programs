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

## Recommended next: compare and combine integral images

Add a degreewise comparison of two selected generator algebras on the same
flag variety. Useful outputs are containment, equality, their sum and
intersection, and the indices when rational spans agree. For unequal spans,
report ranks and the relevant finite torsion separately; do not manufacture
a finite full-ambient index.

These are implementable integer-lattice operations on the existing exact
matrix reports. They are more informative than comparing Smith lists alone:
two sublattices can have the same Smith factors in different directions.
They would expose the interactions needed for the next research question.

## Joint Levi factors and more flag steps

The generic `BorelMap` already accepts arbitrary standard parabolics, and
`GeneratorAlgebra` already analyzes chosen invariant inputs. The current
`two_step_symplectic(n,b,k)` helper builds the middle and tail factors of
`GL(a) x GL(b) x Sp(2k)`, with `a=n-b-k`, and their joint algebra.
The extension is to supply equally explicit classical multi-block input
helpers and systematic integral comparisons, rather than merely allowing
more omitted nodes (which is already supported).

The first proposed experiment is to compare the middle-only, tail-only,
and joint images in two-step symplectic flags; then do the analogous B/D
families. Determine the joint saturation quotient and compare it with the
known single-factor results. General primitive ranges, partition formulas,
and a description of cross-factor gluing are research questions. Rational
kernel presentations alone cannot establish them.

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

1. Integral lattice comparisons and sums/intersections.
2. Two-step joint-image experiments, then classical multi-block helpers.
3. Exceptional integral tables with explicit normalizations.
4. Equivariant computations after choosing the coefficient-ring API.

The first two steps most directly build on the existing code and the
single-factor B/D deductions. Equivariant and quantum work should be
separate projects with their own mathematical specifications and checks.
