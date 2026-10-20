# Extensions after the orthogonal formulas

The full-paper project is on hold. The priorities below concern useful library
features and mathematical experiments. Proposed extensions are not claims of
new theorems or novelty; each needs its own proof and literature check.

The [research ledger](RESEARCH_LEDGER.md) is the maintained index of results,
partial results, counterexamples and open questions. The user's priority is
generalizable statements. New exceptional cases should test a stated general
question rather than become the default endpoint of the project.

## Other cohomology theories and groups

The [research proposal](../articles/research-program-proposal/README.md) now
includes extension track E. The following projects are adopted research
directions, not implemented backends or new general theorems.

| Project | First target | Dependency or limitation |
| --- | --- | --- |
| E1: equivariant cohomology | Integral module presentation and proved ordinary-specialization comparison | Fix the source and integrality conventions in Aim A first; account for base-change kernels |
| E2: K-theory | Restricted associated-bundle image and a filtered comparison with ordinary cohomology | The literal source lies in R(L); Springer polynomial grading need not equal the augmentation/gamma filtration |
| E3: oriented cohomology | A compatible source and specialization theorem, starting with connective K-theory | Use a formal group algebra and explicit coefficients/completion; additive generation is not a theorem for every formal group law |
| E4: quantum cohomology | A multiplicative lift with q=0 recovery, then its kernel/image | Define the source and generator lifts; ordinary quotient relations generally deform |
| E5: affine Grassmannians | Transport an E4 image through the localized Peterson correspondence | The target uses Pontryagin homology; record quotient and localization losses |
| E6: affine/Kac–Moody flags | Degree-bounded nil-Hecke calculations and a separate definition/obstruction for a Springer analogue | Published flag-cohomology operators extend; ordinary finite-dimensional faithful trace projection does not |

Published starting points include
[equivariant oriented characteristic maps](https://arxiv.org/abs/1409.7111),
[equivariant K-Schubert calculus](https://arxiv.org/abs/1101.5911),
[quantum Borel/Toda presentations](https://arxiv.org/abs/alg-geom/9607001),
[Lam–Shimozono's Peterson correspondence](https://arxiv.org/abs/0705.1386),
[Mare–Mihalcea's specified affine quantum deformations](https://arxiv.org/abs/1409.3587)
and [Kostant–Kumar's nil-Hecke construction](https://doi.org/10.1073/pnas.83.6.1543).
The prospective contribution is the representation-dependent image theory,
not the ambient frameworks. Exact literature comparisons remain necessary.

Activate E1 first, followed by E2/E3 and then E4/E5. E6's characteristic-map
calculations can proceed independently of its foundational Springer question.
The ledger preserves OPEN-14/15 and adds OPEN-17 through OPEN-20 for the
new questions. Other isogeny forms and reductive groups remain OPEN-10.

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

The [geometry project](SPRINGER_GEOMETRY.md) now separately computes complete
G2 minimal/adjoint quotient maps and full critical determinants, with generic
degrees 6 and 18 proved and independently checked. E6/E7 minimal/adjoint
quadratic quotient coordinates are checked by Laurent substitution and
tensor Casimirs. The [complete E6 extension](E6_GEOMETRY.md) now gives all
higher minimal/adjoint coordinates and both full determinant products, with
independent matrices and exact unisolvent certificates. The
[E7 extension](E7_GEOMETRY.md) supplies all seven coordinates through degree
18 and both full 133-dimensional products, with independent matrices and
exact certificates. The [generic-degree extension](EXCEPTIONAL_DEGREES.md)
gives E6 degrees 1080/2754 and E7 degrees 7272/33664 for minimal/adjoint maps
on simply connected groups, using the full toric boundary correction.
Target branch loci and minimum Cayley degrees remain further work.

**Completed case studies.** The [first complete G2/F4 tables](EXCEPTIONAL_INTEGRAL.md)
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
restrictions and multiply equivariant Schubert classes. Ambient classes
specialize to ordinary classes when the parameters vanish, but an image
module's faithful specialization requires its own comparison, including
possible kernel/base-change effects. Published inputs include
[double theta polynomials](https://math.umd.edu/~harryt/papers/dbltheta.pdf)
and [double eta polynomials](https://math.umd.edu/~harryt/papers/dbleta.pdf).

This is a substantial backend extension, not an unimplemented proof of
those published formulas. For integral image questions, the coefficient
ring becomes a polynomial ring over Z, which is generally not a PID.
Ordinary Smith normal form therefore does not apply unchanged; Fitting
ideals, carefully specified specializations, or finite graded Z-lattices
would be appropriate objects to investigate.

## Quantum products and affine transport

Quantum multiplication adds curve-degree parameters and Gromov-Witten
corrections to the ordinary cup product. A focused implementation could
start with a Grassmannian family and published
[quantum Pieri rules](https://arxiv.org/abs/0809.4966).
The ordinary Springer-image formulas do not automatically define a quantum
Springer map or prove corresponding integral image results.
Track E4 now specifies a source/lift/specialization milestone before those
image claims. Once a lift exists, E5 asks for its image under the published
localized Peterson correspondence, with the quotient and inverted classes
explicit. This route does not require first inventing a trace Springer
morphism on an infinite-dimensional group.

## Generalization priorities

1. A uniform generic-degree method that accounts for representation lattices,
   degenerate faces and boundary multiplicities; use additional exceptional
   cases to test its hypotheses.
2. Structural criteria for integral Borel-image saturation, prime support and
   exponents, starting from the general source-lattice theorem.
3. A precise integral bridge from Hiller-type chains to the Smith filtration,
   extending the existing incidence and primitive-prefix mechanisms if possible.
4. Representation comparisons beyond the coprime-index lemma, including direct
   sums and overlapping weight spaces.
5. Uniform critical-divisor and harmonic-localization statements, followed by
   target branch/nonproperness calculations where the source data suffice.

The ledger also retains the many-block joint-image questions, closed ambient
B/D joint multiplicities and further exceptional parabolics. Track E gives
equivariant, K-theoretic, oriented, quantum and affine/Kac–Moody extensions
their own ordered milestones and missing definitions. The existing degree-one source and irreducible
lattice classifications are completed general results, not questions to
reopen by merely surveying more representations.
