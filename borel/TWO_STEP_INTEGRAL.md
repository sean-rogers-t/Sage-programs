# First integral joint-image calculations on two-step flags

These are exact, bounded calculations for fixed normalized inputs. They
extend the earlier rational two-step presentation and single-tail integral
calculations. They do not constitute an arbitrary-rank theorem or a claim
of originality. The full-paper project remains on hold.

Subsequent development: the [type C proof](TYPE_C_JOINT_PRIMITIVITY.md)
now establishes joint saturation, primitivity of the separate-image sum,
and the constant-only intersection for all positive block sizes. The
original finite tables below are retained as evidence and examples.
The subsequent [orthogonal joint-image proof](ORTHOGONAL_JOINT_IMAGES.md)
proves the primitive-sum and constant-intersection assertions in B and D,
with integral kernel and rank formulas. Ambient orthogonal joint saturation
is false in general; the proof excludes odd-prime factors, while closed
multiplicities of the powers of 2 remain open here.

## Groups, flags and input lattices

Write n=a+b+k with a,b,k positive. The two omitted nodes are a and a+b;
the flag consists of isotropic subspaces of dimensions a and a+b. The
classical Levi factors are GL(a), GL(b), and the remaining classical
tail of rank k. For type D the experiment uses k=2, so the nodes avoid
the end-node ambiguity; the residual root system is D2.

Let x_1,...,x_n denote the standard Cartesian coordinates used in the
library's normalized Springer/Borel formulas. Fix the middle inputs

$$c_j=e_j(x_{a+1},\ldots,x_{a+b}),\quad 1\leq j\leq b.$$

For B and C fix the tail inputs

$$e_i=e_i(x_{a+b+1}^2,\ldots,x_n^2),\quad 1\leq i\leq k.$$

For D use e_1,...,e_(k-1) and q=product_(r=a+b+1)^n x_r; q^2 is the
omitted last elementary function. The polynomial degrees are j, 2i,
and k, respectively. The coefficient ring for products of these fixed
inputs is Z. These are not silently rescaled literal SO-character inputs.

Write L_d, M_d, J_d for the images of the middle, tail and joint
Z-algebras in H^(2d)(G/P; Z), and S_d=L_d+M_d. The computation records
their ranks and Smith lists, the intersection, both containments, and
J_d/S_d, S_d/L_d and S_d/M_d. See the
[comparison API and algorithms](INTEGRAL_COMPARISONS.md).

## Reproducible range

| Type | (a,b,k) | Degrees tested | Joint saturation indices observed |
| --- | --- | --- | --- |
| C3 | (1,1,1) | 0–9 | 1 |
| C4 | (1,1,2) | 0–12 | 1 |
| C4 | (1,2,1) | 0–14 | 1 |
| C4 | (2,1,1) | 0–12 | 1 |
| C5 | (1,2,2) | 0–8 | 1 |
| B3 | (1,1,1) | 0–9 | 1, 2, 4 |
| B4 | (1,1,2) | 0–12 | 1, 4, 8 |
| D4 | (1,1,2) | 0–12 | 1, 4, 8, 32 |

There are 96 degree rows. An index of 1 includes zero images; the JSON
records ranks so that vanishing is not mistaken for full cohomological
surjectivity. Every type C joint rank agrees with the established rational
Hilbert-function formula from the
[two-step derivation](experiments/TWO_STEP_DERIVATION.md). That agreement
checks rational dimensions; it cannot by itself establish integral saturation.

In every tested positive degree, L_d intersect M_d=0. In degree zero
both are Z times the identity class. Every tested J_d/S_d is torsion-free,
although it often has positive free rank. Thus the separate-image sum is
primitive *inside the joint image* in these examples. The orthogonal joint
image itself can still fail to be primitive in integral cohomology.
These distinct statements are both visible in the next table.

## Representative results

| Type and blocks | d | rank L_d | rank M_d | rank J_d | Joint Smith factors | J_d/S_d |
| --- | --- | --- | --- | --- | --- | --- |
| C3, (1,1,1) | 3 | 1 | 0 | 2 | (1,1) | Z |
| C4, (1,2,1) | 4 | 3 | 1 | 5 | (1,1,1,1,1) | Z |
| C5, (1,2,2) | 6 | 4 | 2 | 9 | nine ones | Z^3 |
| B3, (1,1,1) | 3 | 1 | 0 | 2 | (1,2) | Z |
| B3, (1,1,1) | 7 | 0 | 0 | 1 | (4) | Z |
| B4, (1,1,2) | 6 | 1 | 1 | 3 | (2,2,2) | Z |
| D4, (1,1,2) | 4 | 1 | 3 | 5 | (1,1,2,2,2) | Z |
| D4, (1,1,2) | 8 | 0 | 1 | 3 | (2,4,4) | Z^2 |

For C3 in degree 3 the new direction is beta(c_1 e_1). Neither the
middle algebra nor the tail algebra alone supplies it. In B3 the same
degrees and inputs produce a joint rank-two lattice with a factor 2;
in degree 7 both separate images vanish while a mixed image survives
with factor 4. Rational kernel computations alone cannot detect this.

## Connection with the papers and previous computations

The Borel characteristic map and chosen invariant inputs remain those
used throughout the library's Springer/Borel calculations. Parts II and
III supply the surrounding characteristic-map framework and explicit
generator computations; these tables investigate additional integral
information about products of two chosen Levi-factor input sets.
They do not change the Springer morphism itself.

The earlier B/D formulas describe the specified single-tail sources.
One cannot obtain the present joint Smith lists merely by concatenating
the separate lists: the relative positions of images and the mixed
products matter. The new API makes that distinction exact. The finite
tables motivated the subsequent type C theorem linked above. General
type B/D joint saturation is false, while the later orthogonal proof now
establishes primitivity of the sum inside the joint image. Its remaining
ambient Smith multiplicities are a separate problem. The proof notes
identify their classical inputs separately from the deductions and make
no novelty claim.

## Verification and rerunning

Run inside SageMath, from the `borel` directory:

```sh
PYTHONPATH=src sage -python -m unittest tests.test_lattices -v
PYTHONPATH=src sage -python experiments/two_step_integral_images.py
```

The [experiment](experiments/two_step_integral_images.py) writes
[the complete exact tables](two-step-integral-results.json), including
source hashes and runtime. It verifies all containments and the
rank/intersection identity, and compares the two quotient descriptions
from the second isomorphism theorem.

For each rank-three/four case it independently reconstructs the middle,
tail and joint matrices in degrees 0–4 using Cartesian divided differences,
rather than the production fundamental-weight sparse implementation.
It additionally checks B3 degree 7 and D4 degree 8, including their
factors 4. All 111 independently reconstructed matrices have their Smith
lists checked by gcds of minors. The C5 calculations use the production
matrix engine plus the rational Hilbert check; no independent C5 matrix
verification is claimed.

The comparison tests separately use explicit divisibility conditions,
Bezout and least-common-multiple examples, gcds of minors, and determinant
index identities. They check 1,296 integer points across 16 pairs of skew
rank-two lattices, along with proper subspaces, zero images, nontrivial
free and torsion quotients, and context validation.
