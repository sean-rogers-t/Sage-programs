# Complete integral Springer character images on G2/B

This finite calculation completes the full-flag case alongside the
[two maximal-parabolic images](G2_CHARACTER_FORM.md). Use the complex
group G2, Bourbaki fundamental weights omega1,omega2, and normalized
trace Springer morphisms with identity differential at 1. Polynomial
degree d means cohomology H^(2d). Integral characters mean virtual
characters with integer coefficients.

## The statement in the papers' notation

For the seven-dimensional V(omega1), the earlier source theorem proves

\[
 \mathcal C_7=\mathbb Z[6\omega_1,2\omega_2].
\]

For B the Levi is T, so no invariant-taking remains. Set
x=Theta_7^*(6omega1), y=Theta_7^*(2omega2) on the character side, and
X=xi_(omega1)^B(x)=beta(6omega1), Y=xi_(omega1)^B(y)=beta(2omega2)
in cohomology. The full literal image has the integral presentation

\[
 \boxed{J_{7,T}=\operatorname{im}(\xi_{\omega_1}^{B}|_{\mathbb Z[x,y]})
 \cong\mathbb Z[X,Y]/(X^2-3XY+3Y^2,\;Y^6).}
\]

Both X and Y have polynomial degree one. An integral basis is
Y^j, X Y^j for 0<=j<=5. The ring has rank twelve; it is a proper
full-rank subring of the integral cohomology ring. The complete nonzero
Smith factors in the integral Schubert basis are:

| Polynomial degree d | Cohomology degree | Rank | Minimal-representation factors |
| --- | --- | --- | --- |
| 0 | 0 | 1 | 1 |
| 1 | 2 | 2 | 2, 6 |
| 2 | 4 | 2 | 12, 12 |
| 3 | 6 | 2 | 24, 144 |
| 4 | 8 | 2 | 288, 288 |
| 5 | 10 | 2 | 576, 1728 |
| 6 | 12 | 1 | 3456 |

All higher degree images vanish. These factors describe the absolute
quotients H^(2d)(G2/B;Z)/J_(7,T)^d, not just a root/weight comparison.

For the fourteen-dimensional adjoint V(omega2), the proved full source is
Z[24omega1,8omega2]. Thus J_(14,T)^d=4^d J_(7,T)^d. The same displayed
integral presentation holds in the adjoint generators, and its factors
are 4^d times the minimal list. The source classification extends this
to every nontrivial irreducible G2 representation by the already-proved
Dynkin-index scaling law, with exact rational scale I_lambda/2.

## Why the kernel is complete over Z

The quadratic is twelve times the W-invariant
3omega1^2-3omega1 omega2+omega2^2, so Borel kills it. The weight omega2
descends to G2/P2, which has complex dimension five; its sixth power
therefore vanishes on G2/B as well. Both relations are also checked by
the two computations below.

With lexicographic order X>Y the relations are monic with coprime
leading monomials X^2,Y^6. Polynomial division over Z gives the free
rank-twelve quotient with the displayed basis. Every homogeneous block
of those twelve basis images is rationally independent, as verified
by its integer image matrix. Hence the map from this free quotient to
integral cohomology is injective: a possible extra integer kernel would
remain a rational kernel. This establishes the full integral kernel in
all degrees, rather than inferring it from a bounded rational search.

The full torus source is the previously proved input. The characteristic
map and divided-difference framework are published inputs from Kumar
Part I and standard Schubert calculus; this calculation and the kernel
deduction are recorded here without a literature-priority claim.

## API

```python
from sage_borel import g2_character_algebra

A = g2_character_algebra()  # None means B; 1 and 2 retain the earlier API.
assert A.integral_degree(3, include_ambient=True).smith_factors == (24, 144)
assert A.integral_degree(6).smith_factors == (3456,)
assert A.integral_degree(7).smith_factors == ()
x, y = A.domain.gens()
assert set(A.presentation().relations) == {x*x-3*x*y+3*y*y, y**6}
adj = g2_character_algebra(representation="adjoint")
assert adj.integral_degree(1).smith_factors == (8, 24)
assert g2_character_algebra(1).integral_degree(1).smith_factors == (6,)
```

The returned algebra also provides the existing embedded lattice,
membership, preimage and comparison methods. The `representation` keyword
accepts `minimal` and `adjoint` on B, P1 and P2. The default remains the
minimal representation; the original calls with node 1 or 2 are unchanged.

## Independent verification

The [primary calculation](experiments/g2_full_flag.py) uses library BGG
evaluation and Sage Smith arithmetic. Its
[certificate](experiments/g2-full-flag.json) records every source monomial
column through degree seven, all normal-basis images and both representations.

The [native audit](experiments/audit_g2_full_flag.py) imports neither Sage
nor runtime/primary code. It constructs the twelve Weyl matrices and six
positive roots from the rank-two Cartan datum, derives coroot pairings,
and reconstructs every column by root/coroot Chevalley multiplication.
It determines Smith factors from gcds of all one- and two-row minors.
All 72 monomial columns, the two relations, and the twelve normal-basis
images for each representation agree. See its
[passed report](experiments/g2-full-flag-audit.json).

The Cartan datum, characteristic-map theorem and earlier full-source
theorem are shared mathematical inputs. This is independent arithmetic
verification, not a separate proof of the earlier source theorem.

```sh
PYTHONPATH=src sage -python experiments/g2_full_flag.py
python experiments/audit_g2_full_flag.py
```
