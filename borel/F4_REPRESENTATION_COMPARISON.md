# F4: comparing the minimal and adjoint Springer sources

Changing from the 26-dimensional representation V(omega4) to the
52-dimensional adjoint representation V(omega1) gives a precise degreewise
scaling law for the integral sources and Borel images:

\[
\boxed{\mathcal A_{52,L}^{d}=3^d\mathcal A_{26,L}^{d},
\qquad J_{52,L}^{d}=3^dJ_{26,L}^{d}.}
\]

This holds for every standard Levi L and every polynomial degree d.
Here polynomial degree d corresponds to H^(2d). The source algebras are
identified inside the **common polynomial ring** S_Q(t^*) using the
respective Springer pullbacks. The assertion does not identify the two
Laurent-character functions or multiply one Springer morphism by a constant.

## Definitions and published inputs

For rho in {26,52}, let theta_rho be the normalized trace-projection
Springer morphism, so d(theta_rho)_1=id, and set

\[
\mathcal C_\rho=\{f\in S_{\mathbb Q}(\mathfrak t^*):
                       \theta_\rho^*f\in\mathbb Z[X(T)]\},\qquad
\mathcal A_{\rho,L}=\mathcal C_\rho^{W_L},\qquad
J_{\rho,L}=\beta(\mathcal A_{\rho,L}).
\]

Beta is the Borel map on the same Cartan polynomials for both representations.
The actual virtual characters represented by f are theta_rho*(f), and
these do depend on rho. F4 has trivial center; numbering and weight
coordinates are Sage/Bourbaki.

Trace projection is the construction in
[Kostant--Michor, Section 2.1](https://arxiv.org/pdf/math/0109066).
The polynomial-character/Borel construction is from
[Kumar, Part I](https://kumar.math.unc.edu/papers/kumar76.pdf).
[Kumar--Xie, Part III, Section 6](https://kumar.math.unc.edu/papers/kumar86.pdf)
supplies the minimal-representation generators implemented by the catalogue.
As Cartan polynomials they also generate S_Q(t^*)^(W_L) when we change
representation; their Laurent pullbacks then change.

The preceding [F4 proof](F4_CHARACTER_FORM.md) establishes
C_26=Sym_Z(6 Q_long), where Q_long is the long-root lattice.
The adjoint identification and scaling argument below are deductions
recorded here. No claim is made that they are stated in those papers, nor
that they are new in the literature. A comprehensive priority search has
not been carried out.

## The adjoint torus source

The adjoint weights are all 48 roots, each with multiplicity one, and zero
with multiplicity four. The nonzero weights of V(omega4) are the 24 short
roots, each with multiplicity one. With fundamental-weight coordinates,
the adjoint trace Gram matrix is

\[
M_{52}=\begin{pmatrix}
36&-18&0&0\\
-18&36&-36&0\\
0&-36&72&-36\\
0&0&-36&72
\end{pmatrix}=3M_{26}.
\]

Both weighted weight lattices D_rho span X(T), so the
[linear-lattice formula](CHARACTER_SOURCES.md#the-linear-lattice-directly-from-weights)
gives

\[
\Lambda_{52}=3\Lambda_{26}
=18Q_{\rm long}
=\mathbb Z(18\omega_1)\oplus\mathbb Z(18\omega_2)
 \oplus\mathbb Z(36\omega_3)\oplus\mathbb Z(36\omega_4).
\]

Independently, the Laurent pullbacks B_i of
x=(18omega1,18omega2,36omega3,36omega4) have integral coefficients and a
unit coefficient minor at

\[
(-1,0,0,1),\quad(-1,0,1,-1),\quad
(-1,0,1,0),\quad(-1,1,-1,0).
\]

This certifies the complete linear lattice. It is not sufficient on its
own to determine higher-degree integral characters.

For that step, use the integer unimodular row change

\[
U=\begin{pmatrix}
1&-1&0&1\\
0&1&-1&1\\
0&0&0&1\\
0&0&1&-2
\end{pmatrix}.
\]

In descending lexicographic order, UB has the following leading data:

| Characteristic | Leading exponents | Leading coefficients |
|---|---|---|
| 0 and every odd prime | (2,-1,0,0), (1,1,-2,0), (1,0,0,0), (1,0,0,-1) | 1, 1, 2, 1 |
| 2 | (2,-1,0,0), (1,1,-2,0), (1,0,-1,1), (1,0,0,-1) | 1, 1, 1, 1 |

Both exponent matrices have nonzero determinant over Q. Consequently,
distinct monomials in UB have distinct leading Laurent exponents in each
prime characteristic. The four pullbacks are algebraically independent
modulo every prime. The [least-denominator argument](CHARACTER_SOURCES.md#exact-criterion)
now gives, including inhomogeneous polynomials,

\[
\boxed{\mathcal C_{52}
=\mathbb Z[18\omega_1,18\omega_2,36\omega_3,36\omega_4]
=\operatorname{Sym}_{\mathbb Z}(18Q_{\rm long}).}
\]

The characteristic-two case is essential to this certificate: the third
characteristic-zero leading coefficient is even.

## Why scaling follows for every Levi and every degree

The source identifications imply the equality of homogeneous lattices
C_52^d=3^d C_26^d inside S_Q^d(t^*). Multiplication by 3^d is injective
and commutes with W_L. An element 3^d f is invariant if and only if f is.
Taking invariants therefore proves A_(52,L)^d=3^d A_(26,L)^d.
Applying the same Q-linear Borel map proves the asserted image equality.

More generally, this argument applies whenever two representations have
certified sources Sym_Z(Lambda) and Sym_Z(c Lambda), for a positive rational
c. The conclusion requires those full-source identifications; proportional
trace matrices or linear lattices alone are insufficient.

In F4, c=3 is an integer. If r_d is the rank of J_(26,L)^d, then

\[
J_{26,L}^{d}/J_{52,L}^{d}
\cong(\mathbb Z/3^d\mathbb Z)^{r_d}.
\]

In degree zero this quotient is zero. In the maximal-parabolic cases r_d
is the full Betti number. If s_1,...,s_(r_d) are the previous ambient Smith
factors, the adjoint factors are exactly 3^d s_1,...,3^d s_(r_d).
Thus all primary factors other than the 3-primary ones are unchanged.
After localizing at any prime other than 3, the embedded images agree.

Degreewise multiplication by 3^d respects products, since
3^d 3^e=3^(d+e). It gives an abstract graded Z-algebra isomorphism from
each minimal source/image onto its adjoint counterpart. Their embeddings
in the common polynomial/cohomology rings are different.

## The maps themselves are different

The Laurent monomial t1=e^(omega1) occurs only in the adjoint root support.
Its coefficient vector in (theta_52*(omega_i)) is

\[
(1/9,\ 1/6,\ 1/9,\ 1/18),
\]

whereas all four coefficients are zero for theta_26. Thus these normalized
Springer maps are not constant multiples of one another. This distinction
is also checked directly in the experiment and unit tests.

## Complete image comparison and API

The degree-one adjoint images are 18 Z h for P1/P2 and 36 Z h for P3/P4,
where h=beta(omega_i) is the primitive divisor class. In degree two:

| Parabolic | Minimal Smith factors | Adjoint Smith factors |
|---|---|---|
| P1 | 36 | 324 |
| P2 | 36, 72 | 324, 648 |
| P3 | 72, 144 | 648, 1296 |
| P4 | 144 | 1296 |

Every row of the [74-row minimal table](F4_CHARACTER_TABLES.md) transforms
by the same proved rule. The
[comparison dataset](experiments/f4-adjoint-comparison-results.json) records
all 74 adjoint image lattices, Smith lists and relative quotient factors.
These are derived tables, not 74 new independent Borel computations.

```python
from sage_borel import f4_character_source

minimal = f4_character_source(2)  # existing default remains unchanged
adjoint = f4_character_source(2, representation="adjoint")
old = minimal.degree(2).image()
new = adjoint.degree(2).image()
assert new.smith_factors == (324, 648)
assert new.compare(old).left_in_right
assert old.quotient_by(new).torsion_factors == (9, 9)
```

The accepted representation names are `"minimal"` and `"adjoint"`.
For other highest weights, the general `springer_character_map` and
`IntegralCharacterSource` APIs remain available, but they require their own
successful all-primes certificate. No classification is inferred from
these two F4 representations.

## Verification scope

The [Sage experiment](experiments/f4_adjoint_comparison.py) independently
reconstructs the long-root orbit, compares weights with the Weyl character
calculation, checks the identity differential and linear unit minor, and
verifies both leading-term cases. Fresh API checks cover degrees 0..4 for
all four parabolics. Direct homogeneous and filtered Laurent-intersection
certificates cover degrees 0..3 for all four parabolics (32 lattices).

The [standard-library Python audit](experiments/audit_f4_adjoint_comparison.py)
reconstructs both root orbits, trace projection, the unit minor and both
characteristic cases without Sage or runtime-library imports. It checks
all 74 scaled embedded lattices and their Smith factors by determinantal
divisors, and reconstructs all 32 direct Laurent certificates. Its
[saved result](experiments/f4-adjoint-comparison-audit.json) identifies exactly
which checks are fresh and which earlier Borel tables are reused.

```bash
PYTHONPATH=src sage -python experiments/f4_adjoint_comparison.py
python experiments/audit_f4_adjoint_comparison.py
PYTHONPATH=src sage -python tools/verify_character_release.py
```

Earlier minimal-representation reports retain their historical source
hashes. The new report binds the exact baseline table and audit it uses;
it does not present historical hashes as hashes of the current runtime.
