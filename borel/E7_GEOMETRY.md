# Complete E7 quotient coordinates and full Cayley determinants

Date: 2026-10-20. Version 0.19 computes all seven basic quotient coordinates
for the minimal V(omega_7), dimension 56, and adjoint V(omega_1), dimension
133, with the normalized trace-projection Springer morphisms. Both full
133-dimensional Cayley determinants are explicit exact products. The
calculations have exact unisolvent certificates and independent
representation-matrix checks. These are repository deductions and
computations, without a claim of originality or external review.

## Notation and relation to the papers

G is the simply connected complex group E7, T its maximal torus, W its Weyl
group, and theta_rho:G -> g the Springer morphism. We use Sage/Bourbaki
numbering, simple coroot coordinates w_i on t, and t_i=e^{omega_i} on T.
Write chi_i for the character of V(omega_i). Their dimensions are
(133,912,8645,365750,27664,1539,56).

The map is the normalized trace projection studied by
[Kostant--Michor, Section 2](https://arxiv.org/pdf/math/0109066).
We apply [Rogers's explicit formula, Theorem 1](https://arxiv.org/pdf/1701.01538v1):

\[
 M_\rho=\sum_\mu m_\mu\mu\mu^t,\qquad
 \theta_\rho(t)=M_\rho^{-1}\sum_\mu m_\mu\mu e^\mu(t).
\]

Here M_56=12 C and M_133=36 C, with C the E7 Cartan matrix and det C=2.
Thus d theta at the identity is the identity for both representations.
This extends the [complete E6 calculation](E6_GEOMETRY.md) using the same
normalization. It studies the original Springer map before applying the
Borel map: the positive-degree G-invariants computed here would vanish in
the Borel coinvariant quotient. The earlier integral Levi/Borel image
calculations are a different output of the same explicit Springer formulas.

## Basic invariants on the Lie-algebra quotient

Let Omega be the 56 weights of V(omega_7), all of multiplicity one. Put

\[
 p_d(h)=\operatorname{tr}_{56}(h^d)=\sum_{\mu\in\Omega}\mu(h)^d.
\]

For both source representations, we use the ordered target coordinates

\[
 f=(p_2/24,p_6,p_8,p_{10},p_{12},p_{14},p_{18}).
\]

The first coordinate is (h,h)/2=(1/2)w^t C w, with roots of squared length
two. The remaining coordinates are unscaled power traces, without a d!
denominator. In particular the adjoint-source calculation uses these same
56-dimensional target traces, allowing a direct comparison of the two maps.

These are basic invariants, not just seven invariant functions of the
expected degrees. Invariance follows from the complete Weyl orbit. Their
Jacobian has the exact identity

\[
 \det\left(\frac{\partial f_i}{\partial w_j}\right)
 =2676292792892522496000\prod_{\alpha>0}\alpha(h).
\]

Indeed the determinant is alternating, hence divisible by every positive
root. Its degree is 1+5+7+9+11+13+17=63, the number of positive roots. The
quotient is therefore constant. Exact evaluation at
h=C^{-1}(1,2,3,4,5,6,7)^t=(65,96,129,190,151,107,57)^t gives the displayed
nonzero constant. Algebraic independence follows; equality of the Hilbert
series with Q[t]^W then shows that these invariants generate Q[t]^W.

## All seven quotient coordinates in fundamental characters

The API supplies ordinary, complete Sage polynomials
F_{rho,d}(chi_1,...,chi_7) satisfying

\[
 F_{\rho,d}(\chi(t))=f_d(\theta_\rho(t)).
\]

All coefficients are in the [runtime data](src/sage_borel/data/e7_geometry.json).
The [readable formula report](experiments/e7-geometry-formulas.json) prints
every polynomial. Higher coordinates are fully determined in characters;
their torus pullbacks additionally retain efficient sums of powers.

| degree | 56-dimensional map | adjoint map |
| --- | ---: | ---: |
| 2 | 4 | 5 |
| 6 | 41 | 78 |
| 8 | 105 | 223 |
| 10 | 234 | 535 |
| 12 | 501 | 1231 |
| 14 | 972 | 2515 |
| 18 | 3101 | 8697 |

The table counts nonzero fundamental-character monomials.

The quadratic formulas recover the earlier tensor-Casimir calculation:

\[
 F_{56,2}=\chi_7^2/192-\chi_1/24-\chi_6/144-5/48,
\]
\[
 F_{133,2}=\chi_1^2/1296-5\chi_1/648-\chi_3/1296
                         -5\chi_6/1296-19/1296.
\]

Each raw data record stores an integer polynomial divided by its scale:
24^d for the minimal map and 36^d for the adjoint map. Degree two has the
additional target normalization of 24. The public constructor applies it.
As an additional check, the highest monomial of the raw integer numerator
has coefficient 2*3^d+54 on chi_7^d for the minimal map and coefficient 24
on chi_1^d for the adjoint map. The first follows from the pairings
(mu,omega_7)=+/-3/2 once each and +/-1/2 with multiplicity 27 each; the
second uses the twelve weights with each pairing +/-1 with the highest
root omega_1. All seven coordinates satisfy these leading-coefficient checks.
The character evaluator uses

\[
 \Lambda^2(56)=\chi_6+1,\quad
 \Lambda^3(56)=\chi_5+\chi_7,\quad
 \Lambda^4(56)=\chi_4+\chi_6+1,\quad
 \Lambda^2(133)=\chi_3+\chi_1.
\]

Together with the 56 weights, 126 roots and the 912-dimensional character
chi_2, these identities evaluate all seven characters. The 912 character
has 632 distinct weights; a separate root-string/Freudenthal calculation
checks every multiplicity without using Sage's WeylCharacterRing.

## Why the large character formulas are proved

If lambda is the source highest weight, every exponent of
f_d(theta_lambda(t)) is at most d lambda in root order. A product of d
source weights has this property. Fundamental-character monomials have
their indicated dominant highest weight with coefficient one and only lower
weights beneath it. Highest-weight triangular subtraction proves that the
complete finite basis is

\[
 \mathcal B_{d,\lambda}
 =\{\chi^a:a_i\geq0,\ C^{-1}(d\lambda-a)\in\mathbb Z_{\geq0}^7\}.
\]

This is a proved support bound. In degrees 2,6,8,10,12,14,18 its sizes are
(4,42,106,241,502,973,3102) for lambda=omega_7 and
(5,79,223,553,1231,2515,8697) for lambda=omega_1. All quotient degrees
are even, so the bases through degree 14 nest inside their degree-14
space. We solve degrees through 14 together, then degree 18 separately.

For efficiency, finite-field linear algebra and Chinese remaindering
discover integer numerator coefficients. They are **not** the final proof.
For each basis we save a square set of rational torus evaluation points
whose matrix is invertible modulo 8388593. No denominator vanishes modulo
this prime. Therefore the same evaluation matrix is invertible over Q.
We then check the proposed integer coefficients at **every** saved point
using exact integer-cleared rational arithmetic. The true answer and the
proposed answer belong to the proved space, and an invertible rational
evaluation matrix makes equality on those points a proof of global equality.
No coefficient-height assumption or probabilistic polynomial identity claim
is needed. Failure of an exact equation fails the computation.

The numerator integrality also follows directly: 24 mu(theta_56) and
36 mu(theta_133) are integral Laurent polynomials for every mu in Omega.
An integral invariant Laurent polynomial has integer fundamental-character
coefficients by the same highest-weight subtraction. For exact row clearing,
the highest-root marks (2,2,3,4,3,2,1) bound the coordinates of each
fundamental representation's weights, hence its torus denominators.

All points, coefficient lists, ranks, scales and exact check counts are in
the [minimal low](experiments/e7-minimal-low-interpolation.json),
[minimal degree-18](experiments/e7-minimal-high-interpolation.json),
[adjoint low](experiments/e7-adjoint-low-interpolation.json) and
[adjoint degree-18](experiments/e7-adjoint-high-interpolation.json) certificates.

## Full 133-dimensional Cayley determinants

Let J_T be the logarithmic 7 by 7 derivative of theta on T. Write
T_rho(chi)=det J_T and varrho for the half-sum of positive roots. Then

\[
 \boxed{\Psi_\rho(t)=T_\rho(\chi(t))e^{2\varrho}(t)
                   \prod_{\alpha>0}r_{\rho,\alpha}(t)^2,\qquad
 r_{\rho,\alpha}=\frac{\alpha(\theta_\rho(t))}{e^\alpha(t)-1}.}
\]

The [general differential argument](SPRINGER_GEOMETRY.md#a-general-full-differential-formula)
gives the alpha and -alpha eigenvalues as e^alpha r_alpha and r_alpha.
There are seven Cartan directions and 63 opposite-root pairs: this is the
full 133-dimensional determinant.
In these fundamental-weight torus coordinates, e^{2 varrho}=(t_1...t_7)^2.

The root quotients are explicitly

\[
 \boxed{r_{56,\alpha}=\frac1{12}
       \sum_{\substack{\mu\in\Omega\\(\mu,\alpha)=-1}}e^\mu,\qquad
 r_{133,\alpha}=\frac1{36}\left(2+2e^{-\alpha}
       +\sum_{\substack{\beta\in R\\(\beta,\alpha)=-1}}e^\beta\right).}
\]

The two sums have 12 and 32 terms respectively. Proof: group the trace
moments into alpha strings. The minuscule strings of length one give the
first sum after division by e^alpha-1. The adjoint module has 32 such
strings and the string alpha,0,-alpha; the latter contributes
2(e^alpha-e^{-alpha}), giving 2+2e^{-alpha} after division. The trace
indices give the denominators 12 and 36. Exact Laurent division also
checks each quotient in the implementation.
These Laurent sums extend the quotients across e^alpha=1; evaluation uses
the sums there, so the formula also applies at nonregular torus points.

T_56 and T_133 have 24 and 68 nonzero fundamental-character terms, supplied
in full by the API and readable report. Their exact certificates use
the bounds 7 omega_7 and 7 omega_1, with basis sizes 61 and 135:
[minimal](experiments/e7-minimal-determinant-interpolation.json) and
[adjoint](experiments/e7-adjoint-determinant-interpolation.json).
W-equivariance conjugates J_T, so its determinant is invariant. Expanding
its seven rows gives products of seven source weights, proving these
support bounds. Multiplication by det M_rho=2 I_rho^7 clears the matrix
normalization and leaves an integral invariant Laurent polynomial.
The minimal Cartan determinant has odd central character, consistent with
its degree seven; the quotient coordinates all have even central character.

The product above is retained as an exact product rather than expanded into
one enormous polynomial. It describes the critical locus on T explicitly:
T_rho=0 or one of the displayed r_alpha=0. The monomial e^{2 varrho} never
vanishes on T. Conjugation invariance determines the full regular function
on G from this restriction. At the identity,
T_rho(chi(1))=r_alpha(1)=Psi_rho(1)=1 for both maps.
These explicit functions are the localization factors in Kostant--Michor's
Theorem 3.2 for the two E7 representations; their differential is the
trace-pairing formula in Proposition 2.4(7).

There is also a useful identity connecting the two requested calculations.
Let c=2676292792892522496000 be the basic-invariant Jacobian constant above.
The fundamental-character Jacobian on T is the Weyl denominator
delta=e^{varrho} product_{alpha>0}(1-e^{-alpha}): it is alternating with
highest weight varrho and highest coefficient one. The chain rule gives
det(D_chi F_rho) delta = c product_{alpha>0} alpha(theta_rho) T_rho.
Combining this with the root product proves the exact identity

\[
 \boxed{\Psi_\rho(t)T_\rho(\chi(t))
    =c^{-2}\det(D_\chi F_\rho)(\chi(t))^2.}
\]

It holds across zeros by regularity. A separate regression checks this
relation using the explicit quotient polynomials and determinant products.

## Independent verification

The [matrix oracle](experiments/e7_matrix_oracle.py) starts with a hardcoded
E7 Cartan matrix, builds the 56-dimensional minuscule weight graph and
checks 231 Chevalley/Serre relations. Commutators construct a 133-element
Lie basis. It builds its adjoint action by checking all 17,689 brackets.
This construction uses neither WeylCharacterRing, SpringerCharacterMap
nor interpolation code.

For each representation, the [audit](experiments/audit_e7_geometry.py)
computes the trace-projection differential from those matrices. Its 49
Cartan entries and all 126 root blocks agree coefficientwise with the
primary computation. Weight grading therefore checks the entire
133-dimensional determinant without expanding its full product.

The audit independently evaluates matrix traces at all certificate points,
rechecks support containment and evaluation-matrix rank, and checks all
32,727 exact scalar quotient equations and 196 Cartan-determinant equations.
It reconstructs chi_2 by root strings and Freudenthal, checks the exterior
identities separately with WeylCharacterRing, and records full-determinant
finite-field fixtures. See [the report](experiments/e7-geometry-audit.json).
The independent representation construction is separate; Sage's exact
arithmetic and the optimized monomial-evaluation kernel are shared. The
audit verifies the exponent identities used by that kernel.

## API

```python
from sage.all import QQ
from sage_borel import (e7_adjoint_quotient, e7_cayley_determinant,
                       e7_fundamental_characters)
point = tuple(map(QQ, [2,3,5,7,11,13,17]))
characters = e7_fundamental_characters(point)
for representation in ("minimal", "adjoint"):
    quotient = e7_adjoint_quotient(representation=representation)
    assert quotient.degrees == (2,6,8,10,12,14,18)
    for polynomial, torus_function in zip(quotient.character_polynomials,
                                        quotient.torus_pullbacks):
        assert polynomial(*characters) == torus_function.evaluate(point)
    determinant = e7_cayley_determinant(representation=representation)
    assert len(determinant.root_factors) == 63
    assert determinant.evaluate([1]*7) == 1
    assert determinant.evaluate(point) > 0
```

The quotient exposes `.invariants`, `.character_polynomials` and
`.torus_pullbacks`; the latter have `.evaluate(point)` and guarded
`.expand(max_terms=...)`. The determinant exposes `.torus_jacobian_matrix`,
`.torus_character_polynomial`, `.root_factors`, `.evaluate(point)` and
guarded `.expand(...)`. Rational points and finite-field reductions at
primes with invertible denominators are supported. Bounds apply after
individual Sage operations; they are not peak-memory guarantees.

## Reproduction and remaining scope

From borel/, run each of the six cases with Sage (OMP_NUM_THREADS=1 and
OPENBLAS_NUM_THREADS=1 avoid oversubscription in modular linear algebra):

```sh
PYTHONPATH=src sage -python experiments/e7_quotient_interpolation.py minimal low
PYTHONPATH=src sage -python experiments/e7_quotient_interpolation.py minimal high
PYTHONPATH=src sage -python experiments/e7_quotient_interpolation.py adjoint low
PYTHONPATH=src sage -python experiments/e7_quotient_interpolation.py adjoint high
PYTHONPATH=src sage -python experiments/e7_quotient_interpolation.py minimal determinant
PYTHONPATH=src sage -python experiments/e7_quotient_interpolation.py adjoint determinant
PYTHONPATH=src sage -python experiments/package_e7_geometry.py
PYTHONPATH=src sage -python experiments/audit_e7_geometry.py
PYTHONPATH=src sage -python -m unittest tests.test_e7_geometry -v
```

The interpolation experiment needs NumPy and a Cython/C++ toolchain in
addition to Sage. The installed API only reads the certified data and does
not compile or interpolate anything. A failed or inconclusive reconstruction
raises an error; it does not write a passed certificate.

The higher quotient coordinates and full determinant products are complete
for these two E7 representations. The subsequent [degree calculation](EXCEPTIONAL_DEGREES.md)
gives degrees **7272** (minimal) and **33664** (adjoint) on simply connected E7,
including all boundary contributions. The API returns these through
`generic_degree()`. Target branch loci and single expanded invariant
expressions for the full determinants remain uncomputed. These degrees
do not determine the minimum Cayley degree of E7.
