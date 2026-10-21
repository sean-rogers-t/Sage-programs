# Complete E6 quotient coordinates and full Cayley determinants

Date: 2026-10-20. Version 0.18 computes all basic quotient coordinates for
the minimal V(omega_1), dimension 27, and adjoint V(omega_2), dimension 78,
with the normalized trace-projection Springer morphisms. Both full
78-dimensional Cayley determinants are given as explicit exact products.
The calculations have independent representation-matrix and unisolvent
certificates. These are repository deductions and computations, without
a claim of originality or external review.

## Notation and connection to the papers

G is the simply connected complex group E6, T its maximal torus, W its Weyl
group, and theta_rho:G -> g the Springer morphism. We use Sage/Bourbaki
numbering, simple coroot coordinates w_i on t, and t_i=e^{omega_i} on T.
Put chi_i=character of V(omega_i). Their dimensions are
(27,78,351,2925,351,27). In particular chi_3=Lambda^2 chi_1,
chi_4=Lambda^3 chi_1, chi_5=chi_3^dual and chi_6=chi_1^dual.

The trace-projection map is the generalized Cayley map of
[Kostant--Michor, Section 2](https://arxiv.org/pdf/math/0109066).
The coordinates are computed from
[Rogers's explicit formula, Theorem 1](https://arxiv.org/pdf/1701.01538v1):

\[
 M_\rho=\sum_\mu m_\mu\mu\mu^t,\qquad
 \theta_\rho(t)=M_\rho^{-1}\sum_\mu m_\mu\mu e^\mu(t).
\]

Here M_27=6 C and M_78=24 C, where C is the E6 Cartan matrix. Thus
d theta_1=id in both cases. There is no independent scalar choice for the
two determinants. This continues the [general geometry construction and
G2 calculation](SPRINGER_GEOMETRY.md), retaining the positive-degree
G-invariants killed by the Borel map. It does not change the earlier
integral Levi/Borel image calculations.

## Basic invariants on the Lie-algebra quotient

Let Omega be the 27 weights of V(omega_1), each of multiplicity one, and set

\[
 p_d(h)=\operatorname{tr}_{27}(h^d)=\sum_{\mu\in\Omega}\mu(h)^d.
\]

The ordered target coordinates used for both maps are

\[
 f=(p_2/12,p_5,p_6,p_8,p_9,p_{12}).
\]

The first coordinate is (h,h)/2 with long roots of squared length two:
p_2/12=(1/2) w^t C w. This matches the quadratic normalization in version
0.17. The other five coordinates are unscaled traces, not divided by d!.

These are Weyl invariant because Omega is a whole Weyl orbit. Their degrees
are the fundamental E6 degrees 2,5,6,8,9,12. Their Jacobian satisfies the
exact identity

\[
 \det\left(\frac{\partial f_i}{\partial w_j}\right)
 =423263232000\prod_{\alpha>0}\alpha(h).
\]

Proof: the determinant is alternating, hence divisible by every positive
root. Its degree is 1+4+5+7+8+11=36, the number of positive roots, so it is
a constant times the displayed product. At h=C^{-1}(1,2,3,4,5,6)^t
=(27,39,53,76,56,31)^t, direct rational evaluation gives the nonzero constant
423263232000. This proves algebraic independence. The polynomial subalgebra
they generate has the same Hilbert series as Q[t]^W, so it is all of Q[t]^W.
Thus the computed six pullbacks describe the entire map T/W -> t/W.

## All quotient coordinates in fundamental characters

The API returns exact polynomials F_{rho,d}(chi_1,...,chi_6) satisfying
F_{rho,d}(chi(t))=f_d(theta_rho(t)). Their numbers of nonzero terms are:

| degree | 27-dimensional map | adjoint map |
| --- | ---: | ---: |
| 2 | 3 | 5 |
| 5 | 16 | 20 |
| 6 | 34 | 93 |
| 8 | 82 | 270 |
| 9 | 106 | 272 |
| 12 | 370 | 1459 |

For example, the quadratic formulas recover

\[
 F_{27,2}=\chi_1^2/54-\chi_3/36-5\chi_6/36,
\]
\[
 F_{78,2}=\chi_2^2/576-\chi_1\chi_6/144-\chi_2/192-\chi_4/576-1/64.
\]

All coefficients, including the degree-12 formulas, are in the small
[runtime coefficient file](src/sage_borel/data/e6_geometry.json).
[The readable formula report](experiments/e6-geometry-formulas.json) prints
each complete polynomial. The API constructs ordinary Sage polynomials;
these are not unevaluated sums standing in for missing character formulas.
The raw data represents each coordinate by an integer polynomial divided
by its recorded scale; degree two has the additional target normalization
of 12. `e6_adjoint_quotient()` applies this normalization.

## Why exact interpolation proves the formulas

If lambda is the source highest weight and d the target degree, every
Laurent exponent in f_d(theta_lambda(t)) is at most d lambda in root order.
Indeed every source weight is lambda minus a nonnegative integral sum of
simple roots, and multiplying d terms preserves that condition.

For a dominant weight nu=sum a_i omega_i define the fundamental-character
monomial chi^a. Its highest weight is nu with coefficient one; all its
other weights are lower. Highest-weight triangular subtraction therefore
proves that every invariant Laurent polynomial with the indicated support
bound has a unique expression in the finite basis

\[
 \mathcal B_{d,\lambda}=
 \{\chi^a:a_i\geq0,\ C^{-1}(d\lambda-a)\in\mathbb Z_{\geq0}^6\}.
\]

This is a completeness bound, not a guessed truncation in total polynomial
degree. The basis counts in degrees 2,5,6,8,9,12 are
(3,19,34,82,125,370) for lambda=omega_1 and
(5,49,93,270,425,1459) for lambda=omega_2.

We choose rational torus points with coordinates among +/-1,...,+/-5.
Their evaluation matrix on the relevant basis has full rank modulo
1000003. All denominators are prime to that modulus, so the rational
evaluation matrix is invertible. We solve the resulting integer-cleared
system over QQ, check the equations exactly, and obtain integer numerator
coefficients after scaling by 18^d for the minimal source or 24^d for the
adjoint source. No floating-point interpolation, rational reconstruction
from a few primes, or unproved coefficient-height bound is used.

The minimal calculation shares a basis of size 82 for degrees 2,5,8 and a
basis of size 370 for degrees 6,9,12. The adjoint calculation uses size
1459 for all six degrees. The output's support is checked against the
dominance bound. The exact certificate points and coefficients are stored in
[the minimal report](experiments/e6-minimal-quotient-interpolation.json) and
[the adjoint report](experiments/e6-adjoint-quotient-interpolation.json).

## Full Cayley determinants

Let J_T be the logarithmic 6 by 6 derivative of theta on T, and write
T_rho(chi)=det J_T. The full determinant is

\[
 \boxed{\Psi_\rho(t)=T_\rho(\chi(t))\,e^{2\varrho}(t)
               \prod_{\alpha>0}r_{\rho,\alpha}(t)^2,\qquad
 r_{\rho,\alpha}=\frac{\alpha(\theta_\rho(t))}{e^\alpha(t)-1},}
\]

where varrho is the half-sum of positive roots. This follows from the
root-direction derivative in the [general proof](SPRINGER_GEOMETRY.md#a-general-full-differential-formula):
the alpha direction contributes e^alpha r_alpha and its opposite contributes
r_alpha. There are six Cartan directions and 72 root directions, giving
the full 78-dimensional determinant.

For E6 the root quotients have particularly simple explicit forms. Using
the inner product with long roots of squared length two, put

\[
 \Omega^-_\alpha=\{\mu\in\Omega:(\mu,\alpha)=-1\},\qquad
 R^-_\alpha=\{\beta\in R:(\beta,\alpha)=-1\}.
\]

Then |Omega^-_alpha|=6 and |R^-_alpha|=20, and

\[
 \boxed{r_{27,\alpha}=\frac16\sum_{\mu\in\Omega^-_\alpha}e^\mu,\qquad
 r_{78,\alpha}=\frac1{24}
       \left(2+2e^{-\alpha}+\sum_{\beta\in R^-_\alpha}e^\beta\right).}
\]

Proof: group the trace moments into alpha strings. The minuscule module
has six strings of length one in this direction; their differences factor
by e^alpha-1. The adjoint module has twenty such strings, plus the string
with weights alpha,0,-alpha. Its nonzero moment contribution is
2(e^alpha-e^{-alpha}), yielding 2+2e^{-alpha} after division. The indices
6 and 24 give the stated denominators. Every root quotient is also
checked by exact Laurent division in the implementation.

For the minimal map T_27=P_27/139968, where

\[
\begin{aligned}
 P_{27}={}&9\chi_2^3-4\chi_1\chi_2\chi_6+4\chi_2\chi_5\chi_6
 -16\chi_3\chi_6^2+16\chi_6^3-54\chi_2^2\\
 &-48\chi_1\chi_3-36\chi_2\chi_4+48\chi_3\chi_5
 +72\chi_1\chi_6-72\chi_5\chi_6+216\chi_4.
\end{aligned}
\]

T_78 has 49 terms; its complete expression is in the readable formula
report and `e6_cayley_determinant(representation="adjoint").torus_character_polynomial`.
The same exact interpolation proof uses dominance bounds 6 omega_1 and
6 omega_2, of sizes 34 and 93. The corresponding
[minimal](experiments/e6-minimal-determinant-interpolation.json) and
[adjoint](experiments/e6-adjoint-determinant-interpolation.json) certificates
record all points and coefficients.

The product is the full determinant, not only the torus Jacobian. We retain
the product instead of expanding all 36 factors into a very large single
fundamental-character or Laurent polynomial. It gives an exact description
of the critical locus: on T, T_rho vanishes or one of the explicit r_alpha
vanishes. The monomial e^{2 varrho} never vanishes. As a conjugation-invariant
regular function, Psi is determined everywhere on G by its restriction to T.
At t=1, both normalizations give T_rho(chi(1))=r_alpha(1)=Psi(1)=1.

## Independent checks

The [matrix oracle](experiments/e6_matrix_oracle.py) starts from a hardcoded
E6 Cartan matrix and constructs the minuscule weight graph. It builds
27 by 27 simple raising/lowering matrices and checks 168 Chevalley/Serre
relations. Commutators produce a 78-element Lie basis. It then computes
the 78 by 78 adjoint matrices, checking all 6084 brackets against that basis.
The matrix construction uses neither WeylCharacterRing nor SpringerCharacterMap
nor the interpolation code.

For each representation the audit computes the trace-projection differential
directly from these matrices: B_ij=tr(X_j X_i), A_ij(t)=tr(rho(t)X_j X_i),
and differential B^{-1}A. Weight grading gives a 6 by 6 Cartan block and
72 scalar root blocks. All Cartan entries and all 72 root blocks agree
coefficientwise with the primary calculation, for both representations.
Consequently their product verifies the full determinant identity without
an enormous expansion of that product.

The audit also recomputes each quotient coordinate from direct diagonal
matrix traces at every unisolvent point, checks that the candidate lies in
the proved finite-dimensional space, and independently checks full rank of
each evaluation matrix. This verifies 10110 exact scalar quotient identities
and 127 Cartan-determinant identities. Thus the checks prove the formulas
globally via unisolvence; they are not merely generic-looking spot checks.
Exterior-character identifications are separately checked with Sage's
WeylCharacterRing. Sage rational arithmetic and matrix operations are shared;
the independent representation construction and trace formulas are separate.

See [the audit program](experiments/audit_e6_geometry.py) and
[its report](experiments/e6-geometry-audit.json). The seven regression methods
also cover fresh rational/finite-field points, every root-string formula,
old quadratic normalizations, exact point handling, expansion guards and
independent full-determinant fixtures.

## API

```python
from sage.all import QQ
from sage_borel import (e6_adjoint_quotient, e6_cayley_determinant,
                       e6_fundamental_characters)
point = tuple(map(QQ, [2,3,5,7,11,13]))
chars = e6_fundamental_characters(point)
for representation in ("minimal", "adjoint"):
    quotient = e6_adjoint_quotient(representation=representation)
    assert quotient.degrees == (2,5,6,8,9,12)
    for polynomial, torus_function in zip(quotient.character_polynomials,
                                        quotient.torus_pullbacks):
        assert polynomial(*chars) == torus_function.evaluate(point)
    determinant = e6_cayley_determinant(representation=representation)
    assert len(determinant.root_factors) == 36
    assert determinant.evaluate([1]*6) == 1
    assert determinant.evaluate(point) > 0
```

The quotient's `.invariants` are explicit Cartan polynomials. Its
`.character_polynomials` are the complete fundamental-character formulas.
Its `.torus_pullbacks` are `PowerSumPullback` objects retaining sums of powers;
use `.evaluate(point)` or explicitly bounded `.expand(max_terms=...)`.
`e6_fundamental_characters(point)` evaluates the six source coordinates
without constructing large representation weight tables. The determinant
object exposes `.torus_jacobian_matrix`, `.torus_character_polynomial`,
`.root_factors`, `.evaluate(point)` and explicitly bounded `.expand(...)`.
The general `SpringerGeometry(s).factored_cayley_determinant()` is also
available for other supported types.

The bounds apply after individual Sage operations; they are not peak-memory
or wall-clock guarantees. Rational and number-field points give exact
values. Finite-field evaluations reduce these characteristic-zero identities
only at primes where their denominators are invertible.

## Reproduction and remaining scope

Run from borel/ with Sage:

```sh
PYTHONPATH=src sage -python experiments/e6_quotient_interpolation.py minimal
PYTHONPATH=src sage -python experiments/e6_quotient_interpolation.py adjoint
PYTHONPATH=src sage -python experiments/e6_quotient_interpolation.py minimal determinant
PYTHONPATH=src sage -python experiments/e6_quotient_interpolation.py adjoint determinant
PYTHONPATH=src sage -python experiments/package_e6_geometry.py
PYTHONPATH=src sage -python experiments/audit_e6_geometry.py
PYTHONPATH=src sage -python -m unittest tests.test_e6_geometry -v
```

All requested higher coordinates and full determinant products are complete
for these two E6 representations. The subsequent [degree calculation](EXCEPTIONAL_DEGREES.md)
gives degrees **1080** (minimal) and **2754** (adjoint) on simply connected E6,
including the adjoint boundary-curve correction. The API returns these through
`generic_degree()`. Single expanded invariant polynomials for the full
determinants and target branch loci remain uncomputed. The corresponding
higher E7 calculations are in the [E7 extension](E7_GEOMETRY.md).
