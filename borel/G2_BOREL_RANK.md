# Vertical ranks of the G2 Springer maps

## Scope and result

Work over an algebraically closed field of characteristic zero. Let G be
the simple group of type G2, B=TU a fixed Borel, and theta=theta_rho the
trace-projection Springer map, normalized by d theta at 1 = identity.
We compare rho=V(omega_1), of dimension 7, and rho=V(omega_2), the
14-dimensional adjoint representation. For G2 the root and weight
lattices coincide. The coordinates x=e^{omega_1}(t), y=e^{omega_2}(t)
are characters, not diagonal entries of a defining matrix.

For fixed t, define F_t:U -> u by theta(tu)=h(t)+F_t(u), where
h(t)=theta(t) lies in t. This note determines **rank dF_t at every point**
for these two representations. A centralizer normal form is used on the
root hypertori. This is a vertical derivative: torus directions, negative
root directions, scheme-theoretic fibers and target branch loci are
separate questions.

The useful general result is the centralizer rank formula in Section 2.
It applies in every finite reductive type. The complete G2 table is its
first exceptional-group application here.

## 1. Published inputs and status

[Kostant--Michor](https://arxiv.org/abs/math/0109066), Sections 2 and 4,
provide the trace-projection construction, equivariance, centralizer
compatibility and Jordan-decomposition results. The torus computation
uses the character/trace-form description underlying
[Rogers's explicit determination](https://arxiv.org/abs/1701.01538).
Ordered root subgroups and Chevalley commutators are standard.

The rank reduction below is an elementary deduction from equivariance
and unipotent conjugation. The G2 table is a written deduction with exact
symbolic verification. We make no priority claim for either statement;
the source review in this step checked the relevant Kostant--Michor
results, not the entire literature on equivariant maps or orbit ranks.

The [automatic reconstruction](ROOT_BOREL_API.md) and its earlier
[recurrence proof](BOREL_RECONSTRUCTION.md) supply the computational
input. [The new audit](g2-borel-rank-verification.json) compares all six
symbolic coordinates with full trace calculations in dimensions 7 and
14, and checks the determinantal ideals used for the rank strata.

## 2. General centralizer rank theorem

Let G be a connected reductive group in characteristic zero, B=TU, and
theta:G -> g a regular conjugation-equivariant map satisfying
theta(B) subset b and theta(tU) subset h(t)+u. These properties hold for
the Springer maps under consideration. Put

    U_t = Z_U(t),       u_t = Lie(U_t),       psi_t(z)=theta(tz)-h(t).

Every element of tU is U-conjugate to tz for some z in U_t. At any such
normal form, put X=theta(tz). Then

    rank dF_t at tz = dim( [u,X] + im(d psi_t at z) ).                 (1)

The notation on the left means the derivative of theta restricted to
tU, identified with the derivative of F_t on U. Its rank at any
U-conjugate of tz is the same. Both spaces on the right lie in u.

**Proof.** Eliminate noncentral root coordinates in increasing height
using conjugation by root subgroups. At root alpha, the coefficient
to invert is 1-e^alpha(t), so coordinates with e^alpha(t) != 1 can be
removed without changing lower heights. The remaining element is tz
with z in U_t. Equivalently, one may conjugate the semisimple part into
T inside B and use the commuting Jordan decomposition.

Consider m:U x U_t -> tU, (v,w) -> vtwv^{-1}, at (1,z). Left or right
translation identifies its tangent directions with

    (Ad(tz)-1)u + u_t = u.

Indeed, split u into eigenspaces of Ad(t). Ad(z) commutes with Ad(t)
and is unipotent on each eigenspace. On an eigenspace with eigenvalue
r != 1, r Ad(z)-1 is invertible. The eigenvalue-one space is u_t,
already supplied by varying w. Thus dm is surjective. Differentiate
theta(m(v,w))=Ad(v)theta(tw). Its image is exactly the right side of
(1). Surjectivity of dm proves the equality. Conjugation by U preserves
tU and intertwines the two derivatives, proving rank invariance.

**Regular-torus corollary.** If e^alpha(t) != 1 for every root, U_t=1
and every point of tU is U-conjugate to t. Consequently

    rank dF_t = #{alpha > 0 : alpha(h(t)) != 0}
              = #{alpha > 0 : q_alpha(t) != 0},                      (2)

where

    q_alpha(t) = e^alpha(t) alpha(h(t))/(e^alpha(t)-1).

This rank is constant on tU. In particular, mixed terms cannot raise
the rank above the count of nonzero diagonal factors on this open set.
Rank variation in the unipotent parameters requires a nontrivial U_t.

Formula (1) is stronger than counting diagonal factors: a direction
with zero leading coefficient can contribute through [u,X]. It is also
stronger than a graph of nonzero commutators, since the rank of their
actual linear combinations accounts for dependencies and cancellations.

## 3. Conventions and exhaustion of the singular torus

Use Sage Chevalley vectors E_1,...,E_6 in the root order

| i | beta_i | e^{beta_i}(t) | beta_i(h), for h=h_1 H_1+h_2 H_2 |
| --- | --- | --- | --- |
| 1 | alpha_1 | x^2/y | 2h_1-h_2 |
| 2 | alpha_2 | y^2/x^3 | -3h_1+2h_2 |
| 3 | alpha_1+alpha_2 | y/x | -h_1+h_2 |
| 4 | 2alpha_1+alpha_2 | x | h_1 |
| 5 | 3alpha_1+alpha_2 | x^3/y | 3h_1-h_2 |
| 6 | 3alpha_1+2alpha_2 | y | h_2 |

The short roots have indices 1,3,4. The nonzero brackets for i<j are

    [E1,E2]=-E3, [E1,E3]=-2E4, [E1,E4]=-3E5,
    [E2,E5]=-E6, [E3,E4]=3E6.                                     (3)

The six root hypertori have the following one-to-one parametrizations
by s != 0:

| i | (x,y) |
| --- | --- |
| 1 | (s,s^2) |
| 2 | (s^2,s^3) |
| 3 | (s,s) |
| 4 | (1,s) |
| 5 | (s,s^3) |
| 6 | (s,1) |

Their multiple intersections are exactly:

* (1,1), where U_t=U;
* (-1,1), (1,-1), (-1,-1), with central root pairs (1,6), (4,2), (3,5);
* (zeta,1), zeta^2+zeta+1=0, with central positive roots 2,5,6.

This exhaustion follows by pairwise substitution in the six monomial
equations. Away from these five nonidentity points and the identity,
a point on their union has centralizer root system A1.

## 4. The complete rank rule

### 4.1 Off the root hypertori

Compute h(t) by the existing torus API and use (2). In G2 the possible
ranks there are 6, 5, and 0: a nonzero h can annihilate at most one
positive-root direction, because distinct positive roots in G2 are
nonproportional. If h=0, theta is identically zero on tU.

This is an exact evaluation criterion on the entire open set, without
listing or claiming to resolve the zero scheme of h there.

### 4.2 A single central root and h(t) != 0

Suppose only beta_i centralizes t. Write z=exp(v E_i). Equivariance
under T forces

    theta(tz) = h(t) + q_i(t) v E_i.

The toral part annihilates beta_i: the root reflection fixes t and
therefore fixes h(t). Since h != 0 and G2 has rank two, every other
root has nonzero value on h. Thus the transverse five directions in
(1) have rank five, and

    rank dF_t = 6 if q_i(t) != 0; 5 if q_i(t)=0.                    (4)

The rank is constant on tU in this case as well. Here are the central
factors on the parametrizations of Section 3:

| Representation | short central root | long central root |
| --- | --- | --- |
| minimal | (s^2+4s+1)/(6s) | (s^2+1)/(2s) |
| adjoint | (5s^2+2s+5)/(12s) | (s+1)^2(s^4-2s^3+4s^2-2s+1)/(8s^3) |

These equations give all rank-five exceptions in the scope of (4).
The exclusion h != 0 and the exclusion of multiple hypertori matter.

### 4.3 Single-central-root points with h(t)=0

For the minimal representation there are none: the zeros of h on a
short hypertorus have s=1,-1, and on a long hypertorus have
s=1,-1 or s^2+s+1=0. All are multiple intersections.

For the adjoint representation there are six additional points. In each
of the three long parametrizations, take

    3s^2-2s+3=0.

These points have exactly one central positive root, q_i=10/27, and
theta(t exp(v E_i))=(10/27)v E_i. Formula (1) and (3) give:

| central root i | v=0 | v != 0 |
| --- | --- | --- |
| 2 | 1 | 3 |
| 5 | 1 | 2 |
| 6 | 1 | 1 |

At these same six torus points the minimal map has nonzero toral part
and central factor 1/3, hence vertical rank 6 everywhere on tU.

For example ad(E_2) on u has image spanned by E_3,E_6; ad(E_5)
has image E_6; E_6 is central in u. None of these images contains E_i,
so the rank is one plus that commutator rank when v != 0.

Completeness follows from the exact common zero polynomials of the two
toral coordinates: on a short hypertorus they are proportional to
(s-1)(s+1)^3, and on a long hypertorus to
(s-1)(s+1)(s^2+s+1)(3s^2-2s+3), after removing Laurent units.

### 4.4 The three involutions

For each row below, write the centralizer normal form

    z=exp(A E_short) exp(B E_long).

These two root groups commute. The parameters A,B refer to z **after
U-conjugation**, not generally to the original ordered coordinates of u.
At all three involutions h=0. For the minimal map,

    theta(tz) = (A/3) E_short - B E_long;

for the adjoint map,

    theta(tz) = -(2A/3) E_short.

The full ranks are:

| t; (short,long) | minimal rank | adjoint rank |
| --- | --- | --- |
| (-1,1); (1,6) | 5 if A != 0; 2 if A=0 | 4 if A != 0; 1 if A=0 |
| (1,-1); (4,2) | 4 if (A,B) != (0,0); 2 otherwise | 3 if A != 0; 1 if A=0 |
| (-1,-1); (3,5) | 4 if A != 0; 3 if A=0,B != 0; 2 otherwise | 3 if A != 0; 1 if A=0 |

**Proof of every stratum.** For the minimal map, add the two constant
central derivative columns E_short,E_long to the matrix ad(X) on u.
For the adjoint map only E_short is added. Use (3). In the first row
the three transverse arrows form the E_1 string E_2 -> E_3 -> E_4 -> E_5;
their coefficient is a nonzero multiple of A. In the second row the
transverse matrix has rank two exactly when (A,B) != (0,0). In the third
row its rank is two for A != 0, one for A=0,B != 0, and zero otherwise.
For the adjoint map set B's coefficient in X to zero. The listed
central derivative directions are independent of these transverse
images. This proves the table, including the intersections where
generic-rank computations alone would not suffice.

The verifier additionally checks the radicals of the corresponding
minor ideals: they are (A), (A,B), the unit ideal, or zero exactly as
this description requires.

### 4.5 The two order-three points and the identity

At t=(zeta,1), zeta^2+zeta+1=0, use

    z=exp(b E_2) exp(e E_5) exp(f E_6).

Here h=0. All three central factors equal -1/2 for the minimal map
and 5/8 for the adjoint map. Hence the derivative inside U_t has image
the entire span of E_2,E_5,E_6 at every z. Modulo this space, the only
nonzero transverse commutator is a nonzero multiple of b E_3, arising
from [E_1,E_2]. Therefore **both** representations have rank 4 when
b != 0 and rank 3 when b=0. The coefficient of E_6 in theta(tz)
may include a be term; it does not affect this rank conclusion.

At t=(1,1), all six q_i equal 1. The triangular coordinate map has
Jacobian determinant 1, so its vertical rank is 6 everywhere.

Together Sections 4.1--4.5 cover all torus points and all unipotent
parameters in the stated normal forms.

## 5. Independent checks and what generalizes

[The verifier](tools/verify_g2_borel_rank.py) constructs the minimal
7-dimensional representation from its explicit sl2 weight strings,
checks the Chevalley--Serre relations and the 14-dimensional span, and
computes theta(tu) by matrix exponentials and trace pairings. It checks
the entire symbolic formula, not a list of numerical evaluations. The
adjoint check uses full 14-dimensional Lie-algebra matrices and the
same trace definition. Neither route uses the triangular fixed-vector
recurrence being tested. Root-vector signs are aligned to Sage's basis;
the adjoint route shares Sage's Chevalley constants. This is independent
computational construction, not a separate CAS or external review.

Exact quadratic number fields check both conjugate order-three points
and the six additional adjoint points simultaneously. Polynomial
minor ideals check all the stated normal-form rank strata. No inference
from floating-point ranks or generic samples enters the proof.

The centralizer formula (1), the regular-torus corollary, and the
height-elimination normal form generalize to every finite type and every
regular equivariant map with the stated Borel properties. The six-root
table, the special parameter polynomials and their ranks are G2-specific.
In higher rank a nonzero h may annihilate a larger root subsystem, so
the rank-two simplification in (4) must not be used without modification.

## Next step and purpose

Turn (1) into a centralizer-based rank evaluator that avoids expanding
all Borel coordinates. Its first test should be a singular F4 torus
point with a proper, nonabelian centralizer. Agreement with a full
adjoint trace derivative would establish that the reduction scales
beyond a rank-two table. The goal is reusable rank and collapse
calculations from smaller root systems, not another list of generic
determinants. Full fiber geometry remains a later, separate objective.
