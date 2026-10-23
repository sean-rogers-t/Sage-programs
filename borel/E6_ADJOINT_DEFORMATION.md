# E6 adjoint degree by deformation of the character

2026-10-22. This note reproduces the existing boundary correction 1890
without using the Segre-class formula to derive it. It transfers the
half-spin local method to an exceptional root polytope, with additional
checks for the middle character layer and two-facet crossings.
It is a written characteristic-zero deduction with exact supporting
checks, not external review, formal verification or a priority claim.

## 1. Statement and conventions

Let G be simply connected E6 over C and rho its adjoint representation
V(omega_2), of dimension 78 in the numbering of our Kumar-paper notes.
The Springer/generalized Cayley map is normalized by trace projection:

\[
 \operatorname{tr}(\operatorname{ad}(\theta_\rho(g))\operatorname{ad}(Y))
 =\operatorname{tr}(\operatorname{Ad}(g)\operatorname{ad}(Y)).
\]

Its torus restriction is an invertible constant linear transformation of
the logarithmic gradient of the adjoint character. That character is
6+sum_(alpha in Phi)e^alpha; its constant term has zero derivative.
We first work on the root torus with character lattice M=Q(E6), using
simple roots as an integral basis. The simply connected torus has
character lattice P(E6) and maps to this root torus with degree three.
No parabolic, cohomology ring or integral-image comparison enters here.

**Theorem.** A general independent-coefficient deformation of the root
character has exactly 1890 generic-fiber branches specializing to the
toric boundary at the undeformed character. Their limits are distributed as

| Limiting stratum | Number of strata | Branches per stratum | Total |
| --- | ---: | ---: | ---: |
| Open critical curve in a facet | 810 | 1 | 810 |
| Triple junction inside a facet | 540 | 2 | 1080 |
| Critical point where two facets meet | 810 | 0 | 0 |

Consequently the root-torus degree is 2808-1890=918 and

\[
 \boxed{\deg\theta_{\mathrm{ad},E6}=3\cdot918=2754}
\]

on the simply connected group. The group, torus and adjoint-quotient
degrees agree. This is the degree of the specified adjoint Springer map,
not the minimum Cayley degree of E6. The numbers 1890 and 2754 were already
proved in [the exceptional degree note](EXCEPTIONAL_DEGREES.md); the
contribution of this note is a separate deformation derivation and an
explicit transfer criterion.

## 2. Total intersection and the actual facet lattice

Let Delta be the E6 root polytope in Q(E6), and X_Delta its normal
projective toric variety, with polytope line bundle L. Its normalized
volume is 2808. There are 54 facets, all of primitive height one, with
16 vertices each. They form two Weyl orbits of size 27, exchanged by
negation. Negation intertwines the root-character gradient with its
negative, so the two orbits give equal generic counts.

A facet has D5 half-spin geometry. There are 270 four-dimensional
crosspolytope faces, each lying in exactly two facets; the other 432
four-dimensional faces and all faces of smaller dimension are simplices.
The normal chart at each crosspolytope face is smooth and unimodular in
Q(E6). Its face character has three singular zeros, all with nonsingular
four-dimensional tangent Hessian. These account for the 810 crossings
in the theorem. A simplex face has no singular zero.

These statements are rechecked from Sage's E6 root lattice and the actual
polytope, without importing the previous degree tables or their custom
root-reflection implementation. The facet character requires a degree-two
torus cover to write it in five independent D5 epsilon coordinates.
The degree of this cover is separate from the final central factor three.

Here is an explicit embedding of Q(E6) into Z^6, with coordinates
(normal,x_1,...,x_5), computed from the D5 Cartan inner products:

\[
 B=\begin{pmatrix}
 1&0&0&0&0&0\\
 0&0&0&0&0&1\\
 0&0&0&0&1&-1\\
 0&0&0&1&-1&0\\
 0&1&1&-1&0&0\\
 0&-1&1&0&0&0
 \end{pmatrix}.
\]

Its image is Z times the even-sum sublattice of Z^5, and |det B|=2.
The 72 roots split into layers of sizes 16,40,16. On the covering
torus, with s the inverse normal character, the root character is exactly

\[
 \chi=s^{-1}f(x)+g(x)+s k(x),\qquad
 f=\frac{\prod_{i=1}^5(1+x_i)+\prod_{i=1}^5(1-x_i)}2,
 \quad k=f(x_1^{-1},...,x_5^{-1}),
\]
\[
 g=\sum_{i<j}(x_i+x_i^{-1})(x_j+x_j^{-1}).                 \tag{1}
\]

The deck involution is (s,x)->(s,-x). It acts freely on every facet
torus point. The difference lattice of the 16 facet roots is already
saturated in its hyperplane inside Q(E6); the factor two comes from
expressing that lattice in the larger epsilon lattice. Thus a count on
the displayed cover must be divided by two exactly once.

## 3. Critical strata and the middle-layer identity

The face equations f=df=0 say that at least two x_i equal -1 and at
least two equal +1. There are two sorts of strata on the cover:

* 30 open curves: exactly two minus signs, two plus signs, and one free
  coordinate q in C* minus {1,-1};
* 20 points: sign patterns with two signs of one kind and three of the other.

The first count is 5!/(2!2!)=30; the second is 2*binomial(5,2)=20.
The deck involution pairs the curves and points, giving 15 curves and
10 triple points in the actual facet torus.

Put Z_i=x_i+x_i^-1. On a critical curve the four fixed Z_i are
2,2,-2,-2. Their sum is zero and their pairwise-product sum is -8.
Thus g=-8 independently of q. Moreover

\[
 D_i g=(x_i-x_i^{-1})\sum_{j\ne i}Z_j=0                 \tag{2}
\]

for every i: the first factor vanishes on each fixed coordinate, and
the sum vanishes for the free coordinate. At a triple point every
x_i is a sign, so all derivatives vanish, while g=-8 again.
Hence the entire first derivative of the middle layer vanishes on
the critical locus. The constant six in the adjoint character would
replace g by g+6 and leaves this identity unchanged.

This is the new check needed to transfer the facet calculation.
For a general character the middle layer can have a nonzero derivative
along a critical stratum and change the count; see the limitations in
[the general degree audit](GENERAL_CAYLEY_DEGREE.md).

## 4. Deformation and the count inside a facet

Choose h with independent general coefficients on all 72 roots, and
replace chi by chi+eta h. This auxiliary deformation need not be a
representation character or Weyl invariant. Keep a general target y.
For general nonzero eta every proper face polynomial has nonsingular
zero set, so the generic fiber consists of 2808 reduced torus points.
At eta=0 the original gradient is dominant, since the trace-normalized
Springer map has identity derivative at the identity. Choose y so this
fiber is finite and reduced as well.

Properness of X_Delta makes every generic solution specialize to the
original torus fiber or the boundary. Each original torus point receives
one branch. It remains to count all boundary limits. Genericity conditions
below form finitely many nonempty Zariski-open conditions on (h,y).

On the facet cover write h=s^-1 H+J+s K. With D_i=x_i partial_(x_i),
the cleared gradient equations are

\[
 f+\eta H-s^2(k+\eta K)-y_0s=0,
\]
\[
 D_i f+\eta D_iH+s(D_i g+\eta D_iJ)
       +s^2(D_i k+\eta D_iK)-y_i s=0.                   \tag{3}
\]

Here y_0 pairs the target with the outward normal; its sign has been
chosen so that the first displayed equation has this form.
It is nonzero for a general target. The first equation eliminates s.

Near a critical stratum use Mobius transverse coordinates at x_i=+/-1.
The function f is a unit times A+B, where A and B are products on
disjoint variable sets. On an open curve they have sizes (2,2);
at a triple they have sizes (2,3) or (3,2). Unit coefficients can be
absorbed into one variable in each product with invertible Jacobian.
The radial-derivative identity gives

\[
 (A,B)\subset\mathfrak m_{\rm trans}J_{\rm trans}(f).
                                                               \tag{4}
\]

Consequently along a branch approaching the stratum, the transverse
equations in (3) imply J_trans(f)=O(eta,s), and f and its derivatives
in free parameters are o(eta,s). The normal equation then gives
s=O(eta). By (2), the middle-layer term s D_i g is o(eta) there.
Thus it does not change the leading reduced equations. When contact
coordinates are used, the free equation is still obtained by restricting
the original logarithmic system (3), not by arbitrarily replacing the
nonzero-target gradient by that of a unit multiple of f.

On a critical curve H restricts to a general linear polynomial A_0+B_0q.
The first-order normal equation gives s/eta=H/y_0. The free equation is

\[
 qH'(q)-\kappa H(q)=0,\qquad \kappa=y_q/y_0.             \tag{5}
\]

It has exactly one simple root
q=kappa A_0/((1-kappa)B_0). For general parameters this root is nonzero,
is not +/-1, and has H(q) nonzero. Thus s has order exactly one in eta.
At a triple point H is a general nonzero constant and the same normal
conclusion holds.

For completeness, the coefficients of H and its transverse first jets
really can be chosen independently where needed. For each fixed exponent
of the free q-coordinate, the allowed binary exponents on the other four
coordinates form one parity class. They affinely span Q^4. At a triple,
the even binary exponents in five coordinates affinely span Q^5.
Evaluation at nonzero signs therefore has ranks five and six respectively
for value plus transverse first derivatives. The verifier checks these
small rational matrices directly. The restrictions of a single global
general perturbation satisfy all finitely many required open conditions.

After elimination of s, the transverse equations have leading form

\[
 \prod_{j\ne i}a_j=\eta\alpha_i,\qquad
 \prod_{j\ne i}b_j=\eta\beta_i,                         \tag{6}
\]

with nonzero general leading constants. Equation (4) makes elimination
terms containing f an invertible change of the Jacobian generators.
The middle-layer terms have no new constant contribution, by (2).
More precisely (6) has eta times units on its right sides; it is not
an assertion that all higher-order terms vanish.

A product with k variables contributes k-1 branches. Dividing its
equations forces all variable orders to be 1/(k-1); after that rescaling
the leading solutions are simple. The formal implicit function theorem
lifts precisely those solutions after a finite base change in eta.
The free equation (5), when present, also has nonzero derivative.
There are no additional valuation patterns. Hence each open curve
contributes (2-1)(2-1)=1, and each triple contributes (2-1)(3-1)=2.

The count on the facet cover is therefore 30+2*20=70. The deck action
pairs all branches and is etale at their facet limits, so the actual
facet contribution is 35. The divisor s=0 is distinct from the facet
torus cover, and there is no extra normal ramification factor.

## 5. Two-facet crossings: an all-orders exclusion

It remains to rule out branches limiting to the codimension-two critical
points. This is essential; a sum of the facet counts alone is not a proof.
Choose an actual unimodular chart with normal coordinates s,t and four
tangent logarithmic coordinates u at one such point. Because both facets
have height one, the cleared root character has the form

\[
 F(s,t,u)=f_0(u)+s A(u)+t B(u)+s^2 C(u)+t^2 D(u)+st E(u)+\cdots,
\]

and the constant target monomial is st. Its equations are

\[
 D_uF-y_u st=0,\quad
 F-sF_s-y_1 st=0,\quad F-tF_t-y_2 st=0.                \tag{7}
\]

At the critical point u=0, f_0=df_0=0 and its Hessian H_0 is invertible.
The first four equations, including the eta perturbation, eliminate u
uniquely as a formal function of s,t,eta. Write the remaining equations
as R_1,R_2.

First, at eta=0 both R_i are divisible by st **to all orders**.
Each of the two incident facets has a critical curve approaching the
point. To see this independently of the earlier Segre computation, take
the free coordinate q in (1) to zero or infinity. Its limiting character
is the four-variable even or odd parity polynomial. Its six balanced
sign zeros on the cover are paired into the three actual cross-face
critical points. Thus every point is approached from both incident
facets by a boundary critical curve. The nonsingular tangent Hessian
makes these curves agree with the unique tangent solution on t=0 and
s=0 respectively. Both R_i therefore vanish identically on each axis;
in C[[s,t]] they are divisible by st.

Next their coefficients of st are exactly -y_1 and -y_2. This can be
checked without knowing higher normal terms. Let a=dA(0), b=dB(0).
The linear term of the tangent solution is
u=-H_0^-1(a s+b t). The st coefficient in F after substitution is
E(0)-a^T H_0^-1 b. The st coefficient in sF_s is the same, as is the
one in tF_t. They cancel in both normal equations of (7). Terms in
y_u st affect u only to second order and do not change this computation.
The script separately checks the pure s^2 and t^2 Schur identities from
the actual E6 roots at all three prototype critical points.

If h_0 denotes the value of the perturbing face polynomial at the point,
both reduced normal equations have eta coefficient h_0: the first-order
change in u contributes nothing because df_0=0. A general perturbation
has h_0 nonzero at every one of the finitely many crossing points.
We therefore have the exact formal form

\[
 R_i=st\,U_i(s,t)+\eta V_i(s,t,\eta),\qquad
 U_i(0,0)=-y_i,\quad V_i(0,0,0)=h_0.                  \tag{8}
\]

Take U_2 R_1-U_1 R_2. The st term cancels, leaving eta times a unit,
because its constant coefficient is h_0(y_1-y_2). For a general target
y_1 is different from y_2. Thus eta belongs to the completed local
fiber ideal: there is no branch with eta nonzero specializing to this
point. This is an all-orders statement, not an inference from quadratic
terms alone. The finite list of normal pairs gives finitely many target
hyperplanes to avoid.

All smaller faces are simplices and have no singular zero. We have now
exhausted all boundary limits, including the endpoints of the open curves.

## 6. Conservation of number and comparison with the earlier proof

There are 54 facets. Each receives 35 branches; crossings receive none.
Therefore

\[
 54\cdot35=54\left(15+2\cdot10\right)=810+1080=1890.
\]

Subtracting from 2808 leaves 918 original root-torus points at a general
target. The central triple cover gives 2754 for simply connected E6.

Only after obtaining this result do we compare the old calculation:
7*810-4*540-2*810=1890. That formula assigns global curve contributions
and Segre gluing corrections; the present deformation assigns limiting
branches to disjoint open strata. Their individual summands need not
agree. In particular, the old negative crossing correction does not
mean deformed branches must limit to the crossing.

## 7. What generalizes, with explicit hypotheses

Two local principles transfer beyond this E6 example.

**Height-one facet transfer.** Suppose a cleared character is
f+s g+O(s^2), its critical strata have two-product transverse germs,
and the restriction of dg to each stratum vanishes in tangent directions.
Suppose also that the restricted perturbation has the stated independent
coefficients and transverse first jets. Then elimination of the normal
coordinate gives the free equation dH-kappa H=0, and each of its simple
solutions with H nonzero contributes (k-1)(l-1) transverse branches.
Full stationarity dg=0, as in E6, is a stronger sufficient hypothesis.
If dg has transverse components, first-jet freedom must still ensure
the resulting leading transverse constants are nonzero. The facet
cover degree must be accounted for separately.

This follows by the elimination and order argument of Section 4.
The number of free solutions must be proved for the actual restricted
support; it is one for the E6 critical curves. A general representation
need not have any of these product, stationarity or jet properties.

**Two-normal exclusion.** If tangent elimination at a crossing gives
R_i=st U_i+eta V_i, with the two V_i having the same nonzero constant
and U_1(0) different from U_2(0), no deformed branch limits there.
The proof is the unit combination in Section 5. To use this principle
one must prove the all-orders st divisibility, not merely observe the
leading terms. Three or more boundary divisors, a singular normal chart,
or a degenerate tangent Hessian require a separate argument.

These principles and the independent E6 checks are the methodological
gain. They do not constitute an algorithm for every representation,
an all-exceptional degree classification, or a new value of the E6 degree.
Their published ingredients are standard toric top intersection,
conservation of number, and formal implicit functions, together with
[Kostant--Michor's degree comparison, Corollary 3.3](https://arxiv.org/pdf/math/0109066).
For the toric volume formula see Cox--Little--Schenck, *Toric Varieties*,
Theorem 13.4.3 ([authors' book page](https://dacox.people.amherst.edu/toric.html)).
No exhaustive novelty search for these local formulations was performed.

## 8. Independent verification and next goal

```bash
PYTHONPATH=src sage -python tools/verify_e6_adjoint_deformation.py
```

The [fresh verifier](tools/verify_e6_adjoint_deformation.py) and
[report](e6-adjoint-deformation-verification.json) reconstruct all 72 roots,
the two facet orbits, the index-two cover, the full 16/40/16 character
layers, the codimension-two orbit and smooth normal charts. They check
stationarity on every sign stratum, the Mobius product identity, the
first-jet ranks, the generic free-coordinate solution, simple transverse
product-gradient lengths, actual-root Schur coefficients, and the
universal normal cancellation. They do not import the old Segre verifier
or use a specialized finite-field degree as proof. Both implementations
use Sage arithmetic, so the independence is in constructions and proof
routes, not separate computer-algebra systems.

**Follow-up completed:** the [F4 minimal calculation](F4_MINIMAL_DEGREE.md)
now proves degree 96 with no boundary correction and a finite-flat torus
corollary. Its next target is the F4 minimal/adjoint multiplicity family.

**Original next goal:** inspect the F4 minimal weight polytope and
apply these criteria to its boundary characters. The purpose is to turn
the validated transfer method into a degree calculation not yet supplied
by this E6 benchmark. Success means either a complete degree proof with
every contributing face accounted for, or an explicit boundary germ
showing which of the current hypotheses must be extended. Start with
face and normal-layer classification before a large elimination run.
