# Type A translated Springer maps: character criteria and realization

Date: 2026-10-23. Characteristic zero, over C. These are written deductions
from trace projection and established character identities, with independent
finite matrix checks. No literature-priority claim is made.

## 1. The question and the answer

The [fiber assessment](CAYLEY_FIBER_ASSESSMENT.md) reduced type A translated
maps to matrix polynomials and asked which orders actual representations
can realize. Here is an arbitrary-rank answer for a named family, together
with a character criterion applying to general type A centralizer factors.

Let G=SL_m, V=C^m, and

    rho_a = direct sum over p=0,...,m of a_p copies of exterior^p V,
    a_p nonnegative integers,
    D = sum over p=1,...,m-1 of a_p binomial(m-2,p-1) > 0.

We use the actual simply connected group, its full character lattice, and
the normalized trace Springer map

    B_rho(theta_rho(g),Y) = tr(rho(g) d rho(Y)),
    B_rho(X,Y) = D tr(XY),        d theta_rho at 1 = identity.

No parabolic, integral cohomology lattice, or Borel characteristic map is
involved in this theorem. Upper triangular U is used only when restricting
the resulting nilpotent formulas to a fixed Borel.

**Theorem A (all exterior mixtures and all central points).** For zeta^m=1
and any nilpotent m-by-m matrix N,

    theta_rho(zeta(I+N)) = sum over j=1,...,m-1 of b_j N^j,

    b_j = (-1)^(j-1)/D *
          sum over p=j,...,m-1 of
          a_p zeta^p binomial(m-j-1,p-j).                    (1)

In particular theta_rho(zeta I)=0. The first nonzero b_j is the exact
degeneration order r. Equivalently, the map is polynomially equivalent
on the nilpotent cone, and on the upper nilradical, to N -> b_r N^r.
It is an isomorphism precisely when b_1 is nonzero. It is never identically
zero for a nontrivial exterior mixture: the triangular system (1) is
invertible in a_1,...,a_(m-1).

**Theorem B (closed classification at -I).** Suppose m is even. Let P be
the unique polynomial of degree at most m-2 with

    P(p)=a_p for p=1,...,m-1.

Write s=degree P and lambda for its leading coefficient. Then

    r = m-1-s,
    b_r = (-1)^(s+1) s! lambda / D.                         (2)

Thus every order 1,...,m-1 occurs for an actual faithful representation
of SL_m at this same central point: take a_p=p^s for p=1,...,m-1 and
s=m-1-r. These are honest direct sums with positive integral multiplicities;
they contain the defining representation. The case s=0 gives

    E = exterior^* V,
    theta_E(-(I+N)) = -N^(m-1)/2^(m-2).                    (3)

The trivial summands in E do not affect theta or D.

**Theorem C (the zero map is also realizable).** For even m, the faithful
representation E tensor E has nondegenerate trace form but

    theta_(E tensor E)(-u)=0 for every unipotent u in SL_m. (4)

Together, B and C realize every possible order of a nonzero equivariant
nilpotent map in type A, and also the zero map, for every even m. This is
an existence classification across representations, not a classification
of all representations or all noncentral fibers. The large representations
are explicit witnesses; no minimal-dimension assertion is made.

## 2. Proof by a determinant identity

The published exterior-character identity is

    sum_p z^p chi_(exterior^p V)(g) = det(I+zg).

Differentiate along g exp(epsilon Y). Its trace-pairing gradient is

    det(I+zg) * zg(I+zg)^(-1),                              (5)

followed by traceless projection. This rational expression is a polynomial
matrix: it equals zg adj(I+zg). Thus specialization at z=-1 is legitimate
after cancellation, including where I+zg is singular.

For g=I+N with N nilpotent, det(I+zg)=(1+z)^m. Expanding the inverse as
a finite series in N, the coefficient of N^j for j>=1 in (5) is

    (-1)^(j-1) z^j (1+z)^(m-j-1).                         (6)

The scalar term disappears under traceless projection, while every
positive power of N is traceless. Extract z^p, multiply by a_p zeta^p,
and sum. This proves (1), including its normalization. For completeness,
the quadratic index of exterior^p V is binomial(m-2,p-1): evaluate its
trace form on diag(1,-1,0,...,0), whose nonzero wedge weights consist of
the subsets containing exactly one of the first two positions. Its squared
trace is 2 binomial(m-2,p-1), versus 2 in V.

If p_max is the largest p in 1,...,m-1 with a_p nonzero, (1) at j=p_max
is a nonzero scalar. This proves the nonzero-map assertion without positivity
or genericity assumptions beyond D nonzero.

At zeta=-1, let Delta P(x)=P(x+1)-P(x). Equation (1) becomes

    b_j = (-1)^(m-j) Delta^(m-j-1) P(j) / D.              (7)

Finite differences of order greater than s annihilate P, and
Delta^s P=s! lambda. Thus b_j=0 for j<m-1-s and (2) follows. The endpoints
a_0,a_m are irrelevant; the interpolation polynomial uses only 1,...,m-1.
This proves the classification for every exterior mixture, not only the
particular witnesses p^s. For P=1 the only surviving term is j=m-1 and
D=2^(m-2), proving (3).

To justify polynomial equivalence, write p(x)=b_r x^r v(x), v(0)=1.
The change of variable Y=N v(N)^(1/r) has a finite polynomial expansion
modulo N^m and a polynomial inverse obtained by formal series reversion.
Then p(N)=b_r Y^r. It preserves nilpotence and strict upper triangularity.
For r>1 its differential at zero is zero, so it cannot be an isomorphism;
for r=1 the same coordinate change proves it is one. These are scheme
isomorphisms on the source and consequently transport fiber schemes.

Finally chi_E(g)=det(I+g), so chi_(E tensor E)=chi_E^2. For g=-u with
u unipotent, chi_E(g)=0. The differential of chi_E^2 in every group
tangent direction is consequently zero. Trace projection gives (4).
The trace form of E tensor E is 2 dim(E) B_E (the cross terms vanish on
sl_m), hence remains nondegenerate. Since E contains 1 and V, E tensor E
also contains V and is faithful. This proves C.

## 3. General character criterion on type A centralizer factors

Let G now be connected semisimple, rho have nondegenerate trace form,
t in a maximal torus, and H=Z_G(t)^0. Take a simple factor s=sl_m of the
derived Lie algebra of H, with a fixed defining matrix realization.
Write B_rho restricted to s as D_s tr(XY), D_s nonzero.

The operator rho(t) commutes with the action of s. Decompose into its
eigenspaces W_lambda, viewed as s-modules. For 2<=d<=m define the
primitive trace index I_d(W) by

    tr_W(d rho(X)^d) = I_d(W) tr_V(X^d)
                       modulo products of positive-degree invariants.

This fixes the normalization: I_d(V)=1. It is well-defined because the
invariant algebra for sl_m has basic degrees 2,3,...,m, each once.
Set M_d(t)=sum_lambda lambda I_d(W_lambda). Then, on the nilpotent cone
of this factor,

    theta_rho(t exp X)-theta_rho(t)
        = sum over k=1,...,m-1 of
          M_(k+1)(t)/(D_s k!) X^k.                        (8)

For this displayed formula the output is projected to s. When all derived
factors of H have type A, the entire translated map is the product of these
factor maps with zero central component, as in the previous assessment.
With other factors present, (8) still gives the restriction/projection
on the chosen factor; it asserts no classification for the other factors.

**Proof.** The twisted invariant P_d(X)=tr(rho(t)d rho(X)^d) has primitive
coefficient M_d(t). Its differential in direction Y is
d tr(rho(t)d rho(X)^(d-1)d rho(Y)); cyclicity applies because rho(t)
commutes with s. All positive-degree invariants vanish on nilpotent X,
so the differential of every product of at least two of them vanishes.
Thus for 2<=d<=m the gradient at X is d M_d(t) X^(d-1) in the defining
trace pairing. For d>m every homogeneous invariant is decomposable, so
its gradient vanishes there. Expanding exp(d rho(X)) in the defining
trace equation and dividing by D_s proves (8). Its k=0 term is removed
by subtracting theta_rho(t). This also explains why higher powers of
d rho(X), which need not vanish at exponent m, nevertheless contribute
nothing to the projected map.

Hence the first nonzero twisted primitive index gives the order; if all
M_2,...,M_m vanish the factor map is zero. Formula (8) bridges established
higher-index calculations to the translated map. On its own it is an
extraction criterion. Theorem B supplies the additional closed classification
and realizability statement for a concrete arbitrary-rank family.

For exterior powers the same determinant character identity gives

    I_d(exterior^p V)
       = sum_(l=1)^p (-1)^(l-1) l^(d-1) binomial(m,p-l).   (9)

Indeed expand log det(I+z exp X), then retain the terms linear in
tr(X^d) modulo decomposable invariants. This is a standard higher-index
character calculation, not a newly defined invariant. Equations (8) and
(9) describe exponential coordinates; (1) uses u=I+N. Their higher
coefficients differ, but log(I+N)=N+... preserves the first nonzero order
and its coefficient. Confusing these two coordinates gives incorrect
higher coefficients.

## 4. Relation to previous work and literature status

- Rogers' torus formulas and the character API give the trace Springer map.
  This note describes its values on central translates of unipotent elements.
  It does not change the Springer/Borel cohomology construction or prove a
  new integral-image or generic-degree formula.
- [Kostant--Michor, 4.12--4.13](https://arxiv.org/pdf/math/0109066)
  already define the translated centralizer map and give an invertible-
  differential isomorphism result. The present examples describe what
  can happen when the relevant derivative vanishes, outside that hypothesis.
- [Sobaje, Section 3.1](https://arxiv.org/pdf/1708.04153), together with
  the regular-nilpotent centralizer argument recorded in our assessment,
  supplies the established polynomial framework. We claim no new abstract
  classification of equivariant nilpotent maps.
- [van Ritbergen--Schellekens--Vermaseren, Section 5.1--5.3,
  especially (49)--(54)](https://arxiv.org/pdf/hep-ph/9802376)
  develops characters and generalized indices. The exterior generating
  identity and primitive-index method used here belong to that literature.
  The paper was inspected directly, including its normalization conventions;
  our convention is explicitly I_d(V)=1. Its source attribution points
  further back; no claim to have read every earlier index paper is made.

The finite-difference classification and faithful realization statements
are deductions written and checked here. Targeted searches for this exact
Springer-map formulation did not identify it, but were not an exhaustive
literature review. These elementary proofs may well be known or immediate
to experts. They settle the stated repository question without establishing
publication priority or resolving an explicitly posed Kostant--Michor problem.

## 5. Independent verification and next research gate

The [verifier](tools/verify_type_a_character_degeneration.py) builds exterior
Lie-algebra matrices by signed wedge substitution, independently of (1),
(8), (9), the character API and the Borel recurrence. It compares their
trace moments with every primitive coefficient in ranks m=2,...,8; checks
central points in exact cyclotomic fields; and verifies all claimed orders
for even m=2,4,6,8. A separate symbolic SL4 calculation uses compound
matrices of a general upper unipotent matrix and the full sl4 trace pairing.
It checks all exterior summands before checking the mixtures. The tensor
square zero-map claim is also checked through explicit tensor matrices
for SL2. These finite tests corroborate the written arbitrary-rank proofs;
they are not an external audit or a separate-CAS verification.

The next substantive gate is **irreducible ambient representations at
noncentral semisimple points**. At a central point any irreducible rho has
rho(t)=lambda I, so its translated map is lambda times the usual unipotent
Springer isomorphism and has order one. Therefore the arbitrary orders
above rely essentially on reducible representations with distinct central
characters. Restriction of a single irreducible ambient representation to
a proper centralizer can have several rho(t)-eigenspaces. Formula (8) is
the appropriate criterion there. A uniform branching-family classification,
or an obstruction proving its orders remain bounded, would advance beyond
the now-solved exterior-mixture realization problem. More central examples
would not do so.
