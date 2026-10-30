# The final E8 adjoint case: a degree count at infinity

26 October 2026. This completes the A2+A1^3 case (seed 12357), which
represents eight standard Levi subsets. It supersedes the outstanding
gap in Section 7 of [the preceding four-case note](E8_FINAL_CASES.md).
The argument uses the original critical ideal, including its embedded
and nonreduced structure. Its characteristic-zero step is a homogeneous
Hilbert-function comparison, not an affine finite-field dimension claim.

## 1. The resulting theorem

Let G be a simple, simply connected complex algebraic group, with Killing
form B, and let its adjoint Springer--Cayley map be defined by

    B(theta_Ad(g),X) = tr(Ad(g) ad(X)).

For every proper Levi subgroup K, put

    Z_K^o = {xi in z(k) : Z_G(xi)=K}.

**Theorem, by the recorded case proofs.** There is a dense Zariski open
U_K in Z_K^o such that the actual scheme theta_Ad^{-1}(xi) intersected
with G_reg is smooth for every xi in U_K. Empty fibers are allowed.
Regular nonsemisimple points are retained. Positive-dimensional transverse
fibers in other cases are allowed; universal quasifiniteness is not asserted.

The new part proved below is the E8 A2+A1^3 case. Its generic nonempty
regular-source fibers, when present, are pure of dimension 12. Generic
nonemptiness is not proved here. With the existing
[Weyl-equivalence certificates](E8_WEYL_TRANSFERS.md), this raises E8
coverage from 247 to all 255 proper standard subsets. The earlier
classical, G2, F4, E6 and E7 results are indexed in
[the proof status](PROOF_STATUS.md) and [the ledger](RESEARCH_LEDGER.md).
This is a completion by cases, not a new uniform root-theoretic proof or
a literature-priority claim.

## 2. The original residual scheme

Use the full adjoint character and the constant-unit A2 elimination
already reconstructed from all E8 roots in
[the character data](e8-final-character-models.json). The remaining
derived A1 characters are x,y,w; a is the commuting A1 character. Put

    A=U+U^-1, B=V+V^-1, C=UV+(UV)^-1,
    H0=A^2+B^2+C^2-ABC-4,
    p=xy+wa, q=xw+ya, r=xa+yw.

Here U=z0^3/z1 and V=z0^6/z1; the exponent determinant is 3, so this
change of central torus coordinates is a finite etale isogeny. The
complete residual scalar potential has the compact expression

    f=x^2+y^2+w^2+a^2-p^2-q^2-r^2+xywa
      -Aqr-Bpr-Cpq-BCp-ACq-ABr-ABC-4.

This is exactly the previously recorded expanded potential. In
P=Q[x,y,w,a,A,B,C], set

    Y=V(f_x,f_y,f_w,H0),
    I=(f_x,f_y,f_w,H0, all 4-by-4 minors of their full 4-by-7 Jacobian).

There are 38 distinct nonzero generators. No radical, generic-parameter
replacement, or derived root-discriminant localization is made. We prove

    dim_C V(I) <= 2.                                    (1)

The central target has dimension three, so (1) excludes a dominant
Jacobian-deficiency image, including every component supporting a
persistent infinitesimal thickening.

## 3. Why the extra components at infinity can be used

Homogenize every one of the 38 generators separately with a new variable
tau (called h in the machine certificates). Let J be the resulting
homogeneous ideal in Q[x,y,w,a,A,B,C,tau]. This construction has extra
components at tau=0. We do not discard or saturate them.

The following standard degree argument is useful. Suppose a homogeneous
scheme has dimension at most d and degree at most e. If known distinct
d-dimensional components lying in tau=0 contribute at least e to its
degree, counted with their generic scheme multiplicities, every
d-dimensional component lies in tau=0. Its tau!=0 chart therefore has
dimension at most d-1. For a homogenized affine ideal, that chart is
V(I) times G_m, giving dim V(I)<=d-2.

The degree statement follows by a graded prime filtration: Hilbert
polynomials add; the top-dimensional contribution of a prime is its
degree multiplied by the length at its generic point. Lower-dimensional
and embedded contributions cannot cancel a positive top-dimensional
contribution. This is the usual multiplicity/cycle formalism; see
[Stacks, cycles associated to coherent sheaves](https://stacks.math.columbia.edu/tag/02QV).
We supply the required components, multiplicities, and upper bound here.

## 4. Four components already account for degree 21

The following homogeneous subschemes are contained in V(J):

    L   = V(A,B,C,tau),
    X_A = V(tau,A,p,Br+Cq),
    X_B = V(tau,B,q,Ar+Cp),
    X_C = V(tau,C,r,Aq+Bp).

All four have dimension four and are distinct. Containment of every
original homogenized generator is checked exactly over Q. It can also
be seen directly for X_A: at tau=A=0 the highest potential is
-p(Br+Cq). On p=Br+Cq=0, its three derived derivatives vanish. After
the H0-Jacobian row removes the A column, the remaining three rows span
at most the two differentials dp and d(Br+Cq), so every maximal minor
vanishes. The other two cases follow by the checked permutation symmetry.

Each X_A, X_B, X_C is an irreducible, generically reduced complete
intersection of degrees 1,1,2,3, hence has degree six. Here are the
irreducibility and boundary details. The quadric p=0 is irreducible.
On it, the common zero set q=r=0 has dimension one; an exact Groebner
basis of (p,q,r) in four variables is

    w^3 a-wa^3, x^2 a-w^2 a, y^2 a-w^2 a,
    xy+wa, xw+ya, yw+xa.

Off that one-dimensional set, Br+Cq=0 is a line bundle over the
irreducible three-dimensional quadric. Its complement has dimension
at most 1+2=3. The hypersurface is a complete intersection of pure
dimension four, so this complement contains no additional component.
It is generically reduced, and Cohen--Macaulayness excludes embedded
components. Thus it is integral. This argument also handles B=0 or
C=0 and the intersections with L; no boundary chart was deleted.

The linear component L has degree one and generic multiplicity at least
three. To see this, work at its generic point, with coefficient field
Q(x,y,w,a) and normal variables A,B,C,tau. Write the three homogenized
derived equations as

    g = M v + tau L0 + E,  v=(A,B,C)^T,  E in (A,B,C,tau)^3.

The matrix M is the derivative, in x,y,w, of (-qr,-pr,-pq). Its
determinant is

    D=-2pqr[-a^3+a(x^2+y^2+w^2)-2xyw],

which is nonzero (its value at (x,y,w,a)=(1,2,3,4) is 61600).
The three equations therefore eliminate A,B,C in the completed local
ring. The fourth equation and every homogenized Jacobian minor have
normal order at least three. Their images in the remaining one-variable
power-series ring lie in (tau^3). Once the dimension upper bound below
is applied, the local quotient is Artinian and its length is at least
three. Consequently these four components contribute at least

    3*1 + 1*6 + 1*6 + 1*6 = 21.                         (2)

## 5. An explicit member supplies the missing leading term

Let H=tau(A^2+B^2+C^2)-ABC-4tau^3 be the homogenized Markov relation,
let adj(M) be the adjugate, and put

    l=adj(M)L0, s=-tau*l, T=adj(M)g,
    Q_A=D^2 tau A+D tau s_A-D^2 BC,
    Q_B=D^2 tau B+D tau s_B-D s_A C,
    Q_C=D^2 tau C+D tau s_C-s_A s_B,
    R=D^3 H-Q_A T_A-Q_B T_B-Q_C T_C.                    (3)

Equation (3) is an explicit characteristic-zero membership R in J:
it is a polynomial combination of H and the original three g's, with
no denominator. It is homogeneous of total degree 30. Its normal-order
three part is

    tau^3 [D(l_A^2+l_B^2+l_C^2)+l_A*l_B*l_C-4D^3];

all other terms have normal order at least five. This follows also by
the three successive divided differences of H in A,B,C. With weighted
degree reverse lexicographic order, weights

    (2,2,2,2,1,1,1,1) on (x,y,w,a,A,B,C,tau),

the leading monomial of R is x^10 y^10 w^6 a tau^3. The coefficient is
nonzero over Q and at both certificate primes. Weighting here chooses
a term order; the Hilbert functions below use the ordinary total grading.
All ideals and members are homogeneous for that ordinary grading.

## 6. The finite certificate and the characteristic-zero conclusion

At primes 32003 and 32009, exact sparse F4 row operations on the original
homogenized generators are stopped after degree 20. The output is used
only as a list of original-ideal members, not asserted to be a complete
Groebner basis. The compact certificates retain 59 selected members
at each prime. Their weighted leading monomials, together with that
of (3), generate a monomial ideal of dimension four and degree 21.

Membership is replayed from the original gradients and all Jacobian
minors, using a pinned msolve source revision and a driver whose only
engine change is the degree stopping condition. The degree is also
computed independently without a computer-algebra Hilbert polynomial:
enumerate all coordinate supports, localize the monomial ideal at its
four-dimensional coordinate primes, and count the remaining standard
monomials. The lengths sum to 21.

For every ordinary degree n, the original homogeneous generators give
an integer coefficient matrix. Reduction modulo a prime cannot increase
its rank, so

    Hilb_Q(Q[x,...,tau]/J,n)
      <= Hilb_Fp(Fp[x,...,tau]/J_p,n)
      <= Hilb_Fp(Fp[x,...,tau]/M_p,n),

where M_p is the leading-monomial ideal just certified. The second
inequality holds because its generators are leading monomials of
actual J_p members. It does not require the selected list to be a
Groebner basis. Thus dim V(J)<=4 and its degree is at most 21. Together
with (2), the entire four-dimensional degree is accounted for at
infinity. Section 3 gives (1) in characteristic zero.

This use of multiplicities is why the previous extra components at
infinity are harmless. It neither assumes that the original affine
scheme is reduced nor replaces it by its radical.

## 7. Returning to the actual group fibers

Off V(I), the original four equations have Jacobian rank four. Their
actual scheme is smooth of dimension three there. The image of V(I)
under any three-dimensional central target map has closure of dimension
at most two. Remove that image and the images of nondominant components,
then apply characteristic-zero generic smoothness to the remaining
smooth source. Its generic transverse fibers are smooth and
zero-dimensional, possibly empty. No regular-sequence hypothesis was
assumed in advance; the Jacobian rank supplies the local codimension.

The central trace changes also require their ramification boundaries.
For a=s+s^-1, the wall s=+/-1 gives the target coordinate
(s-s^-1)f_a=0. The finite quotient (U,V)->(A,B,C) by simultaneous
inversion is etale away from the four points U,V=+/-1. At each such
point both logarithmic derivatives of A,B,C vanish, so both corresponding
central target coordinates vanish. These are proper target subsets.
The full unit-elimination chain rule identifies these derivatives with
the original adjoint-character target coordinates. The central summand
-ABC-4 has been retained for this argument.

On the remaining central cover the changes are etale. Smoothness
therefore transfers to the original central torus. For centralizer-exact
targets, the original group fiber lies in K scheme-theoretically, and
K-regular equals G-regular on that fiber. The smooth regular adjoint
quotient and the unit A2 elimination then lift the transverse statement
to the actual group fiber. Its relative conjugacy dimension is
dim K-rank K=12. Explicit Weyl conjugations transfer it to all eight
standard subsets in this class. See the
[scheme reduction](GENERAL_CAYLEY_FIBER_CRITICAL.md) and
[Weyl-transfer audit](E8_WEYL_TRANSFERS.md).

## 8. Verification, scope, and the next step

- [Sage verifier](tools/verify_e8_last_case.py): exact original generators,
  the rational member (3), all infinity containments and normal orders,
  and full member replays at both primes.
- [Native audit](tools/audit_e8_last_case.py): an independent residual
  reconstruction, rational identities, normal-order argument, and
  coordinate-localization counts of the degree.
- [32003 report](e8-last-32003-verification.json),
  [32009 report](e8-last-32009-verification.json), and
  [native report](e8-last-case-audit.json).
- [Pinned engine builder](tools/build_e8_last_certificate_engine.sh) and
  [small driver](tools/e8_last_f4_driver.c). The builder requires a C
  compiler, GMP/FLINT/MPFR development files, autotools and libtool. It
  installs no software globally. Invoke the Sage verifier with --engine
  pointing to the printed driver path and --prime 32003 or 32009.

The Sp4 exceptional-target control is compatible with the generic
conclusion. The E6 excess example remains allowed. Altered-character
nilpotents remain a counterexample to extending this statement to
arbitrary representations. The new deduction concerns the full adjoint
character; the homogeneous specialization and degree/multiplicity
principles are standard inputs. No new literature-priority claim is made.

The next step is an end-to-end theorem audit and consolidation of the
case dependencies into one statement with exact hypotheses. Its purpose
is to make the completed all-simple result independently reviewable;
success means every reduction and exceptional certificate is linked to
the same actual-scheme theorem. Further Levi examples are unnecessary.
