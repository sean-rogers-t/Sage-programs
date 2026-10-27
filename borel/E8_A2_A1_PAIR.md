# E8 adjoint fibers for the A2+A1^2 Levi class

26 October 2026. GEO-88 closes the class represented by retained Bourbaki
nodes 1235. Its 28 standard subsets have nonempty smooth generic actual
regular-source fibers of pure dimension ten. Current E8 coverage is
**218/255**, leaving **five classes and 37 standard subsets**. The earlier
[A1^3 result](E8_A1_TRIPLE.md) and [Weyl-transfer table](E8_WEYL_TRANSFERS.md)
are preceding checkpoints; their certificates are retained unchanged.

## 1. The theorem and what it does not assert

Let G be the simple simply connected complex group of type E8. Use its
full adjoint representation and Killing-form normalization

    B_Ad(theta(g),H)=tr(Ad(g) ad(H)).

Let K be conjugate to K_{1235}, with derived group SL3 x SL2 x SL2,
and let Z_K^circ={h in z(k): Z_G(h)=K}.

**Theorem.** There is a dense open U_K in Z_K^circ such that, for every
h in U_K, the actual scheme theta^(-1)(h) intersect G_reg is nonempty,
smooth and pure ten-dimensional. The transverse critical quotient is
a finite nonempty reduced scheme. All 28 standard subsets in the
certified Weyl class are covered.

We do not determine its generic length, the full Springer-map degree,
or the singularities at every special central target. This is not the
remaining all-E8 theorem and is not a literature-priority claim.

## 2. A compact full-character formula with justified coordinates

Let a_0,b_0 be the two fundamental SL3 characters and x,y the SL2
characters. The original central torus variables p,q,r,v correspond
to these simple-coroot columns:

    (10,15,20,30,24,18,12,6), (6,9,12,18,15,12,8,4),
    (4,6,8,12,10,8,6,3), (2,3,4,6,5,4,3,2).

Together with e1,e2,e3,e5 these give the central etale cover of K of
degree 12. Pass to the finite etale common cover p=s^2 and set

    A=p*a_0, B=b_0/p, u=s^3,
    w1=s^3*q*r*v, w2=s^3*q*r, w3=s^3*q,
    w4=s^-9*q^-3*r^-2*v^-1,       product(wi)=1,
    S=sum(wi), T=sum_{i<j}(wi*wj), V=sum(1/wi).

The change (a_0,b_0)->(A,B) is an affine coordinate isomorphism over
the whole central torus. The map (s,q,r,v)->(u,w1,w2,w3) is a finite
etale isogeny of degree three, with exponent matrix

    M=[3 0 0 0]
      [3 1 1 1]
      [3 1 1 0]
      [3 1 0 0],       det M=-3.                        (1)

Thus the two torus models have a common finite etale cover, of degrees
two and three respectively. No fractional character is silently
introduced on the original group.

In these coordinates the FULL adjoint character is

    Phi=AB+A*F+B*G+H,
    F=xy+u^-1*S*x+u^-1*V*y+T+u^-2,
    G=xy+u*V*x+u*S*y+T+u^2,
    H=x^2+y^2+Txy+(uS+u^-1 V)x+(u^-1 S+uV)y+SV-3.      (2)

The 248-dimensional root character is independently reconstructed
before testing this identity. In particular all linear terms in A,B,x,y
and the constant term are retained. S,T,V are independent invariant
coordinates for an SL4 torus; they are not independent torus coordinates.
Their ramification is addressed separately in Section 5.

The paired equations Phi_A=B+F, Phi_B=A+G give exact elimination
A=-G, B=-F, with residual potential psi=H-FG. The only inverted
quantities are Laurent units. All remaining derivatives agree by
the chain rule. In the initial normalization,

    D_p(chi) at fixed a_0,b_0
      =D_p(Phi) at fixed A,B + A*Phi_A-B*Phi_B.

Consequently the extra terms vanish on the actual relative critical
scheme. The target coordinates on the common cover obey

    diag(2,1,1,1)*c_original = M^tr*c_model.              (3)

This is an invertible linear change over C. It preserves genericity,
target dimensions and the actual fiber schemes through the etale covers.

## 3. Binary potential and squarefreeness

Work first over the smooth parameter space Z=C*_u x A^3_(S,T,V), of
dimension four. Put

    k=uV+u^-1 S, m=uS+u^-1 V,
    c0=1-SV, d0=u^2+u^-2.

Then psi=alpha*x^2+beta*x+gamma, where

    alpha=c0-y^2-ky,
    beta=-my^2-(T+d0+S^2+V^2)y-Tk,
    gamma=c0*y^2-Tmy+SV-4-T^2-d0*T.                    (4)

We verify all three hypotheses of the
[binary criterion](BINARY_BRANCHING_FIBERS.md#1-a-common-criterion),
with coefficient ring R=C[u,u^-1,S,T,V].

**Coprimality.** Alpha has leading y coefficient -1. At
(u,S,T,V)=(1,0,2,0), alpha=1-y^2 and beta=-4y. These are coprime,
and the degree of alpha is preserved. A common positive-y-degree
factor would survive that specialization; a parameter-only factor is
impossible because alpha is monic up to a unit. Thus gcd(alpha,beta)=1.

**Primitivity.** Define

    P=4alpha^2 gamma_y-2alpha beta beta_y+beta^2 alpha_y.

Its degree in y is five, with leading coefficient

    [y^5]P=2*((uS-u^-1 V)^2+4).                        (5)

This is independent of T. Any irreducible parameter content must
therefore be independent of T as well. Two coefficients give

    [T^2 y^2]P=-k,    [T^2 y]P=-2*(k^2+c0).             (6)

The content must divide k and c0. But k is a unit multiple of the
primitive linear polynomial S+u^2 V, while c0 modulo that polynomial
is 1+u^2 V^2, which is nonzero. Hence gcd(k,c0)=1 and P is primitive.
This argument is over C, not just a rational-factorization test.

**Squarefreeness.** At the same specialization (1,0,2,0),

    P=8y(y^2-3)(y^2+1).                               (7)

Its five roots are distinct, and degree five is preserved. This proves
separability over Frac(R); primitivity then proves squarefreeness in R[y].
The specialization is a polynomial certificate, not an assertion about
a generic Springer fiber's degree or about torus ramification there.
The leading coefficient (5) is never inverted: its zero locus and all
its intersections remain in the reduced hypersurface calculation.

## 4. The full elimination boundary, including its intersections

The exceptional base in the binary criterion is

    E=V(alpha,beta,alpha_y,beta_y,gamma_y) in A^1_y x Z.

Since alpha_y=-2y-k, set y0=-k/2. Exact substitution gives

    4alpha(y0)=Delta=(uV-u^-1 S)^2+4,
    beta(y0)=k*(S^2+V^2+d0*(2-SV)-2T)/4.               (8)

The two irreducible components of Delta=0 over C have equations

    S=u^2 V-2epsilon*i*u,       epsilon=+/-1.

Each is a smooth three-dimensional space with coordinates u!=0,V,T,
and y is determined. The coefficient of T in beta(y0) on it is

    -k/2=-uV+epsilon*i,                                (9)

which is not identically zero on either component. Thus beta=0 cuts
each to dimension at most two, and dim E<=2 as required.

In particular k=0 has NOT been divided out. On that sublocus
V=epsilon*i/u and S=-epsilon*i*u; beta(y0)=0 automatically, and u,T
remain free. This is a two-dimensional sublocus already included in
the estimate. Further equations beta_y=gamma_y=0 can only lower its
dimension. Both components and all additional parameter intersections
are retained; the two Delta components themselves are disjoint.

For clarity, the dimension argument controls the WHOLE actual boundary.
On alpha=0 the first critical equation also gives beta=0. Coprimality
implies dim V(alpha,beta)<=3 in the five-dimensional (y,parameter)
space. Outside E the remaining equation in x is a nonzero polynomial,
so fibers are finite or empty. Over E at most one free x remains.
Therefore the full critical boundary has dimension at most three.

On alpha!=0, eliminate x=-beta/(2alpha). The actual scheme is V(P),
a reduced hypersurface of pure dimension four with nonsmooth locus
of dimension at most three. Images of these two bad loci cannot
dominate a four-dimensional central target. This excludes every
embedded or nonreduced structure supported on the boundary from the
generic target fiber; it does not take the radical of the original ideal.
No assumption that the original equations form a regular sequence
is needed.

These source estimates are unchanged by adding any central Laurent
polynomial Q(u,S,T,V) to psi. Such a change need not preserve dominance;
the nonemptiness argument below uses the actual full adjoint character.

## 5. Every torus-quotient ramification wall is excluded correctly

The map (w1,w2,w3,w4), product(wi)=1, to (S,T,V) is the finite SL4
Weyl quotient. Its logarithmic Jacobian determinant is, up to sign,

    product_{i<j}(wi-wj).                              (10)

It is therefore etale away from exactly six root walls wi=wj. Finite
pullback preserves the dimension estimates in Section 4, but smoothness
can only be pulled back etale-locally; the walls must be controlled.

The twelve ratios wi/wj, i!=j, are exactly the twelve E8 roots whose
pairings with alpha_1,alpha_2,alpha_3,alpha_5 vanish. They form the
commuting A3 root subsystem. This identification is verified in the
original cocharacter basis, including each root's coroot direction.
Its coroot acts by wi*d/dwi-wj*d/dwj, fixing u,A,B,x,y.

For EVERY point, before imposing the critical equations, symmetry gives
the exact Laurent identity

    D_(epsilon_i-epsilon_j) Phi
      =(wi-wj)*(Phi_S + sum_{l not in {i,j}} wl*Phi_T
                          +Phi_V/(wi*wj)).              (11)

Thus on wi=wj the target component B_Ad(h,alpha^vee) is zero. These
are roots outside the Levi. For h in Z_K^circ, alpha(h) cannot be zero,
so B_Ad(h,alpha^vee) is nonzero. The actual fiber has no point on any
of these walls. This also excludes all their intersections and any
nilpotent structure supported there. Equivalently, the wall equation
and the prescribed nonzero target component generate the unit ideal
by (11). No ramified divisor is merely deleted from the source.

On the remaining etale open, the source is smooth of dimension four
away from loci of dimension at most three. After excluding their
nondense target images, characteristic-zero generic smoothness gives
finite smooth generic target quotients. This is the standard generic
smoothness theorem, applied only after establishing the full boundary
and ramification estimates.

## 6. Dominance, nonemptiness and passage back to G

Use the point

    A=B=0, x=y=-1, u=w1=w2=w3=w4=1,
    hence S=V=4, T=6.

All four derived derivatives of Phi vanish. In (A,B,x,y) their Hessian
is

    [0 1 3 3]
    [1 0 3 3]
    [3 3 2 6]
    [3 3 6 2],       determinant -112.

The residual (x,y) Hessian is [[-16,-12],[-12,-16]], of determinant 112.
For the logarithmic u-direction, the full mixed column is (6,-6,0,0)
and the full second derivative is -16. The Schur complement is 56.
The invariant derivatives are (Phi_S,Phi_T,Phi_V)=(2,1,2).
At wi=1 their first logarithmic derivatives vanish, while their Hessians
in (log w1,log w2,log w3) are respectively C,2C,C, with C=I_3+J_3 and
J_3 the all-ones matrix. Thus the central target derivative on the local
critical branch is

    diag(56,6C),       determinant 48384.               (12)

It is invertible, proving dominance. The witness is over target zero
and on the ramification locus, and is used ONLY for dominance. The
dominant component then meets the complement of the excluded target
hyperplanes and bad images, and its image contains a dense open.
Intersecting with Z_K^circ gives nonempty smooth finite quotients.

The [regular-quotient critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md)
now applies to the actual Springer scheme. For centralizer-exact h,
equivariance places the fiber in K and makes K-regularity equivalent
to G-regularity on that fiber. The quotient map on K_reg is smooth
and surjective, with relative dimension dim K-rank G=18-8=10. The
common torus isogenies and the central group cover are finite etale;
smoothness descends, with the actual scheme structure retained.
This proves the theorem, including regular nonsemisimple source points.
No derived torus discriminant is inverted: A,B,x,y range over the whole
derived character quotient. Only the separately controlled central
ramification walls are excluded by centralizer-exact targets.

## 7. Independent checks, known inputs, and what remains open

The [Sage verifier](tools/verify_e8_a2_a1_pair.py) constructs the full
restricted character from E8 roots and records the
[exact certificate](e8-a2-a1-pair-verification.json). The
[native audit](tools/audit_e8_a2_a1_pair.py) reconstructs all weights
using independent Cartan reflections and subset models of fundamental
characters, then checks the algebra in SymPy. Its
[report](e8-a2-a1-pair-audit.json) covers both complex elimination
boundary components, the k=0 subloci, all six actual root walls, their
coroot directions, the Vandermonde identity and the dominance Hessian.
All computations are exact in characteristic zero. The written argument,
not a dimension count at sample parameters, supplies the scheme geometry.

Published inputs are the equivariant trace-projection construction of
[Kostant and Michor](https://www.mat.univie.ac.at/~michor/kostant.pdf),
standard SL4 invariant theory, regular adjoint quotients, and
[generic smoothness in characteristic zero](https://virtualmath1.stanford.edu/~vakil/0506-216/216class5152.pdf).
The binary criterion and regular-quotient reduction were already in this
repository. The advance here is verification of their hypotheses for
the FULL E8 character, including both kinds of boundary. No claim of
literature novelty or resolution of a published open problem follows.

The controls remain consistent: special central targets may be singular
(Sp4); positive excess can occur in other families (E6); and altered
characters cannot inherit this conclusion without the full coefficient
and boundary identities. No universal quasifiniteness claim is made.

The [certified Weyl words](e8-weyl-transfers-verification.json) transfer
the proof to all 28 standard subsets in this class. E8 coverage becomes
190+28=218/255. Still unresolved are A1^4 (7 subsets), A2^2+A1 (8),
A2+A1^3 (8), A3+A1^2 (10), and A3+A2+A1 (4): five classes and 37 subsets.
The generic degree and all-special-target behavior of the present class
also remain undetermined.

**Next target:** A1^4, represented by 1257. Its residual problem has
three affine variables, so the binary argument does not immediately
apply. The goal is to use the already reconstructed root character to
seek a symmetric quadratic-block criterion with all rank strata
controlled, thereby testing how much of the method extends beyond
binary elimination. A complete proof would cover its seven subsets.
