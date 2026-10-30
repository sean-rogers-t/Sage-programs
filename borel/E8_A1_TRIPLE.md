# E8 adjoint fibers for the cubic A1^3 Levi class

26 October 2026. Characteristic zero. GEO-87 closes the A1^3 class,
represented by retained Bourbaki nodes 125. All 21 standard subsets in
its certified Weyl class have nonempty smooth generic regular-source
fibers of pure dimension six. Current E8 coverage is **190/255**, with
**six unresolved classes representing 65 standard subsets**. The earlier
[Weyl-transfer count](E8_WEYL_TRANSFERS.md) is the preceding checkpoint;
its orbit certificates remain unchanged and are used here.

## 1. The theorem

Let G be the simple simply connected complex group of type E8 and let
theta=theta_Ad be defined by

    B(theta(g),H)=tr(Ad(g) ad(H)),

where B is the Killing form. Let K be a Levi conjugate to the standard
Levi K_{125}, whose derived group is SL2^3. Put

    Z_K^circ={h in z(k) : Z_G(h)=K}.

**Theorem.** There is a dense open U_K in Z_K^circ such that for every
h in U_K the actual scheme theta^(-1)(h) intersect G_reg is nonempty,
smooth, and pure six-dimensional. Its transverse critical quotient is
finite, nonempty and reduced. No generic degree, number of quotient
points or all-special-target smoothness is asserted.

The argument is organized around a polynomial-family proposition, so
its mechanism can be reused without assuming an arbitrary character has
the same form. The full adjoint root decomposition verifies that E8 does.

## 2. The full root-derived potential

Complete alpha_1,alpha_2,alpha_5 to eight mutually orthogonal E8 roots.
In simple-coroot coordinates the columns, ordered X,x,y,t,a,b,c,d, are

    (1,0,0,0,0,0,0,0), (0,1,0,0,0,0,0,0),
    (0,0,0,0,1,0,0,0), (1,1,2,2,1,0,0,0),
    (0,0,0,0,0,0,0,1), (1,1,2,2,2,2,2,1),
    (1,2,2,4,3,2,2,1), (2,3,4,6,5,4,2,1).

Their Gram matrix is 2I_8 and their coroot-lattice index is 16. The last
five span the Lie center of K and define a central torus isogeny. Together
with the first three root SL2 groups this gives a finite central etale
cover SL2^3 x (C*)^5 -> K of degree 16. This cover retains the E8 group
and its character lattice; it is not a change of isogeny type.

Let X,x,y be the fundamental characters of the three derived SL2
factors. For independent central torus variables z_0,...,z_4 put

    (t,a,b,c,d)=(z_0+z_0^-1,...,z_4+z_4^-1),
    B=ab+cd, C=ac+bd, D=ad+bc,
    Q=t^2+a^2+b^2+c^2+d^2+abcd-8.

The full 248-dimensional adjoint character is

    Phi = X^2+x^2+y^2+tXxy+BXx+CXy+Dxy
          +tDX+tCx+tBy+Q.                              (1)

This identity includes all linear and constant terms. One direct root
explanation is that the restriction to the eight commuting A1 factors
has eight adjoint summands and fourteen tensor products of four defining
representations. The eight adjoint characters sum to the eight squares
minus 8. Each of the fourteen root supports contributes one fourfold
product. Their dimensions give 8*3+14*16=248. The verifier reconstructs
every weight, including zero with multiplicity eight; the independent
audit derives the fourteen supports from the Euclidean E8 roots.

## 3. A five-parameter binary proposition

For this section t,a,b,c,d are independent affine parameters, and Q
may be ANY polynomial in those parameters. Let Phi_Q be (1), and let

    Y_Q=V(Phi_X,Phi_x,Phi_y) in A^3 x A^5.

**Proposition.** Outside a closed subset of dimension at most four,
Y_Q is smooth of pure dimension five. In particular any morphism
Y_Q -> A^5 has finite smooth generic fibers, possibly empty. The relative
critical scheme is independent of Q. Nonemptiness requires a separate
dominance argument.

Complete the globally nondegenerate X-square. Set

    L=ty+B, M=Cy+tD, X=-(Lx+M)/2,
    alpha=1-L^2/4,
    beta=Dy+tC-LM/2,
    gamma=y^2+tBy+Q-M^2/4.

The actual critical ideal eliminates X over the whole parameter space;
the residual potential is f=alpha*x^2+beta*x+gamma. Every remaining
derived derivative and central derivative agrees by the chain rule.
The only division so far is by the constant 2.

Apply the [binary criterion](BINARY_BRANCHING_FIBERS.md#1-a-common-criterion)
with s=5 and

    P=4alpha^2 gamma_y-2alpha beta beta_y+beta^2 alpha_y.  (2)

Here are full certificates for its three hypotheses.

### Coprimality of alpha and beta

As a polynomial in y, alpha has coefficients -t^2/4, -tB/2, 1-B^2/4.
Their gcd is one: B is independent of t and 1-B^2/4 is nonzero. Thus
alpha is primitive over C[t,a,b,c,d]. At

    (t,a,b,c,d)=(1,1,0,0,1),

we have B=C=0, D=1 and

    alpha=1-y^2/4, beta=y/2.

The degree of alpha is preserved and the specialized pair is coprime.
Any common positive-y-degree factor would retain its degree, because
the leading coefficient of alpha remains nonzero. A parameter-only
common factor is excluded by primitivity. Hence gcd(alpha,beta)=1.

### A primitive squarefree critical polynomial

The degree of P in y is five, with leading coefficient t^4/2. At
t=a=b=c=d=0, P=8y. Therefore the coefficients of P have gcd one: any
common irreducible parameter factor would have to divide t, which the
second identity excludes. At the same specialization (1,1,0,0,1),

    P=y(y^2-3)(y^2-5)/2.                               (3)

This polynomial is squarefree and its degree is still five. Its nonzero
discriminant proves separability over the parameter function field;
primitivity then proves squarefreeness in the full polynomial ring.
These are polynomial identities over Q, and the argument works over C.

### The entire singular elimination boundary

No part of alpha=0 is discarded without a dimension bound. In fact

    alpha_y=-tL/2,  t*alpha-(L/2)*alpha_y=t,
    (alpha,alpha_y)=(t,B^2-4).                         (4)

This is an equality of ideals, not merely an assertion about reduced
supports. The two underlying boundary components are t=0, B=2epsilon,
where epsilon=+/-1. On either component

    beta_y=D-epsilon*C=(a-epsilon*b)(d-epsilon*c).      (5)

The hypersurface ab+cd-2epsilon is irreducible: it is a primitive linear
polynomial in a, since gcd(b,cd-2epsilon)=1. It divides neither factor
on the right of (5). Thus beta_y is nonzero on each component.

The locus

    E=V(alpha,beta,alpha_y,beta_y,gamma_y) in A^1_y x A^5

therefore has dimension at most three. Indeed (4) has two irreducible
components of dimension four, and beta_y cuts each properly; the
remaining equations can only lower dimension. This includes t=0,
parameter intersections and all possible degeneracies of the polynomial
in x. No division by t, B-2, B+2, or beta_y was used.

### Scheme-theoretic conclusion

On alpha!=0, eliminate x=-beta/(2alpha). The critical scheme becomes
the actual hypersurface V(P), which is reduced of pure dimension five.
Its nonsmooth locus has dimension at most four. On alpha=0, beta=0
as well. Coprimality gives dim V(alpha,beta)<=4 in the six-dimensional
(y,parameter) space. Outside E the second critical equation is a nonzero
polynomial in x, so its fibers are finite or empty; over E there is at
most one free x. Hence the entire boundary has dimension at most four.

Remove the closures of the target images of this boundary and the
nonsmooth locus. Their dimensions are at most four in the target A^5.
Over a further dense target open, characteristic-zero generic smoothness
applies to the smooth five-dimensional remainder, giving zero-dimensional
smooth fibers. A finite-type zero-dimensional fiber is finite. This
argument excludes every bad boundary support from generic target fibers,
including any embedded or nonreduced structure supported there; it never
replaces the original critical scheme by its reduction.

This proves the proposition. It does not assume the original equations
form a regular sequence, nor does it claim all central-parameter fibers
are finite or reduced.

## 4. Torus ramification and actual Springer targets

Return to the actual Q in Section 2. Pull the proposition back under the
finite map (C*)^5 -> A^5 given by the five traces. It is etale away from
the five divisors z_i=+/-1. On the entire relative critical scheme,

    mu_i=z_i Phi_{z_i}=(z_i-z_i^-1) Phi_{v_i},
    (v_0,...,v_4)=(t,a,b,c,d).                          (6)

Consequently every ramification divisor has target image in mu_i=0.
Exclude those five target hyperplanes. No source divisor is silently
deleted; its whole target image is controlled, regardless of its scheme
structure or dimension. Finite pullback preserves the other bad-locus
dimension bounds. The smooth locus away from these boundaries is still
five-dimensional, so generic target critical quotients are finite smooth
schemes, possibly empty at this stage.

The actual target coordinates are c_i=B(h,H_i), where H_i are the five
central coroot columns above. They are independent linear coordinates
on z(k). The finite-index basis causes no loss over C. The original
central derivatives, not newly chosen independent coefficients B,C,D,
are the functions mu_i in (6).

The [regular-quotient critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md)
identifies the actual Springer fiber on K_reg with these equations,
plus dim K-rank G=14-8=6 smooth orbit directions. For h with centralizer
K, equivariance gives scheme-theoretic containment theta^(-1)(h) in K,
and regularity in K is equivalent to regularity in G on that fiber.
The quotient map on K_reg is smooth and surjective. Smoothness descends
through the finite central etale cover, with all nilpotents retained.

## 5. Dominance and nonemptiness

At the point

    X=x=y=-1, z_0=...=z_4=1,
    equivalently t=a=b=c=d=2,

the three relative derivatives vanish. Their Hessian is

    [2 6 6]
    [6 2 6],     determinant 224.
    [6 6 2]

The relative incidence is therefore smooth locally, with central torus
parameters as local coordinates. At this point

    (Phi_t,Phi_a,Phi_b,Phi_c,Phi_d)=(-21,6,6,6,6).

All first derivatives of the five traces with respect to logarithmic
torus coordinates vanish, so the mixed derived/central Hessian block
vanishes. The derivative of the central target map on the local critical
branch is exactly

    diag(-42,12,12,12,12),   determinant -870912.

It is invertible, proving dominance. The witness lies over target zero
and on trace ramification; it is used ONLY to prove dominance. The
resulting dominant component meets the complement of those divisors
and of the previously bounded bad loci. Its image contains a dense
target open, which intersects Z_K^circ. Thus the generic quotient is
nonempty, and the regular-quotient presentation proves the theorem.

## 6. Verification, known inputs, and the remaining work

The [Sage certificate](e8-a1-triple-verification.json), generated by
[this verifier](tools/verify_e8_a1_triple.py), checks the full character,
the coroot cover, all chain rules, the binary identities and both Hessians.
The [independent audit](e8-a1-triple-audit.json) uses
[native Euclidean roots and SymPy](tools/audit_e8_a1_triple.py), without
Sage or imports from the verifier. It independently recovers the
orthogonal completion and fourteen supports, proves the exact boundary
ideal identity, and checks the squarefree specialization and dominance.
There are no modular or numerical tests in these certificates.

Published inputs are generic smoothness in characteristic zero (for
example [Vakil, Theorem 3.3](https://virtualmath1.stanford.edu/~vakil/0506-216/216class5152.pdf)),
standard root-system and adjoint-quotient facts, and the equivariant
trace-projection construction. The binary criterion and regular-quotient
presentation were already established in this repository from those
inputs. This continuation verifies their hypotheses for the full E8
character; it is a repository deduction, not a claim of literature
novelty. The proof of the geometric implications is written above;
computer identities alone are not a smoothness theorem.

The standard controls remain intact: Sp4 shows why special targets cannot
be included automatically; the E6 excess example shows this finite-quotient
mechanism is not universal; altered characters must satisfy the actual
coprimality, squarefreeness and boundary hypotheses, which Weyl symmetry
or nonlinear coefficients alone do not imply.

The existing 21 [Weyl words](e8-weyl-transfers-verification.json) transfer
this theorem to the entire A1^3 class. Thus 169+21=190 standard E8 subsets
are covered. The remaining classes are A2+A1^2 (28), A1^4 (7), A2^2+A1 (8),
A2+A1^3 (8), A3+A1^2 (10), and A3+A2+A1 (4): six classes and 65 subsets.

**Next target:** A2+A1^2, represented by 1235. Its paired A2 elimination
again leaves two affine variables. The goal is to reuse the binary
mechanism with its full E8 coefficients and prove every boundary-image
bound. Success would cover 28 subsets; a failed hypothesis must be
recorded explicitly rather than inferred away from a generic sample.
