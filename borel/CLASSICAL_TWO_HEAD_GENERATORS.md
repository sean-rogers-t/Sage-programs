# Three additional generators when the smaller block has size two

**29 October 2026: written integral deduction.** This proves the
short-presentation conjecture of
[the general denominator note](CLASSICAL_SATURATION_GENERATORS.md)
for ell=2 in every rank. It concerns the normalized B image and the
D square image, not the full Euler source. It uses the earlier uniform
prefix and geometric comparison theorems; it is not an external review
or a claim of literature priority.

The smaller block is placed at the head using the integral
[block duality](TYPE_B_SMITH.md#4-block-duality-including-integral-lattices).
This transports the source subring and its saturation through full-flag
pullback and an integral Weyl-group automorphism; it does not claim
that the individual signed-cell prefixes are unchanged by block exchange.

## 1. Statement in the source ring

Let A>=2, n=A+2 and c=ceil(n/2). Write

    R=Z[u,v]/(h_(n-1),h_n),
    u=y_1+y_2, v=y_1 y_2,
    p_j=y_1^j+y_2^j.

The h_j are complete symmetric functions in y_1,y_2, with h_0=1,
h_1=u and h_j=u h_(j-1)-v h_(j-2). The power sums satisfy the same
recurrence, starting p_0=2,p_1=u. Source weight is |y_i|=1;
phi(y_i) is a squared root, so source weight m has cohomological
degree 4m. Interpret fractions below inside R_Q using the injective
normalized map phi.

Let K_e=ker(R->H*(W_(e-1);Z)), as in the general note. Then

    K_1=(p_c,p_(c+1),v^c),
    K_2=(v^c).                                               (1)

The generator v^c in K_1 is redundant after the Grassmannian relations,
but retaining it makes the geometric high-monomial description transparent.
Consequently the full saturation is the following R-module:

    S = R + R z + R w + R t,
    z=p_c/2, w=p_(c+1)/2, t=v^c/4.                           (2)

These are actual integral classes after phi, not merely rational
expressions. Thus three additional classes suffice to generate S as
an algebra over the original image, in every A>=2. This set is minimal
among homogeneous algebra generating sets over R: the three necessary
source weights are c, c+1 and 2c, as proved below.

## 2. Reduce the high-monomial ideal to three elements

All high one-part monomials are p_j for j>=c. Their recurrence expresses
them in the ideal (p_c,p_(c+1)). Every high two-part monomial is v^c
times an ordinary symmetric polynomial: for a>=b>=c it is

    m_(a,b)=v^b p_(a-b) if a>b, and m_(b,b)=v^b.

Thus J_1=(p_c,p_(c+1),v^c) and J_2=(v^c). The general inclusion theorem
gives J_e subset K_e. The top-layer theorem already gives J_2=K_2.
It remains to prove J_1=K_1 integrally.

## 3. Eliminate the Grassmannian relations exactly

Put S0=Z[u,v]. The following polynomial identities hold before any quotient:

    h_(2c-1)=p_c h_(c-1),
    h_(2c)=p_c h_c-v^c,
    h_(2c-2)=p_c h_(c-2)+v^(c-1).                           (3)

For example, multiply h_j=(y_1^(j+1)-y_2^(j+1))/(y_1-y_2)
by p_c and cancel the two cross terms. This is a polynomial identity;
no localization at y_1-y_2 is required.

Using p_(c+1)=u p_c-v p_(c-1), (3) gives

    R/J_1 = S0/(p_c, v p_(c-1), v^c)       if n=2c,
    R/J_1 = S0/(p_c, v p_(c-1), v^(c-1))   if n=2c-1.       (4)

In the odd case c>=3 under A>=2. These are equalities of integral
ideals, not just equalities of radicals.

## 4. Monic integral straightening

Set g_j=v^j p_(c-j), 0<=j<c. In the even case a monic Groebner basis
for (4), in lexicographic order u>v, is

    g_0,g_1,...,g_(c-1), v^c.

In the odd case it is

    g_0,g_1,...,g_(c-2), v^(c-1).

Here is the uniform verification over Z. For j<=c-3,

    v g_j-u g_(j+1)=-g_(j+2).

This recurrence both generates every listed g_j from the first two
and reduces the adjacent S-polynomials. At the even boundary the
next remainder is -2v^c, and the final adjacent S-polynomial is zero.
At the odd boundary, -g_(c-1)=-u v^(c-1) is already divisible by the
last generator, and the S-polynomial of g_(c-2) and v^(c-1) is
-2v^c. All reduce integrally. Nonadjacent pairs satisfy the usual
chain criterion: intermediate leading monomials divide the pair's
least common multiple. Equivalently their S-polynomials telescope
through the adjacent ones. All leading coefficients are one, so no
prime is inverted in this argument.

The standard monomials are

    u^a v^b with a,b>=0 and a+b<c,

with just v^(c-1) removed in the odd case. In particular R/J_1 is
free over Z. Its weight Hilbert polynomial is

    [c+1 choose 2]_q                      if n=2c,
    [c+1 choose 2]_q-q^(2c-2)            if n=2c-1.          (5)

Indeed u^a v^b corresponds to the partition (a+b,b) in the
(c-1)-by-2 rectangle, with weight a+2b.

## 5. Identify the exact kernel, then multiply the new classes

The established defect theorem gives

    rank_q(R/K_1)=[c choose 2]_q + q^(c-1)[floor(n/2) choose 1]_q.

The q-Pascal identity identifies this with (5), in both parity cases.
The surjection R/J_1 -> R/K_1 therefore has a torsion kernel in each
degree. But R/J_1 is free, so that kernel is zero. This proves (1).
The denominator theorem S=R+(1/2)K_1+(1/4)K_2 proves (2).

The additional classes satisfy the explicit product identities

    z^2=2t,       zw=u t,       w^2=2v t,
    zt=wt=t^2=0,
    2z=p_c,       2w=p_(c+1),   4t=v^c.                    (6)

To check these, use y_i^n=0 in the ordinary flag splitting ring and
2c>=n. For instance p_c^2=2v^c and p_c p_(c+1)=u v^c,
since terms assigning two high exponents to one root vanish. The
remaining identities follow in the same way, over Z before taking
the justified fractions.

Equations (6) are useful relations, but are not asserted to be a
complete defining ideal over R: there are also R-linear relations
among 1,z,w,t. Formula (2) is a complete module-generation theorem.
It supplies multiplication by reducing source polynomials in R_Q.

### Minimality of the three homogeneous generators

There are no additional saturation classes below weight c. The element
p_c=s_(c)-s_(c-1,1) is primitive in R, since c<=A, so z is not in R.
At least one new homogeneous generator is therefore needed in weight c.

In weight c+1, the coefficient of s_(c,1) in p_(c+1) is -1, whereas
its coefficient in u p_c is zero by Pieri cancellation. This partition
fits the rectangle, including n=4 and n=5. Hence
p_(c+1) is not in 2R+Z u p_c. It follows that w cannot be generated
by R and a weight-c generator: a second generator is necessary in
weight c+1. Changing z by an integral source class does not change
this parity obstruction.

Finally, v^c=s_(c,c) is primitive in R, so t is not in (1/2)R.
Every element of weight 2c generated by R,z,w lies in (1/2)R:
the only possible product of two new positive-weight generators in
that weight is z^2=2t=v^c/2. All other terms are R-multiples of z or w.
Formula (2) shows that every class of weight below 2c is already in
R[z,w]. Thus t cannot be generated by any lower-weight classes.
A third generator is necessary in weight 2c. Together with (2), this
proves the claimed minimal number and the three required degrees.

## 6. Known inputs, checks and remaining work

Published algebraic inputs are the ordinary Grassmannian presentation,
two-variable symmetric-function recurrences, the monic Groebner-basis
criterion and q-Pascal. The identification with the normalized Springer
image additionally uses our earlier integral prefix theorem. The proof
does not infer integral equality from the Hilbert polynomial alone:
the monic calculation in Section 4 is the essential torsion argument.

The independent finite ideal tests already include two smaller-head-two
cases. A separate symbolic/straightening check is recorded in
[the two-head audit](classical-two-head-generator-audit.json).
Those checks corroborate the identities and reduction rules; Sections
2--5 are the arbitrary-rank proof. Literature novelty is unestablished.

The remaining short-presentation problem has ell>=3 and 1<=e<ell.
The exterior-power Toeplitz formulation in the general note isolates
its integral image question. The next useful step is a higher-variable
monic straightening rule, rather than another list of Smith factors.
