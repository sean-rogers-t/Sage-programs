# The next elementary layer, including the root boundary

28 October 2026. IMG-43. Written deduction, with internal independent
checks; no external-review or literature-priority claim.

## 1. Statement and scope

Retain the full literal defining SL_m/P1 root-polynomial source and
K-theory forget-equivariance defect Delta_m of
[IMG-41](TYPE_A_K_FIRST_CORRECTION.md). The grading is transported from
the homogeneous model. Let p>=5 divide m and put d=p+3,
s=v_p(m), and a_d=sum_(k>=2) P_k(d-2^k+1), where P_k counts
partitions into parts at most k. Then

    (Delta_m)_(p+3),(p) = (Z/p^s)^a_(p+3),  m>p,          (1)
    (Delta_p)_(p+3),(p) = (Z/p)^(a_(p+3)-2),  m=p.        (1r)

For m>p, the first shortened factor in weight p+2 is followed by
restoration of the partition multiplicities in weight p+3. This is a statement about
the actual elementary source, not just its moment subring. When m>p,
we have m>=2p>d, so no root equation or coefficient truncation enters
this weight.

The same proof applies to the auxiliary unrooted elementary model with
independent coefficients through d and mass m=p. It does NOT identify
that auxiliary model with SL_p/P1. The separate root calculation in
section 6 proves (1r). The prime
p=3 is excluded: this weight reaches 2p, where two Witt factors occur.

The proof uses the complete relative moment theorem of
[IMG-42](TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md). In particular, source base
change is not injective in this weight, and the argument below retains
its common source/ideal torsion rather than assuming it away.

## 2. Exhaustive reduction of the first Witt enlargement

Work over Z_(p). Use the moments a_i,q_i, b=(a_p-a_1^p)/p,
g=(q_p-q_1^p)/p=b+H, and alpha=(m-m^p)/p from IMG-41.
Let J_0,C_0 be the actual mass-m moment source and coefficient ideal,
after quotienting out scalar coefficient polynomials. All remaining
ambient coordinates have positive L-degree.

Because d<2p, every new word contains exactly one b or g, multiplied by
a word of weight three in a_1,a_2,a_3,q_1,q_2,q_3. Other coordinate
changes have p-unit diagonal entries by Newton identities. The only
new directions outside the old integral ambient lattice are

    V_3=b(q_3-a_3),       a_1 V_2=a_1 b(q_2-a_2).

They are independent of order p modulo that lattice: their coefficients
on a_p a_2 L and a_p a_1^2 L are respectively 3/p and 2/p, with no
cross contribution. Their p-multiples are in C_0, and both directions
occur in the enlarged coefficient ideal. Cancelling them therefore
changes neither the old relative quotient nor its embedded classes.

Put

    v_A=m b a_2 L,       v_B=m b a_1^2 L,    v_C=m b L q_2,
    w_3=a_3 H,          w_12=a_1 a_2 H,     w_111=a_1^3 H,
    w_1q2=a_1 q_2 H,    w_2q1=a_2 q_1 H,   w_11q1=a_1^2 q_1 H,
    E_3=q_3 H-alpha q_d,
    E_12=q_1 q_2 H-m alpha q_d.

The new source is J_0 plus these eleven directions; the new ideal is
C_0 plus the first nine. This enumeration is complete: the weight-three
words are a_3,q_3; the four products of a_1 or q_1 with a_2 or q_2;
and the four cubics in a_1,q_1. Apply

    q_1^2=m q_2+a_1^2-m a_2

to reduce the words containing q_1^2. Scalar products disappear. This
gives exactly the displayed list and the two cancelled directions.

## 3. Reduce to five independent new directions and one old class

Set Z=q_1 H+m b L. IMG-41's even quadratic identity proves Z is in
J_0 in weight p+1. Its weight-p+2 comparison also proves that
m b a_1 L, a_1^2 H, and a_1 q_1 H are in J_0. Consequently

    v_B, w_111, w_11q1 are in C_0,
    v_A+w_2q1=a_2 Z is in C_0,
    E_12+v_C=q_2 Z-m alpha q_d is in J_0.                  (2)

Write R_3=q_3 q_p-m q_d and K=w_1q2+w_2q1. The integral difference-
trace certificate, now for the pair (3,p), gives

    (m/p) R_3 = -3 K                         modulo C_0. (3)

Here is the division audit. In the formal trace identity for M R_3,
reduce after division by p and specialization M=m. Terms containing M
are integral coefficient-ideal words since p|m. Interior coefficients
binomial(p,j) are divisible by p. The only surviving terms are

    -3(a_1 q_2+a_2 q_1)(q_p-a_p)/p.

Replace (q_p-a_p)/p by H+(q_1^p-a_1^p)/p. Every interior binomial term
in the latter expansion is in C_0 using m L=q_1-a_1. Its final term is
(m/p) times the specialization of

    M^(p-1) L^p (a_1 q_2+a_2 q_1).

This has coefficient degree p+1, weight p+3, and a positive coefficient
factor. IMG-42's filling recurrence, with excess two, puts it in IJ
over Z[1/6], hence over Z_(p). This proves (3) without dividing an
unknown torsion class. In particular K is in J_0 because 3 is a p-unit.

It follows that the new lattices reduce to

    J'=J_0+<v_A,v_C,w_3,w_12,E_3>,
    C'=C_0+<v_A,v_C,w_3,w_12,K>.                           (4)

All five added source classes have order at most p modulo J_0. For the
first four, their p-multiples are already in C_0. For the fifth,

    p E_3=R_3 modulo C_0.                                 (5)

Indeed the other term is q_3 q_1^p-m^p q_d, which vanishes in the old
defect using q_1^2=m q_2 and the proved m-annihilator.

## 4. Independence despite noninjective source base change

Let C=(A/J)_d be the VARIABLE-mass moment ambient/source quotient over
Z_(p)[M]. IMG-42 proves that its entire p-torsion T is one Z/p in
coefficient degree three, killed by M. The quotient C/T is free over
Z_(p), with nilpotent M. Thus ker(M-m:C -> C)=T.

The following elementary version of the Bockstein argument suffices.
For classes c_i in C/(M-m)C, choose lifts with p c_i=(M-m)b_i in C.
If sum r_i c_i=0, then sum r_i b_i is in pC+T. Therefore independence
of the b_i modulo pC+T proves independence of the c_i of order p.
This statement does not assume M-m injective on C.

For the five directions in (4), the lowest coefficient-degree pieces
of b_i modulo p, up to signs, are respectively

    a_p a_2 L,
    M a_p L^3+2a_1 a_p L^2+a_2 a_p L,
    a_3 L^p,
    a_1 a_2 L^p,
    a_p L^3+3a_(p+1)L^2+3a_(p+2)L.                       (6)

These formulas follow directly by subtracting the variable-M and mass-m
expressions. Modulo p the divided difference
((M-M^p)-(m-m^p))/(M-m) is 1-M^(p-1). All other such terms carry an
interior binomial(p,i), so disappear modulo p. The terms involving
a_1^p or M^(p-1) have higher coefficient degree and do not affect (6).
One can also remove them by the high-coefficient recurrence. For E_3,
Lucas' binomial identity gives the last line after discarding a_d,
which is in J.

The third and fifth lines are independent in coefficient degree one:
after removing a_d and M L^d, the a_i L^(d-i), 1<=i<d, form a basis.
For the other three use the explicit quadratic quotient of IMG-42.
Polarize in x_1,x_2, put y=(x_1+x_2)/2, v=(x_1-x_2)/2,
z=L+y, w=v^2, and r=(p+1)/2. Remove the coordinate subspaces in
Z_(p)[y,w] and Z_(p)[z,w], then the distinct unit-pivot polynomials

    P_i=even_v((z+v)^i(y-v)^(d-i)),  1<=i<d.

The remaining quotient has the monomial basis

    y^a z^b w^c,   a,b,c>=1,  a+b+2c=d.

Using (y+v)^p=y^p+v^p in characteristic p, the first, second and fourth
lines of (6) reduce respectively to

    -4 y z w^r,       -4 y^p z w,       -4 y z^p w.       (7)

For example the reductions subtract P_1, P_3-P_2, and P_p,
respectively; any pure y,w or z,w terms have already been removed.
These are three distinct basis monomials for every p>=5. None of the
five low-degree classes can meet T, which is in coefficient degree
three. This proves the required independence modulo pC+T and hence
all five order-p assertions in the actual mass-m ambient quotient.

## 5. The one killed factor is restored

The first four independent directions in (4) are shared by source and
ideal. As their p-multiples are in C_0, they create no further intersection
with J_0. Thus the old defect loses exactly the subgroup generated by K,
and the remaining source extension is E_3 with relation (5).

R_3 is a primitive direction in the old moment defect. Here is a check
which avoids assuming a factorial pairing is perfect. In coefficient
degree two the defect is free over Z_(p), and is generated by

    R_i=q_i q_(d-i)-M q_d,  1<=i<=d/2.

It has rank d/2-1. Its one relation is the alternating binomial relation

    sum_(1<=i<d/2) (-1)^i binomial(d,i)R_i
                  +(-1)^(d/2)binomial(d,d/2)R_(d/2)/2=0.

The R_1 coefficient is -d, a p-unit. Eliminating R_1 leaves an integral
basis containing R_3. Saturation and the rank calculation show this
primitive relation is the entire kernel. The other coefficient degrees
are direct summands of the universal defect and remain so after M=m.

Equation (3) therefore removes exactly the order-p subgroup of the
cyclic factor generated by R_3. Equation (5) adjoins its p-th root E_3,
and

    m E_3=(m/p)R_3=-3K=0

in the new quotient. This restores that factor to order p^s. All other
factors remain unchanged, proving (1). The ordinary image is unchanged:
all added directions in (4) have ordinary value zero.

## 6. The root boundary: one coefficient class and one root class

For m=p, first delete c_(p+1),c_(p+2),c_(p+3) from the auxiliary model.
The only non-scalar extra jets in this weight are

    f_(p+2)=-c_(p+1)L-c_(p+2),
    f_(p+3)=c_(p+1)L^2+2c_(p+2)L+c_(p+3).

Put u=c_(p+1)L^2, v=c_1 c_(p+1)L, w=c_(p+2)L. The coefficient ideal
contains v, p u and p w: use c_1 f_(p+2), c_(p+1)f_2 and c_(p+2)f_1,
discarding scalar coefficient terms. Every source word involving these
extra coefficients or jets is, modulo the ideal, a multiple of u+2w.
Indeed f_(p+2)f_1 gives only -p(u+w)+v modulo scalars, while f_(p+1)
is itself an input coefficient; all other new words are in the ideal.
Lower jets use none of the deleted coefficients. Thus this enumeration
also accounts for the entire intersection with the coefficient-deletion
kernel. The functional extracting the u coefficient modulo p vanishes
on the coefficient ideal and takes value one on u+2w. Deletion removes
exactly one Z/p from the auxiliary defect, in every p>=5.

It remains to impose F=F(L)=f_p=0 over Z_(p)[c_1,...,c_p,L]. Write
J,C for the unrooted source and its coefficient ideal. The needed source membership follows from a low-degree identity:

    2a_1^2 L=q_1q_2-mq_3+2a_2q_1-3a_1a_2+m a_3.        (8a)

This is an exact polynomial identity in variable mass m, obtained by
expanding q_1,q_2,q_3. Every term on the right lies in the source.
Since a_1=-c_1 and 2 is a p-unit, c_1^2 L is already in the unrooted
source. Multiplying by F, itself a source generator, proves

    F c_1^2 L belongs to J.                              (8)

The remaining argument computes the entire intersection and supplies
a nonvanishing detector; membership alone would not determine the loss.
Modulo F times scalar coefficient polynomials, every element of (F)_d
is uniquely of the form

    A F L^3+B F c_1 L^2+C' F c_2 L+D F c_1^2 L.

If this is in J, its ordinary, linear c_1 and linear c_2 coefficients
force

    A=p a,       B=-3a modulo p,       C'=-6a modulo p.   (9)

To see this, the ordinary image in weight d is p Z_(p)L^d. Among pure
jet words of this weight, only F f_3 has ordinary coefficient of
p-valuation one. Modulo p its two linear coefficients are -1 and -2,
while its ordinary coefficient divided by p is 1/3. All other words
have the necessary extra p factors in these linear coefficients.
Words with positive input coefficients give zero to these two
coefficients modulo p. This proves (9); it imposes no condition on D.

Conversely the conditions (9) are sufficient, using (8). Use F f_3,
F c_1 f_2, F c_2 f_1 and (8). The first supplies the free ordinary
generator and the required residues in (9); the others supply
p F c_1 L^2, p F c_2 L, and F c_1^2 L. In the quotient by C,

    p F c_2 L=0,   p F c_1^2 L=0,
    p F c_1 L^2=2 F c_1^2 L.

Thus the torsion killed by F=0 is at most the cyclic group generated
by F c_1^2 L. This class has exact order p, not a smaller order.
For a polynomial of weight p+3, define a rational functional by

    phi(P)=( [c_1^2 L^(p+1)]P
                 +((p+2)/p)[c_1 L^(p+2)]P )/p.

It is integral on C. For a pure jet word W of weight p+2 with ordinary
coefficient A_0, its linear c_1 coefficient is exactly
-(p+2)A_0/p, since each f_j has relative linear coefficient -j/p.
Hence phi(c_1 W)=0. With another c_1 input, the remaining ordinary
coefficient has weight p+1 and is divisible by p. All other coefficient-
ideal words have zero contribution to phi. But

    phi(F c_1^2 L)=1/p.

Consequently the full root intersection removes one free ordinary
summand and exactly this one Z/p. Together with coefficient deletion,
this proves

    (Delta_p)_(p+3),(p)=(Z/p)^(a_(p+3)-2).                (10)

The identity (8a) proves membership before imposing F=0. The full
four-parameter intersection (9) excludes additional losses. In
particular, the argument does not replace the actual source by its
rational span or infer source membership from its vanishing after
F=0. Formula (10) is (1r).

## 7. Checks and remaining scope

The written argument is uniform in p,m; finite checks corroborate it.
The [Sage verifier](tools/verify_type_a_k_next_elementary.py) checks the
complete Witt generator comparison, the ideal identities, independence
of all five classes, and the actual elementary image. The
[native audit](tools/audit_type_a_k_next_elementary.py) independently
reconstructs elementary lattices, the formal trace identity, and the
quadratic reductions. See the [Sage report](type-a-k-next-elementary-verification.json)
and [native report](type-a-k-next-elementary-audit.json).

For m=p, in weights p+3 at p=5,7, the checked stages lose
one factor at coefficient deletion and another at the root equation.
Unlike weight p+2, the root equation now kills a torsion class as well
as the free ordinary image. Thus the previous root argument cannot be
copied. Section 6 gives the separate uniform proof, using (8a), the
full intersection and the order detector. No claim is made at p=3
or in subsequent weights.
