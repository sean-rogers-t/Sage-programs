# A uniform next-layer obstruction, and its relative cancellation

28 October 2026. IMG-42. Written deduction after
[the first corrected multiplicities](TYPE_A_K_FIRST_CORRECTION.md).
This note concerns the UNIVERSAL MOMENT model. The subsequent
[IMG-43 comparison](TYPE_A_K_NEXT_ELEMENTARY.md) transfers the result to
the actual elementary Springer source for every p>=5 dividing m, with
a separate root formula at m=p.

## 1. Precise statement

Let p>=5 be prime, d=p+3, and work over Z_(p). Use

    S=Z_(p)[M,a_1,...,a_d], A=S[L], I=(a_1,...,a_d),
    q_j=M L^j+sum_(i=1)^j binomial(j,i)a_i L^(j-i),
    J=S[q_1,...,q_d], C=A/J, B=A/IJ.

M has weight zero; a_i has weight i and L has weight one. M and every
a_i have coefficient degree one, while L has coefficient degree zero.
Subscripts (k,d) specify coefficient degree and weight, respectively.

The p-primary torsion of C_(3,p+3) is exactly Z/p. The corresponding
torsion of B_(3,p+3) is also Z/p, and B -> C identifies these two
torsion groups. Multiplication by M kills both groups, because the
coefficient-degree-four targets are p-torsion-free. In particular,
(J/IJ)_(3,p+3) is p-torsion-free. These are the ONLY p-torsion blocks
of C_d and B_d: all other coefficient degrees are p-torsion-free.

Consequently the source-image map

    J/(M-m)J -> A/(M-m)A

has a nonzero Z/p kernel in this weight for every integer m divisible
by p. The hypothesis used to identify actual moment images by injective
base change in the preceding proof is therefore false in the next
degree. HOWEVER, this particular kernel contribution disappears after
passing to J/IJ: the matching coefficient-ideal torsion cancels it.
This distinction is essential. A failure of source base change alone
does not establish a new specialization-defect factor. In fact, for the
actual mass-m MOMENT image in the unrooted ambient polynomial ring,

    (moment specialization defect)_(p+3),(p)
           = (Z/p^v_p(m))^(a_(p+3))                          (1)

for every p>=5 dividing m. The actual elementary Springer source has
additional Witt generators and coefficient changes. Formula (1) is not
a claim about its degree-p+3 defect.

## 2. A general norm bound on torsion exponents

Fix k<p and p<=d<2p, so polarization in k slots is integral p-locally.
The norm from symmetric to divided powers is diagonal on monomials,
with diagonal entries alpha!. Each has p-valuation at most one.
It is injective over Z_(p), and its cokernel is killed by p.

Let J_div be the sum of the divided-power source spaces on the binary
hyperplanes. The dual ordinary-polynomial annihilator is the product
of their distinct primitive linear equations, with S_k-invariants.
Its dimension is the same in characteristic zero and characteristic p.
Thus J_div is p-saturated. The commuting norm square gives

    p J_div subset Norm(J),

since p times every local divided monomial has an integral symmetric
preimage. For x in the p-saturation of J, Norm(x) belongs to J_div;
therefore p Norm(x) belongs to Norm(J), and injectivity gives p x in J.
Thus all p-primary torsion of A/J in this band is killed by p.

The same argument applies to IJ. The divided-power coefficient ideal
has dual polynomial annihilator the pure ordinary line plus the
proper-subset product, by IMG-40's coefficient-ideal calculation.
Its characteristic-independent rank proves saturation. Norm preserves
the condition of a positive unshifted coefficient, so p times this
divided ideal lies in Norm(IJ). Therefore all p-primary torsion of
A/IJ in this band is likewise killed by p.

These are general lattice deductions; neither one presumes the ranks
of the ordinary symmetric source or its coefficient ideal are constant.

### Quadratic coefficient degree in every weight

The coefficient-degree-two source is saturated over Z[1/2] in EVERY
weight, without a bound on the weight or prime. Polarize in two slots
x_1,x_2 and put

    y=(x_1+x_2)/2, v=(x_1-x_2)/2, z=L+y, w=v^2.

The invariant ambient ring is Z[1/2][y,z,w], with weights 1,1,2.
The empty and full subset source spaces are respectively the degree-d
parts of Z[1/2][y,w] and Z[1/2][z,w]. They span a coordinate-monomial
submodule. The remaining source space is spanned by the even-in-v
parts of

    (z+v)^i(y-v)^(d-i), 0<=i<=d.

The two endpoint polynomials are already in the empty/full source.
For 1<=i<d, the w=0 term is exactly z^i y^(d-i). These are distinct
mixed monomials, absent from the empty/full source, with unit coefficients.
Thus all the remaining generators have distinct unit pivots, proving
saturation over Z[1/2]. This also gives source rank 2d for d>0.

The coefficient ideal is saturated integrally: its scalar coefficient
monomials together with a_i q_(d-i), 1<=i<d, have distinct unit pivots
M a_i L^(d-i). This proves that coefficient degree two creates no odd
primary torsion in either ambient quotient, in any weight.

### A uniform recurrence for high coefficient degree

There is a separate parameter-free filling lemma. Put d=k+e with e>=2
and k>=e+1. Over Z[1/e!], the entire coefficient-degree-k, weight-d
ambient module lies in J, and IJ is exactly its positive-coefficient
submodule. For e=1 the same conclusion holds when k>=3 over Z[1/2].
No mass specialization or polarization is needed.

Write a_lambda for a product of positive-index coefficients, and consider
X=M^u a_lambda L^t. Its coefficient degree is k=u+length(lambda), and
its deficit is h=t-u. The weight condition gives

    h=e-sum_(i in lambda)(i-1) <= e.

For h<=0 use X=a_lambda M^(u-t)(q_1-a_1)^t. Induct on h>0. Suppose
u>=1, a_1 divides a_lambda, and length(lambda)>=2. Pascal's identity gives

    (a_lambda/a_1)(M L)^(u-1)
                    [M q_(h+2)-(M L)q_(h+1)]
       = X + sum_(i=2)^(h+2) binomial(h+1,i-1)
                   M^u (a_lambda/a_1)a_i L^(t+1-i).          (3a)

Every summand on the right after X has smaller deficit. The left is in
IJ, since M L=q_1-a_1 is a source element and a_lambda/a_1 is a positive
input coefficient. Thus X is in IJ without any division.

If u>=2 and lambda is nonempty, expansion of
a_lambda (M L)^(u-1) q_(h+1) has leading term X. Its next term is
(h+1)M^(u-1)a_lambda a_1 L^(t-1), which is covered by (3a); all subsequent
terms have smaller deficit. This settles every such X. When u=1,
length(lambda)=k-1>=e exceeds the waste sum e-h, so lambda contains
a_1 and has at least two entries: (3a) applies directly.

When u=0, lambda again contains a_1. Expand
(a_lambda/a_1)q_(h+1). Its leading term with one M is already handled,
its next term is (h+1)X, and the rest has smaller deficit. If h<e,
divide by h+1<=e. When h=e>=2, the only such target is X=a_1^k L^e.
Expand instead

    a_1^(k-2) q_2 q_e.

Its all-a_1 term is 2e X. Every other term either has positive M-power,
already handled, or has a coefficient index greater than one, so its
deficit is smaller than e. Since k-2>=1, this is a coefficient-ideal
source word. Division by 2e is allowed over Z[1/e!]. For e=1 the
original single-q step divides only by 2, a unit in that case's ring.
This finishes all positive-coefficient monomials.
Finally, if lambda is empty, u=k>=3; expansion of
(M L)^(u-1)q_(h+1) expresses X as a source element plus the already
handled positive-coefficient terms. This proves J=A in that bidegree.
The reverse inclusion for IJ follows from its definition.

This is a terminating recurrence: (3a) lowers deficit; the other steps
either reach the already handled positive-M case at the same deficit
or lower the deficit. It does
not divide by M, by k, or by a representation rank. Repeated divisions
can occur, so the assertion is localization at e! (at 2 when e=1), not a bound of
one factorial on the common denominator of every expanded witness.

In particular, for EVERY odd prime p and weight d<2p, ALL
nonpolarizable coefficient degrees k>=p are handled by this lemma
(and by M L=q_1-a_1 if k>=d). Indeed e=d-k<=p-1 and k>=e+1;
the e=1 case has k>=3. Both ambient quotients are p-torsion-free in
all these degrees. Thus a first-Frobenius-band proof only needs the
polarizable coefficient degrees k<p. For r=3 this supplies the high
degrees used here, but the recurrence is uniform throughout that band.

## 3. The fixed obstruction in three coordinates

Write W_B for z_0=sum_(i in B) z_i, B subset {1,2,3}. On each W_B
consider tensors in u,v of bidegree (1,3). Define a family T_B as follows:

    T_empty=T_{1}=T_{2}=T_{3}=0.

For B={j,l}, with outside index i, put

    T_B=(u_j v_l-u_l v_j)(v_j+v_l-v_i)(v_j-v_l).

For the full subset put

    T_full=(v_1+v_2+v_3) sum_i u_i v_i^2
                     -(u_1+u_2+u_3) sum_i v_i^3.

This family is equivariant under S_3. It is compatible on every pairwise
intersection: substitution of the corresponding equality of subset
sums in both u and v proves this directly with integer coefficients.
It is NOT the restriction of a global invariant tensor, even if only
the proper subsets are retained.

Here are the fixed integer certificates. For the full family the local
invariant coordinates have dimension 44; the compatibility matrix has
rank 23, with nonunit Smith factors 2,2. Global invariant tensors have
dimension and restriction rank 20, with all restriction Smith factors
one. Thus the gluing obstruction has dimension one over every field of
characteristic other than two. On proper subsets the dimensions are
38,17,20, respectively, with the only compatibility nonunit equal to 2.
Appending the displayed T_B to the global restriction columns gives
rank 21 with all nonzero Smith factors one, both for full and proper
families. It therefore generates both obstruction spaces, and restriction
from the full to the proper family identifies those spaces.

For four nonzero coordinates the full-family dimensions are 80,59,21,
and the proper-family dimensions are 73,52,21. The compatibility and
restriction matrices have only unit nonzero Smith factors. Both
obstructions vanish. The invariant coefficient-condition kernel is
the pure u_0 v_0^3 line for k=3 and k=4; its constraint matrix is primitive.

All these matrices involve fixed dimensions, independent of p and m.
Polynomial substitutions and separate coefficient-dictionary expansions
provide independent exact constructions. The finite-support lemma in
IMG-41 additionally extends the four-coordinate full-family gluing
statement to every k>=4. No inference from a finite prime list is used.

## 4. From the tensor obstruction to source torsion

In characteristic p>=5 and degree p+3<2p, norm has kernel and cokernel
V^(1) tensor Sym^3(V), with the kernel embedded by u_i v_j v_l v_h
mapping to z_i^p z_j z_l z_h. These monomials are distinct, and the
description commutes with restriction.

A source annihilator's norm-cokernel tensor vanishes on every W_B.
The overlap injectivity lemma from IMG-41, with r=3, makes it zero.
Thus the annihilator is in the norm image. A norm preimage restricts
to a compatible kernel-tensor family. The obstruction to subtracting
a global kernel tensor is precisely the one-dimensional space in
section 3. If that obstruction vanishes, the preimage is divisible by
the full product A_3 of the eight binary equations, as in characteristic
zero. Norm is injective on this product by the same tensor injectivity.

Every compatible tensor family actually occurs here. Regard it as
ordinary degree-(p+3) polynomials on the hyperplanes. No three distinct
binary normals are dependent: an affine line through two different
0/1 vectors has no third 0/1 vector. Pairwise-compatible polynomials
therefore glue. Inductively the residual on the next hyperplane is
divisible by the product of the distinct previous intersection equations;
lift the quotient and correct by that product. Average over S_3, since
6 is invertible, to retain invariance.

It follows that the source annihilator dimension increases by exactly
one modulo p. Thus the source rank drops by exactly one. The exponent
bound in section 2 identifies its torsion as one Z/p, not an unknown
higher p-power factor.

For k=4 the norm-cokernel tensor vanishes using five independent normals
and total tensor degree four. The four-coordinate extension certificate
then removes all kernel families. The source annihilator has its
characteristic-zero dimension, so C_(4,p+3) has no p-torsion. This proves
that M kills the torsion in C_(3,p+3).

The same proof works for every 4<=k<p by the finite-support gluing
lemma. Coefficient degree one is directly saturated, degree two was
settled above in every weight, and all k>=p are settled by the filling
recurrence. Thus C_d has no other p-torsion blocks.

## 5. Why the coefficient ideal has the SAME torsion

The norm-cokernel class of an IJ annihilator is on the pure ordinary
line by the coefficient-condition certificate. Subtract that functional.
The norm preimage has kernel-tensor restrictions at the empty subset
and at each singleton, since a polynomial in one variable of degree
p+3 is itself a kernel monomial.

These restrictions can be removed by a global kernel tensor. The empty
and singleton normals are independent coordinates; their paired-coordinate
ideals form a coordinate-monomial arrangement, on which compatible
tensors glue coefficient by coefficient. Average over S_k (k=3 or 4).
The preimage now vanishes on those hyperplanes. On every other proper
B its output-only part is also a kernel tensor: set an outside coordinate
equal to the sum over B and the other outside coordinates to zero.
The already established singleton vanishing then proves the assertion.

We are left with precisely a compatible kernel-tensor family on proper
subsets. Modulo global tensors its obstruction is the proper-family
space of section 3. A zero obstruction leaves the proper-subset product
A_3^-; all obstruction families occur by the polynomial gluing argument.
Hence the coefficient-ideal annihilator dimension increases by exactly
one modulo p. Section 2 again makes the corresponding torsion exactly Z/p.

The integral annihilator lattice has no hidden p-index in this argument.
Norm on A_3^- times polynomials is injective modulo p, because a kernel
tensor vanishing on all proper binary hyperplanes is zero. Its reduction
is disjoint from the primitive pure divided-power ordinary functional.
Those lattices are therefore primitive and already give the whole
characteristic-zero annihilator. The same statement holds for A_3.

Reduction of annihilators modulo p, modulo reductions of their integral
annihilators, is naturally dual to the p-torsion of the corresponding
ambient quotient, since that torsion is killed by p. The full-to-proper
obstruction map is nonzero on the displayed T_B and is an isomorphism.
Consequently the natural map B_(3,p+3) -> C_(3,p+3) identifies their
torsion groups. The proper four-coordinate certificate proves the
coefficient-degree-four target of B is also p-torsion-free. Thus M
kills this shared torsion, and J/IJ has no p-torsion in this block.

For every 4<=k<p the proper-family obstruction vanishes as well.
For k=4 this is a certificate. For k>4 restrict a proper family to
four coordinate labels: the full four-label subset is still proper
globally, so this gives a FULL four-coordinate family. The finite-support
proof assembles its extensions, using three-coordinate overlap injectivity.
The ambient coefficient condition reduces to four labels in the same
way. The preceding norm argument therefore proves saturation of IJ
for these k. Together with the quadratic and high-degree lemmas, this
accounts for every coefficient degree of B_d.

## 6. Actual images and the cancellation criterion

The ambient A_d is free over Z_(p)[M]. Tensoring 0 -> J -> A -> C -> 0
with M=m therefore identifies the kernel of source-image base change
with ker(M-m:C -> C). The shared Z/p in coefficient degree three is
killed by M and by m, so supplies the asserted nonzero kernel for every
p|m. This is the first weight at which the preceding all-coefficient
saturation theorem fails.

For the specialization quotient, use instead

    0 -> J/IJ -> A/IJ -> A/J -> 0.

The kernel contributed to (J/IJ)/(M-m) by actual-image base change is
the image of

    coker(ker(M-m:A/IJ -> A/IJ)
                  -> ker(M-m:A/J -> A/J)).

This is an exact-sequence statement, with the actual coefficient ideal
retained. The shared torsion block maps isomorphically. It is the entire
p-torsion of each ambient quotient, and the torsion-free part of C is
M-nilpotent, so ker(M-m:C -> C) consists of exactly that block. The
torsion-free kernel J/IJ has the usual decomposition into a free ordinary
Z_(p)[M] summand and a finite Z_(p)-free defect killed by M. Its rational
rank is a_d by the characteristic-zero partition calculation. Thus M-m
is injective on J/IJ. Subtracting the matching torsion lift in B shows
that ker(M-m:B -> B) is the same block too. The displayed cokernel is zero.

Consequently actual-image specialization of J/IJ agrees with abstract
base change, even though this is false for J alone. Its defect is
(Z_(p)/m)^(a_d), proving (1). This argument retains the full coefficient
ideal and does not assume a flat source.

## 7. A checked limit: relative saturation fails in higher weight

The cancellation just proved is NOT an all-weight relative-saturation
theorem. In the universal integral moment model, coefficient degree
three and weight twelve, put

    R_(i,j,l)=q_i q_j q_l-M^2 q_12,
    T=R_(1,2,9)+12R_(2,3,7)+21R_(2,4,6)-63R_(2,5,5)
                      -84R_(3,3,6)+252R_(3,4,5)-140R_(4,4,4).

The class of T in J/IJ has EXACT order five. The
[saved certificate](experiments/type-a-k-relative-five-obstruction.json)
gives an explicit integral linear combination of 78 literal coefficient-
ideal words equal to 5T. It also supplies a rational linear functional
which is integral on every one of the 90 coefficient-ideal generators
in this bidegree and has value 4889/5 on T. Thus T is not in IJ.
The certificate includes the full M-homogenization, not a specialization
of M to one. An independent polynomial replay verifies the identity
with variable M and every functional evaluation, without Smith form.

In particular, even a prime greater than the coefficient degree can
occur in the relative torsion. This is an obstruction to extending the
method by claiming relative saturation in every weight. It is not a
claim about the actual elementary Springer source: adding its Witt
generators and coefficient relations may change this class. No general
formula for the first higher failure is inferred from this example.

## 8. Verification and status

The [Sage check](tools/verify_type_a_k_next_layer.py) and the
[native check](tools/audit_type_a_k_next_layer.py) independently construct
the fixed integer tensor arrangements. They check the moment source,
coefficient ideal, and natural torsion map at p=5,7. The Sage check also
tests rational/modular ranks at p=13, including the first positive-
codimension four-coordinate source; that control does not compute
Smith factors. Three actual mass-m moment images check the relative
base-change conclusion. See the
[Sage report](type-a-k-next-layer-verification.json) and
[native report](type-a-k-next-layer-audit.json).

The high-coefficient recurrence has explicitly constructed formal
source witnesses, including literal coefficient-ideal witnesses for
every positive-coefficient target tested. The
[symbolic builder](tools/build_type_a_k_high_coefficient.py) and
[independent Sage replay](tools/verify_type_a_k_high_coefficient.py)
check the substitutions, grading, and denominator prime bound. Their
[build report](type-a-k-high-coefficient-build.json) and
[replay report](type-a-k-high-coefficient-verification.json) bind the
[witness archive](experiments/type-a-k-high-coefficient-witnesses.json.gz).
These examples verify implementations of the written recurrence;
the induction in section 2 is the all-rank argument.

The negative control has a separate
[builder](tools/build_type_a_k_relative_obstruction.py),
[independent replay](tools/audit_type_a_k_relative_obstruction.py),
[build report](type-a-k-relative-obstruction-build.json), and
[audit report](type-a-k-relative-obstruction-audit.json).

## 9. What is known and what remains open

The norm and exact-sequence tools are standard. The fixed arrangement
certificates and their uniform consequences are written deductions in
this repository, with internal independent checks rather than external
review or a priority claim. The p+2 theorem covers the actual Springer source for all odd p|m,
including m=p; IMG-43 extends the actual-source formulas to p+3 for
p>=5, including its distinct root correction.

The relative moment comparison is now proved in this first nonflat
weight. The subsequent IMG-43 proof completes the exhaustive elementary-source
comparison with weight-three multipliers for m>p. At m=p it proves the
coefficient-deletion comparison and the full root intersection, using
a degree-three identity that already puts c_1^2 L in the source.
At p=3 this weight is 2p and two Witt factors can occur; it is explicitly
outside the theorem above.
