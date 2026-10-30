# The first correction to the type A specialization multiplicities

28 October 2026. IMG-41. Written deduction with internal independent checks.
Continue the notation and exact source conventions of
[IMG-40](TYPE_A_K_STABLE_MULTIPLICITIES.md). The representation is the
defining representation of SL_m, the parabolic is P1, and the source is
the full literal integral root-polynomial source, normalized by
theta*(x_i-x_j)=t_i-t_j. Delta_m is the actual K-specialization kernel.
Its grading is transported from the homogeneous model; it is not an
intrinsic grading of K-theory. Work locally at an odd prime p dividing m.

## 1. Target statement and proof dependencies

Put s=v_p(m), d=p+2, and

    a_d=sum_(k>=2) P_k(d-2^k+1),

where P_k(n) counts partitions into parts at most k. The first-correction
theorem is

    (Delta_m)_d,(p) = (Z/p^s)^(a_d-1) + Z/p^(s-1),   m>p,       (1)
    (Delta_p)_(p+2),(p) = (Z/p)^(a_(p+2)-2).                    (2)

A factor Z/p^0 is zero. In particular (2) is zero for p=3.
The root boundary m=p must not be silently included in (1).

The factors have explicit representatives. Put r=(p+1)/2 and
R_i=q_i q_(p+2-i)-m q_(p+2). For m>p use R_1 for the shortened factor,
the integral class

    E'=q_2 (q_p-q_1^p)/p-((m-m^p)/p)q_(p+2)

in place of R_2, and R_3,...,R_r. Add the earlier partition representatives
B_(k,lambda) for k>=3, in weight p+2. All factors except R_1 have full
order p^s. At m=p use only R_3,...,R_r and these k>=3 representatives.
They form the basis in (2).

The proof below has three logically separate parts: moment saturation,
the full elementary-source/coefficient comparison, and the root boundary.
It uses the earlier integral m-annihilator and actual-image descent,
and the exact three-coordinate tensor certificate of IMG-40 section 14.
Standard Newton identities, divided powers and elementary lattice exact
sequences are inputs. The deductions here are not a claim of literature
priority or external review. Degrees beyond p+2 and p=2 are not covered.

## 2. A tensor gluing lemma in arbitrary dimension

Let k>=3. Write W_B for z_0=sum_(i in B) z_i, B a subset of {1,...,k}.
In two sets of coordinates u,v, consider tensors of bidegree (1,2).
An S_k-equivariant family T_B on W_B, compatible on every pairwise
intersection, extends to an ambient S_k-invariant tensor over Z[1/2].
The extension is unique. Also, an ambient invariant tensor restricting
to tensors in the B-coordinates alone for every proper B is a multiple
of u_0 v_0^2.

Here is the reduction to the existing three-coordinate certificate.
For every three-element set S, set all coordinates outside S equal to
zero. Compatibility identifies the restrictions of T_B for all B with
the same B intersect S. Thus one obtains the certified S_3-equivariant
family, with its unique ambient extension on {0} union S.

The extensions agree on overlaps of size two. Indeed restriction to the
four binary hyperplanes in two coordinates is injective on ALL tensors
of bidegree (1,2), without an invariance assumption. To see this, take
the first three normals as coordinates x,y,z. A tensor vanishing on
their hyperplanes is a combination of

    u_x v_y v_z, u_y v_x v_z, u_z v_x v_y.

The fourth normal is y+z-x. Substitution first kills the second and third
coefficients (look at u_y v_z^2 and u_z v_y^2), then the first.
This argument uses no division.

Every tensor monomial involves at most three nonzero coordinate labels.
Its coefficient can therefore be read on a triple containing those labels.
For fewer labels, the triples containing them are connected by overlaps
of size two, so the coefficient is independent of the choice. Assemble
these coefficients to obtain the global tensor. Its local restrictions
agree on every triple and hence agree everywhere. Uniqueness also gives
S_k-invariance. The second assertion follows by restriction to triples
and the second part of the certified three-coordinate lemma.

For k=2 the full-family extension assertion is the separate two-coordinate
integer certificate. No assertion about its proper-subset coefficient
condition is needed: that condition has a larger kernel in dimension two.

This is a coordinate-support proof, not an extrapolation from several k.

### A general finite-support reduction

The same argument gives the following separate structural lemma. Fix
r>=2 and any commutative coefficient ring. If equivariant compatible
tensors of bidegree (1,r) on the binary hyperplanes extend in r+1
nonzero coordinates, then they extend in EVERY k>=r+1 coordinates.
An ambient invariant proper-subset coefficient condition likewise reduces
to r+1 coordinates. This lemma does not assert that the required base
extension or coefficient condition is true for every r.

For completeness, the overlap injectivity is uniform. In r nonzero
coordinates, the first r+1 normals z_0,z_0-z_1,...,z_0-z_r can be taken
as coordinates n_0,...,n_r. A tensor in every paired-coordinate ideal
has the form sum_i b_i u_i product_(j!=i) v_j. Restrict also to
n_0=n_1+n_2, the hyperplane belonging to {1,2}. The squared v_2 term
forces b_1=0 and the squared v_1 term forces b_2=0. Each u_i for i>=3
then forces b_i=0, and the remaining terms force b_0=0. All coefficients
used are units, so this holds over the original ring. Assemble coefficients
on (r+1)-element coordinate sets exactly as above.

Thus a higher Frobenius-layer proof can seek a fixed tensor certificate
depending on the offset r, rather than a growing list of masses or ranks.
Small coefficient degrees, the nonpolarizable high degrees, and the
actual elementary-source comparison still require their own arguments.

## 3. Moment saturation in weight p+2

Use the universal moment image over Z_(p)[M,a_1,...,a_d],

    q_j=M L^j+sum_(i=1)^j binomial(j,i) a_i L^(j-i),
    J=S[q_1,...,q_d], I=(a_1,...,a_d).

Give M and every a_i coefficient degree one. We claim that J and IJ are
p-saturated in every coefficient degree in weight d=p+2.

### 3.1 Coefficient degrees 2<=k<p: source

Polarization is valid since k! is invertible. In degree d<2p, the kernel
of the norm Sym^d -> Gamma^d consists of the distinct monomials z_i^p z_j z_l.
It is V^(1) tensor Sym^2(V); the cokernel has the corresponding divided
power description. These identifications commute with our restrictions.

A source annihilator's norm-cokernel tensor vanishes on every W_B.
For k>=3, four independent normals force it to vanish, since the
intersection of four paired-coordinate ideals has no total-degree-three
element. For k=2 use the four-hyperplane injectivity in section 2.
The restrictions of a norm preimage are a compatible family of kernel
tensors. Glue them by section 2, subtract the global kernel, and obtain
a polynomial divisible by

    A_k=product_(B subset {1,...,k}) (z_0-sum_B z_i).

The same injectivity test shows that norm is injective on this product.
The source annihilator has exactly its characteristic-zero dimension.

### 3.2 Coefficient ideal for k>=3

The second tensor assertion leaves only the pure ordinary norm-cokernel
functional. Subtract it. At the empty subset, subtract a kernel polynomial
so that z_0 divides the norm preimage P. Its restriction at a singleton
is a kernel tensor divisible by that coordinate. All such invariant
restrictions are supplied by these seven global polynomials, where
s_1=sum z_i and s_2=sum z_i^2:

    z_0^(p+2), z_0^(p+1)s_1, z_0^p s_1^2, z_0^p s_2,
    z_0^2 sum z_i^p, z_0 sum z_i^(p+1),
    z_0 s_1 sum z_i^p.

Expansion in inside/outside coordinates proves spanning; only 2 needs
inverting. Subtract them to kill every singleton restriction.

On any other proper B, the output-only term is itself a kernel tensor:
set an outside coordinate equal to sum_B z_i and all other outside
coordinates to zero, using the already established singleton factor.
For |B|>=3 choose two different inside indices j,l and one outside i.
The four independent normals

    sum_B z_i, sum_B z_i-z_j, sum_B z_i-z_l, sum_B z_i-z_i(outside)

force the tensor to vanish. For |B|=2 the four normals instead leave
only the span

    (u_j v_l-u_l v_j)(v_j+v_l-v_i).

It is antisymmetric in j,l; invariance and invertibility of 2 kill it.
Thus every proper-subset factor divides P. Norm is injective on their
product by four independent normals. This gives exactly the earlier
characteristic-zero coefficient-ideal annihilator dimension.

### 3.3 The small and high coefficient degrees

For k=1 the two source rows have unit pivots at a_d and M L^d; the ideal
is generated by a_d. For k=2 the ideal has a direct unit-pivot basis:
the scalar monomials M a_d and a_i a_(d-i), together with

    a_i q_(d-i), 1<=i<d.

Each latter row has its own monomial M a_i L^(d-i) with coefficient one.
Thus its saturation needs no norm-cokernel argument.

The only coefficient degrees not covered by polarization or J=A are
k=d-2=p and k=d-1=p+1. Here is an integral argument covering them.
For k=d-2, monomials with sum(i-1) over their positive coefficient factors
at least two are already in J, by M L=q_1-a_1. The same assertion is in
IJ if a positive coefficient is present. The remaining monomials are

    x_r=M^(r-2) a_1^(d-r) L^r,             2<=r<=d,
    y_r=M^(r-1) a_2 a_1^(d-r-2) L^r,       1<=r<=d-2.

Multiplying q_j by the appropriate powers of M L and the indicated
coefficients gives the following relations in A/J:

    y_r+j y_(r-1)=0                         (2<=j<=r),
    x_r+j x_(r-1)+binomial(j,2)y_(r-2)=0     (3<=j<=r).

The j=2,3 relations kill all y_r after inverting 2. The j=3,4 relations
then kill x_3,...,x_d and give 3x_2=0. The q_2^2 relation at r=4 is

    x_4+4x_3+4x_2+2y_2+4y_1=0,

so 4x_2=0 too, and x_2=0. For IJ use r<=d-1 in the x-relations;
the multiplier then contains a_1. All y-relations already contain a_2.
Since d>=5, this proves IJ is the entire positive-coefficient submodule.
Its complement is the pure ordinary monomial.

For k=d-1 the analogous relations are x_r+jx_(r-1)=0 for j=2,3;
they kill the positive-coefficient monomials after inverting 2.
For k>=d the identity M L=q_1-a_1 handles every monomial directly.

Consequently both J and IJ are p-saturated in every coefficient degree.
The earlier nilpotent-M base-change argument applies: C=(A/J)_d is
Z_(p)-free with nilpotent M, so M-m is injective on C. In the actual
mass-m MOMENT model its specialization defect is

    D_m=(Z_(p)/m)^(a_d).                                      (3)

The coefficient-degree-two summand has the basis

    R_i=q_i q_(d-i)-m q_d, 1<=i<d/2.                         (4)

Indeed these are all quadratic source words after removing the ordinary
lift and coefficient words, and their number equals the computed rank.
They are independent integrally by the saturation just proved.

## 4. Retain every new elementary-source and coefficient direction

Initially use an unrooted elementary model with independent coefficients
through weight d. This equals the actual model in the weight in question
when m>p, since m>=2p>d. The auxiliary m=p model is treated in section 6.

Newton's identities give the coordinates

    b=(a_p-a_1^p)/p, g=(q_p-q_1^p)/p=b+H,
    alpha=(m-m^p)/p,
    H=alpha L^p+sum_(i=1)^(p-1) binomial(p,i)/p
                       (a_i-m^(p-i)a_1^i)L^(p-i).

All these are integral in their specified rings. Other coefficient and
source coordinate changes through weight d have unit diagonal entries.
Because d<2p, a word contains at most one b or g.

Quotient out all scalar coefficient polynomials, which are in IJ. Write
J_0,C_0 for the old moment source and coefficient ideal in the remaining
positive-L ambient module. Define

    V=b(q_2-a_2),                 v=m b a_1 L,
    w_1=a_2 H, w_2=a_1^2 H, w_3=a_1 q_1 H,
    E=q_2 H-alpha q_d.

Only V is outside the old integral ambient module: its coefficient on
a_p a_1 L is 2/p. It has exact order p modulo that ambient module and
pV is in C_0. It occurs in both the new source and ideal, so it can be
cancelled without changing the old quotient or its embedded classes.

After this cancellation the new lattices are precisely

    J'=J_0+<v,w_1,w_2,w_3,E>,
    C'=C_0+<v,w_1,w_2,w_3>.                                  (5)

To verify completeness, multiply b and g by the five weight-two words
a_2,q_2,a_1^2,a_1q_1,q_1^2. Scalar products vanish. The identity
q_1^2=mq_2+a_1^2-ma_2 reduces bq_1^2 and gq_1^2; b a_1q_1 gives v,
g a_1q_1 gives v+w_3, and q_2g gives V+E+alpha q_d.
Thus (5) is an exhaustive word enumeration, not a list of some relations.

## 5. The two killed directions and the one added direction

Let C=(A/J)_d in the variable-mass moment model. The existing Bockstein
isomorphism identifies C_m[p] with ker(M:C/pC -> C/pC). Its use is
justified by the freeness and nilpotence established in section 3.

The classes of v,w_2,w_3 in C_m are zero. Their Bockstein representatives
are combinations of a_1 a_p L, a_1^2 L^p, and M a_1 L^(p+1), with
terms in high coefficient degree already zero by section 3.3. All three
low-degree monomials have a polarized exponent at least p. The source
annihilator in coefficient degree two lies in the norm image by section
3.1, so it annihilates these monomials. They are therefore in J modulo p.
In contrast, the representatives for w_1 and E are respectively

    -a_2 L^p,     a_p L^2+2a_(p+1)L.                         (6)

In coefficient degree one, eliminate a_d and M L^d with the two source
rows. The remaining a_i L^(d-i), 1<=i<d, are a free basis. The two
vectors in (6) are independent for every odd p. Hence w_1 and E add
two independent order-p directions modulo J_0, and only w_1 is in C'.

Here are the actual classes of the ideal directions lying in J_0:

    v=(m/(2p)) R_2,
    w_2=(m/p)(R_1-R_2/2),
    w_3=-(m/(2p)) R_2                  modulo C_0.            (7)

The fractions 1/2 are p-units; m/p is an integer. We prove these
identities rather than infer them from lattice experiments.

First put Z=q_1 H+m bL. The earlier integral even quadratic identity
of weight p+1, with p+1 inverted, places Z in J_0 and gives

    (p+1)(Z-alpha q_(p+1))=sum c_i R_(i,p+1-i) mod C_0.

Its cubic consequence is zero by the previously certified universal
cubic lemma. Thus q_1(Z-alpha q_(p+1))=0 in D_m. Multiplying Z by a_1
gives w_3+v in C_0, proving the last equality once the first is proved.

For the first equality use the integral moment trace certificate for
mR_2. In a formal coefficient/source ring write a_0=Q_0=M. That
certificate is

    S=sum_(i=0)^2 sum_(j=0)^p binomial(2,i)binomial(p,j)
                  a_i a_j Q_(d-i-j)-M^2 Q_d
      -sum_(r=1)^d a_r sum_(i=0)^2
                  binomial(2,i)binomial(p,r-i)Q_(2-i)Q_(p-r+i).

An out-of-range binomial term is zero. Substituting actual q's gives
S=M R_2 by the Hasse trace identity. Every term has an input coefficient.
At M=m, divide by p. Terms with M are in C_0 since p|m; terms with an
interior binomial(p,j) are also in C_0. The only terms left give

    (m/p)R_2 = -2 a_1q_1(q_p-a_p)/p mod C_0.

Replace (q_p-a_p)/p by H+(q_1^p-a_1^p)/p. In the latter summand,
interior binomial terms are visibly in C_0 using mL=q_1-a_1. The final
term is (m/p) times the specialization of
M^(p-1)a_1q_1 L^p, in coefficient degree p+1 and weight p+2; section
3.3 puts it in IJ. Thus (m/p)R_2=-2w_3=2v, as required.

Finally q_1Z, the identity q_1^2=mq_2+a_1^2-ma_2, and pE=R_2 mod C_0
give

    w_2-v=alpha R_1-(m/p)R_2.

Here mV and mw_1 are in C_0. Also alpha=m/p modulo m, and mD_m=0.
This proves the middle equality of (7). The relation pE=R_2 follows
by expanding H; the remaining term q_2q_1^p-m^p q_d is zero by
q_1^2=mq_2 and the m-annihilator.

The R_1,R_2 are independent basis directions in (4). Hence (7) kills
exactly their two order-p subgroups in D_m. Equation (6) proves that
there are no further killed directions: w_1 is independent modulo J_0,
while w_3=-v modulo C_0. It also proves that E supplies exactly one
new order-p quotient direction. The identities already proved give
mE=(m/p)R_2=2v=0 in the new quotient, including in the auxiliary model.
Together with pE=R_2 this restores the second cyclic factor to
order p^s, while the R_1 factor has order p^(s-1). All other summands
are unchanged. This proves (1).

The higher-coefficient partition representatives remain a basis of the
moment summands: each individual exponent is at most p-1 when k>=3
and the total weight is p+2. The integral triangular pairing therefore
still has unit diagonal factorials. Equation (7) changes only R_1,R_2,
so their images, together with the listed new quadratic representatives,
give the asserted elementary basis.

In particular the earlier bound on R_1 is sharp here: its p-primary
order is p^(s-1). This concerns m>p and weight p+2 only.

## 6. The auxiliary coefficients and the root boundary m=p

At m=p, first take independent elementary coefficients c_1,...,c_(p+2).
Define jets through degree p+2 by

    f_j=[u^j] sum_i (-1)^i c_i u^i(1+Lu)^(p-i), c_0=1.

The negative binomial exponents are used only in a finite truncation.
Newton translation gives exactly the elementary model of sections 4--5.
Here

    f_p=F(L), f_(p+1)=c_(p+1),
    f_(p+2)=-c_(p+1)L-c_(p+2).

Set c_(p+1)=c_(p+2)=0. In weight p+2 the kernel on source quotients is
generated by c_(p+1)L: it is in the source by the displayed jet, and its
coefficient-ideal intersection is precisely p times that generator,
from c_(p+1)f_1. The other removed directions are scalar and already in
the ideal. The ordinary image is unchanged. Thus this deletion removes
one copy of Z/p from the specialization defect. It is a genuine torsion
direction; ignoring it would give the wrong root-boundary answer.

It remains to impose F=0 in the resulting unrooted model over
Z_(p)[c_1,...,c_p]. Modulo coefficient-ideal elements F c_2 and F c_1^2,
the ambient F-ideal in this weight consists of

    A F L^2+B F c_1 L.

Its intersection with the source is exactly the lattice

    A=p a, B=-2a mod p.                                    (8)

For necessity, the ordinary image in this weight is p Z_(p)L^(p+2).
In a source word without F every positive jet has ordinary coefficient
divisible by p, so its ordinary coefficient has at least two such factors.
The only source word that can have ordinary coefficient of valuation one
is F f_2. Its ordinary coefficient divided by p is (p-1)/2, while its
c_1 L^(p+1) coefficient modulo p is 1. Other source words contribute
zero to this linear coefficient modulo p after the ordinary term has
been accounted for. The proposed element has linear coefficient B-A,
which is B modulo p. This proves B=-2a modulo p.

Conversely the lattice in (8) is generated by

    p F L^2-2F c_1 L,       p F c_1 L.

The first is 2/(p-1) times F f_2 modulo F c_2, and the second is in
IJ using c_1 F f_1 and F c_1^2. The first maps to a generator of the
free ordinary image p Z_(p)L^(p+2). Imposing F=0 therefore removes
only that free summand and preserves the entire remaining defect.

Formula (1) for the auxiliary model with s=1 had a_d-1 copies of Z/p.
Coefficient deletion removes one and the root relation removes no
further torsion. This gives (2).

To see the asserted root-boundary basis, take the coefficient of
c_(p+1)L modulo p on the auxiliary defect. It vanishes on the coefficient
ideal and on every listed R_i and higher B_(k,lambda). On E' it is
2 alpha(p+1), a unit modulo p, since a_j has leading coefficient -j c_j.
It takes value 1 on the deleted class c_(p+1)L. Therefore this deleted
class has a nonzero E' coordinate in the auxiliary basis. Quotienting
leaves precisely the displayed remaining basis; the root relation does
not change it.

## 7. Checks, status, and next obstruction

The [Sage verifier](tools/verify_type_a_k_first_correction.py) checks ten
actual elementary-jet pieces, including their explicit generators, three
separate root-boundary stage comparisons, and five moment/Witt lattice
comparisons including all ideal class identities and Bockstein independence.
The [native audit](tools/audit_type_a_k_first_correction.py) independently
builds five elementary image lattices, checks symbolic trace certificates
with a variable mass, and verifies the fixed high-coefficient elimination.
The two-coordinate and three-coordinate tensor matrices are separately
built by polynomial substitution and integer coefficient dictionaries.
Their [Sage report](type-a-k-first-correction-verification.json) and
[native report](type-a-k-first-correction-audit.json) bind the proof sources.
These finite checks support the written uniform deduction; they are not
an external mathematical review. The root comparison retains the
auxiliary coefficient torsion and the full F-ideal.

The next degree p+3 has Frobenius tensors of bidegree (1,3), and the
elementary word comparison has weight-three multipliers.
[IMG-42](TYPE_A_K_NEXT_LAYER_OBSTRUCTION.md) proves the complete relative
MOMENT comparison in this degree for p>=5: a common source/ideal torsion
block cancels, despite failure of source base change alone. It does not
by itself complete the elementary word comparison or the root boundary.
[IMG-43](TYPE_A_K_NEXT_ELEMENTARY.md) subsequently completes that word
comparison and both actual-source formulas for p>=5, including a
separate root-intersection argument. For p=3
this degree also reaches 2p, so words with two Witt generators appear.
No extension at p=3 is assumed.
