# A partition formula for the specialization defect through the residue characteristic

**Subsequent completion:** [IMG-41](TYPE_A_K_FIRST_CORRECTION.md) gives the
full corrected formula in weight p+2 for odd p, including its distinct
root boundary m=p. The unresolved p+2 statements below describe the
historical IMG-40 stage and are superseded by that proof.

28 October 2026. IMG-40. This is a continuation of
[IMG-39](TYPE_A_K_INTEGRAL_RECURRENCE.md), for the defining representation
of SL_m, P=P1, and the full literal integral root-polynomial Springer
source. It is a formula in a specified range, not a claim to have computed
the entire defect in every rank. The grading below is transported from
the universal homogeneous model; it is not an intrinsic grading of K-theory.

## 1. Statement and scope

Use the actual image and specialization kernel

    R=Z[t_1^+-1,...,t_m^+-1]/(product t_i-1), I_T=(t_i-1),
    A_T=K_T^0(SL_m/P1)=R[L]/product_i(L-t_i),
    f_j=[u^(m-j)] product_i(L+u-t_i),
    J_T=R[f_1,...,f_(m-1)],
    Delta_m=ker(J_T/I_T J_T -> ordinary image).

The normalized Springer pullback is theta^*(x_i-x_j)=t_i-t_j.
The [homogeneous descent](TYPE_A_K_GRADED_LIFTS.md) identifies this kernel
with the graded kernel for the universal polynomial F and augmentation
c_1=...=c_m=0. That descent retains the actual image and integral torsion.

For k>=1 let P_k(n) be the number of partitions of n into parts at most k,
with P_k(0)=1 and P_k(n)=0 for n<0. Define

    a_d=sum_(k>=2) P_k(d-(2^k-1)).                         (1)

Only finitely many terms contribute. Equivalently,

    sum_(d>=0) a_d q^d
      =sum_(k>=2) q^(2^k-1)/product_(i=1)^k(1-q^i).       (2)

**Theorem (stable range, including its prime boundary).** Let p divide m,
s=v_p(m), and 0<=d<=p. Then

    (Delta_m)_d localized at p = (Z/p^s)^(a_d).           (3)

Thus every nonzero factor in this range has the full order p^s. There
are no factors of smaller p-power order here. For the layer polynomials
H_(p,r) defined in IMG-39, every 1<=r<=s has the same truncation

    H_(p,r)(q) = sum_(d=0)^p a_d q^d modulo q^(p+1).

**First extension.** If p is odd, the same factor formula holds
also in degree d=p+1. Section 13 proves this with an integral source-and-
coefficient comparison. In that degree replace the single representative
R_(1,p)=q_1 q_p-m q_(p+1) in (3b) by E_p of section 12. Keeping the old
representative would give a nonprimitive class. This extension does not
claim the case p=2, d=3. Section 13.4 includes the root case m=p.

The first values for d=0,...,10 are

    0, 0, 0, 1, 1, 2, 2, 4, 4, 6, 7.

The factors have explicit integral representatives. Write q_j for the
Newton power sum of the elementary Springer generators f_1,...,f_j.
For k>=2 and a weakly increasing k-tuple lambda of nonnegative integers
with sum(lambda)=d-(2^k-1), put e_i=2^(i-1)+lambda_i. Then

    B_(k,lambda)=product_(i=1)^k q_(e_i)-m^(k-1) q_d.    (3b)

Their classes form a basis of the group in (3) as a Z/p^s-module.
The representatives involve no division. The assertion that they are
a basis is p-local and has the same d<=p scope as (3).

This is uniform in m, including mixed m and arbitrary v_p(m). It does
not determine degrees beyond these ranges, in particular the remaining
small-prime and overflow parts. Sections 2--6 prove d<p using Newton coordinates,
polarization and the perfect differential pairing. Sections 10--11 handle
d=p separately, retaining the first Witt correction and divided powers.

Equivalently, using IMG-39's m-annihilator, for EVERY m>=2 and d>=0,

    (Delta_m)_d tensor Z[1/d!]
      = (Z[1/d!]/m Z[1/d!])^(a_d).                     (3a)

For m<=d both sides vanish, since every prime dividing m has been
inverted. Formula (3a) does not recover the primes removed by localization.
Weight zero has zero defect directly.

The written deduction below combines standard tools with a subset-indexed
hyperplane calculation. Computational checks corroborate separate parts of the partition proof.
The further annihilator in section 14 explicitly uses a fixed exact
integer tensor certificate as one lemma of its all-prime proof. External review and priority are
not established.

## 2. A universal moment image with a variable mass

Introduce a variable M, coefficient variables a_i of weight i, and L of
weight one. M has weight zero. For work through weight d it suffices to
use a_1,...,a_d and the moments

    q_0=M,
    q_j(L)=M L^j+sum_(i=1)^j binomial(j,i) a_i L^(j-i).

Put S=Z[M,a_1,...,a_d], A=S[L], J=S[q_1,...,q_d], and I=(a_1,...,a_d).
All assertions below concern weights at most d, so omitting higher
moments does not change them. Give M and every a_i coefficient degree
one and L coefficient degree zero. Every q_j has coefficient degree one.
The ordinary weight-d image, for d>0, is

    (J_0)_d=M Z[M] L^d,

a free rank-one Z[M]-module. Let D_d be the weight-d kernel of J/IJ -> J_0.
This D_d is an auxiliary universal moment defect, not yet Delta_m.

**Integral moment trace lemma.** M D_d=0.

Define the S-linear operator T on S[Y] by T(Y^r)=q_r(L), including
T(1)=M. The binomial identity gives

    T(H)=M H(L)+sum_(r>=1) a_r H^{[r]}(L),               (4)

where H^{[r]} is the r-th Hasse derivative, and only r<=d occurs in weight d.
Also

    q_j(Y+t)=sum_(r=0)^j binomial(j,r) q_(j-r)(Y)t^r.

Therefore the Hasse coefficients of a polynomial in the q_j evaluate
into J. Reduce an element of the kernel modulo IJ to R(q), where R has
coefficients in Z[M]. Its ordinary image is zero as a polynomial in M,L,
so H(Y)=R(q(Y)) belongs to I S[Y]. Hence T(H) is in IJ, while (4) differs
from M R(q) by an element of IJ. This proves the lemma integrally.
There is no root equation, rational division, or assumption that T is
multiplicative.

The quotient C_d=(A/J)_d is a finite Z-module in the sense of being
finitely generated, not necessarily finite as a set. Indeed ML=q_1-a_1,
so M^d A_d is contained in J_d. It has only coefficient degrees below d:
for a monomial M^u a_lambda L^r of coefficient degree k,
k=u+length(lambda)>=d implies u>=r, since r+length(lambda)<=d.
It can then be written using (ML)^r in J. These observations justify
the finite coefficient-degree arguments used below.

## 3. Polarization identifies the source spaces

Fix a field K of characteristic zero or characteristic p>d. In coefficient
degree k<=d, polarization identifies the weight-d ambient space with

    V_(k,d)=K[L,x_1,...,x_k]_(d)^(S_k).

Explicitly a product a_(i_1)...a_(i_k), with a_0=M, maps to the average
of x_1^i_1...x_k^i_k over S_k. The factor k! is invertible. A q_j in a
slot instead gives (x_i+L)^j in that slot. Multiplying by M inserts a
zero-exponent slot; it does not disappear from the coefficient grading.

For B subset {1,...,k}, let W_B be the span of the k linear forms

    x_i+L for i in B;   x_i for i not in B.

The polarized J_(k,d) is the invariant part of the sum of Sym^d(W_B).
The polarized (IJ)_(k,d) is the invariant part of the sum, over B, of
the degree-d part of the ideal generated by x_i, i not in B, in Sym(W_B).
The latter condition is exactly that at least one unshifted input
coefficient has positive index. A factor M alone is NOT such an input.
These descriptions enumerate all source words, not only quadratic words.

Taking invariant parts commutes with these sums because k! is invertible.
The symmetric ambient Hilbert series is

    H_k(q)=1/((1-q) product_(i=1)^k(1-q^i)).              (5)

## 4. Compute the two annihilators

Use dual variables z_0,z_1,...,z_k and the differential pairing
<z^alpha,x^beta>=alpha! if alpha=beta and zero otherwise, with x_0=L.
It is perfect in degree d because d! is invertible. Its restriction to
invariants is perfect as well, by the self-adjoint averaging operator.

Define

    ell_B=z_0-sum_(i in B) z_i,
    A_k^-=product_(B proper subset {1,...,k}) ell_B,
    A_k=A_k^- ell_{1,...,k}.

All these linear factors are distinct over every field. Their products
are S_k-invariant, with degrees 2^k-1 and 2^k, respectively.

The annihilator of Sym^d(W_B) is the degree-d part of (ell_B).
Intersecting these principal ideals gives their product, since the
distinct linear forms are pairwise nonassociate irreducibles. Thus

    J_(k,d)^perp = (A_k K[z_0,...,z_k])_d^(S_k).         (6)

For the coefficient ideal, a dual polynomial P must satisfy, for every
proper B,

    P(sum_(i in B) z_i,z_1,...,z_k)
        belongs to K[z_i : i in B].                    (7)

This follows directly by pairing against monomials having a positive
power in an unshifted slot. For k>=2 and d>0, (7) has exactly the solutions

    (IJ)_(k,d)^perp
      = K z_0^d direct-sum (A_k^- K[z_0,...,z_k])_d^(S_k).  (8)

Here is the completeness argument. For B empty, P(0,z)=0. For singleton
B={i}, symmetry and homogeneity give P(z_i,z)=c z_i^d with the SAME
constant c for every i. Subtract c z_0^d. The result is divisible by
z_0 and each z_0-z_i. For any other proper B, choose i outside B.
Its restriction in (7) is independent of z_i but is divisible by
sum_(j in B) z_j-z_i. Over a polynomial domain this forces that restriction
to vanish. Therefore every proper-subset factor divides the result,
and their product divides it. The converse is immediate. The sum is
direct because z_0^d is not divisible by z_0-z_1.

For k=1 the separate answer is (IJ)_(1,d)^perp=z_0 K[z_0,z_1]_(d-1).
Together with (6) it gives dim(J/IJ)_(1,d)=1, hence no defect in that
coefficient degree. Coefficient degree zero likewise contributes no defect.

For k>=2, subtracting (6) from (8) and using (5) gives

    dim (J/IJ)_(k,d) = 1+P_k(d-(2^k-1)).                (9)

The ordinary image has dimension one in every coefficient degree k>=1.
Consequently D_(k,d) has dimension P_k(d-(2^k-1)). This also explains
why successive new sectors begin at weights 3,7,15,31,... .

## 5. Why the field computation gives integral lattices

Work over Z_(p), with p>d. For each k<=d the source and coefficient-ideal
spaces are images of explicit finite integer matrices: expand all
products of a_i and q_j of the two specified degrees. Their dimensions
over Q and over F_p are equal by sections 3--4, including the corresponding
formula for J from (6). Equality of these ranks means that both images
are p-saturated in their ambient free Z_(p)-module. Thus the quotients
A/J and J/IJ in these bidegrees are Z_(p)-free. This is an integral
saturation argument; rational dimensions alone would not suffice.

In coefficient degrees k>=d, J=A by section 2. It follows that C_d is
Z_(p)-free and M acts nilpotently on it. Also D_d is Z_(p)-free of rank
a_d, with M acting as zero by (4). Since the ordinary image M Z_(p)[M]
is free over Z_(p)[M], its surjection has a module splitting, and

    (J/IJ)_d = Z_(p)[M] direct-sum D_d,
    D_d = (Z_(p)[M]/(M))^(a_d).                         (10)

There is an explicit dual partition basis for each coefficient degree k:
take A_k^- s_lambda(z_1,...,z_k), where lambda ranges over partitions
of d-(2^k-1) of length at most k. Quotienting by ell_{1,...,k} sets
z_0=sum z_i, so the usual integral Schur basis supplies the indicated
basis. These functionals vanish on the ordinary lift M^(k-1)q_d, since
A_k^- restricted to the full hyperplane contains product_i z_i.
The perfect pairing and saturation above identify them with the dual
of D_(k,d). This statement is p-local; no unsupported global integral
Schur basis is asserted at primes <=d.

For explicit primal representatives one may instead use the dual monomial
symmetric basis. Put r=d-(2^k-1), and index its partitions by weakly
increasing tuples lambda of length k and sum r. On the full hyperplane,
A_k^- has leading monomial z_1 z_2^2 ... z_k^(2^(k-1)), with coefficient
one, in lexicographic order z_k>...>z_1. The leading monomial of m_mu
is z_1^mu_1...z_k^mu_k. Therefore the coefficient matrix

    [z_1^(1+lambda_1)...z_k^(2^(k-1)+lambda_k)]
                   (A_k^- restricted to ell_full=0)*m_mu

is triangular with diagonal entries one. Pairing the formal version of
(3b), with M in place of m, against A_k^- m_mu multiplies the column
by product_i e_i!, a unit at p>d. Its subtracted ordinary term pairs
to zero. Thus these relations are a Z_(p)-basis of D_(k,d), proving
the explicit assertion (3b) after mass specialization.

## 6. Specializing M=m preserves the actual image

On the finite free Z_(p)-module C_d, multiplication by M is nilpotent.
For a nonzero integer m, M-m is injective: if (M-m)v=0, repeated
application gives m^N v=M^N v=0 for large N, and freeness forces v=0.
The same holds in every lower weight used in the source construction.

The sequence 0 -> J -> A -> C -> 0 therefore stays injective on the left
after imposing M=m, in the weights in question. This proves that the
specialized module is the ACTUAL moment image J_m in A_m. It is not an
assumption that an image commutes with nonflat base change.

After additionally setting a_i=0, (10) becomes

    (J_m/I J_m)_d = Z_(p) direct-sum (Z_(p)/m)^(a_d).

Its ordinary image is the free rank-one group m Z_(p)L^d. The torsion
subgroup is exactly the specialization kernel.

Finally, for the original homogeneous Springer source in weight d<p|m,
one has d<m, so F(L)=0 introduces no relation in this weight. Use the
Newton coordinate changes

    a_i=(-1)^i p_i(c),
    q_j=p_j(f_1,...,f_j).

Each change is triangular with diagonal coefficient plus or minus i or j,
which is a unit for indices <=d. The power-sum formula for the shifted
roots is exactly the q_j formula in section 2 with M=m. Thus both the
coefficient ring and the FULL source in weights <=d are identified.
This proves (3) for the actual K-specialization defect.

The existing full-flag and good-prime parabolic transfer also carries the
new multiplicity formula to those sources; see section 15 for its exact
hypotheses. This is a corollary of the earlier integral Artin-basis theorem.

## 7. Primitive classes and a failed simplification

For k>=2 put D=2^k-1. If p>D and p|m, the moment relation

    R_k=q_1 q_2 q_4 ... q_(2^(k-1))-m^(k-1) q_D          (11)

has order exactly p^v_p(m) in (Delta_m)_D localized at p. It is the
primitive class in the newly appearing coefficient-degree-k sector.
Indeed, on the full hyperplane,

    A_k^- = product_(empty != B subset {1,...,k})
                          (sum_(i in B) z_i).

The coefficient of z_1 z_2^2 ... z_k^(2^(k-1)) is one: successively
the prescribed exponents force every subset to choose its largest index.
Its differential pairing with the first term of (11) is therefore
product_(i=0)^(k-1) (2^i)!, a p-local unit. The second term pairs to zero.

In particular, multiplication by q_1=f_1 does NOT annihilate the entire
defect. For p>7 dividing m, set R=q_2 q_4-m q_6, an ordinary relation of
weight six. Then

    q_1 R = R_3 modulo m Delta_m,

because their difference is m times q_1 q_6-m q_7. By IMG-39 and (11),
q_1 R has exact order p^v_p(m). This is a uniform counterexample, not
just a failed rank-seven test. Low-degree annihilation cannot be
promoted to a presentation of the whole defect.

There is, however, a sharp uniform replacement:

    f_1^2 Delta_m=0 for every m>=2.                    (11b)

The elementary Newton identity, written in the moment coordinates but
integral in the original generators, is

    q_1^2-m q_2=a_1^2-m a_2.

The right side is in the coefficient augmentation ideal, which acts
as zero on J/IJ. Thus q_1^2 acts as m q_2 on the defect. IMG-39's
m-annihilator proves (11b). No localization or bound on the weight is
needed. The degree-seven family above shows that replacing f_1^2 by f_1
is false. The same reasoning gives q_1^2 D=0 in the universal moment
model, using M D=0.

## 8. The complete degree-three piece, including the small primes

The earlier order detector can also be upgraded to a complete statement:

    (Delta_m)_0=(Delta_m)_1=(Delta_m)_2=0,
    (Delta_m)_3 = Z/t_m, for m>=3,
    t_m=m/2 if m=2 modulo 4, and t_m=m otherwise.        (12)

Delta_2=0 in every degree. This uses IMG-39's independent order detector;
the extra argument needed here is generation of the entire degree-three
kernel, not another lower bound.

Write b_j=binomial(m,j), h=gcd(m-1,2). The primitive ordinary relation
in degree two is

    A z_1^2-B z_2,  A=(m-1)/h, B=2m/h.

It specializes in the source to A c_1^2-B c_2, in IJ. Thus degree two
has no defect. For m>=4 the ordinary degree-three relation lattice is
the kernel of the row (m^3, m*b_2, b_3). It contains the primitive vector
(A,-B,0), coming from the degree-two relation times z_1. The quotient
by that vector is rank one and torsion-free. Therefore the full defect
in degree three is cyclic.

Let u=2 when m is even and u=1 otherwise. The gcd of the first two
entries of that row is m^2/u; the gcd of all three is m/gcd(m,3).
For a primitive generator of the remaining ordinary relation, the absolute
z_3 coefficient is m*gcd(m,3)/u. The class R_m from IMG-39 has z_3
coefficient -3m/gcd(m-2,6). It is therefore 3 times that generator if
m=1 modulo 3, and once that generator otherwise. In either case this
multiple is a unit modulo m. Since m annihilates the whole defect,
R_m generates it. Its independently proved order is exactly t_m.
For m=3 the root relation removes f_3, and f_1^3-3f_1 f_2 leaves the
single generator f_1 f_2=R_3, already known to have order three.

Multiplication by f_1 kills this degree-three piece, despite the higher
counterexample in section 7. For an integral identity put
q=(m-2)/gcd(m-2,6) and S=-q c_1 c_2+(3m/gcd(m-2,6))c_3. Expansion gives

    (f_1-c_1)R_m
      =2q(c_1^2 f_2-c_2 f_1^2)+S(f_1-c_1).

Every term on the right is in IJ; so is c_1 R_m. No division by a
nonunit is involved. This proves f_1(Delta_m)_3=0 in every rank.

## 9. Independent checks and remaining work

The [elementary-jet verifier](tools/verify_type_a_k_stable_multiplicities.py)
builds complete homogeneous image and coefficient-ideal lattices over Z
in its checked degrees, at prime-power and mixed m, without using moments
or apolarity. The [native audit](tools/audit_type_a_k_stable_multiplicities.py)
independently constructs the dual restriction constraints, compares exact
rational ranks with modular ranks, checks symbolic moment-trace identities,
and checks primitive subset-product coefficients. It imports neither Sage
nor the jet builder. The reports are
[jet verification](type-a-k-stable-multiplicities-verification.json) and
[native audit](type-a-k-stable-multiplicities-audit.json).
No finite test establishes the all-rank theorem.

A finite negative control is m=10, p=5, d=7<m: the actual 5-primary
piece is (Z/5)^3, whereas blindly extending (1) would predict (Z/5)^4.
This failure occurs before the root equation in degree m, so the
arithmetic restriction is not merely an overflow issue. The same verifier
records this as an out-of-range counterexample, not as evidence for an
unstated small-prime formula.

The remaining range is d>p outside the first extension above. Newton changes and the differential
pairing cease to be unimodular starting at d=p, but sections 10--11 resolve
that boundary. At d>=m the root relation also enters. Divided powers and
integral elementary generators must be retained.
For example, IMG-39's degree-three detector already distinguishes the
2-primary behavior when m is 2 modulo 4. Extending (3) unchanged is not
justified. The next task is to replace the perfect-pairing step by an
integral divided-power lattice comparison and determine its Smith factors.
That would address the remaining small-prime multiplicities structurally.

Published inputs are Newton identities, polarization, the differential
pairing, integral symmetric-polynomial/Schur bases, elementary unique
factorization, and standard exact-sequence base change. The moment trace,
the paired source/ideal calculation (6)--(9), and their application (3)
are the deductions written here. A targeted search is not a priority audit.

Related background, not a proof of this integral theorem:

- [Dufresne--Kohls, The separating variety for the basic representations of
  the additive group](https://arxiv.org/abs/1201.5933) studies graph closures
  and separating varieties in characteristic zero. Our moment translation
  is a basic additive action, but a set-theoretic graph description does
  not identify this integral specialization kernel or its scheme structure.
- [Ottaviani, Apolarity theory](https://people.dimai.unifi.it/ottavian/apolarity.pdf)
  gives the standard differential-pairing framework. The needed good-
  characteristic and saturation arguments are supplied explicitly above.

## 10. The prime boundary: compare the actual elementary and moment defects

There is a further uniform reduction at weight d=p. Let p|m, and let
D_(m,p)^mom denote the ACTUAL specialization kernel in weight p of

    Z_(p)[a_1,...,a_p][q_1,...,q_p]
        subset Z_(p)[a_1,...,a_p,L],
    q_j=m L^j+sum_(i=1)^j binomial(j,i)a_i L^(j-i),

with coefficient augmentation a_i=0 and no root equation. Then

    (Delta_m)_p localized at p = D_(m,p)^mom.            (13)

Comparison (13) alone does not compute its right-hand side. Section 11
supplies the separate saturation argument needed to finish d=p.

Here is a proof of the comparison. Below weight p all Newton changes
are integral isomorphisms over Z_(p). At weight p use the integral
coordinates

    b=(a_p-a_1^p)/p,       g=(q_p-q_1^p)/p.              (14)

Newton identities and the characteristic-p power-sum congruence show
that these are integral polynomials in the elementary coefficients and
elementary source generators, respectively. Their highest-coordinate
coefficients are units: b=-c_p plus lower terms, while g has coefficient
(-1)^(p-1) on f_p. These are the usual first Witt-coordinate corrections,
not a new divisibility theorem. Explicitly,

    g=b+(m-m^p)L^p/p
       +sum_(i=1)^(p-1) (binomial(p,i)/p)
                           (a_i-m^(p-i)a_1^i)L^(p-i).   (15)

All coefficients are integral when p|m. Kill the entire scalar coefficient
space in weight p, which is already in IJ. The remaining ambient space
has basis L^r times coefficient monomials of weight p-r, r>=1; it uses
only coefficient coordinates of index <p. The coefficient-ideal lattice
in this space is also the SAME for the moment and elementary models:
its generators are coefficients of positive weight times source pieces
of weights <p. The new coefficient a_p or b contributes only a scalar.

In this common ambient quotient the elementary image is obtained from
the moment image by adjoining g, with p g=q_p-q_1^p. Their quotient is
cyclic of order p: the ordinary moment image is m Z_(p)L^p, whereas g
has ordinary coefficient (m-m^p)/p, of valuation v_p(m)-1. Therefore the
map of these two image quotients onto their ordinary-image quotient is
an isomorphism. The induced specialization kernels are isomorphic.
This argument uses equality of the actual coefficient-ideal lattices;
it does not assume flatness of a ghost-coordinate change.

If m=p, perform this comparison first without imposing F(L)=0, including
f_p=F(L) in the auxiliary source. In weight p, quotienting by F removes
exactly Z_(p)f_p. Its ordinary image is the primitive generator L^p, so
it removes the free ordinary summand and preserves the defect. In lower
weights nothing changes. This proves (13) also at the root boundary.

For reference, one potentially troublesome Frobenius relation already
has an integral lift in the variable-mass moment source:

    [M^(p-1)(q_p-a_p)-(q_1-a_1)^p]/p
      =sum_(i=1)^(p-1) (binomial(p,i)/p)
                     M^(i-1) a_i (q_1-a_1)^(p-i) in IJ. (16)

This follows by direct binomial expansion. It handles this relation only;
it does not prove that all modular relations are generated by it. The
completeness argument below uses the entire hyperplane arrangement,
rather than an assertion that this one lift settles the boundary.

The standard ghost-coordinate background is in
[Hazewinkel, Witt vectors, Part 1, section 9](https://arxiv.org/abs/0804.3888).
Its identities explain (14); they do not compute our specialization defect.

## 11. Divided-power saturation in degree p

We prove that both full moment lattices J and IJ are p-saturated in
weight p. This finishes the theorem, by the mass-specialization argument
of section 6 and the actual-image comparison (13). In coefficient degrees
k>=p, section 2 already gives J=A; it also gives IJ equal to the span of
ambient monomials containing some a_i of positive index. Indeed, extract
one such a_i and write the remaining L powers using ML=q_1-a_1. This is
a saturated coordinate submodule. It remains to treat 1<=k<p, where
averaging over S_k is still legitimate.

Work over F_p and put V=span(L,x_1,...,x_k). The correct dual of Sym^p V
is the divided-power space Gamma^p(V^*), with the coefficient pairing.
The natural norm map

    N: Sym^p(V^*) -> Gamma^p(V^*),
    N(z^alpha)=product_i(alpha_i!) z^[alpha]

has kernel the span of the pure p-th powers, and cokernel another copy
of that Frobenius-twisted linear space. These statements follow on the
monomial basis: the only vanishing factorials in total degree p occur
at alpha=p e_i. The construction commutes with restriction to any W_B.
Thus the usual polynomial restriction arguments apply modulo the pure
p-th powers; they cannot be applied without this qualification.

### 11.1 The full source

Let gamma annihilate every Sym^p W_B. Its class in the norm cokernel
restricts to zero on every W_B, and hence belongs to the intersection
of the normal lines F_p ell_B. Two distinct normals already have zero
intersection. Thus gamma=N(P). For every B the restriction P|W_B is
the p-th power of a linear form h_B on W_B. On pairwise intersections
these linear forms agree, because taking p-th powers is injective on
a polynomial ring over F_p.

The h_B extend to a single ambient linear form h. Here the specific
subset arrangement matters: no three distinct normals
ell_B are linearly dependent, over ANY field. After normalizing their
z_0 coefficients to one, dependence would make three distinct vertices
of the binary cube collinear. On a coordinate where two vertices differ,
the line parameter must be zero or one, so the third is an endpoint.
For completeness, fix ell_0 and choose lifts of h_B. Pairwise agreement
allows the lift for B to be written h_0+t_B ell_0 modulo ell_B.
Compatibility of B and C gives

    (t_B-t_C) ell_0 in span(ell_B,ell_C).

The three-normal independence forces t_B=t_C. This gives a common h;
with just two hyperplanes, their pairwise agreement already suffices.

Consequently P-h^p vanishes on all W_B, so P-h^p=A_k Q. Conversely all
N(A_k Q) annihilate the full source. The norm is injective on A_k Q:
a nonzero p-th power of a linear form cannot be divisible by two
distinct linear factors. Averaging, or uniqueness of the decomposition,
retains the S_k-invariant subspaces. Therefore the dimension of the
source annihilator is exactly that in (6), also in characteristic p.

### 11.2 The full coefficient ideal

The class of an annihilator of IJ in the norm cokernel must, on every
proper W_B, involve only the output coordinates indexed by B. For B
empty this forces that class to be a multiple of the class of z_0^[p].
The functional z_0^[p] itself annihilates IJ: every source word in IJ
has a positive unshifted x_i factor and hence no pure L^p term. Subtract
that multiple. The remaining annihilator is N(P), with P invariant after
averaging. Its restriction to each proper W_B is an output-only
polynomial plus a linear combination of the pure p-th powers of the
outside variables.

At B empty, subtract a global pure p-th power so that P(0,z)=0.
Invariance makes this subtraction possible invariantly. Thus z_0 divides
P. At B={i}, divisibility by z_i eliminates all outside pure powers,
leaving c z_i^p, with the same c for every i. Subtract c z_0^p. Now
z_0 and all z_0-z_i divide P; neither subtraction changes N(P).

For any nonempty proper B choose i outside B and put s_B=sum_(j in B)z_j.
The restriction R=P(s_B,z) has the form

    R=Q(z_j:j in B)+sum_(j outside B) d_j z_j^p,

and is divisible both by s_B and by s_B-z_i. Setting z_i=s_B first
forces d_j=0 for outside j!=i and Q=-d_i s_B^p. Hence
R=d_i(z_i-s_B)^p. Since s_B and z_i-s_B are relatively prime, divisibility
by s_B forces d_i=0. Thus R=0 for every proper B, and A_k^- divides P.
For k>=2 the norm is injective on A_k^- Q, again because the product
has at least two distinct factors. We have obtained the precise direct sum

    (IJ)_(k,p)^perp = F_p z_0^[p] direct-sum
                       N((A_k^- F_p[z])_p^(S_k)), k>=2.  (17)

Its dimension equals the characteristic-zero answer (8). For k=1,
N(z_0 F_p[z_0,z_1]_(p-1)) loses exactly the line z_0^p, and the extra
line z_0^[p] restores it. Thus the dimension is again the answer of
section 4 and this coefficient degree has no defect.

### 11.3 Integral conclusion and representatives

The ranks of BOTH original integer image matrices now agree over Q and
F_p in weight p. This proves p-saturation, not merely equality of reduced
supports or a dimension calculation for a different source. In combination
with the lower-weight result, sections 5--6 apply unchanged: the actual
moment defect in weight p is (Z/p^v_p(m))^a_p. Comparison (13) gives the
same answer for the elementary Springer source, including m=p.

The basis (3b) also survives. For k>=2 every exponent e_i is at most p-1,
so its factorial remains a p-unit. The same triangular pairings against
N(A_k^- m_lambda) prove primitivity and independence. These annihilators
vanish on the ordinary lift M^(k-1)q_p. The extra divided-power line in
(17) accounts for that ordinary direction and does not add a defect.
Under (13) the moment relations are the same integral expressions in
the original f_j. At m=p, interpret q_p by Newton's identity with f_p=0.

The native audit checks the actual moment-word matrices at d=p, for
p=2,3,5,7 and every 1<=k<p, using coefficient pairing and separate
rational/modular elimination. It does not reuse the differential-pairing
matrix, which would be invalid here. These finite checks corroborate the
argument; the no-three-normal-dependence and divided-power arguments are
the general proof. Nothing here establishes saturation in degree p+1 or
higher by this argument alone. Section 13 treats the first mixed
Frobenius degree with an additional bilinear gluing argument.

## 12. Uniform relations for the first Witt correction beyond the range

The following are identities in the ACTUAL elementary-source defect, not
a multiplicity formula. They are valid for every odd prime p dividing m.
Use Newton moments q_j for every j, including j>=m, where they are
polynomials in f_1,...,f_(m-1) with f_m=0. Put

    n=p+1, r=(p+1)/2, alpha=(m-m^p)/p,
    g=(q_p-q_1^p)/p,
    R_(i,j)=q_i q_j-m q_(i+j),
    E_p=q_1 g-alpha q_n.

The R_(i,j) and E_p have zero ordinary image, so define defect classes.
For 2<=i<r put c_i=(-1)^i binomial(n,i)/p and put
c_r=(-1)^r binomial(n,r)/(2p). These are integers, and sum c_i=1.
Then

    p E_p=R_(1,p),
    (p+1) E_p=sum_(i=2)^r c_i R_(i,n-i),               (18)

and consequently

    (m/p) R_(1,p+1)+sum_(i=2)^r c_i q_1 R_(i,p+1-i)=0. (19)

Thus the degree-p+1 class contributed by q_1 g is already generated by
moment relations p-locally: p+1 is a unit. This does NOT prove that the
whole degree-p+1 defect is unchanged. The accompanying enlargement of
the coefficient ideal must also be controlled.

Here is an integral proof. For even n, the standard alternating-binomial
quadratic expression

    Q_n(v)=m v_n+sum_(i=1)^(r-1) (-1)^i binomial(n,i)v_i v_(n-i)
                   +(-1)^r binomial(n,r)v_r^2/2

is unchanged by simultaneous translation of all roots. Therefore
Q_n(q)=Q_n(a), where a_j=(-1)^j p_j(c) and a_0=q_0=m. This also follows
by expanding sum_(u,v)((L-t_u)-(L-t_v))^n and dividing the symmetric
sum by two; the displayed coefficients are integral. Let

    b=(a_p-a_1^p)/p,
    S=sum_(h=0)^(r-1) q_1^(2(r-1-h)) a_1^(2h).

The first Witt correction b is integral in the original coefficients.
Using q_1^2-a_1^2=m(q_2-a_2), the quadratic identity gives the EXACT
polynomial identity

    (m/p)(q_n-a_n) - n(m/p)(q_2-a_2) S
       - n(q_1 g-a_1 b)
       + sum_(i=2)^r c_i(q_i q_(n-i)-a_i a_(n-i))=0.   (20)

All divisions displayed here have already been justified integrally.
Modulo IJ, S becomes q_1^(p-1). The relation q_1^2=m q_2 in J/IJ,
together with m Delta_m=0, gives

    q_2 q_1^(p-1)=m^(p-1) q_n,
    q_1^(p+1)=m^p q_n.

For example, in the first equality write q_1^(p-1)=m^(r-1)q_2^(r-1)
and use that q_2^r-m^(r-1)q_n is a defect class; r-1>=1. Substitution
in (20), using sum c_i=1, proves the second equality of (18). The first
follows from p g=q_p-q_1^p. Multiplying by q_1, using (11b) and
m(q_2 g-alpha q_(p+2))=0, proves (19), since

    (p+1)alpha = m/p modulo m.

These reductions remain valid when the root equation makes the ordinary
image in the indicated degree zero. They use the full source and the
actual augmentation ideal, rather than a rational ghost-coordinate model.

The cubic combination in (19) is p-torsion in the universal moment
defect, by the quadratic identity and M D=0. Section 14 proves that it
vanishes, using a tensor lemma checked by exact integer certificates.
Completeness of the resulting relations in degree p+2 remains open;
neither (19) nor that vanishing is a full multiplicity formula.

## 13. The first mixed degree, with the coefficient enlargement retained

**Theorem.** For every odd p|m,

    (Delta_m)_(p+1) localized at p = (Z/p^v_p(m))^(a_(p+1)). (21)

The proof has two independent steps: saturation of the full variable-mass
moment lattices, and comparison of the actual elementary and moment
images. Merely adding the divided source generator g would not suffice.

### 13.1 Moment saturation in weight p+1

First take 2<=k<p. In degree p+1 the norm map has kernel and cokernel
isomorphic to V^*^(1) tensor V^*. Its kernel consists of the polynomials
sum b_ij z_i^p z_j. Distinct pairs (i,j) give distinct monomials, so this
description is injective. Its cokernel description follows from divided-
power comultiplication to Gamma^p tensor Gamma^1 followed by the pure-
power quotient. Both descriptions commute with restriction. All matrices
here have entries in F_p, so the Frobenius twist can be identified with
the original vector space for the following linear-algebra argument.

A bilinear form which restricts to zero on each of three hyperplanes
with independent normals is zero. Indeed its polynomial in two sets
of variables lies in each ideal (ell_i(z),ell_i(w)); after a coordinate
change, their intersection has no element of bidegree (1,1). Applying
this to ell_empty, ell_{1}, ell_{2} shows that a source annihilator has
zero norm-cokernel class. Write it N(P). On every W_B, P restricts to
a Frobenius-linear bilinear form T_B, and the T_B agree on intersections.

These invariant T_B glue to one global bilinear form. Split them into
symmetric and alternating parts, since 2 is invertible. Symmetric parts
are ordinary quadratic polynomials. Compatible polynomials on this
hyperplane arrangement always glue: inductively, a residual on the last
hyperplane vanishes on each previous intersection, hence is divisible
by the product of their distinct linear equations. Lift the quotient
and correct by that product. The equations are distinct because no
three normals are dependent, as proved in section 11.

For alternating parts, invariance under S_B times S_(B-complement)
forces the local form to be

    t_b (sum_(i in B) dz_i) wedge (sum_(i outside B) dz_i), b=|B|.

Within either block a transposition changes the sign of a coefficient,
so there are no additional invariant terms. Empty and full B give zero.
For nonempty B properly contained in a proper C, compatibility forces
t_b=t_c: on the intersection, the wedge of the sum over B and the sum
outside C is nonzero. Thus all t_b agree (there is only one for k=2).
The ambient form t dz_0 wedge sum_i dz_i has exactly these restrictions.
The symmetric part can be averaged to be invariant. Subtracting the
global Frobenius-linear polynomial now leaves P divisible by A_k.
The norm is injective on A_k times polynomials, by the same three-normal
bilinear test. Hence the source-annihilator dimension is again (6).

For IJ, its invariant norm-cokernel class is a multiple of the pure
z_0^[p+1] class. To check this, the empty subset first restricts its
bilinear form to a combination of z_0 tensor z_0, z_0 tensor sum z_i,
and sum z_i tensor z_0. Singleton subsets eliminate the last two terms.
The pure class itself annihilates IJ. Subtract it and write the rest N(P).
On each proper W_B, P is an output-only polynomial plus a Frobenius-
linear polynomial. At the empty subset subtract a global kernel element
so that z_0 divides P. At a singleton, divisibility by z_i leaves only

    a z_i^(p+1)+b z_i^p sum_(j!=i)z_j+c z_i sum_(j!=i)z_j^p.

Invariance makes a,b,c independent of i. Subtracting a combination of
z_0^(p+1), z_0^p sum z_i, and z_0 sum z_i^p kills all singleton
restrictions, without changing N(P). Thus z_0 and every z_0-z_i divide P.

Now let 2<=|B|<k and restrict to z_0=s_B=sum_(i in B)z_i. Choose an
outside i. Setting z_i=s_B and the other outside variables to zero
shows that the output-only part is itself Frobenius-linear. The entire
restriction is therefore a bilinear kernel polynomial. It vanishes on
the three hyperplanes with normals s_B, s_B-z_j, s_B-z_i, where j is
inside B. They are independent since |B|>=2 and i is outside. The
bilinear test forces the restriction to vanish. Hence A_k^- divides P;
the norm is injective on that product as well. The dimensions of both
annihilators equal the characteristic-zero answers.

Coefficient degree one is directly saturated: its source is spanned by
a_(p+1) and q_(p+1), and the coefficient of M L^(p+1) in the latter is
one; its ideal is spanned by a_(p+1). In coefficient degree k=p=d-1,
polarization is unavailable but an integral argument works. The only
ambient monomials not already visibly in J are

    x_r=M^(r-1) a_1^(d-r) L^r, 1<=r<=d.

Using q_2 gives x_r+2x_(r-1)=0 in A/J for r>=2. Using q_3 gives
x_r+3x_(r-1)=0 for r>=3. Thus x_2=0 and 2x_1=0, and all the x_r vanish
p-locally. For r<=d-1 these same relations are in IJ, because the
multiplier contains a positive coefficient; d>=4 allows r=3 here.
Consequently IJ is exactly the positive-coefficient coordinate submodule
p-locally. Degrees k>=d were already handled in section 2. This proves
both saturations in weight p+1, including the averaging boundary k=p.

Sections 5--6 now give the actual MOMENT defect (Z/p^s)^a_(p+1).

### 13.2 The two new coefficient-ideal directions do not kill old classes

Here m>p ensures p+1<m; no root equation occurs. Work over Z_(p) and
kill all scalar coefficient polynomials, already in the respective
augmentation ideals. The moment coefficient a_p is replaced by
a_1^p+p b in the elementary model, and g=b+H, where H is the positive-L
part in (15). All lower Newton coordinates agree. The common positive-L
moment ambient module embeds in the elementary one with index p,
in the direction bL. The only new ideal generators in weight p+1 are

    v=m bL=(m/p)(a_p-a_1^p)L,       w=a_1 H.

They both lie in the embedded moment ambient module, and p v,p w are
in the OLD moment coefficient ideal. The only additional source word
to consider is q_1 g. Identity (20), with p+1 inverted, places its
positive-L part in the old moment image. Thus, writing J_0,C_0 for
the old moment image and coefficient ideal in this ambient module,

    J_elementary=J_0+Z_(p)v+Z_(p)w,
    C_elementary=C_0+Z_(p)v+Z_(p)w.                    (22)

It remains to prove that v,w are independent modulo J_0. This is the
step that prevents silently losing old defect classes.

Let C=(A/J)_(p+1) in the VARIABLE-MASS moment model. It is Z_(p)-free
by 13.1, and M is nilpotent. Actual-image base change gives
C_m=C/(M-m)C for the mass-m ambient/image quotient. The exact sequence
for M-m and reduction modulo p give a Bockstein isomorphism

    C_m[p] -> ker(M:C/pC -> C/pC).                     (23)

For v its image is -(a_p-a_1^p)L modulo p, since M(a_p-a_1^p)L=0 in C.
The class a_1^p L vanishes p-locally by the k=p argument above, so this
is -a_p L. Put y=a_1 L^p. Expansion of the integral old source word
a_1(q_p-q_1^p-a_p+a_1^p) gives in C

    (M-M^p)y+p W(M)=0,
    W(M)=a_1 sum_(i=1)^(p-1) (binomial(p,i)/p)
                         (a_i-M^(p-i)a_1^i)L^(p-i).

The new vector w is represented at mass m by alpha y+W(M).
Its image under (23) is -(1-M^(p-1))y=-y modulo p. The last equality
follows from M(1-M^(p-1))y=0 modulo p and nilpotence of M.
The two classes a_p L and a_1 L^p are independent in coefficient degree
one of C/pC: eliminate a_(p+1) and M L^(p+1) using the two source rows,
leaving all a_i L^(p+1-i), 1<=i<=p, as a free basis. Therefore v,w give
independent order-p classes modulo J_0. Equation (22) induces an
isomorphism of the old and new specialization defects. This proves (21).

### 13.3 Basis correction and precise boundary

Only one factorial in the partition-basis pairing fails to be a unit:
the coefficient-degree-two representative R_(1,p). Replace it by
E_p=R_(1,p)/p as a CLASS, using the integral formula (18), whose factor
p+1 is a p-unit. Its diagonal pairing becomes (p-1)!, a unit; all other
diagonal pairings remain units. In the variable-mass moment module use
(p+1)^(-1) sum c_i R_(i,p+1-i); this avoids introducing M/p there. This proves the stated corrected basis
without dividing a polynomial arbitrarily in the actual source.

The root case m=p requires the separate comparison below. The proof
does not extend to p=2: the alternating/symmetric
split and the high-coefficient saturation use that 2 is invertible.
Section 14 treats one degree-p+2 cubic obstruction. Complete multiplicities
in that degree remain unresolved.

### 13.4 Retaining the root boundary m=p

First use the unrooted auxiliary elementary model through weight p+1,
with coefficient variables c_1,...,c_(p+1). Its elementary source
coordinates can be defined integrally by

    f_j=[u^j] sum_(i=0)^(p+1) (-1)^i c_i u^i (1+Lu)^(p-i), c_0=1.

Only j<=p+1 is used, so the negative exponent in the final term poses
no problem. Newton's identity or the logarithmic derivative of this
series gives exactly the moment translation formula with mass p.
In particular f_(p+1)=(-1)^(p+1)c_(p+1), a scalar, and f_p is the
usual monic F(L). Thus the proof of 13.2 applies to this auxiliary
model as well: it has the same first Witt coefficient and source
corrections, with independent a_(p+1).

Set c_(p+1)=0. In weight p+1 this removes only a scalar coefficient
direction, already in the augmentation ideal, and removes the redundant
scalar f_(p+1). It leaves J/IJ and its defect unchanged. The remaining
unrooted source is generated by f_1,...,f_p over Z_(p)[c_1,...,c_p].

Now impose F=f_p=0. In weight p+1 the ambient F-ideal is spanned by
F L and F c_1. An element a F L+b F c_1 belongs to the unrooted image
only if p divides a: its ordinary image is a L^(p+1), while the entire
ordinary source image in that weight is p Z_(p)L^(p+1): each positive
jet below p has coefficient divisible by p, and f_p f_1 attains p. Conversely,
p F L=F(f_1+c_1) is in the source and F c_1 is in IJ. Hence

    (F A intersect J)_(p+1)=Z_(p) F f_1+Z_(p) F c_1.

The second generator is already in IJ. The first maps to the primitive
generator p L^(p+1) of the unrooted ordinary image. Quotienting removes
that free ordinary summand and preserves the entire defect, including
its nilpotents and torsion. This proves (21) and the corrected basis
also when m=p. No analogous higher-weight intersection calculation is
being assumed.

## 14. Vanishing of the cubic obstruction and a first reduced-order class

**Theorem (with an exact finite tensor certificate).** For every odd
prime p|m, integrally in the actual elementary-source defect,

    (m/p) [q_1 q_(p+1)-m q_(p+2)] = 0.                (24)

This is an upper bound for the order of this specific class. It is not
a determination of its exact order or of every Smith factor in degree
p+2. In particular, no assertion that it generates a direct summand is
being made.

Let T_p=sum_(i=2)^((p+1)/2) c_i q_1 R_(i,p+1-i) in the universal
variable-mass moment defect, with the c_i of section 12. The quadratic
identity gives p T_p=(p+1)q_1 R_(1,p)=0. The last equality follows from
q_1^2=M q_2 modulo IJ and M D=0. We prove that the coefficient-degree-
three part of D in weight p+2 has no p-torsion. Thus T_p=0, without
dividing a torsion class by p. Its image in the actual elementary source
vanishes under every mass and coefficient specialization, regardless of
whether that specialization is flat. Equation (19) then proves (24).

### 14.1 The fixed tensor lemma and its certificates

For k=3 use the eight hyperplanes W_B of section 3. Consider tensors of
bidegree (1,2), written u_i v_j v_l with j<=l. All coordinate changes
use the same binary-subset equations in u and v.

The following two assertions hold in every characteristic other than two:

- Every S_3-equivariant family of local tensors on W_B, compatible on
  every pairwise intersection, is the restriction of an ambient invariant
  tensor.
- An ambient invariant tensor whose restriction to each proper W_B uses
  only indices in B is a multiple of u_0 v_0^2.

Here is the exact finite certificate, with its construction specified so
that it can be replayed. For each size b=0,1,2,3 take orbit sums of local
monomials u_i v_j v_l, 1<=i,j,l<=3, under S_b times S_(3-b). Transport
them to all subsets of that size. These give 28 family coordinates.
Ambient S_3 orbit sums, allowing index zero fixed by S_3, give 11
coordinates. Pairwise compatibility is imposed by substituting the
linear equation sum_B z_i=sum_C z_i, with a coefficient +1 or -1
as pivot, in both tensor factors. Discarding zero constraints gives
a 110 by 28 integer matrix. Its nonzero Smith factors are sixteen ones
and one two. The ambient restriction matrix has rank 11 and all eleven
nonzero Smith factors one, and its product with the compatibility
matrix is zero. Thus it identifies the full compatibility kernel over
Z[1/2] and over every field of odd characteristic.

For the second assertion, retain from each proper-subset restriction
all monomials with an outside index. This matrix on the 11 ambient
coordinates has ten nonzero Smith factors, all one, and kernel exactly
the coordinate u_0 v_0^2. These are integer lattice statements, not
floating-point ranks or tests at a few residue characteristics.

The [Sage verifier](tools/verify_type_a_k_stable_multiplicities.py) constructs
the matrices by polynomial substitution. The [native audit](tools/audit_type_a_k_stable_multiplicities.py)
constructs tensor coefficient dictionaries and performs separate integer
Smith calculations. Both reports contain the ranks, nonunit Smith
factors and identified kernel. This fixed tensor lemma is the computational
certificate used in the following all-prime deduction; it is not an
exceptional-group or rank census.

### 14.2 Apply the lemma to the actual full moment lattices

Take p>=5 and k=3<p. In degree p+2 the norm kernel is
V^*^(1) tensor Sym^2(V^*), embedded by u_i v_j v_l -> z_i^p z_j z_l.
It is injective since p+2<2p. The cokernel is the corresponding divided-
power tensor; Sym^2 and Gamma^2 are naturally isomorphic here since
2 is a unit. These descriptions commute with restriction.

A source-annihilator's norm-cokernel class vanishes on W_B for all B.
Four independent normals ell_empty,ell_{1},ell_{2},ell_{3} force this
tensor to vanish: the intersection of their paired-coordinate ideals
has no bidegree-(1,2) element. Its norm preimage therefore restricts
on each W_B to a Frobenius-quadratic tensor. The first tensor lemma
glues these tensors globally. Subtract the resulting norm-kernel
polynomial; the remainder is divisible by A_3. The norm is injective
on that product by the same four-normal test. This recovers precisely
the characteristic-zero source-annihilator dimension.

For IJ the second tensor lemma leaves only the pure ordinary functional
z_0^[p+2] in the norm cokernel. Subtract it. As before, subtract a kernel
polynomial at the empty subset so that z_0 divides the norm preimage P.
At a singleton, its Frobenius-quadratic restriction is invariant under
the two outside indices and divisible by the singleton variable. All
such terms are supplied by restrictions of the following seven global
kernel polynomials, where s=sum z_i and t_2=sum z_i^2:

    z_0^(p+2), z_0^(p+1)s, z_0^p s^2, z_0^p t_2,
    z_0^2 sum z_i^p, z_0 sum z_i^(p+1), z_0 s sum z_i^p.

Their restrictions span by expanding s=z_i+sum_(j!=i)z_j; the only
division needed is by two. Subtract them to kill all singleton
restrictions. Now z_0 and all three z_0-z_i divide P.

For a proper two-element B={j,l} and its outside index i, set
s_B=z_j+z_l. The output-only part of P|W_B is itself Frobenius-quadratic,
by setting z_i=s_B. Thus the whole restriction is such a tensor and
is divisible by z_j,z_l,z_j+z_l,z_j+z_l-z_i. In tensor coordinates,
the bidegree-(1,2) intersection of these four paired ideals is the
one-dimensional span

    (u_j v_l-u_l v_j) (v_j+v_l-v_i).

To verify this, first impose z_j,z_l,z_j+z_l-z_i as three coordinate
normals. The only monomials are u_1 v_2 v_3,u_2 v_1 v_3,u_3 v_1 v_2;
the fourth equation leaves their first-minus-second combination. It
changes sign when j and l are interchanged. The required invariance
therefore forces it to vanish. All proper-subset restrictions of P
are zero, so A_3^- divides P. The norm is injective on this product
as well, since it includes the four independent normals above.

Both annihilator dimensions now agree with characteristic zero. The
original integer source and coefficient-ideal images are p-saturated
in coefficient degree three and weight p+2. Therefore that part of
J/IJ, and hence D, is p-torsion-free.

For p=3 polarization in coefficient degree three is unavailable. The
two implementations instead expand ALL moment words in weight five
and coefficient degree three directly. The source has rank 16 and
nonunit Smith factor 2; the coefficient ideal has rank 15 and nonunit
Smith factor 2. Both are 3-saturated. This is a finite universal-mass
certificate, so covers every mass m; it is not a check at one rank.
It finishes T_3=0 and thus (24) for all odd primes.

The next unresolved step is completeness of the first Witt relations
in degree p+2, including every coefficient-ideal direction. Formula
(24) explains a reduced-order direction, but it does not by itself
count the surviving factors or prove the suggested first correction.

## 15. Formal full-flag and good-prime parabolic corollaries

Retain the SAME defining representation, literal integral root-polynomial
source and transported homogeneous grading. Write [j]_q=1+q+...+q^(j-1)
and A(q)=sum a_d q^d from (2). By the already proved
[integral Artin transfer, section 7](TYPE_A_K_INTEGRAL_RECURRENCE.md), the
SL_m/B defect is a direct sum of projective defects shifted by the
lengths in S_(m-1). Therefore in the proved degree range its p-primary
piece has only factors p^v_p(m), with multiplicity

    [q^d] A(q) product_(j=1)^(m-1) [j]_q.             (25)

For a proper parabolic with Levi block sizes n_1,...,n_s summing to m,
put g=gcd(n_1,...,n_s). If p|m and p does NOT divide g, the parabolic
transfer in the same earlier note gives multiplicity

    [q^d] A(q) product_(j=1)^(m-1)[j]_q
                     / product_(i=1)^s product_(j=1)^(n_i)[j]_q. (26)

All factors again have order p^v_p(m). To see that no smaller factors
can be hidden in this statement, the full-flag decomposition contains
the unshifted parabolic defect as a summand. A smaller factor in a
degree in the range would already occur in that full-flag degree.
The multiplicities then satisfy convolution with the Levi length
polynomial; its constant term is one, giving (26) by induction on d.
Interpret the quotient as a formal power series, and use it ONLY in
the proved range d<=p (also d=p+1 for odd p). It is not a claim about
the complete Hilbert series in higher degrees.

These corollaries require no averaging or factorial inversion beyond
the precise hypotheses of the established transfer. At primes dividing
every Levi block size, that transfer is unavailable and (26) is not
claimed. The full-flag formula (25) has no such additional exclusion.
