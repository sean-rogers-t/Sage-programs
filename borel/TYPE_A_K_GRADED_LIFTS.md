# Full-source overflow products: homogeneous descent and exact certificates

**Historical note, superseded on 28 October:** the uniform identities
left open below are now proved by [an integral trace recurrence](TYPE_A_K_INTEGRAL_RECURRENCE.md).
The homogeneous descent and m=7 certificates remain valid. Their finite
audit's negative uniform-proof flag records the scope of that method.

27 October 2026. Continuation of IMG-37/38. The defining SL_m/P1 Springer
K-image and normalization remain those of [the specialization theorem](TYPE_A_K_SPECIALIZATION.md).
The Cayley-fiber program stays shelved. No new API or PDF is involved.

## 1. Outcome and limits

The requested uniform conclusion is

    m f_i f_j in I J,  1<=i<=j<m, i+j>=m.                 (1)

For prime-power m these are precisely the overflow relations in the
[complete annihilator criterion](TYPE_A_K_ANNIHILATOR.md), up to p-local
units. For mixed m, that criterion uses the monomial lifts g_i instead;
(1) alone is then not the whole overflow test.

**At this historical stage, (1) was not proved in arbitrary rank.** We proved a uniform descent
to a homogeneous full-source problem and used it to produce exact
certificates covering every quadratic relation for m=3,4,7. The first
two are controls. The m=7 result is new finite coverage: it proves
7 Delta_7=0 without computing the whole image module, whose earlier
standard-basis calculation exceeded the available practical memory.
Combined with IMG-38, the full Delta_7 has exponent exactly 7. Its full
multiplicities have not been calculated here.

This is not an all-rank overflow theorem. The structural advance is that
bounded homogeneous identities suffice and can be checked without
trusting a normal-form computation or a source-degree truncation guess.

## 2. Allowing the constant coefficient is justified

Temporarily use the full coefficient ring

    S^u=Z[c_1,...,c_m],
    A^u=S^u[L]/F(L),
    J^u=S^u[f_1,...,f_(m-1)] subset A^u,
    C^u=A^u/J^u.

Here F(X)=sum_(k=0)^m (-1)^k c_k X^(m-k), c_0=1, and

    f_j=sum_(k=0)^j (-1)^k binomial(m-k,j-k)c_k L^(j-k).

Give L weight 1 and c_k weight k. All these rings and C^u are graded,
concentrated in nonnegative weights. The element c_m-1 is a
nonzerodivisor on C^u: if (c_m-1)x=0 and x is nonzero, the lowest
nonzero homogeneous component of x contributes its negative to the
lowest component of the product, with nothing to cancel it.

The short exact sequence 0 -> J^u -> A^u -> C^u -> 0 therefore remains
injective on the left after imposing c_m=1. Indeed the possible kernel
is Tor_1^(S^u)(C^u,S^u/(c_m-1)), which is exactly the kernel of
c_m-1 on C^u. The resulting image is precisely the original J, since
the original generators and coefficient ring are recovered.

Consequently J/IJ is J^u/(c_k-binomial(m,k):1<=k<=m)J^u.
This uses an injectivity proof, not an assumption that images commute
with arbitrary coefficient specialization. The same assertion for the
ambient ring is immediate from its free S^u-module presentation.

## 3. Integral translation moves the target to the origin

Simultaneously translate all roots and L by 1. Algebraically this is
the integral automorphism

    L -> L+1,
    c_k -> sum_(r=0)^k binomial(m-r,k-r)c_r.

Its inverse is translation by -1. The identity
F_(translated coefficients)(L+1+z)=F_c(L+z) proves that it preserves
F(L)=0 and fixes every f_j. Thus it preserves J^u. It takes the ideal
(c_k-binomial(m,k)) to (c_1,...,c_m): the displayed coefficient change
is triangular with diagonal entries one and the indicated constant terms.

We obtain an isomorphism of the ACTUAL specialized image and ordinary
image maps, hence of their kernels:

    Delta_m = ker(J^u/(c_1,...,c_m)J^u -> Z[L]/L^m).       (2)

The map on the right sends f_j to binomial(m,j)L^j. Formula (2) is
graded. It retains torsion and the full scheme algebra; it is not a
statement after tensoring with Q or replacing the image by its saturation.

This is also not the trace-zero SL_m cohomological specialization c_1=0.
The nonzerodivisor argument above uses c_m-1, which has constant term -1;
it does not justify imposing c_1=0 on C^u. Earlier cohomological and
K-theoretic defect calculations therefore need not agree.

## 4. A uniform bounded-degree proof target

Put S^u_+=(c_1,...,c_m). By (2), (1) is equivalent to

    m f_i f_j in S^u_+ J^u.                              (3)

This is a homogeneous membership question in degree i+j. In particular,
all overflow products lie in degrees m,...,2m-2, in every rank.
For prime-power m the lower-degree quadratic relations in the earlier
criterion have degree i+j as well. Therefore checking their integral
lifts through degree 2m-2 settles the ENTIRE annihilator in that rank:
the ordinary kernel is generated as an ideal by those quadratic
relations, and multiplication preserves their annihilation.

This does not assert that Delta_m itself vanishes above degree 2m-2.
It is a bound on degrees of ideal generators relevant to the annihilator.
No truncation of an unspecified infinite source has been assumed.

For exact computation eliminate c_m, rather than any source generator:

    c_m=(-1)^(m+1) sum_(k=0)^(m-1) (-1)^k c_k L^(m-k).

Thus A^u is the polynomial ring Z[L,c_1,...,c_(m-1)]. Use its weighted
monomials in each degree as an integral basis. The degree-d source image
is generated inductively by

    c_k (J^u)_(d-k),  1<=k<=m,
    f_k (J^u)_(d-k),  1<=k<m,

with (J^u)_0=Z. The coefficient ideal in degree d is generated by
the first set only. This is an exhaustive induction on weight.

Crucially, a POSITIVE membership certificate needs even less: each
auxiliary polynomial only needs an explicit integral construction from
these products in smaller degree. It is not necessary to certify that
the auxiliary polynomials form a complete basis. This observation makes
the independent replay independent of Hermite normal form algorithms.

## 5. What the finite certificates actually prove

The [builder](experiments/type_a_k_graded_lifts.py) stores forward integer
linear combinations constructing each auxiliary row from c_k and f_k
times earlier rows. For each pair i<=j it stores an integer identity
expressing the target using c_k times earlier source rows.

For i+j>=m the target is m f_i f_j. For i+j<m it is

    m (D f_i f_j-N f_(i+j)),
    N/D=binomial(m,i)binomial(m,j)/binomial(m,i+j)

in lowest terms. For m=p^s the denominator D is prime to p, by the
valuation inequalities in the preceding note. Dividing by D and the
normalizing source units is therefore legitimate over Z_(p). Together
with the known absence of torsion at other primes, these identities
prove m Delta_m=0 for each certified prime-power rank.

The [independent auditor](tools/audit_type_a_k_graded_lifts.py) uses only
Python integer polynomial arithmetic. It imports neither Sage, Singular,
SymPy nor the builder. It reconstructs the coefficients and f_j from
their defining formulas, verifies every auxiliary construction, verifies
every final identity, and checks that every pair was included exactly
once. It separately verifies translation invariance and the inverse
coefficient substitution by exact binomial sums through m=24. The
all-rank translation argument is the polynomial identity in section 3. A final identity is in S^u_+J^u because every coefficient candidate
has an explicit factor c_k and an inductively verified source factor.

The [audit report](type-a-k-graded-lifts-audit.json) covers m=3,4,7.
The compressed certificates retain all integer identities. None is a
one-coefficient slice, numerical evaluation, rational membership test,
or assertion based only on a Smith-factor list.

## 6. What F(L)=0 does not establish by itself

Monic division reduces f_i f_j to an ambient polynomial of degree <m
in L. At ordinary specialization its coefficients vanish for i+j>=m,
so it lies in I A. This proves that it represents a defect class, not
that m times that class is zero. Even expressing the remainder in the
linear span of 1,f_1,...,f_(m-1) can require p-denominators; the FULL
multiplicative source matters.

The homogeneous reduction eliminates this ambiguity from the finite
certificates: their coefficients are explicitly integers and their
source factors are explicitly constructed. What was still missing was a
formula for those identities valid for ALL m,i,j, or a structural
integral argument that ensures their existence. A successful Hermite
calculation in several ranks is not such an argument.

Status of inputs: standard grading, Tor, translation and integral linear
algebra are general tools. The ordinary source relation presentation is
the preceding repository deduction. The homogeneous descent and finite
full-source identities are established here, with no external review or
literature-priority claim. No literature classification is asserted.

**Historical next target, now settled by IMG-39:** seek an induction in the homogeneous model which constructs
the coefficient-ideal identities integrally. It must control divisibility
at p, not just solve over Q. IMG-39 supplies such a construction using
difference trace, including the mixed-prime relations. The finite
certificates in this note remain independent controls of that theorem.
