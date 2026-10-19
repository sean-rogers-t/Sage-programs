# Type D units from an ordinary Grassmannian over F2

**Status: complete repository derivation of every normalized type-D
unit multiplicity, 19 October 2026.** The proof supplies the integral
mod-2 comparison, a closed formula for odd tails, and an exact Gaussian
intersection count for all even tails. It does not determine the higher
type-D Smith multiplicities or prove the full square-sector Smith
comparison. These are unpublished deductions from the identified
Schubert inputs, without an originality claim or external expert review.

## 1. Setup and the exact result

Put
\[
X_D=\operatorname{OG}(a,2n),\qquad n=a+k,\qquad a\ge1,\quad k\ge2.
\]
The fixed normalized generators are
\[
q_i=e_i(x_{a+1}^2,\ldots,x_n^2),\qquad
p=x_{a+1}\cdots x_n,\qquad p^2=q_k.
\]
Their polynomial degrees are \(2i,k\); polynomial degree \(d\)
means cohomological degree \(2d\). Let \(u_D(a,k;d)\) count the
nonzero Smith factors equal to 1 in the degree-\(d\) integral image.
Here \(q_k\) may be included as a redundant square-only generator.

Use the ordinary Grassmann ring
\[
R_{a,k}=\mathbb F_2[z_1,\ldots,z_k]^{S_k}/
(h_{a+1},\ldots,h_{a+k}),\qquad \deg z_j=1.
\tag{1}
\]
It has the Schur basis \(s_\mu\), \(\mu\subseteq a^k\).
For \(a=1\), interpret \(R_{0,k}=\mathbb F_2\).
Write \(F(f)=f^2\) for Frobenius, and set a homogeneous component
with a negative or nonintegral degree equal to zero.

**Exact rank theorem.** The number of type-D unit factors is
\[
\boxed{
u_D(a,k;d)=\dim_{\mathbb F_2}\left(
F(R_{a,k}^{d/2})+
e_k F(R_{a-1,k}^{(d-k)/2})\right).
}\tag{2}
\]
The two subspaces are calculated in the same ring \(R_{a,k}\).
The second expression means: compute the square in \(R_{a-1,k}\)
and apply the injective linear map
\[
T:R_{a-1,k}\longrightarrow R_{a,k},\qquad
T(s_\nu)=s_{\nu+(1^k)}.
\tag{3}
\]
Thus \(T\) is multiplication by \(e_k\) on a lift, not a
claimed ring embedding of the smaller Grassmannian.

In particular, the square-only D sector and the normalized type-B
image on \(\operatorname{OG}(a,2n+1)\) have **the same unit
multiplicities in every rank and degree**. Equality of their higher
Smith factors does not follow from this mod-2 argument.

## 2. Why the Chern subring really is the Grassmann ring

Let \(c_i=c_i(Q)\), where \(Q=\mathbb C^{2n}/S\), and reduce
cohomology modulo 2. The published special Schubert formula is
\[
c_i=\tau_i\ (i<k),\qquad
c_k=\tau_k+\tau'_k,\qquad
c_i=2\tau_i\ (i>k).
\tag{4}
\]
See [Buch--Kresch--Tamvakis, *A Giambelli formula for even orthogonal
Grassmannians*, equation (7)](https://math.umd.edu/~harryt/papers/oggiam.pdf).
Consequently \(c_i=0\) modulo 2 for \(i>k\).
The identity \(c(S)c(Q)=1\), with \(\operatorname{rank}S=a\),
then gives the relations \(h_j(c_1,\ldots,c_k)=0\) for \(j>a\),
where the arguments denote **elementary** generators: the complete
functions are the coefficients of
\((1+c_1t+\cdots+c_kt^k)^{-1}\) over \(\mathbb F_2\).
Thus there is a well-defined map
\[
\Phi:R_{a,k}\longrightarrow H^{2*}(X_D;\mathbb F_2),\qquad
e_i(z)\longmapsto c_i.
\tag{5}
\]

This map is injective, integrally justified by the Schubert basis.
For \(\mu\subseteq a^k\), put \(\lambda=\mu'\subseteq k^a\).
The dual Jacobi--Trudi formula gives
\[
\Phi(s_\mu)=\det(c_{\lambda_i-i+j}).
\tag{6}
\]
Every part of \(\lambda\) is at most \(k\), so the extra raising
factors in BKT Theorem 1 are absent:
\(\lambda_i+\lambda_j\le2k<2k+j-i\) whenever \(i<j\).
Its normalization exponent is also zero. Therefore (6) is the
reduction of the integral class
\[
[Y_\lambda]=
\begin{cases}
\tau_\lambda&\text{if \(k\) is not a part of \(\lambda\)},\\
\tau_\lambda+\tau'_\lambda&\text{if \(k\) is a part of \(\lambda\)}.
\end{cases}
\tag{7}
\]
These vectors have disjoint, nonempty supports in the integral typed
Schubert basis. Their reductions are linearly independent over
\(\mathbb F_2\). This proves the asserted injection and identifies
the entire Chern-generated subring modulo 2 with (1). The transpose
\(\lambda=\mu'\) is essential: the BKT determinant uses the
\(c_i\) as complete-function indices, whereas (1) sends the
elementary symmetric functions to the \(c_i\).

## 3. Proof of the rank theorem

The normalized Chern formula and Euler formula give
\[
q_i=c_i^2+2\sum_{j=1}^i(-1)^j c_{i+j}c_{i-j}
\equiv c_i^2\pmod2,
\qquad
p=\tau_k-\tau'_k\equiv c_k\pmod2.
\tag{8}
\]
The Euler sign has no effect modulo 2. Every input monomial is a
square-only monomial or \(p\) times one, because \(p^2=q_k\).
Under (5), their reductions consequently span
\(F(R_{a,k})+e_kF(R_{a,k})\), with the degree conventions in (2).

For clarity, the smaller ring in (2) introduces no assumption.
The ordinary Grassmann projection
\(R_{a,k}\to R_{a-1,k}\) kills exactly the span of the Schur
classes with first part \(a\). By the vertical Pieri rule,
\[
e_k s_\nu=
\begin{cases}
s_{\nu+(1^k)}&\nu_1<a,\\
0&\nu_1=a
\end{cases}
\quad\text{in }R_{a,k}.
\tag{9}
\]
Thus its kernel equals the kernel of multiplication by \(e_k\).
The latter multiplication factors as the projection followed by (3).
The projection is a ring map and hence commutes with Frobenius, so
\(e_kF(R_{a,k})=T(F(R_{a-1,k}))\), degree by degree.

Finally, reducing an integer coefficient matrix modulo 2 gives rank
equal to its number of odd nonzero Smith factors. Every nonzero
normalized D factor is a power of 2 by the localization argument in
[TYPE_D_SMITH.md](TYPE_D_SMITH.md). Its odd factors are therefore
precisely its unit factors. The injection (5) shows that no additional
relations are introduced when computing this rank in (1). This proves
(2).

The type-B unit theorem computes the rank of the same Frobenius map
on \(R_{a,k}\). This proves the square-sector assertion as well.
It is an equality of unit counts, not an assertion that the integral
square-sector lattices have already been identified.

## 4. A closed result when the tail dimension is odd

For \(b\ge0\), write \(\mathcal U_{b,k}(q)\) for the type-B unit
polynomial, taking \(\mathcal U_{0,k}=1\). Explicitly, for \(b>0\),
put \(N=b+k\), \(L=\min(b,k)\), \(r=\lfloor N/2\rfloor\),
\(c=\lceil N/2\rceil\), and set
\[
\mathcal U_{b,k}(q)=
\sum_{j=0}^{\lfloor L/2\rfloor}
q^{j(c-L+j)}{r\brack j}_q{c-j\brack L-2j}_q.
\tag{10}
\]
Its variable records half the even polynomial degree.

If \(k\) is odd, the two terms in (2) have opposite polynomial
degree parities. Hence
\[
\boxed{
\sum_{d\ge0}u_D(a,k;d)z^d
=\mathcal U_{a,k}(z^2)+z^k\mathcal U_{a-1,k}(z^2).
}\tag{11}
\]
This is an unconditional closed all-rank type-D unit formula. It
does not determine the higher factors in either degree parity.

## 5. Explicit intersection classes for an even tail

Suppose \(k=2h\), and put \(u=\lceil a/2\rceil\),
\(v=\lfloor a/2\rfloor\). Define
\[
I_{a,2h}(q)=\sum_m q^m\dim_{\mathbb F_2}
\left(F(R_{a,2h}^m)\cap
e_{2h}F(R_{a-1,2h}^{m-h})\right).
\tag{12}
\]
The exact rank theorem says
\[
\sum_m u_D(a,2h;2m)q^m
=\mathcal U_{a,2h}(q)+q^h\mathcal U_{a-1,2h}(q)-I_{a,2h}(q).
\tag{13}
\]

**Proved intersection classes.** For every partition
\(\mu\subseteq a^h\) with \(\mu_h\ge v+1\),
\[
\boxed{
s_\mu^2=e_{2h}s_{\mu-(1^h)}^2\ne0
\quad\text{in }R_{a,2h}.
}\tag{14}
\]
The classes on the left are linearly independent as \(\mu\) varies.
They give the coefficientwise lower bound
\[
\boxed{
I_{a,2h}(q)\ \ge\
q^{h(v+1)}{h+u-1\brack h}_q.
}\tag{15}
\]

### Proof using the two-quotient

The Schur expansion of \(s_\mu(z)^2=s_\mu(z_1^2,\ldots,z_{2h}^2)\)
has support on partitions \(\nu\) with empty 2-core. Its coefficient
on such a partition is, modulo 2,
\(c_{\alpha\beta}^{\mu}\), where \((\alpha,\beta)\) is the
2-quotient of \(\nu\). This is the classical SXP rule; see
[Wildon, *A generalized SXP rule proved by bijections and involutions*, equation (1)](https://arxiv.org/pdf/1508.07030).
Signs disappear over \(\mathbb F_2\).

Pad \(\nu\) to length \(2h\), and use beta numbers
\(\nu_i+2h-i\). For empty 2-core there are \(h\) even and
\(h\) odd beads. Dividing the even beads by 2, and subtracting 1
then dividing the odd beads by 2, gives the beta sets for
\(\alpha,\beta\). The bound \(\nu\subseteq a^{2h}\) is
equivalent to
\[
\alpha\subseteq u^h,\qquad\beta\subseteq v^h.
\tag{16}
\]
Also \(\nu_{2h}>0\) is equivalent to \(\alpha_h>0\).
In that case the quotient of \(\nu-(1^{2h})\) is
\((\beta,\alpha-(1^h))\). Subtracting 1 swaps the two bead
parities; since their counts are equal, empty 2-core is preserved.

If \(\nu_{2h}=0\), then \(\alpha_h=0\).
A necessary Littlewood--Richardson bound is
\(\mu_h\le\alpha_h+\beta_1\le v\), contradicting our
hypothesis. One direct justification of the bound is to realize
\(s_\beta\) as a constituent of the product of elementary functions
indexed by its columns. This is a product of \(\beta_1\)
vertical-strip operations, each adding at most one cell to a given
row. Schur positivity ensures that every constituent of
\(s_\alpha s_\beta\) satisfies that bound. Thus the left side
of (14) has no support on partitions without a full \(2h\)-column.

If \(\nu_{2h}>0\), the determinant twist in \(h\) variables gives
\[
c_{\alpha\beta}^{\mu}
=c_{\alpha-(1^h),\beta}^{\mu-(1^h)}.
\tag{17}
\]
All partitions in this equality have at most \(h\) parts. Apply
the SXP rule to the right side of (14), and use the vertical Pieri
rule (9). Equations (16)--(17) identify every surviving coefficient.
Partitions with nonempty 2-core have coefficient zero on both sides.
This proves the equality.

To prove independence, for each candidate \(\mu\) choose the target
partition whose 2-quotient is
\[
\alpha=\mu-(v^h),\qquad\beta=(v^h).
\]
It fits (16). For another candidate \(\widetilde\mu\), its
coefficient in this row is
\(c_{\alpha,(v^h)}^{\widetilde\mu}
=\delta_{\mu,\widetilde\mu}\): in \(h\) variables multiplication
by \(s_{(v^h)}\) is the determinant twist. The corresponding
coefficient submatrix is the identity. Finally, subtracting
\(((v+1)^h)\) from the candidate partitions leaves the rectangle
\((u-1)^h\), giving (15).

### The general upper bound closes the intersection formula

**Theorem.** For every \(a\ge1\) and \(h\ge1\),
\[
\boxed{
I_{a,2h}(q)=q^{h(v+1)}{h+u-1\brack h}_q.
}\tag{18}
\]
Combining this with (13) gives the closed all-rank formula
\[
\boxed{
\sum_m u_D(a,2h;2m)q^m
=\mathcal U_{a,2h}(q)+q^h\mathcal U_{a-1,2h}(q)
-q^{h(v+1)}{h+u-1\brack h}_q.
}\tag{19}
\]
The explicit classes (14) form a basis of the entire intersection.
Together (11) and (19) determine every normalized D unit multiplicity.

Here is the missing upper bound. On root polynomials in characteristic
2, use the derivation \(D(z_i)=z_i^2\). It preserves symmetric
polynomials, kills every square, and satisfies
\[
D(e_k)=e_1e_k,\qquad
D(h_j)=e_1h_j+(j+1)h_{j+1}.
\tag{20}
\]
The first identity uses exactly \(k\) roots. The second follows from
the complete-function generating series. Since every \(h_j\) with
\(j>a\) vanishes in \(R_{a,k}\), its defining ideal is stable
under \(D\), which therefore descends to that ring.

If \(x=F(f)=e_kF(g)\) lies in the intersection (12), applying
\(D\) yields \(e_1e_kF(g)=0\). By the injection (3), this implies
\(e_1F(g)=0\) in \(R'=R_{a-1,2h}\). Consequently multiplication
by \(e_k\) identifies the intersection with a subspace of
\[
\ker\left(e_1:\operatorname{im}(F:R'\to R')\to R'\right),
\tag{21}
\]
followed by a degree shift \(h\) in (12).

We now calculate the kernel in (21), including its rank over
\(\mathbb F_2\), grading \(F(R'^s)\) by \(s\).
The case \(a=1\) is addressed below; first assume \(a\ge2\). Put
\[
n'=a+2h-1,\quad r'=h+u-1=\lfloor n'/2\rfloor,
\quad c'=h+v=\lceil n'/2\rceil.
\]
Let \(N\) be a maximal-rank square-zero operator on
\(V=\mathbb C^{n'}\), and define
\[
W=\{U\in\operatorname{Gr}(2h,V):NU\subset U\},\qquad
Z=\{U\in W:\operatorname{rank}(N|_U)\le h-1\}.
\tag{22}
\]
The rank of \(F\) is the degreewise Betti rank of \(W\), by
[the unit-layer proof](TYPE_B_UNIT_SMITH.md). This statement holds
also when \(2h>n'/2\): annihilator duality replaces the Grassmann
plane dimension by its complement and preserves the square-zero
Jordan type.

The rank of \(e_1F\) is the degreewise Betti rank of \(Z\).
We justify this additional identification instead of assuming that a
rational rank remains unchanged modulo 2. Apply
[Wildon's Theorem 1.1](https://arxiv.org/pdf/1508.07030) with
\(r=2\) and \(\tau=(1)\). The Schur expansion of
\(s_{(1)}(s_\lambda\circ p_2)\) is supported on partitions
with 2-core \((1)\). Over \(\mathbb Z\), its coefficient on a
partition \(\nu\) with 2-quotient \((\alpha,\beta)\) is exactly
\[
\varepsilon_\nu c_{\alpha\beta}^{\lambda},\qquad
\varepsilon_\nu\in\{1,-1\},
\]
where the sign depends on the target row \(\nu\), not on the
source column \(\lambda\). The coefficients are zero on all other
cores. Thus the signs can be removed by integral row operations;
modulo 2 they disappear. With \(2h\) beads, the core \((1)\) has
\(h+1\) even beads and \(h-1\) odd beads. Truncation to
\((a-1)^{2h}\) therefore gives, up to these row signs, precisely
the integral row vectors of
\[
s_\alpha s_\beta,\qquad
\alpha\subseteq(v-1)^{h+1},\quad
\beta\subseteq u^{h-1},
\tag{23}
\]
in the ordinary Grassmann Schur coefficient lattice for
\(\operatorname{Gr}(2h,n')\). The quotient sizes total the
source degree: the extra box is in the 2-core and contributes no
extra shift to (23). Here the integral matrix is defined columnwise
on the rectangle Schur basis by
\(s_\lambda\mapsto s_{(1)}(s_\lambda\circ p_2)\), followed by
rectangle truncation. It is not asserted to be an integral ring
endomorphism of the finite Grassmann quotient. Its reduction is
exactly the linear map \(e_1F\).

For \(a\ge2\), put \(I=\operatorname{im}N\),
\(K=\ker N\), and consider the incidence variety
\[
Y=\{(F_0,U):F_0\subset I,\ \dim F_0=h-1;
\ F_0\subset U\subset N^{-1}F_0,\ \dim U=2h\}.
\tag{24}
\]
It is a Grassmann bundle with base \(\operatorname{Gr}(h-1,r')\)
and fiber \(\operatorname{Gr}(h+1,c')\). Its middle ambient
bundle \(N^{-1}F_0/F_0\) has total Chern class 1, as follows from
\(0\to K/F_0\to N^{-1}F_0/F_0\to F_0\to0\).
Thus its integral cohomology is the tensor product of the two
Grassmann rings. Schur branching identifies the rows of ambient
restriction to \(Y\) exactly with (23).

The proper map \(f:Y\to Z\) is a Chow envelope over every extension
field of \(\mathbb C\). For an invariant \(U\) with
\(j=\operatorname{rank}(N|_U)\le h-1\),
\[
NU\subset U\cap I,\qquad
\dim(U\cap I)\ge 2h-j-(c'-r')\ge h-1.
\]
One can choose an \((h-1)\)-plane \(F_0\) between those spaces
by linear algebra over the same field. Taking this lift at the
generic point of any irreducible subvariety and closing its graph
proves surjectivity of \(f_*\) on integral Chow groups.

Both \(W\) and \(Z\) are finite closed unions of Iwahori cells
in the same affine Grassmannian lattice model used in
[TYPE_B_UNIT_SMITH.md, Sections 3--4](TYPE_B_UNIT_SMITH.md#3-the-invariant-subspace-variety-and-its-affine-schubert-cells).
For \(Z\), the added condition is closed and invariant: every
Iwahori automorphism commutes with multiplication by \(t=N\),
so it preserves \(\operatorname{rank}(N|_U)\).
Their integral homology is consequently free and even, with
algebraic Schubert-cycle bases. The envelope therefore makes
\(f^*\) a primitive injection on cohomology.

Ambient integral restriction
\(H^*(\operatorname{Gr}(2h,n');\mathbb Z)\to H^*(Z;\mathbb Z)\)
is surjective. This is the same application of
[Anderson's integral affine-Grassmannian Chern-generation theorem](https://doi.org/10.4153/S0008439524000122)
as in the unit-layer proof: closed even-cell unions admit surjective
restriction from affine-Grassmannian cohomology, and the Sato index
bundle restricts to \([V]-[U]\). Its Chern classes come from the
ordinary ambient Grassmannian. Thus the restriction to \(Y\) has
primitive image of rank equal to the Betti rank of \(Z\). Its
coefficient matrix has source columns indexed by \(\lambda\) and
target rows indexed by \((\alpha,\beta)\); those rows are (23).
Every nonzero Smith factor of this matrix is 1. Removing the row
signs in the integer SXP matrix gives this same matrix, so its rank
is unchanged modulo 2. This proves the asserted rank formula for
\(e_1F\).

The complement \(W\setminus Z\) has maximal rank \(h\).
Choose \(NU\in\operatorname{Gr}(h,r')\); then
\(U\cap K=NU\). The remaining choices form an affine bundle of
rank \(h(c'-h)=hv\). Its compact-support polynomial is
\[
q^{hv}{r'\brack h}_q.
\tag{25}
\]
The closed rank filtration has free even successive pieces, so
\(P_W(q)-P_Z(q)\) is exactly (25). This difference is the
kernel polynomial in (21). For \(a=1\), the same statement is
immediate: \(W\) is a point, \(Z\) is empty, and (25) equals 1.

The shift by \(e_k\) in (12) multiplies this upper bound by
\(q^h\). It therefore gives
\(I_{a,2h}(q)\le q^{h(v+1)}{h+u-1\brack h}_q\).
The proved independent classes (14)--(15) give the reverse
coefficientwise inequality. Equality proves (18), and those classes
form the full intersection basis.

### An independent elementary proof for the two-dimensional tail

For \(h=1\), equality in (15) has the following independent proof:
\[
\boxed{
I_{a,2}(q)=q^{v+1}(1+q+\cdots+q^{u-1}),\qquad
u=\lceil a/2\rceil,\quad v=\lfloor a/2\rfloor.
}\tag{26}
\]
Consequently (13), with (26), is a closed type-D unit formula for
every \(a\) when \(k=2\).

Here is the upper bound, in addition to the independent classes already
constructed. In the two-quotient coefficient coordinates, the
Frobenius image is the image of symmetric polynomials in
\[
T=\mathbb F_2[A,B]/(A^{u+1},B^{v+1}).
\]
The Euler image is \(A\) times the image of symmetric polynomials in
\[
T'=\mathbb F_2[A,B]/(A^u,B^{v+1}).
\]
Indeed \(c_{(i),(j)}^\lambda\) is the coefficient of \(A^iB^j\)
in \(s_\lambda(A,B)\); subtracting the full two-column exchanges
the quotient components and lowers the even component by one. The
source Grassmann width bounds introduce no further restriction:
the maximum surviving degrees in \(T,T'\) are \(a,a-1\),
respectively. This identifies the two image spaces exactly.

The image of symmetric polynomials in either truncated rectangle is
easy to describe. Coefficients agree on any monomial and its swap
when both survive. If only one survives, its coefficient is free.
This follows directly from the basis of monomial symmetric functions.

Let an element of the intersection have degree \(m\) in these
coordinates, and write it as \(A g\), where \(g\in T'\) has
degree \(m-1\). Denote its coefficients by
\[
x_i=[A^iB^{m-1-i}]g,
\quad
\max(0,m-1-v)\le i\le\min(u-1,m-1).
\]
Symmetry of \(g\) equates \(x_i\) with
\(x_{m-1-i}\) when both monomials survive. Symmetry of \(Ag\)
equates \(x_i\) with \(x_{m-2-i}\) whenever both corresponding
monomials lie in the common square \(0\le\deg_A,\deg_B\le v\).
Together these reflections, whose centers differ by one, equate all
consecutive surviving coefficients. More explicitly, put
\(L=m-1-v\) and \(H=u-1\) when \(m>v\). If \(u=v\),
the first reflection is available for \(i>L\) and the second on
the whole interval \([L,H]\); their composition sends every
\(i>L\) to \(i-1\). If \(u=v+1\), the first is available on
the whole interval and the second on \([L,H-1]\); the same
composition again sends \(i>L\) to \(i-1\).

For \(1\le m\le v\), the analogous interval is
\([0,m-1]\), and the same reflections equate all its entries.
The coefficient of \(A^m\) in \(Ag\) must then be zero because
the coefficient of \(B^m\) is zero. This forces all entries to
vanish. For \(v+1\le m\le a\), their equality leaves at most
one free coefficient. There are no surviving monomials for \(m>a\).
The case \(m=0\) is zero because the Euler image has positive
degree. The classes (14) with \(\mu=(m)\) supply a nonzero vector
for every \(v+1\le m\le a\). This proves (26).


## 6. What remains open and the literal integral form

The calculation above determines the number of unit factors. It does
not establish the full square-sector D/B Smith comparison or identify
all gluing flags at positive 2-adic levels. Higher Smith exponents
still require additional integral comparison and gluing information
from [ORTHOGONAL_SMITH.md](ORTHOGONAL_SMITH.md).

For the literal defining/vector polynomial-character lattice, every
positive-degree normalized factor is multiplied by \(2^d\) in degree
\(d\). The unit statements here refer to the normalized lattice;
see [TYPE_D_SMITH.md](TYPE_D_SMITH.md) for that normalization and the
already proved exact Euler-sector comparison.
