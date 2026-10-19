# Type D square Smith factors when one block has size two

**Status: internally audited repository deduction, 19 October 2026.**
This note proves equality with the normalized type-B square image in every
degree when \(\min(a,k)=2\), and gives a general upper bound on every
square-sector exponent. Its two-variable argument has been independently
reviewed within this session. These are unpublished deductions from the
published Schubert inputs below, without an originality claim or external
expert certification. The general B/D square comparison for larger blocks
is not proved here.

## 1. Statement and normalization

Let
\[
X_D=\operatorname{OG}(a,2n),\qquad
X_B=\operatorname{OG}(a,2n+1),\qquad n=a+k,
\]
and use the same normalized square inputs
\[
q_i=e_i(x_{a+1}^2,\ldots,x_n^2),\qquad 1\le i\le k.
\]
In D, \(q_k=p^2\), where \(p=x_{a+1}\cdots x_n\).
Only polynomials in the \(q_i\) are included; the separate Euler sector
\(p\mathbb Z[q]\) is not part of this theorem. Polynomial degree \(d\)
means cohomological degree \(2d\). Smith lists retain unit factors and
refer to saturation in the selected rational image, not necessarily all
cohomology in that degree.

**Theorem.** If \(\min(a,k)=2\), the normalized D square-sector and B
Smith lists agree in every degree. Put \(A=\max(a,k)\),
\(u=\lceil A/2\rceil\), and \(v=\lfloor A/2\rfloor\). In degree
\(d=2m\), their common list is
\[
\boxed{\{2^{\epsilon(r,s)}:0\le s\le r\le A,\quad r+s=m\},}
\tag{1}
\]
where
\[
\epsilon(r,s)=
\begin{cases}
0,&r\le u\text{ and }s\le v,\\
2,&s>u,\\
1,&\text{otherwise}.
\end{cases}
\tag{2}
\]
Odd degrees have zero square image. The partition labels count factors;
they do not identify each original Schur input as an individual Smith
vector. The formula follows from the comparison proved below and the
[independent type-B two-row theorem](TYPE_B_TWO_ROW_SMITH.md).

For a multiplicity form, put
\[
G_A(z)={A+2\brack2}_z,\qquad
E_0(z)={u+2\brack2}_z-
\begin{cases}0&A\text{ even},\\z^{2u}&A\text{ odd},\end{cases}
\]
\[
E_2(z)=z^{2u+2}{v+1\brack2}_z,
\qquad E_1(z)=G_A(z)-E_0(z)-E_2(z).
\tag{3}
\]
The coefficient of \(z^m\) in \(E_j\) counts the factors \(2^j\)
in degree \(2m\). A Gaussian binomial with upper entry smaller than two
is zero.

## 2. The published inputs

The following inputs, rather than the Smith conclusion, are published:

- The integral type-C theta Giambelli basis and B/C row rescaling are in
  [Buch--Kresch--Tamvakis, A Giambelli formula for isotropic Grassmannians](https://math.umd.edu/~harryt/papers/isogiamfinal.pdf),
  Definition 1 and Theorem 1, and Section 2.1. If
  \(\nu(\lambda)=\#\{i:\lambda_i>k\}\), the Chern-class-preserving
  comparison sends the C class \(\sigma_\lambda\) to
  \(2^{\nu(\lambda)}\tau_\lambda\) in B.
- The even Giambelli formula, primitive invariant basis
  \(Y_\lambda=\tau_\lambda\) or \(\tau_\lambda+\tau'_\lambda\), and
  finite specialization are in
  [Buch--Kresch--Tamvakis, A Giambelli formula for even orthogonal Grassmannians](https://math.umd.edu/~harryt/papers/oggiam.pdf),
  Theorems 1 and 4 and equation (16).
- The quaternionic Grassmannian Schur basis and the compact fibration
  used below are in
  [Pragacz--Ratajski, Formulas for Lagrangian and orthogonal degeneracy loci](https://www.impan.pl/~pragacz/download/comp3868.pdf),
  Appendix A, Sections 10.3 and 10.5. The integral comparison is also
  explained in [TYPE_B_SMITH.md, Section 1](TYPE_B_SMITH.md#1-a-primitive-symplectic-comparison-lattice).

The integral block dualities and common square-input kernel are already
proved in [the square comparison note](TYPE_D_SQUARE_COMPARISON.md).
They allow us to put the two-element block in the head. We therefore
prove the theorem for \(a=2\), \(k=K\), and \(M=2K+2\).

## 3. A two-variable integral Chern model

Let \(x,y\) be the two Chern roots of the dual tautological bundle, and
set
\[
u_1=x^2,\quad u_2=y^2,\quad h=x+y,\quad e=xy,
\qquad s=u_1+u_2,\qquad t=u_1u_2.
\]
Write \(H_j=h_j(u_1,u_2)\), with \(H_0=1\) and \(H_j=0\) for
\(j<0\), and use the ordinary Grassmann ring
\[
R=\mathbb Z[u_1,u_2]^{S_2}/(H_{K+1},H_{K+2}).
\tag{4}
\]
Its integral Schur basis consists of \(s_{(r,s)}(u_1,u_2)\) with
\(0\le s\le r\le K\). In particular, \(H_j=0\) for every \(j>K\).

The integral type-C cohomology ring for \(\operatorname{IG}(2,2K+4)\)
has the presentation
\[
C=R[e,h]/(e^2-t,\ h^2-s-2e).
\tag{5}
\]
To see its integrality explicitly, use the compact fibration
\(\operatorname{LG}(2,4)\to\operatorname{IG}(2,2K+4)
\to\operatorname{HGr}(2,K+2)\). The base has the integral Schur ring
\(R\), with its characteristic classes pulled back to the head squares.
The fiber \(\operatorname{LG}(2,4)=Q^3\) has integral basis
\(1,h,e,he\), where \(h^2=2e\) and \(e^2=0\). The global classes
\(h=c_1(S^*)\) and \(e=c_2(S^*)\) restrict to these fiber classes.
Integral Leray--Hirsch therefore supplies the same \(R\)-basis globally.
The identities \(h^2=s+2e\), \(e^2=t\) prove (5), since their monic
quotient is already free on exactly this basis. Equivalently, eliminate
\(s=h^2-2e\) and \(t=e^2\) before imposing (4).

Here \(\deg h=1\), \(\deg e=2\), and a root \(u_i\) has polynomial
degree two. Hence the even-degree part of (5) is precisely
\[
C^{\mathrm{even}}=R\oplus eR.
\tag{6}
\]
The squared head and tail elementary series are inverse over
\(\mathbb Z\). Thus the selected square-input image in C is exactly
the summand \(R\). Head Schur classes and the original tail Schur
classes differ by conjugation and the sign \((-1)^{|\lambda|}\), an
integral change of source basis.

For the universal quotient Chern classes \(c_j=h_j(x,y)\), the
generating series gives
\[
c_{2j}=H_j+eH_{j-1},\qquad c_{2j+1}=hH_j.
\tag{7}
\]
All \(c_i\) with \(i>M\) vanish, since the quotient bundle has rank
\(M\).

There is a Chern-preserving integral homomorphism from (5) to the D
Chern subring: the two relations \(H_{K+1}=H_{K+2}=0\) hold there
by the same squared elementary-series identity. In D we also have
\(e p=0\). This is the vanishing of the global Euler invariant on the
complete flag; rational vanishing is integral because cohomology is
torsion-free and partial-to-full pullback is injective. Since
\(p^2=q_K=(-1)^KH_K\), we obtain \(eH_K=0\) in D. In particular,
\[
c_M=H_{K+1}+eH_K=0\qquad\text{in D}.
\tag{8}
\]
Equation (8) uses the Euler relation; the rank of the quotient bundle
alone would not imply it.

## 4. The even-degree Schubert rows that survive

For a two-row partition \((r,s)\), the type-C theta polynomial is
\[
\Theta_{(r,s)}=
\begin{cases}
c_rc_s-c_{r+1}c_{s-1},&r+s\le2K+1,\\
c_rc_s+2\displaystyle\sum_{j=1}^s(-1)^jc_{r+j}c_{s-j},
&r+s\ge2K+2.
\end{cases}
\tag{9}
\]
The even-orthogonal unnormalized Giambelli polynomial has the same
two branches, with its second branch beginning at \(2K+1\).
Consequently the two polynomials are **equal whenever \(r+s\) is
even**. The single-row formula is the same in both types.

The finite C partitions have length at most two and first part at most
\(M\); the D partitions have first part at most \(M-1\). For every
retained partition \(\lambda\) of even size, (9) therefore maps to
\[
\Theta_\lambda\longmapsto2^{\nu(\lambda)}Y_\lambda
\qquad\text{in D}.
\tag{10}
\]
The missing C partitions have the form \((M,2s)\), \(0\le s\le K\),
where \(s=0\) means the single row \((M)\). Every raising correction
involving the first row of (9) has a Chern index greater than \(M\),
so
\[
\Theta_{(M,2s)}=c_Mc_{2s}=e\,s_{(K,s)}(u_1,u_2)
\qquad\text{in C}.
\tag{11}
\]
Indeed \(c_M=eH_K\),
\(H_KH_s=s_{(K,s)}\) in the rectangle ring, and
\(tH_KH_{s-1}=0\): multiplying \(s_{(K,s-1)}\) by the determinant
adds one box to each row and exceeds the width \(K\). Formula (8)
shows that all classes in (11) map to zero in D.

Use any fixed integral Schur source basis of \(R\), and let
\(C_\lambda\) denote its coefficient row on
\(\Theta_\lambda\). By the B/C comparison, the corresponding B row
is \(2^{\nu(\lambda)}C_\lambda\). By (10), the D row lattice is
generated by precisely these weighted rows for the retained partitions.
If \(Y_\lambda\) is a sum of two typed classes, its two coefficient
rows are identical; this introduces no extra factor of two. Thus it
remains to show that each missing weighted B row is redundant over
\(\mathbb Z\).

## 5. An integral relation eliminating every missing row

Let \(\ell_s:C^{\mathrm{even}}\to\mathbb Z\) extract the coefficient
of \(e s_{(K,s)}\) in the direct sum (6). It vanishes on every
source element \(f\in R\), while (11) gives the identity matrix on
the missing Schubert classes. Expanding \(f\) in the C Schubert
basis consequently gives the exact row relation
\[
C_{(M,2s)}=-\sum_{\lambda\ \mathrm{retained}}
\ell_s(\Theta_\lambda)C_\lambda.
\tag{12}
\]
Only rows of the same degree occur.

For a missing row,
\(\nu(M,2s)=1+\mathbf 1_{2s>K}\). If this exponent is two, multiplying
(12) by four expresses its weighted row as an integral combination of
the retained weighted rows, whose exponents are at most two.

Suppose instead that \(2s\le K\), so the missing exponent is one.
The sole remaining divisibility issue in (12) is
\[
\ell_s(\Theta_\lambda)\equiv0\pmod2
\quad\text{when }\nu(\lambda)=2.
\tag{13}
\]
Write such a retained partition as \(\lambda=(r,t)\), where
\(r>t>K\) and \(r\le2K+1\). Its total degree is even, so \(r,t\)
have the same parity. All nonleading corrections in its second branch
of (9) have a factor two. If both parts are odd, (7) gives leading
product \(h^2H_iH_j=(s+2e)H_iH_j\), whose \(e\)-coefficient is even.
If \(r=2i\), \(t=2j\), the leading \(e\)-coefficient is
\[
H_iH_{j-1}+H_{i-1}H_j.
\tag{14}
\]
In the degree of (12),
\(i+j=K+s+1\) and
\(s\le K/2<i,j\le K\). The two-variable Pieri rule therefore
gives coefficient one on \(s_{(K,s)}\) in each summand of (14):
all of \(i,j,i-1,j-1\) are at least \(s\), and each is at most
\(K\). Their sum has coefficient two. This proves (13).

Multiplying (12) by two is now an integral combination of retained
weighted rows as well. Every missing B row belongs to the D row
lattice; the reverse containment was already explicit in (10).
The two row lattices are equal. Transposition preserves nonzero Smith
factors, proving equality of the selected image Smith lists for
\(a=2\). Integral block duality proves (1)--(3) for either block order.

## 6. A separate bound in every rank

**Proposition.** For arbitrary \(a,k\ge1\), every nonzero normalized
D square-sector factor divides \(2^{\min(a,k)}\).

This proposition does not identify the multiplicities of those factors.
It is useful beyond the two-row theorem but is not a proof of the general
B/D equality.

Let \(C_D\) be the integral Chern subring, represented by symmetric
polynomials in the \(a\) head roots. On the integral polynomial ring,
delete every monomial having an odd exponent in a head root. This
coefficient projection commutes with head permutations and fixes exactly
the head-square symmetric polynomials. It descends to \(C_D\): over
\(\mathbb Q\), the full-flag D kernel
\[
(e_1(x^2),\ldots,e_{n-1}(x^2),x_1\cdots x_n)
\]
is stable under every individual head sign change. Its parity projections
therefore preserve the kernel. If two integral representatives agree in
cohomology, their projected difference vanishes rationally and then
integrally by torsion-freeness. Integral partial-to-full injectivity
finishes the descent. This argument uses coefficient deletion on the
polynomial representatives; it does not assert that the same projection
preserves every integral cohomology class outside the Chern subring.

The image of the descended projection is precisely the selected square
subring \(L\), by inversion of the head and tail squared elementary
series. Thus \(L\) is an additive direct summand of \(C_D\), degree by
degree. In particular,
\[
C_D\cap(L\otimes\mathbb Q)=L.
\tag{15}
\]

The integral invariant cohomology lattice has basis \(Y_\lambda\):
each is one typed Schubert class or the sum of a swapped pair. The
published even Giambelli formula gives
\(2^{\nu(\lambda)}Y_\lambda=R^\lambda c_\lambda\in C_D\), with
\(\nu(\lambda)\le a\). Hence \(2^aH^*(X_D;\mathbb Z)^+\subset C_D\).
The saturation of \(L\) lies in this invariant lattice, since it is
the kernel of \(1-\iota\) and hence primitive. For any saturated class
\(z\), (15) now gives
\[
2^az\in C_D\cap(L\otimes\mathbb Q)=L.
\]
Every factor divides \(2^a\). Apply integral square-sector block
duality to obtain the asserted bound by \(2^{\min(a,k)}\).

For \(k=1\), the Grassmannian denotes the parabolic omitting both fork
nodes, rather than a single maximal-parabolic node; the same geometric
argument and low-rank interpretation apply.

## 7. Verification and the remaining scope

The independent
[Chern/Littlewood--Richardson verifier](experiments/verify_type_d_two_row_smith.py)
uses no Borel, Weyl-group, or divided-difference engine. For
\(1\le K\le12\), in every square degree \(0\le d\le4K\), it checks
the unimodular C Schubert change of basis, the exact relations (12),
the parity condition (13), equality of the two weighted integer row
modules, and the factors in (1)--(2). Its
[report](experiments/type_d_two_row_smith_verification.json)
records the script and proof hashes. These finite checks corroborate
the proof; they are not substituted for it.

This result completes every square-sector Smith layer for the two-row
family. The corresponding general comparison for \(\min(a,k)\ge3\)
and the positive-level gluing needed for the full D image remain
separate questions; see [the candidate note](TYPE_D_CANDIDATE.md) and
[the integral gluing theorem](ORTHOGONAL_SMITH.md).

For the literal defining/vector polynomial-character lattice, multiply
every nonzero normalized degree-\(d\) factor by \(2^d\), as explained
in [the type-D normalization note](TYPE_D_SMITH.md#literal-polynomial-character-form).
No claim about a spin-representation Springer map is made here.
