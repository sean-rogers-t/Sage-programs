# Integral Hankel balancing for symplectic Springer-image ideals

22 October 2026. The preceding [Schur-minor theorem](SYMPLECTIC_K_MINORS.md)
provides uniform generators. Here the doubled-block convolution is
simplified in every rank. The proof applies published Hankel ideal
theory, with the integral lifting and the specific Springer-image
consequences written explicitly. No priority claim is made.
The PDF and LaTeX article are not revised.

## 1. Statement over an arbitrary commutative ring

Let A be a commutative ring, p(u) monic of degree r, and H multiplication
by q(u) on V=A[u]/p in the basis 1,u,...,u^(r-1).
Put D_j=I_j(H), D_0=A, and D_j=0 outside 0<=j<=r.

**Theorem.** If 0<=a<=b<=r, a+b=c+d and
a<=c<=d<=b, then

\[
 \boxed{D_aD_b\subseteq D_cD_d.} \tag{1}
\]

In particular D_(j-1)D_(j+1) is contained in D_j^2. Hence, for
0<=ell<=2r,

\[
 \boxed{\sum_{a+b=\ell}D_aD_b
       =D_{\lfloor\ell/2\rfloor}D_{\lceil\ell/2\rceil}.} \tag{2}
\]

These are integral ideal equalities/containments, including at two and
after arbitrary coefficient specialization. They are not valid for
arbitrary unstructured matrices. The generic-matrix characteristic-two
failure is recorded in [Bruns-Vetter, Remarks 10.14(g)-(h)](https://www.home.uni-osnabrueck.de/wbruns/brunsw/detrings.pdf).
The additional multiplication/Hankel structure is essential here.

## 2. An integral Hankel representative in every degree

Let epsilon:V->A extract the coefficient of u^(r-1) in the unique
representative of degree less than r. Define

\[
 T_{ij}=\epsilon(u^{i+j}\bmod p),\qquad0\le i,j<r.
\]

Its entries are zero for i+j<r-1 and one for i+j=r-1.
Reversing its rows makes it triangular with diagonal one, so
det(T)=(-1)^(r(r-1)/2). It is unimodular over A, including in
inseparable or characteristic-two situations. This is a coefficient
pairing, not the field trace pairing.

For M=TH one has

\[
 M_{ij}=\epsilon(q u^{i+j}\bmod p)=\mu_{i+j}. \tag{3}
\]

Thus M is Hankel, with 2r-1 moments mu_0,...,mu_(2r-2), and
I_j(M)=I_j(H) for every j. It suffices to prove (1) for a generic
square Hankel matrix over Z and then substitute its moments.
The substitution need not be flat: a polynomial ideal containment
over Z remains a containment under every ring homomorphism.

## 3. Known field results and the passage to integers

Write X=(x_(i+j))_(0<=i,j<r), with independent variables
x_0,...,x_(2r-2), and let L_j be its minor ideals.
For c=1, [Le Dinh Nam, Theorem 3.25 and Corollary 3.26](https://arxiv.org/abs/1008.3843)
give the primary decomposition of products of generic Hankel minor
ideals and a diagonal-term Groebner basis consisting of natural
products of minors over a field. Remark 3.3(a) identifies the ideals
from the different Hankel shapes. These results build on Conca's
Hankel straightening theory; they are published inputs.

Over Q, the product decomposition uses the quantities
gamma_t(b,a)=max(b-t+1,0)+max(a-t+1,0).
Convexity gives gamma_t(b,a)>=gamma_t(d,c) when (b,a) is the more
unequal pair with the same sum. The decomposition therefore gives
L_a L_b contained in L_c L_d over Q. We use this field containment
only as the starting point for an integral argument.

Here is why no integer denominator is needed. The natural product
Groebner generators can be chosen to be products of minors of the
square X itself, with leading coefficients one:

1. For N=2r-1 moments, the initial generators of a j-minor ideal
   are products x_(v1)...x_(vj), with successive indices separated
   by at least two.
2. Each such product occurs as the diagonal term of a minor of the
   square matrix X: take row indices floor(v_k/2) and column indices
   ceil(v_k/2). Both lists are strictly increasing and bounded by r-1.
3. Products of these square minors therefore cover the leading
   generators of the corresponding product ideal. The published
   product theorem implies that they are a Groebner basis over Q.
   Their coefficients are integers and their leading coefficients
   are one.

Now divide any integer generator of L_a L_b by this monic Groebner
basis of L_c L_d. Multivariate division uses only subtraction of
integer monomial multiples, so every quotient and remainder has
integer coefficients. The remainder is zero because the polynomial
belongs to L_c L_d over Q. This exhibits an integral ideal membership,
proving L_a L_b contained in L_c L_d over Z.

If a=0, the needed containment follows directly by Laplace expansion
of each (c+d)-minor into c-minors times d-minors; no field theorem is
needed. Iterating the elementary balancing step proves (1). Substituting
(3) proves it for every monic multiplication matrix. In (2), every
summand is contained in the balanced product, and that product is
itself one of the summands. This proves the equality.

This argument uses rational primary decomposition plus a monic
integer Groebner lift. It does not infer integrality merely from
agreement over fields or from a rational ideal calculation.

## 4. A single balanced-product formula for all Springer Fitting ideals

Keep the defining Sp_(2n)(C) representation, P=P1, the full literal
source and R=Z[t_i^+-1], with s_i=t_i+t_i^-1. Equivalently use
Z[s_1,...,s_n] and then base change. Put

\[
 C=K_T^0(Sp_{2n}/P_1)/R[L-L^{-1}],\quad n=2r+\varepsilon,
 \qquad p(s)=\prod_i(s-s_i)=E(s^2)+sO(s^2).
\]

Let H multiply O modulo E in even rank, or E modulo O in odd rank;
the latter multiplier may first be reduced modulo its monic modulus.
For 0<=i<n set

\[
 k=n-i,\qquad \ell=k-\varepsilon,\qquad
 a=\lfloor\ell/2\rfloor,\quad b=\lceil\ell/2\rceil.
\]

The doubled-module presentation and (2) give

\[
 \boxed{\operatorname{Fitt}_i(C)=2^k I_a(H)I_b(H).} \tag{4}
\]

For i>=n the ideal is R. The previous convolution sum has disappeared:
only one pair of half-block minor ideals is required in every rank.
When ell is even it is a square; when ell is odd it is the product of
two adjacent minor ideals. Each half-block ideal has the initial-column
Schur generators proved in the companion note.

This is a general structural formula. It does not identify the
specialized ideals' primary decompositions or give closed numerical
Smith multiplicities in all ordinary ranks. Generic Hankel primary
decompositions need not survive the moment substitution as primary
decompositions; intersections and primary components do not generally
commute with arbitrary base change.

The separate [ordinary Hankel specialization](SYMPLECTIC_K_ORDINARY_SMITH.md)
now completes all ordinary multiplicities, using Euler's finite binomial
continued fraction and the known orthogonal Hankel diagonalization.
That arithmetic step is additional to the ideal-balancing theorem here.

More generally, for m copies of a monic multiplication presentation,
write ell=md+s, 0<=s<m. Repeated pair balancing gives

\[
 I_\ell(\operatorname{diag}(H,\ldots,H))
       = I_d(H)^{m-s}I_{d+1}(H)^s. \tag{5}
\]

This corollary describes repeated copies; it does not assert such a
formula for unrelated summands.

## 5. Explicit quartic blocks: ranks eight and nine

Let p=u^4+P*u^3+Q*u^2+R0*u+S and q=a*u^3+b*u^2+c*u+d.
Define A=(a,b,c,d), and define B by its six polynomial generators

\[
 \begin{aligned}
 B_{01}&=d^2-aR0d+acS,&
 B_{02}&=cd-aQd+abS,\\
 B_{03}&=bd-aPd+a^2S,&
 B_{12}&=c^2-aQc-bd+abR0,\\
 B_{13}&=bc-aPc-ad+a^2R0,&
 B_{23}&=b^2-aPb-ac+a^2Q .
 \end{aligned} \tag{6}
\]

Put v=(-S(b-aP),-aS-R0(b-aP),
d-aR0-Q(b-aP),c-aQ-P(b-aP)).
The four generators of D are

\[
 \begin{aligned}
 D_{012}&=v_0B_{12}-v_1B_{02}+v_2B_{01},\\
 D_{013}&=v_0B_{13}-v_1B_{03}+v_3B_{01},\\
 D_{023}&=v_0B_{23}-v_2B_{03}+v_3B_{02},\\
 D_{123}&=v_1B_{23}-v_2B_{13}+v_3B_{12}.
 \end{aligned} \tag{7}
\]

These are the initial-column two- and three-minors, respectively.
Thus I_1(H)=A, I_2(H)=B, I_3(H)=D and I_4(H)=(f), f=det(H).
The minor ideals of H+H, in increasing size, are now proved to be

\[
 R,\ A,\ A^2,\ AB,\ B^2,\ BD,\ D^2,\ fD,\ (f^2). \tag{8}
\]

For n=8 take (P,Q,R0,S)=(e_2,e_4,e_6,e_8) and
(a,b,c,d)=(e_1,e_3,e_5,e_7); the actual multiplier has an irrelevant
overall minus sign. For n=9 use the same P,Q,R0,S and

\[
 (a,b,c,d)=(e_1e_2-e_3,\ e_1e_4-e_5,\
            e_1e_6-e_7,\ e_1e_8-e_9).
\]

Write F=product_(i<j)(s_i+s_j), so (f)=(F). The complete lists are

| Fitting index i | n=8 | n=9 |
| --- | --- | --- |
| 0 | (256F^2) | (512F^2) |
| 1 | 128FD | 256FD |
| 2 | 64D^2 | 128D^2 |
| 3 | 32BD | 64BD |
| 4 | 16B^2 | 32B^2 |
| 5 | 8AB | 16AB |
| 6 | 4A^2 | 8A^2 |
| 7 | 2A | 4A |
| 8 | R | (2) |
| i>=9 | R | R |

All ideals on the right have explicit polynomial generators (6)-(7).
The table holds integrally, including characteristic two, without a
transversality or distinct-root assumption.

## 6. Independent certificates and remaining questions

[The verification script](tools/verify_symplectic_balancing.py) checks
the coefficient-pairing Hankel equivalence over independent integer
polynomial variables, verifies the quartic polynomial generator lists,
and reconstructs universal integer Hankel identities for the critical
containment I_1 I_3 contained in I_2^2. These identities are polynomial
identities in seven independent moments, not finite point samples.
Finite comparisons with the original 2n projective inclusion provide
a separate check of the rank-eight/nine tables.
[Report](symplectic-balancing-verification.json).

The all-rank proof rests on the stated Hankel product theorems and
the integral lift in section 3; the finite checks do not supply its
arbitrary-rank quantifier.

What is resolved: the balanced-product formula in every rank, the
quartic integral containment, and explicit rank-eight/nine ideal lists.
What remains: arithmetic primary structure after the moment substitution,
minimal or simpler generators of larger half-block ideals, and closed
ordinary Smith multiplicities. Larger isotropic flags and B/D quadric
targets remain separate work. The new API optimization follows (4);
large expansions in individual spectral variables may still be costly.
