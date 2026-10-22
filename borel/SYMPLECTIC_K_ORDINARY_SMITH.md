# Closed ordinary symplectic projective K-image Smith factors

22 October 2026. This completes the ordinary Smith list for the defining
symplectic projective family. The argument identifies the specialized
Hankel moments with Euler's finite binomial continued fraction, then
applies the established orthogonal-polynomial diagonalization of Hankel
matrices over the 2-local integers. It is a deduction from known results,
not a claim of a new continued fraction or Hankel diagonalization theorem.
Literature priority of this particular Springer-image application is not
established. The PDF and LaTeX exposition are unchanged.

## 1. What ordinary means and which image is being classified

Ordinary means **nonequivariant**, here still in integral complex K-theory.
The map from torus-equivariant K-theory forgets the equivariant structure
and augments every torus character t_i to 1, so s_i=t_i+t_i^(-1) becomes 2.
It does not change K-theory to ordinary singular cohomology.

Let G=Sp_(2n,C), rho its defining representation and P=P1. The variety is
CP^(2n-1). Use the literal character-integral Springer source, with the
Part I normalization

\[
 \theta_\rho^*(\epsilon_i)=(t_i-t_i^{-1})/2,
 \qquad \Lambda_\rho=2\mathbb Z\epsilon_1+\cdots+2\mathbb Z\epsilon_n.
\]

This is the [previously specified full Levi source](OTHER_THEORIES_HEIGHT.md#1-the-source-and-the-space),
not the unscaled generator ring. Put

\[
 K=K^0(\mathbb{CP}^{2n-1})=\mathbb Z[v]/v^{2n},\quad
 L=1+v,\quad b=L-L^{-1},\quad J=\mathbb Z[b]\subset K.
 \tag{1}
\]

The higher Levi generators add no independent image classes. The object
classified is the additive finite group C_n=K/J, equivalently the
integer Smith factors of the inclusion J into K. Both lattices have
rank 2n. The n unit factors contribute no summands to C_n.

## 2. Complete formula

Write n=2r+epsilon, epsilon in {0,1}, and let s_2(a) be the number of
ones in the binary expansion of the nonnegative integer a (s_2(0)=0).
For 0<=j<r define

\[
 a_{n,j}=\begin{cases}
 4+8j+s_2(r-j-1)-s_2(r+j),&n=2r,\\
 8+8j+s_2(r-j-1)-s_2(r+j+1),&n=2r+1.
 \end{cases} \tag{2}
\]

**Theorem.** The complete Smith list of J into K, in divisibility order,
is

\[
 \boxed{\underbrace{1,\ldots,1}_{n},\quad
 (2\text{ if }\epsilon=1),\quad
 2^{a_{n,0}},2^{a_{n,0}},\ldots,
 2^{a_{n,r-1}},2^{a_{n,r-1}}.} \tag{3}
\]

For n=1 the paired list is empty, giving (1,2). Thus

\[
 C_n\cong(\mathbb Z/2)^\epsilon\oplus
 \bigoplus_{j=0}^{r-1}(\mathbb Z/2^{a_{n,j}})^2. \tag{4}
\]

The a_(n,j) are strictly increasing. Consequently the multiplicity of
2^e is epsilon for e=1, two when e=a_(n,j) for some j, and zero otherwise;
the unit multiplicity is n. This is a closed formula for every
multiplicity, without minors, Smith elimination or a recursive rank table.

Equivalent factorial-valuation expressions, useful for checking (2), are

\[
 a_{2r,j}=3+6j+\nu_2\left(\frac{(r+j)!}{(r-j-1)!}\right),\qquad
 a_{2r+1,j}=6+6j+\nu_2\left(\frac{(r+j+1)!}{(r-j-1)!}\right).
 \tag{5}
\]

Examples, omitting the n unit factors:

| n | Exponents of 2 in the nonunit factors |
| --- | --- |
| 1 | 1 |
| 2 | 3, 3 |
| 3 | 1, 7, 7 |
| 4 | 4, 4, 10, 10 |
| 5 | 1, 7, 7, 15, 15 |
| 6 | 3, 3, 12, 12, 18, 18 |
| 8 | 5, 5, 11, 11, 19, 19, 25, 25 |

## 3. Reduction to an ordinary Hankel matrix

The [integral doubled-block theorem](SYMPLECTIC_K_FITTING.md) is already
proved over the full spectral coefficient ring, hence remains valid
under s_i=2. Write

\[
 (s-2)^n=E_n(s^2)+sO_n(s^2),\qquad
 (p,q)=\begin{cases}(E_n,O_n),&n\text{ even},\\
                    (O_n,E_n),&n\text{ odd}.
       \end{cases} \tag{6}
\]

Here p(u) is monic of degree r. Let H be multiplication by q on
Z[u]/p. The previous theorem says, as abelian groups,

\[
 C_n\cong\operatorname{coker}(2H)^2\oplus(\mathbb Z/2)^\epsilon.
 \tag{7}
\]

The [unimodular coefficient pairing](SYMPLECTIC_K_BALANCING.md#2-an-integral-hankel-representative-in-every-degree)
identifies H with the r-by-r Hankel matrix (mu_(i+j)), where

\[
 \mu_k=\operatorname{coeff}_{u^{r-1}}(q u^k\bmod p).
 \tag{8}
\]

Equivalently mu_k is the coefficient of u^(-k-1) in the expansion at
infinity of q(u)/p(u). This follows by monic division: for a remainder
of degree below r its leading coefficient is the coefficient of u^(-1)
after division by p; polynomial parts have no negative powers.
Thus the moment generating function is obtained directly from (6).

Define

\[
 A_n(t)=\sum_{k\ge0}\binom n{2k}t^k,\quad
 B_n(t)=\sum_{k\ge0}\binom n{2k+1}t^k,\quad
 F_n(t)=B_n(t)/(nA_n(t)). \tag{9}
\]

Binomial expansion of (6), reversing the coefficients, gives

\[
 \sum_{k\ge0}\mu_k t^k=
 \begin{cases}
 -2n F_n(4t),&n\text{ even},\\
 -2n\bigl(F_n(4t)-1\bigr)/t,&n\text{ odd}.
 \end{cases} \tag{10}
\]

For odd n the constant polynomial part of q/p is -2n; its removal
explains the shifted sequence in (10). Missing this shift would give
the wrong odd-rank factors.

## 4. Published continued fraction and its 2-local coefficients

[Euler, E750, Sections VI--VII](https://arxiv.org/pdf/math/0507459)
gives the finite binomial continued fraction

\[
 F_n(t)=\cfrac{1}{1+\cfrac{c_1t}{1+\cfrac{c_2t}{1+\cdots}}},\qquad
 c_k=\frac{n^2-k^2}{(2k-1)(2k+1)},\quad c_n=0.
 \tag{11}
\]

This is a rational-function identity for positive integral n: the
fraction terminates. Analytic convergence is irrelevant. Euler writes
the reciprocal using odd denominators 1,3,5,...; normalizing those
denominators gives (11). In analytic notation it is also
F_n(x^2)=tanh(n arctanh(x))/(nx).

Work over Z_(2), the rationals with odd denominator. Set

\[
 \alpha_k=-4c_k\quad(1\le k<n),\qquad\alpha_k=0\quad(k\ge n).
 \tag{12}
\]

Every alpha_k belongs to Z_(2), since its displayed denominator is
odd. Then F_n(4t) is the standard S-fraction with successive
numerators alpha_k*t and denominator signs minus.

The usual S-to-J contraction (see
[DLMF 3.10.11--3.10.12](https://dlmf.nist.gov/3.10#E11)) gives the monic
orthogonal-polynomial recurrence coefficients below. They can also
be obtained by pairing consecutive levels in the finite fraction:

| Moment sequence normalized to constant 1 | beta_0 | beta_j, j>=1 | lambda_j, j>=1 |
| --- | --- | --- | --- |
| F_n(4t) | alpha_1 | alpha_(2j)+alpha_(2j+1) | alpha_(2j-1) alpha_(2j) |
| (F_n(4t)-1)/(alpha_1*t) | alpha_1+alpha_2 | alpha_(2j+1)+alpha_(2j+2) | alpha_(2j) alpha_(2j+1) |

For clarity, the second row follows by writing the first S-tail as
T_2. Since F=1/(1-alpha_1*t*T_2),
(F-1)/(alpha_1*t)=1/(1-alpha_1*t-alpha_2*t*T_3).
Pairing the remaining levels gives the indicated coefficients.
For odd n>=3 alpha_1 is nonzero; n=1 is handled separately by (7).
All beta and lambda coefficients are 2-local integers.

## 5. An actual unimodular diagonalization, not just determinants

For either row above define monic polynomials

\[
 P_{-1}=0,\ P_0=1,\qquad
 P_{j+1}(u)=(u-\beta_j)P_j(u)-\lambda_jP_{j-1}(u).
 \tag{13}
\]

Their coefficient matrix U, using P_0,...,P_(r-1), is lower
unitriangular over Z_(2). It and its inverse are therefore integral
over Z_(2). The classical moment identity, stated as
[Miller--Stanton, Theorem 1 and its proof](https://arxiv.org/pdf/1704.03539),
says that the normalized moment Hankel matrix becomes
diag(1,lambda_1,lambda_1 lambda_2,...) under U(-)U^t.
This identity holds over the universal coefficient ring and specializes
to the 2-local coefficients here. Scaling back by 2*mu_0 gives

\[
 U(2\mu_{i+j})U^t=\operatorname{diag}(d_0,\ldots,d_{r-1}),\qquad
 d_j=2\mu_0\prod_{k=1}^j\lambda_k. \tag{14}
\]

The two cases of (10) give

\[
 2\mu_0=\begin{cases}-4n,&n\text{ even},\\
                    16n(n^2-1)/3,&n\text{ odd}.
       \end{cases} \tag{15}
\]

All denominators in U and (14) are odd. This is the essential
integrality check: a diagonalization over Q alone would not determine
Smith factors at two.

For n even, each n^2-(2k-1)^2 is odd, so

\[
 \nu_2(d_j)=2+\nu_2(n)+4j+
                    \sum_{k=1}^j\nu_2(n^2-4k^2). \tag{16}
\]

For n odd, each n^2-4k^2 is odd, so

\[
 \nu_2(d_j)=4+4j+
                    \sum_{k=0}^j\nu_2(n^2-(2k+1)^2). \tag{17}
\]

In even rank, the successive increments in (16) are at least six;
in odd rank those in (17) are at least seven. Thus these diagonal
entries already occur in increasing 2-adic divisibility order.
Writing n=2r or 2r+1 in the products proves (5). Legendre's classical
identity nu_2(a!)=a-s_2(a) proves (2).

Finally the original inclusion (1), in the v-power basis, has columns
1,b,...,b^(2n-1). Since b=2v+O(v^2), it is triangular with diagonal
1,2,...,2^(2n-1). Its determinant is a power of two. Hence there are
no odd-primary Smith factors: the 2-local computation determines
the entire integer Smith list. Adding the n unit factors to (7)
proves (3). QED.

The last exponent reduces to 4(n-1)-s_2(n-1) for n>=2, recovering
the [previous largest-factor theorem](OTHER_THEORIES_HEIGHT.md#5-an-arbitrary-rank-integral-k-theory-consequence).
The total exponent is n(2n-1), also recovered from the triangular
inclusion. Those earlier invariants are now consequences of the full list.

## 6. Independent checks, attribution and remaining scope

The [verification script](tools/verify_symplectic_ordinary.py) compares
the closed formula with exact integer Smith elimination on the binomial
half block, and independently with powers of L-L^(-1) in the original
truncated projective ring. Projective Euler-characteristic coordinates
provide a third integral-basis check. It also verifies the rational
continued fraction exactly and reconstructs the explicit 2-local
orthogonal coefficient matrix and diagonal identity (14).
The [report](symplectic-ordinary-verification.json) states the finite
ranges and source hashes. These are corroboration; the arbitrary-rank
proof is Sections 3--5 with the cited published inputs.

**Known:** the projective K-ring, the literal Springer coordinates from
Part I, Euler's finite binomial continued fraction, S-to-J contraction,
orthogonal Hankel diagonalization, Smith theory and Legendre's identity.
**Deduction here:** their identification with this specific ordinary
Springer-image Hankel matrix and the complete list (2)--(4).
**Still separate:** arithmetic primary structure of the equivariant
image after moment substitution, global equivariant cyclic factors in
higher ranks, larger isotropic flags and other representations.
This theorem does not settle any of those questions, or every ordinary
symplectic flag variety. A targeted literature search located the
classical mechanism, but did not establish priority of the application.
