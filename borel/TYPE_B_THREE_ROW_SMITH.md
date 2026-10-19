# A closed formula for all type-B Smith factors when the smaller block is three

This note determines every normalized Smith multiplicity when
`min(a,k)=3`. The formulas count bounded partitions and do not require
divided differences, a coefficient matrix, or Smith reduction. They follow
from integral geometry of three filtered coefficient lattices, together
with Schur multiplication. No novelty claim is made.

The general normalization and block duality are proved in
[TYPE_B_SMITH.md](TYPE_B_SMITH.md). The complete answers for smaller block
sizes one and two are there and in
[TYPE_B_TWO_ROW_SMITH.md](TYPE_B_TWO_ROW_SMITH.md), respectively.

## Formula

Put the two block sizes in the form `(A,3)`, where `A>=3`. Write `d=2m`
for polynomial degree, so the corresponding cohomological degree is `4m`.
All odd polynomial-degree images vanish. Define

$$
U=\left\lceil\frac{A+1}{2}\right\rceil,\qquad
V=\left\lfloor\frac{A+1}{2}\right\rfloor,\qquad
T=\left\lfloor\frac{A+2}{2}\right\rfloor.
$$

Let the Gaussian coefficient count rectangle partitions, and set

$$
\begin{aligned}
G_A(q)&={A+3\brack3}_q,\\
F_A(q)&=\sum_{t=0}^{V}\sum_{s=t}^{U}\sum_{r=s}^{A}q^{r+s+t},\\
H_A(q)&=q^{3(T+1)}{A-T+2\brack3}_q.
\end{aligned}\tag{1}
$$

The polynomial `F_A` therefore counts partitions with three padded parts
`r>=s>=t>=0`, subject to `r<=A`, `s<=U`, and `t<=V`. Define one more
polynomial by the following closed formulas:

$$
\begin{aligned}
J_{2b}(q)
&={b+3\brack3}_q
 +q^{b+1}\left({b+2\brack2}_q-q^{2b}\right)
 -q^{2b+1}(1+q+\cdots+q^{b-1}),\\
J_{2b+1}(q)
&={b+3\brack3}_q+q^{b+1}{b+2\brack2}_q.
\end{aligned}\tag{2}
$$

If `M_e(m)` denotes the number of Smith factors `2^e` in degree `2m`,
the complete answer is

$$
\boxed{
\begin{aligned}
\sum_m M_0(m)q^m&=J_A(q),\\
\sum_m M_1(m)q^m&=F_A(q)-J_A(q),\\
\sum_m M_2(m)q^m&=G_A(q)-H_A(q)-F_A(q),\\
\sum_m M_3(m)q^m&=H_A(q).
\end{aligned}}
\tag{3}
$$

In particular all four displayed differences have nonnegative integer
coefficients. Their sum is `G_A`, the rank generating polynomial. The
formulas determine multiplicities; they do not assign an exponent to each
original Schur generator individually.

For example, with `A=3`, the Smith exponent lists as `m` runs from 0 to 9
are

$$
(0),\ (0),\ (0,0),\ (0,0,1),\ (0,1,1),\ (1,1,1),\
(1,1,2),\ (1,2),\ (2),\ (3).
$$

For the defining/vector Springer integral form, every degree-`d` image is
multiplied by `2^d`. Thus the four factors in degree `d=2m` become
`2^(2m)`, `2^(2m+1)`, `2^(2m+2)`, and `2^(2m+3)`, with exactly the
multiplicities in (3). This separate normalization is proved in
[INTEGRAL_LITERATURE.md](INTEGRAL_LITERATURE.md).

## 1. Filtered symplectic row lattices determine the answer

By integral block duality, use a head block of size three and a tail block
of size `A`. Work in `X_C=IG(3,2n)`, where `n=A+3`. In degree `2m`, use
the head squared-root Schur source basis

$$
f_\lambda=s_\lambda(y_1^2,y_2^2,y_3^2),
\qquad \lambda\subseteq A^3,\quad |\lambda|=m.
$$

The total squared-root relation changes the original tail Schur basis to
this basis by conjugation and signs, an integral unimodular change. Let
`C` be its integral symplectic Schubert coefficient matrix and let `r` be
its number of columns. Its full row lattice is `Z^r`, by the primitive
quaternionic Grassmann comparison. For a minimal right signed-permutation
representative `w`, let `nu(w)` be its number of negative entries. Define

$$
P_t=\operatorname{span}_{\mathbb Z}
\{C_w:\nu(w)\le t\}\subseteq\mathbb Z^r,
\qquad 0\le t\le3,
$$

and put `R_t=rank(P_t)` and `R_(-1)=0`.

We will prove that all the `P_t` are primitive. For a chain of primitive
sublattices of a free abelian group, there is a common adapted integral
basis: successively choose bases of the free quotients. In this basis
`P_t` consists of the first `R_t` coordinates.

The type-B matrix has row `2^nu(w) C_w`. Its row lattice is

$$
\sum_{t=0}^3 2^t P_t.
\tag{4}
$$

Indeed, replacing exact-level rows by cumulative prefixes does not change
the sum, since any extra contribution from an earlier level is an integer
multiple of that earlier weighted row. In the adapted basis, (4) is a
diagonal lattice with exactly `R_t-R_(t-1)` diagonal entries equal to
`2^t`. The nonzero Smith factors of a matrix and its transpose agree.
Thus it remains to compute the four prefix ranks and prove their
primitivity.

Two of these calculations are already available:

* `P_0` is primitive by the integral graph-incidence argument in
  [TYPE_B_THREE_ROW_UNITS.md](TYPE_B_THREE_ROW_UNITS.md), and its rank
  generating polynomial is `J_A` in (2).
* `P_2` is primitive by
  [TYPE_B_HIGHEST_SMITH.md, equation (3)](TYPE_B_HIGHEST_SMITH.md#1-a-stronger-integral-coefficient-statement).
  In the present head Schur basis it is exactly the coordinate lattice
  with `lambda_3<=T`. Its rank generating polynomial is `G_A-H_A`.

Also `P_3=Z^r`, whose rank generating polynomial is `G_A`. The new step is
the middle prefix `P_1`.

## 2. The incidence geometry for the middle prefix

Fix a complex Lagrangian `L_0` in the ambient symplectic space. As in the
Schubert convention checked in
[TYPE_B_HIGHEST_SMITH.md, Section 2](TYPE_B_HIGHEST_SMITH.md#2-the-schubert-union-selected-by-nu--ell-1),
the cell indexed by `w` has

$$
\dim(\Sigma\cap L_0)=3-\nu(w).
$$

Consequently the Schubert union whose restriction records `P_1` is

$$
Z_1=\{\Sigma\in X_C:\dim(\Sigma\cap L_0)\ge2\}.
$$

Resolve this incidence condition by remembering a two-plane:

$$
Y=\{(S,\Sigma):S\subseteq\Sigma\cap L_0,\ \dim S=2\}
\longrightarrow Z_1.
$$

Over `Gr(2,L_0)=Gr(2,n)`, this is the bundle of isotropic lines in
`S^perp/S`, with fiber `IG(1,2n-4)=CP^(2A+1)`.

The pullback from `Z_1` to `Y` is a primitive injection. On every Schubert
cell, its canonical echelon coordinates give a regular two-plane section:
take the span of the first two positive pivot vectors in
`Sigma cap L_0`. The closure of this section maps with degree one onto
the corresponding Schubert variety. Hence homology pushforward surjects
onto the integral Schubert homology basis of `Z_1`. These homology groups
are free, so universal coefficients make cohomology pullback primitive.
This is the same degree-one-cycle argument used for the line incidence
resolution in the highest-factor proof; it does not infer primitivity
merely from birationality.

Choose a compatible quaternionic structure. The quaternionic span `H(S)`
has rank two, and its orthogonal complement `E` has quaternionic rank
`n-2=A+1`. Let `Q=HGr(1,E)` be the relative quaternionic projective bundle
over `Gr(2,n)`. Taking the quaternionic span of the residual isotropic line
gives a map `Y->Q` with fiber `LG(1,2)=CP^1`. Base and fiber have free
integral even cohomology. The integral Serre spectral sequence therefore
shows that its cohomology pullback is primitive as well.

Primitive target inclusions preserve coefficient row lattices, since an
integral left inverse gives the reverse row containment. The source row
lattice on `Q` is consequently exactly `P_1`.

## 3. The middle prefix is a complete-symmetric monomial lattice

Let `U_Q` be the quaternionic tautological line on `Q`, with first Borel
class `u`. The integral quaternionic projective bundle theorem gives a
Leray--Hirsch basis

$$
\sigma_\theta(S^*)\,u^r,
\qquad \theta\subseteq(A+1)^2,\quad 0\le r\le A.
\tag{5}
$$

Here `sigma_theta` is the ordinary Schur class on `Gr(2,n)`. The integral
quaternionic Schur bundle basis follows, for example, from
[Panin--Walter, Theorem 11.2, equations (11.7)--(11.8)](https://arxiv.org/pdf/1011.0649).
If `z_1,z_2` are the Chern roots of `S^*`, the source quaternionic bundle
splits on `Q` into `H(S)` and `U_Q`. Its Schur addition formula is

$$
f_\lambda\longmapsto
\sum_{\rho,r}c^{\lambda}_{\rho,(r)}
s_\rho(z_1^2,z_2^2)\,u^r.
\tag{6}
$$

Every nonzero fiber term already has `r<=A`, because `lambda_1<=A`.
There is no reduction of a higher power of `u` and hence no hidden
base-dependent projective bundle term.

The SXP rule says that the ordinary Schur coefficient of
`s_rho(z_1^2,z_2^2)` at `theta` vanishes unless `theta` has empty 2-core;
otherwise it is

$$
\epsilon_\theta c^\rho_{\alpha,\beta},\qquad
\epsilon_\theta\in\{1,-1\},
\tag{7}
$$

where `(alpha,beta)` is its 2-quotient. See
[Wildon, *A generalized SXP rule proved by bijections and involutions*, equation (1)](https://arxiv.org/abs/1508.07030).
The sign in (7) depends on the target row, so it does not affect the
integer row lattice.

For completeness, pad `theta` to two parts. Its two beta numbers are
distinct integers in `[0,A+2]`. Empty 2-core means that one is even and
one is odd, as for the empty partition. Dividing the runners by two
shows that the possible quotients are exactly

$$
\alpha=(s),\quad 0\le s\le U,
\qquad
\beta=(t),\quad 0\le t\le V.
\tag{8}
$$

Interchanging runners merely interchanges the two factors. Combining
(6)--(8) and associativity of Littlewood--Richardson multiplication, the
row `(theta,r)` is the Schur coefficient vector of

$$
\pm h_r h_s h_t,\qquad
0\le r\le A,\quad 0\le s\le U,\quad 0\le t\le V,
\tag{9}
$$

restricted to the source Schur rectangle `A^3`. All such triples occur.

Products of three complete symmetric functions have no Schur terms with
more than three rows. Passing from the full symmetric function group to
the rectangle in (9) therefore discards exactly the terms of first part
greater than `A`. This is precisely the usual cohomological quotient for
`Gr(3,A+3)`.

Inside that quotient, the `h_gamma=h_(gamma_1)h_(gamma_2)h_(gamma_3)`
indexed by `gamma subseteq A^3` are an integral basis. Indeed, their
Schur expansions have Kostka coefficients, with diagonal one and with
all other partitions dominating `gamma`. Discarding the terms outside
the rectangle leaves a unitriangular integer matrix on its Schur basis.

Sort the three indices in (9) into `r>=s>=t`. Their allowed multisets
are exactly

$$
r\le A,\qquad s\le U,\qquad t\le V.
\tag{10}
$$

The necessity follows since at most one original index can exceed `U`
and at most two can exceed `V`. Sufficiency follows by assigning the
sorted indices to the caps `A,U,V`. Thus `P_1` is the span of a subset
of an integral basis, indexed by (10). It is primitive, and its rank
generating polynomial is `F_A` in (1).

## 4. The four rank increments

We have now proved that the entire filtered flag is primitive, with rank
generating polynomials

$$
R_0:J_A,\qquad R_1:F_A,\qquad
R_2:G_A-H_A,\qquad R_3:G_A.
$$

Taking the consecutive rank differences in (4) gives (3). This proves
all multiplicities uniformly for every `A>=3` and every degree, rather
than only the examples used to suggest the formulas.

## 5. Equivalence with the compact Gaussian formulas

The formulas above also agree symbolically, in every width, with the
compact expression suggested by the general defect-stratum calculation.
Here is an algebraic verification, separate from finite checks.

Write `n=A+3`, `c=ceil(n/2)`, and `r=floor(n/2)`. For this section only,
write `[v]=(1-q^v)/(1-q)` for a q-integer. Put

$$
\mathcal U(n,l)=
\sum_j q^{j(c-l+j)}{r\brack j}_q{c-j\brack l-2j}_q,
\tag{11}
$$

where `max(0,l-c)<=j<=min(r,floor(l/2))`. In particular,

$$
\begin{aligned}
\mathcal U(n,3)&={c\brack3}_q+q^{c-2}[r][c-1],\\
W:=\mathcal U(n-2,2)&={c-1\brack2}_q+q^{c-2}[r-1],\\
\mathcal U(n-4,1)&=[c-2],\qquad
\mathcal U(n-6,0)=1.
\end{aligned}\tag{12}
$$

We claim that (3) is precisely

$$
\boxed{\sum_m M_e(m)q^m
=q^{ec}{r\brack e}_q\mathcal U(n-2e,3-e),
\qquad 0\le e\le3.}\tag{13}
$$

First, Gaussian Pascal immediately gives
`J_A=mathcal U(n,3)`. In the odd-width case `A=2b+1`, use `c=r=b+2`
and

$$
J_A={b+2\brack3}_q+q^b(1+q){b+2\brack2}_q
={b+2\brack3}_q+q^b[b+2][b+1].
$$

For `A=2b`, use `c=b+2`, `r=b+1`; the same Pascal identity simplifies
(2) to

$$
J_A={b+2\brack3}_q+q^b[b+1]^2.
$$

Next we verify the only additional middle-prefix identity needed:

$$
F_A-J_A=q^c[r]W.\tag{14}
$$

If `A=2b+1`, sorting the capped triples in (1) according to whether their
largest part exceeds `b+1` gives

$$
F_A={b+4\brack3}_q+q^{b+2}[b]{b+3\brack2}_q,
\qquad W={b+2\brack2}_q.
$$

Gaussian Pascal and the elementary q-integer identity

$$
[b+2][b+1]-[b][b+3]=q^b[2]
$$

then give

$$
F_A-J_A=q^{b+2}[b+2]{b+2\brack2}_q,
$$

which is (14). If `A=2b`, the analogous capped-triple decomposition is

$$
\begin{aligned}
F_A&={b+4\brack3}_q-q^{3b+3}
 +q^{b+2}[b-1]\left({b+3\brack2}_q-q^{2b+2}\right),\\
W&={b+2\brack2}_q-q^{2b}.
\end{aligned}
$$

The excluded triples and pairs are respectively `(b+1,b+1,b+1)` and
`(b+1,b+1)`. Substitution of (2), followed by Gaussian Pascal and

$$
[b+1]^2-[b-1][b+3]=q^{b-1}[2]^2,
$$

gives

$$
F_A-J_A=q^{b+2}[b+1]
\left({b+2\brack2}_q-q^{2b}\right).
$$

This proves (14) in the other parity as well. Both q-integer identities
follow directly on substituting `[v]=(1-q^v)/(1-q)`.

The highest factor in (3) already has the compact form
`H_A=q^(3c)*Gaussian(r,3)`, since `T=c-1`. Finally, the q-Vandermonde
identity gives

$$
{n\brack3}_q
={c\brack3}_q
 +q^{c-2}[r]{c\brack2}_q
 +q^{2c-2}{r\brack2}_q[c]
 +q^{3c}{r\brack3}_q.
$$

Use

$$
{c\brack2}_q=[c-1]+q^2{c-1\brack2}_q,
\qquad [r][r-1]=[2]{r\brack2}_q,
\qquad [c]=[2]+q^2[c-2].
$$

Together with (12), these rewrite the same total as

$$
{n\brack3}_q
=\mathcal U(n,3)+q^c[r]W
 +q^{2c}{r\brack2}_q[c-2]
 +q^{3c}{r\brack3}_q.
$$

Thus the remaining factor-four polynomial in (3) is exactly the `e=2`
case of (13). This completes the symbolic equivalence for all four
factors. It does not prove (13) for larger smaller-block sizes.

## Verification and scope

The independent verification in
[verify_type_b_partition_formulas.py](experiments/verify_type_b_partition_formulas.py)
checks these formulas against the recorded BGG and Cartesian symplectic
coefficient calculations. Its
[report](experiments/type_b_partition_formula_verification.json)
records the dataset hashes and exact finite scope. Those finite checks
corroborate the proof; they are not its replacement. The general case
`min(a,k)>=4` still requires intermediate-prefix arguments beyond this
note, although the largest-factor formula is proved for every rank.
