# Closed type-B Smith multiplicities when one block has size two

This note completes the multiplicity calculation for `min(a,k)=2` in the
normalization of [TYPE_B_SMITH.md](TYPE_B_SMITH.md). It is a deduction from
the classical two-row theta Giambelli formula and the type-B/type-C
comparison. No claim of a new theorem in the literature is made.

## Statement

Let `a,k >= 1`, `n=a+k`, and use the normalized tail generators

$$
q_i=e_i(x_{a+1}^2,\ldots,x_n^2),\qquad
X_B=\operatorname{OG}(a,2n+1).
$$

Let `L_B(a,k;d)` be the integral image of degree-`d` polynomials over Z in
these generators; `q_i` has polynomial degree `2i`. Its Smith factors are
computed in the integral Schubert lattice, or equivalently in the
saturation of its rational image. For odd `d` the selected image is zero.

Suppose `min(a,k)=2`, and set

$$
A=\max(a,k),\qquad u=\lceil A/2\rceil,\qquad v=\lfloor A/2\rfloor.
$$

For `d=2m`, the nonzero Smith multiset is

$$
\boxed{\{\,2^{e(r,s)}:
             0\le s\le r\le A,\quad r+s=m\,\},}\tag{1}
$$

where

$$
\boxed{
e(r,s)=
\begin{cases}
0,&r\le u\text{ and }s\le v,\\
2,&s>u,\\
1,&\text{otherwise}.
\end{cases}}\tag{2}
$$

Thus only `1`, `2` and `4` occur. The partition label in (1) labels the
diagonal form obtained below; the theorem does not assert that each
unmodified Schur input already has a single Schubert term.

When `A=2v` is even, the statistic is particularly simple:
`e(r,s)` is the number of parts of `(r,s)` exceeding `v`. When `A=2v+1`
is odd, it is the number of parts exceeding `v+1`, with the additional
value `e(v+1,v+1)=1` at the central square.

## Multiplicities without enumerating partitions

Write `M_j(m)` for the number of factors `2^j` in degree `d=2m`, and let
`[N choose 2]_t` denote the Gaussian binomial coefficient, taken to be
zero when `N<2`. Put

$$
G_A(t)={A+2\brack2}_t,
$$

$$
E_0(t)={u+2\brack2}_t
       -\begin{cases}0,&A\text{ even},\\t^{2u},&A\text{ odd},\end{cases}
\qquad
E_2(t)=t^{2u+2}{v+1\brack2}_t,
$$

$$
E_1(t)=G_A(t)-E_0(t)-E_2(t).
\tag{3}
$$

Then `M_j(m)=[t^m]E_j(t)`. In particular,

$$
\operatorname{rank}L_B(a,k;2m)=[t^m]G_A(t),\qquad
[\operatorname{sat}(L_B):L_B]=2^{M_1(m)+2M_2(m)}.
$$

Here a rank-zero degree has no nonzero factors and index one in its zero
saturation. Formula (3) follows directly from (2): the unit partitions fill
the `u`-by-`2` rectangle, with its top square removed when `A` is odd;
subtracting `u+1` from both parts identifies the exponent-two partitions
with the `(v-1)`-by-`2` rectangle.

For example, `A=4` gives the following complete calculation:

| `m=d/2` | Smith factors |
| --- | --- |
| 0 | `(1)` |
| 1 | `(1)` |
| 2 | `(1,1)` |
| 3 | `(1,2)` |
| 4 | `(1,2,2)` |
| 5 | `(2,2)` |
| 6 | `(2,4)` |
| 7 | `(4)` |
| 8 | `(4)` |

## Reduction to a two-variable source

The integral block duality proved in
[TYPE_B_SMITH.md, Section 4](TYPE_B_SMITH.md#4-block-duality-including-integral-lattices)
allows the two-element block to be the head. We may therefore take
`a=2`, `k=A`, and use the two Chern roots `y_1,y_2` of the dual
tautological bundle. Put

$$
f_s=s_{(m-s,s)}(y_1^2,y_2^2),\qquad
\max(0,m-A)\le s\le\lfloor m/2\rfloor.
\tag{4}
$$

These give an integral source basis for the selected image. Indeed, the
total squared-root relation makes the head and tail elementary generating
series inverse. The involution interchanging elementary and complete
symmetric functions consequently sends a tail Schur class indexed by
`lambda` to the head Schur class indexed by `lambda'`, up to the unit
`(-1)^|lambda|`. Thus (4) is the same image lattice as the original tail
generators, with an integral change of source basis.

## The Schubert identities used in the proof

Let `c_j=h_j(y_1,y_2)`, with `c_0=1` and `c_j=0` for `j<0`.
These are the Chern classes of the universal quotient bundle: its total
Chern class is the inverse of that of the tautological bundle.
For a partition `mu`, write `Theta_mu(c)` for its type-C theta Giambelli
class, and `tau_mu` for the integral type-B Schubert class. The standard
rational comparison preserves the `c_j` and sends

$$
\Theta_\mu(c)\longmapsto 2^{\ell_A(\mu)}\tau_\mu,
\qquad \ell_A(\mu)=\#\{i:\mu_i>A\}.
\tag{5}
$$

The required primary reference is
[Buch--Kresch--Tamvakis, *A Giambelli formula for isotropic Grassmannians*](https://math.umd.edu/~harryt/papers/isogiamfinal.pdf):
Definition 1 and Theorem 1 on PDF p.3 give the raising operator; Section 2.1
on PDF p.8 gives (5) and equation (12). The finite Schubert partitions fit
the `2`-by-`(2A+2)` rectangle. The equal-row vanishing used below is also
given in Section 1.3, equation (9), PDF p.6.

For two rows `(b,c)`, that raising operator specializes to

$$
\Theta_{(b,c)}=
\begin{cases}
c_b c_c-c_{b+1}c_{c-1},&b+c\le2A+1,\\[2pt]
c_b c_c+2\displaystyle\sum_{j=1}^{c}(-1)^j
                     c_{b+j}c_{c-j},&b+c>2A+1.
\end{cases}
\tag{6}
$$

The second expression may also be applied formally at `(b,c)=(m,m)`.
For `m>A` it is zero in cohomology. Equivalently, the squared-root
relation gives `h_m(y_1^2,y_2^2)=0` for `m>A`, since the complementary
tail has rank `A`.

## Proof in degrees `m <= A`

Set `S_j=s_(2m-j,j)(y_1,y_2)` for `0<=j<=m`. The elementary two-variable
plethysm identity is

$$
f_s=\sum_{j=2s}^{m}(-1)^j S_j.
\tag{7}
$$

It follows either by cancelling adjacent monomials or by multiplying
`s_(m-2s)(y_1^2,y_2^2)` by `(y_1y_2)^(2s)`.
In the present degree every partition `(2m-j,j)` is `A`-strict and fits
the finite Schubert rectangle. Since `2m<=2A`, the first branch of (6)
identifies `S_j=Theta_(2m-j,j)`.

Let `t=floor(m/2)`. Make the unimodular column change

$$
g_s=f_s-f_{s+1}\quad(0\le s<t),\qquad g_t=f_t.
$$

By (7), `g_s=S_(2s)-S_(2s+1)` for `s<t`. The last column is `S_m` when
`m` is even, and `S_(m-1)-S_m` when `m` is odd. The Schubert supports of
these columns are disjoint. In (5), their row exponents are

$$
\ell_A(2m-j,j)=\mathbf1_{2m-j>A},
$$

because `j<=m<=A`. Each column therefore has gcd

$$
2^{\varepsilon_s},\qquad
\varepsilon_s=
\mathbf1_{\,2m-\min(2s+1,m)>A}.
\tag{8}
$$

Disjoint supports mean that integer row operations reduce each column to
its gcd independently. Thus these are all the nonzero Smith factors.
For `r=m-s`, (8) is zero exactly when `r<=u` and `s<=v`.
No exponent two occurs in this degree range. This proves (1)--(2) for
`m<=A`, including the unit in degree zero.

## Proof in degrees `m > A`

Continue to use `S_j` from the previous section. Define the formal
Pfaffian expressions

$$
T_j=c_{2m-j}c_j+
          2\sum_{h=1}^{j}(-1)^h c_{2m-j+h}c_{j-h}.
$$

Expanding the two-row Jacobi--Trudi determinant gives the exact
telescoping identities

$$
T_j=\sum_{i=0}^{j}(-1)^{j-i}S_i,
\qquad S_j=T_j+T_{j-1},\qquad T_{-1}=0.
$$

Substitution into (7) therefore gives

$$
f_s=T_{2s-1}+(-1)^m T_m.
$$

The last term vanishes for `m>A`. For every allowed source pair
`r=m-s>=s>=m-A>=1`, the first term is a genuine theta class:

$$
\boxed{f_s=\Theta_{(2r+1,\,2s-1)}.}\tag{9}
$$

Its parts are distinct positive integers, its first part is at most
`2A+1`, and its degree is `2m>2A+1`. Thus it is an allowed finite
`A`-strict Schubert partition and the second branch of (6) applies.
Different `s` give distinct Schubert rows. Formula (5) turns (9) into

$$
f_s\longmapsto
2^{1+\mathbf1_{2s-1>A}}\tau_{(2r+1,\,2s-1)}.
\tag{10}
$$

The first part always exceeds `A` because `r>=m/2>A/2`.
Consequently (10) is already diagonal: its exponent is two precisely
when `s>u`, and otherwise one. There are no unit factors for `m>A`.
Together with (8), this proves the theorem in every degree.

## Normalization and verification

The theorem concerns the normalized generators `q_i`, without rescaling.
For the defining/vector Springer polynomial-character intersection
described in [INTEGRAL_LITERATURE.md](INTEGRAL_LITERATURE.md), the actual
degree-`d` lattice is `2^d L_B(a,k;d)`. Its factors are consequently
`2^(d+e(r,s))`. This follows from the separate integral-form calculation;
it is not implicit in (1).

The [verification script](experiments/verify_type_b_two_row.py) checks the
two-variable identities symbolically and compares (1) with both recorded
finite Smith datasets, including either orientation of the two blocks.
The proof above applies to all `A>=2`; the data are a check of its
implementation and conventions. The general case `min(a,k)>=3` still
requires a higher-row multiplicity argument.
