# Full type-D gluing for a two-dimensional head and an even tail

**Status: internally audited complete repository deduction, 19 October 2026.**
The argument below has received a fresh adversarial audit within this
session. It is unpublished, without an originality claim or external
expert certification. It determines the full normalized image when the
head has dimension two and the tail dimension is even. It does not use a
full head/tail block duality: the Euler generator belongs to the specified
tail.

## 1. Statement

Let \(k=2h\geq2\), and set
\[
X=\operatorname{OG}(2,2(k+2)).
\]
Use the normalized inputs
\[
q_i=e_i(x_{\mathrm{tail}}^2),\qquad
p=\prod x_{\mathrm{tail}},\qquad p^2=q_k.
\]
Polynomial degree \(d\) means cohomological degree \(2d\).
Write \(L_+\) for the square-input image and \(L_-\) for the image of
\(p\mathbb Z[q]\) in a fixed degree. Let \(S_\pm\) be their
saturations in integral cohomology and let
\[
S=((L_++L_-)\otimes\mathbb Q)\cap H^{2d}(X;\mathbb Z),
\qquad W=S/(S_+\oplus S_-).
\]

**Theorem.** For a two-dimensional head and every even \(k\),
\[
\dim_{\mathbb F_2}W=
\begin{cases}1&d=2k,\\0&d\ne2k.\end{cases}
\tag{1}
\]
In the one nonzero case both gluing directions are represented by
original image classes in \(L_+\) and \(L_-\). Thus all gluing is unit
gluing, and the full Smith list is the union of the two known sector
lists, except that one pair \((1,1)\) in degree \(2k\) is replaced by
\((1,2)\).

Equivalently, with the type-B multiplicity polynomial
\(\mathcal B_{a,k,e}\) defined in
[TYPE_D_CANDIDATE.md](TYPE_D_CANDIDATE.md),
\[
\sum_d M_D(2,k;e,d)z^d
=\mathcal B_{2,k,e}(z^2)+z^k\mathcal B_{1,k,e}(z^2)
 +(\delta_{e,1}-\delta_{e,0})z^{2k}.
\tag{2}
\]
The variable in \(\mathcal B\) records half the square polynomial
degree. Units are retained and all indices concern the selected rational
image, rather than the full cohomology of \(X\).

The sector comparisons used in (2) are separately proved in
[TYPE_D_TWO_ROW_SMITH.md](TYPE_D_TWO_ROW_SMITH.md) and
[TYPE_D_SMITH.md](TYPE_D_SMITH.md#exact-euler-sector-comparison-with-type-b).
The new step is (1).

## 2. The two-root Chern model and the quotient by square classes

Let \(x,y\) be the roots of the dual tautological two-plane, and put
\[
u=x^2,\quad v=y^2,\quad s=u+v,\quad t=uv,\quad e=xy,\quad H_j=h_j(u,v).
\]
Take \(H_0=1\) and \(H_j=0\) for \(j<0\). Define
\[
R=\mathbb Z[u,v]^{S_2}/(H_{k+1},H_{k+2}),\qquad
R_0=\mathbb Z[u,v]^{S_2}/(H_k,H_{k+1}).
\tag{3}
\]
Their Schur rectangles have widths \(k\) and \(k-1\), respectively,
and height two.

The integral symplectic comparison ring is
\[
C=H^*(\operatorname{IG}(2,2k+4);\mathbb Z)
  =R[e,x+y]/(e^2-t,(x+y)^2-s-2e).
\]
Its even polynomial-degree part is \(R\oplus eR\). The integral
Leray--Hirsch justification, using fiber
\(\operatorname{LG}(2,4)=Q^3\), is given in
[TYPE_D_TWO_ROW_SMITH.md](TYPE_D_TWO_ROW_SMITH.md#3-a-two-variable-integral-chern-model).
In particular these are integral free modules, not only rational
presentations.

Write \(c_j=h_j(x,y)\) for the universal quotient Chern classes.
The generating series gives
\[
c_{2j}=H_j+eH_{j-1},\qquad c_{2j+1}=(x+y)H_j.
\tag{4}
\]
Put \(M=2k+2\). For a two-row \(k\)-strict partition \(\lambda=(r,l)\),
the symplectic Schubert polynomial is
\[
\Theta_{(r,l)}=
\begin{cases}
c_rc_l-c_{r+1}c_{l-1},&r+l\leq2k+1,\\
c_rc_l+2\sum_{j=1}^l(-1)^jc_{r+j}c_{l-j},&r+l\geq2k+2.
\end{cases}
\tag{5}
\]
For one row use \(\Theta_{(r)}=c_r\).
For even size, the even-orthogonal unnormalized raising polynomial is
the same: its threshold differs only at the odd value \(2k+1\).
Thus for retained partitions of length at most two and first part at
most \(M-1\),
\[
                         Y_\lambda=2^{-\nu(\lambda)}\Theta_\lambda,
\qquad \nu(\lambda)=\#\{i:\lambda_i>k\}.
\tag{6}
\]
These \(Y_\lambda\) form the integral plus-eigenlattice basis:
one Schubert class if \(\lambda\) has no part \(k\), or the sum
\(\tau_\lambda+\tau'_\lambda\) if it has a part \(k\).

The published inputs are the symplectic Giambelli formula and the even
Giambelli formula and typed basis:
[Buch--Kresch--Tamvakis, Theorem 1](https://math.umd.edu/~harryt/papers/isogiamfinal.pdf)
and
[Buch--Kresch--Tamvakis, Theorems 1 and 4, equation (16)](https://math.umd.edu/~harryt/papers/oggiam.pdf).

The Chern-preserving map from \(C\) to the D Chern subring kills
\(c_M=eH_k\). Its even kernel is exactly the span of the missing
Schubert classes
\[
\Theta_{(M,2j)}=e\,s_{(k,j)}(u,v),\qquad 0\leq j\leq k.
\tag{7}
\]
For completeness, \(c_M=eH_k\), \(H_kH_j=s_{(k,j)}\) in \(R\),
and \(tH_kH_{j-1}=0\). All raising terms with first Chern index
greater than \(M\) vanish. The retained classes map to the independent
classes \(2^{\nu(\lambda)}Y_\lambda\), so there are no other even
kernel vectors. Consequently the even D Chern subring has the
integral additive decomposition
\[
                              C_D^{\mathrm{even}}=R\oplus eR_0.
\tag{8}
\]

Over \(\mathbb Q\), the square sector is exactly \(R_\mathbb Q\).
Define the extraction map
\[
E:H^{\mathrm{even\ polynomial\ degree}}(X;\mathbb Z)^+
       \longrightarrow (R_0)_\mathbb Q
\tag{9}
\]
by taking the coefficient of \(e\) in (8), extended rationally.
Its kernel is \(S_+\), the integral square-sector saturation.
Let \(T\) be its image. Therefore
\[
                    H^+(X;\mathbb Z)/S_+\cong T.
\tag{10}
\]
All statements below may be read degree by degree.

## 3. The exact integral extraction lattice

For \(\gamma=(r,l)\subseteq(k-1)^2\), define
\[
\tau(\gamma)=
\begin{cases}
0&r\leq h-1,\\
2&l\geq h+1,\\
1&\text{otherwise}.
\end{cases}
\tag{11}
\]

**Lemma.** The extraction lattice is the following diagonal lattice
in the specified Schur coordinates:
\[
                 T=\bigoplus_{\gamma\subseteq(k-1)^2}
                         2^{-\tau(\gamma)}\mathbb Z\,s_\gamma.
\tag{12}
\]

Here are the coefficient calculations, including all required
divisibility.

First, \(T\) contains the entire integer lattice \(R_0\). Indeed every
\(e\,s_\gamma\) is an integral Chern polynomial, so its extraction
\(s_\gamma\) lies in \(T\).

For two indices \(i>j\geq0\), (4) gives
\[
\begin{aligned}
E(c_{2i}c_{2j})&=H_iH_{j-1}+H_{i-1}H_j,\\
E(c_{2i+1}c_{2j+1})&=2H_iH_j.
\end{aligned}
\tag{13}
\]
When the first branch of (5) applies, these become
\[
\begin{aligned}
E(\Theta_{(2i,2j)})&=H_{i-1}H_j-H_iH_{j-1}=s_{(i-1,j)},\\
E(\Theta_{(2i+1,2j+1)})&=H_iH_j-H_{i+1}H_{j-1}=s_{(i,j)}.
\end{aligned}
\tag{14}
\]
The formulas include \(j=0\) with \(H_{-1}=0\); the one-row formula
has \(E(c_{2i})=H_{i-1}\).
For the second branch, all nonleading raising terms in (5) are
twice integral polynomials. The two-variable Pieri rule gives
\[
H_iH_{j-1}+H_{i-1}H_j\equiv s_{(i-1,j)}\pmod2.
\tag{15}
\]
Indeed the two products have identical Schur summands with second
part \(0,\ldots,j-1\); the second product has the one additional
summand with second part \(j\). This also holds after rectangle
truncation.

We prove containment in the right side of (12). For \(\nu=0\),
(6) has integral extraction. For \(\nu=1\), write the partition as
\((r,l)\) with \(r>2h\), \(l\leq2h\), and even total size.
If both parts are even, \(r=2i,l=2j\) with \(i>h\), so both (14)
and (15) show that every odd extraction coefficient has first Schur
part at least \(h\). If both are odd, \(r=2i+1,l=2j+1\) with
\(i\geq h\), \(j\leq h-1\). In the first branch (14) its odd
coefficient again has first part at least \(h\), and in the second
branch (13) has no odd extraction coefficients at all. Thus dividing
by \(2\) never puts a half into a coordinate with first part at most
\(h-1\).

For \(\nu=2\), both parts exceed \(2h\) and are distinct.
The extraction degree in \(u,v\) is at least \(2h+1\); a partition
with first part at most \(h-1\) has size at most \(2h-2\).
Those low coordinates are therefore absent. Formula (13), (15),
and the factor two on all raising corrections show that odd
extraction coefficients can occur only at
\((i-1,j)\) with \(j>h\). Thus division by \(4\) produces at most
halves in the remaining coordinates with second part at most \(h\),
and quarters only in coordinates with second part at least \(h+1\).
This proves containment.

For the reverse containment, take any \(\gamma=(r,l)\) with
\(r\geq h\) and \(0\leq l\leq r\leq k-1\), and choose
\[
                         \lambda=(2r+2,2l).
\tag{16}
\]
This is retained and \(k\)-strict, since \(2r+2\leq2k<M\)
and \(2r+2>2l\). Equations (14)--(15) give
\[
                         E(\Theta_\lambda)\equiv s_\gamma\pmod2.
\tag{17}
\]
If \(l\leq h\), it has \(\nu=1\), so \(E(Y_\lambda)\) is
\(s_\gamma/2\) plus an integer vector. If \(l>h\), it has
\(\nu=2\), so \(E(2Y_\lambda)\) has that form. Since \(R_0\subset T\),
this proves that \(s_\gamma/2\in T\) for every coordinate with
first part at least \(h\).

For \(l>h\), now use \(Y_\lambda\) itself in (16).
It has extraction \(s_\gamma/4\) plus a half-integral vector.
Every other coordinate in that vector has the same degree
\(|\gamma|\geq2h+2\); hence none has first part at most \(h-1\).
All those halves have just been proved to belong to \(T\).
Subtracting them proves \(s_\gamma/4\in T\).
Together with the integer coordinates this proves (12).

## 4. Gluing is a divisibility question in \(T\)

The Euler sector in head dimension two is especially simple.
By the integral Euler comparison with the smaller B Grassmannian,
it has rank one in each degree
\[
                           d=k+2j,\qquad 0\leq j\leq k,
\tag{18}
\]
and zero otherwise. Its saturation is the primitive line on
\[
Z_\lambda=\tau_\lambda-\tau'_\lambda,\qquad
\lambda=\text{the sorted pair }(k,2j).
\tag{19}
\]
The actual Euler image has coefficient \(1\) on that line if
\(2j\leq k\), and coefficient \(2\) if \(2j>k\).
This is the one-dimensional-head type-B quadric formula transported
by the Euler comparison.

The diagram involution swaps \(\tau_\lambda,\tau'_\lambda\).
Modulo 2, \(Z_\lambda\) and \(Y_\lambda\) have the same class.
All other invariant-basis classes have disjoint Schubert supports.

**Criterion.** The space \(W\) in degree (18) is nonzero precisely
when
\[
                              E(Y_\lambda)\in2T.
\tag{20}
\]
It then has dimension one.

To prove the criterion, \(S_-\) has rank one and its mod-2 image is
the line on \(Z_\lambda\). For primitive sector lattices in a free
ambient lattice, the gluing space identifies with the intersection
of their mod-2 images; this is the integral gluing lemma in
[ORTHOGONAL_SMITH.md](ORTHOGONAL_SMITH.md#lattice-lemma).
Therefore \(W\ne0\) if and only if
\[
                          Y_\lambda\in S_++2H^+(X;\mathbb Z).
\]
It is harmless to require the correcting class to be invariant:
if an invariant difference equals \(2h\), torsion-freeness forces
\(\iota h=h\). Taking the quotient (10) gives exactly (20).
Equivalently, the entire criterion follows by lifting the equality
of the two sector reductions in the typed Schubert basis.
No statement about a rational eigenspace alone is being used.

## 5. Testing every possible Euler degree

Write \(k=2h\) as before.

If \(0\leq j<h\), the typed partition is \((k,2j)\) and has
\(\nu=0\). Its first-branch formula (14) gives exactly
\[
                            E(Y_{(k,2j)})=s_{(h-1,j)}.
\tag{21}
\]
This coordinate has \(\tau=0\) in (11). Every vector in \(2T\)
has an even integer coefficient on it. Hence (21) is not in \(2T\).

If \(h<j\leq k\), the partition is \((2j,k)\) and has \(\nu=1\).
Its second-branch formula, together with (15), gives
\[
                E(Y_{(2j,k)})=\tfrac12s_{(j-1,h)}
                                    +\text{an integer vector}.
\tag{22}
\]
The coefficient at the displayed coordinate is an odd half-integer.
That coordinate has \(\tau=1\) in (11); its coefficient in every
vector of \(2T\) is integral. Thus (22) is not in \(2T\), regardless
of its other coordinates.

For \(j=h\), the partition is \((k,k)\) and its first branch has
zero extraction. In fact it is the original square-input class
\[
\begin{aligned}
Y_{(k,k)}
 &=c_k^2-c_{k+1}c_{k-1}\\
 &=H_h^2-sH_hH_{h-1}+tH_{h-1}^2
 =t^h=s_{(h,h)}(u,v).
\end{aligned}
\tag{23}
\]
The last identity follows by induction from
\(H_j=sH_{j-1}-tH_{j-2}\): if
\(F_j=H_j^2-sH_jH_{j-1}+tH_{j-1}^2\), substitution gives
\(F_j=tF_{j-1}\), and \(F_1=t\).
Thus (20) holds exactly in this case. Its degree is \(d=2k\).

Moreover \(Y_{(k,k)}\in L_+\), and
\(Z_{(k,k)}\in L_-\): in (18) this is \(2j=k\), the last
primitive degree of the one-head B image. Their half-sum is the
individual typed Schubert class \(\tau_{(k,k)}\).
This is a unit-unit gluing pair. The lattice lemma changes precisely
\((1,1)\) to \((1,2)\), proving (1)--(2).

## 6. A separate full-image exponent bound

The following bound is valid for every head size \(a\), independently
of the preceding two-head calculation:
\[
               \text{every full normalized D factor divides }2^a.
\tag{24}
\]
It is not a head/tail symmetric bound.

Let \(A\) be the polynomial subring generated by the head Chern classes
and the tail Euler class \(p\) on the original \(X_D\).
The ordinary eta formula expresses each typed Schubert class with
denominator at most \(2^a\): a type-zero class uses \(2^{-\nu}\),
and a positive-type class has one additional possible half but
\(\nu+1\leq a\). The remaining factors are integral Chern polynomials
and \(p\). Therefore \(2^aH^*(X_D;\mathbb Z)\subset A\).

There is an integral additive retraction from \(A\) onto the full
input image \(L=\mathbb Z[q,p]\). Represent a class by
\(F_0(x_{\rm head})+pF_1(x_{\rm head})\), reduce \(p^2=q_k\),
and delete all monomials with an odd exponent in any head root in
both \(F_i\). This fixes precisely the square-polynomial and
Euler-times-square polynomial images. It descends: on the complete
flag the rational D characteristic ideal is generated by the squared
invariants and the global Euler invariant, which is zero. Individual
head sign changes, leaving \(p\) fixed, preserve that ideal.
Coefficient deletion preserves its rational kernel; torsion-free
cohomology and injective complete-flag pullback give integral descent.
Unlike the relative incidence case with nonzero base Euler, the
relation \(c_a(S^*)p=0\) is preserved by this projection.

Hence \(L\) is primitive inside \(A\). If \(z\) is an integral class
in its rational span, then \(2^az\in A\) and the retraction fixes
\(2^az\), so \(2^az\in L\). This proves (24).

This particular retraction gives the bound \(2^a\). The separate
[full exponent-bound proof](TYPE_D_FULL_EXPONENT_BOUND.md) also proves
the bound \(2^k\), by moving the specified roots to a new head and
projecting onto its even-sign invariants. Together the two arguments
give \(2^{\min(a,k)}\). The relocated source is a head invariant
algebra, rather than the full tail algebra for the exchanged blocks.

## 7. Scope and normalization

This completes the even-tail gluing problem for \(a=2\), and the
one-head case is the already proved quadric formula.
The separate [general gluing theorem](TYPE_D_EVEN_GLUING_THEOREM.md)
now establishes exhaustion in every head size; the present proof is an
independent calculation for the two-head family.

All statements use the normalized fixed input lattice. For the literal
defining/vector character form, multiply every degree-\(d\) Smith
factor by \(2^d\), as in
[TYPE_D_SMITH.md](TYPE_D_SMITH.md#literal-polynomial-character-form).

## 8. Concrete finite verification

The Sage script
[verify_type_d_two_head_gluing.py](experiments/verify_type_d_two_head_gluing.py)
and its
[report](experiments/type_d_two_head_gluing_verification.json)
check every even polynomial degree for \(k=2,4,\ldots,24\).
They use the independently audited two-root Chern/Littlewood--Richardson
model from the two-row square verifier. No Borel, Weyl-group,
divided-difference, or candidate-verifier code is imported.

The checks establish exact equality of the scaled rational extraction
lattice with the diagonal lattice (12) in all 324 degree components.
All 168 Euler-line tests satisfy (20) exactly in the 12 middle
components. Separately, the script constructs the full coefficient
matrix in the individual typed Schubert basis, using equal square
coefficients and opposite Euler coefficients on a typed pair.
All 324 actual integer Smith lists agree with (2).

These finite checks corroborate the written all-even-tail proof.
They do not replace the integral coordinate argument or the cited
Schubert formulas. The report records the exact script, proof, and
shared Chern-model byte hashes.
