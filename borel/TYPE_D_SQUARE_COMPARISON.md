# Type D square images: an integral comparison and its remaining gap

**Status: repository deductions, 19 October 2026.** This note proves the
common square-input kernel, block duality, a finite primitive range, and
a one-sided comparison with type B. It does **not** prove equality of the
general B and D square-sector Smith lists. The finite checks below are
corroboration. No originality or external mathematical certification is
claimed.

## 1. Fixed inputs and lattices

Let \(a,k\geq1\), \(n=a+k\), and set
\[
X_B=\operatorname{OG}(a,2n+1),\qquad
X_D=\operatorname{OG}(a,2n).
\]
In D the usual maximal-parabolic notation applies when \(k\geq2\).
The boundary \(k=1\) means the parabolic omitting both fork nodes; it is
not a single-node maximal parabolic. Low-rank cases are interpreted by
their Grassmannian geometry.

Use the normalized square inputs
\[
q_i=e_i(x_{a+1}^2,\ldots,x_n^2),\qquad 1\leq i\leq k.
\]
In D the last input is \(q_k=p^2\), where
\(p=x_{a+1}\cdots x_n\). The present note concerns only the square
subalgebra, not the additional odd powers of \(p\).
Polynomial degree \(d\) means cohomological degree \(2d\).
Let \(L^B_d,L^{D,+}_d\) be the integer spans of all degree-\(d\)
polynomials in these inputs, inside integral cohomology. Smith factors
are taken in the saturation of each selected rational image, with unit
factors retained. Odd polynomial degrees have zero square image.

The integral Schubert bases and the Chern-class realization used here
are those of Buch--Kresch--Tamvakis:
[the odd/symplectic paper, Section 2.1](https://math.umd.edu/~harryt/papers/isogiamfinal.pdf)
and [the even paper, equations (5)--(7), Proposition 1, and Theorem 4](https://math.umd.edu/~harryt/papers/oggiam.pdf).
Proposition 1 is a rational invariant-subring statement; it does not
assert the square-sector integral comparison sought below.

## 2. The square-input kernels agree integrally

**Deduction.** The B and D maps from the fixed integer square-input
algebra have the same kernel. Their degree-\(2m\) images are free on
the images of
\[
s_\lambda(x_{a+1}^2,\ldots,x_n^2),\qquad
\lambda\subseteq a^k,\quad |\lambda|=m.
\tag{1}
\]
In particular their common rank is the number of these rectangle
partitions. This assertion concerns a basis of the image lattice;
the displayed classes need not individually be Smith basis vectors.

### Proof

First work over \(\mathbb Q\) on the complete flags. Put
\(u_i=x_i^2\), \(E=x_1\cdots x_n\), and
\(A=\mathbb Q[u_1,\ldots,u_n]\). The rational characteristic-map
kernels are
\[
I_B=(e_1(u),\ldots,e_n(u)),\qquad
I_D=(e_1(u),\ldots,e_{n-1}(u),E).
\tag{2}
\]
These are the ordinary additive Borel presentations; for the
characteristic-map framework and invariant kernel after localization
see [Calmes--Petrov--Zainoulline, Theorems 6.4 and 6.9 and Sections
13.6--13.7](https://arxiv.org/pdf/0905.1341).

There is a parity decomposition of polynomial rings
\[
\mathbb Q[x_1,\ldots,x_n]
=\bigoplus_{J\subseteq\{1,\ldots,n\}}x_J A,
\qquad x_J=\prod_{j\in J}x_j.
\]
The generators \(e_i(u)\) preserve every parity component. Multiplying
by \(E\) sends the \(J\)-component to the complementary parity
component. Its parity-zero component can therefore come only from
\(E\cdot E A=e_n(u)A\). Consequently
\[
I_D\cap A=(e_1(u),\ldots,e_n(u))A=I_B\cap A.
\tag{3}
\]
This uses the decomposition by parity, rather than an averaging
argument about integral cohomology.

Intersect (3) with the ring of symmetric polynomials in the tail
\(u_{a+1},\ldots,u_n\). The type-B kernel and integral rectangle
basis have already been proved by the primitive quaternionic
comparison in [TYPE_B_SMITH.md, Sections 1--2](TYPE_B_SMITH.md#1-a-primitive-symplectic-comparison-lattice).
Thus the rational D square-input kernel is the same rectangle ideal.
Both integral flag cohomology groups are torsion-free. An integer input
vanishes integrally precisely when its rational image vanishes.
Hence their fixed integral kernels agree as well, and (1) is an
integral image basis in D as it is in B.

Passing from a partial flag to its complete flag does not change this
kernel. Indeed its pullback is an integral primitive injection, as
explained in the next section.

## 3. Integral block duality

**Deduction.** For the square sector,
\[
\operatorname{Smith}_{D,+}(a,k;d)
=\operatorname{Smith}_{D,+}(k,a;d).
\tag{4}
\]

### Proof

The compact complete-to-partial flag projection has simply connected
base and free even integral cohomology in its base and flag fiber.
Its integral Serre spectral sequence collapses: the possible source
and target bidegrees of every differential have incompatible parities.
The base edge map is injective and its cokernel has a filtration with
free associated groups. It is therefore primitive. Composing an image
lattice with this pullback preserves its nonzero Smith factors.

On the compact complete flag \(SO(2n)/T\), right multiplication by
normalizer representatives gives integral Weyl-group ring
automorphisms. Coordinate permutations belong to \(W(D_n)\); no
individual sign change is required. Thus a permutation moves the
chosen tail alphabet to a head alphabet of the same size, preserving
its image-lattice Smith factors.

The global squared elementary invariants vanish in flag cohomology:
they vanish rationally by (2), hence integrally by torsion-freeness.
Writing \(E_H(t),E_T(t)\) for the two squared elementary series gives
\[
E_H(t)E_T(t)=1.
\]
Formal inversion expresses each coefficient of one series as an
integer complete-symmetric polynomial in the coefficients of the
other. Their generated integral image subrings coincide. Moving the
complementary block to the tail and using primitive flag pullback
again proves (4). This argument is valid for the \(k=1\) boundary
parabolic as well.

## 4. An improved finite primitive range and a closed boundary family

**Deduction.** All nonzero normalized square-sector Smith factors are
\(1\) whenever
\[
d\leq\max(a,k).
\tag{5}
\]

The stable full D image
\(\mathbb Z[q_1,\ldots,q_{k-1},p]\) is primitive by the integral
eta-polynomial specialization in
[TYPE_D_SMITH.md](TYPE_D_SMITH.md#general-deduction-stable-primitivity-and-a-finite-stable-range).
Its square subalgebra, with \(q_k=p^2\), is an additive direct summand:
retain precisely the even powers of \(p\) in the polynomial monomial
basis. Thus its stable image is primitive as well. The integral finite
Schubert specialization is an isomorphism in degree \(d\) when
\(a\geq d\), so every square factor is \(1\) in that range.
Apply (4) to obtain the range \(k\geq d\) too. This proves (5).
It is a sufficient range, not a claim that it is optimal.

When \(\min(a,k)=1\), write \(A=\max(a,k)\). The entire square-sector
list agrees with type B and has just one factor in each even degree
\(d=2m\), \(0\leq m\leq A\):
\[
\operatorname{Smith}_{D,+}(a,k;2m)=
\begin{cases}
(1),&2m\leq A,\\
(2),&A<2m\leq2A.
\end{cases}
\tag{6}
\]
All other degrees have empty lists. For \(a=1\), use
\(X_D=Q^{2k}\) and \(q_i=(-1)^ih^{2i}\). The powers \(h^d\) are
primitive through degree \(k\) and twice a primitive class above it;
in middle degree \(h^k\) is the primitive sum of the two ruling
classes. The integral quadric formulas are stated in
[Dimitrov--Roth, Section 2.8, equation (2.8.1)](https://arxiv.org/pdf/1707.06840)
and applied in [TYPE_D_SMITH.md](TYPE_D_SMITH.md#quadric-deduction-for-a1).
Equation (4) covers the other block order. Comparing with the
one-block B formula proves the equality claimed in this boundary
family. Equation (6) concerns the square sector alone.

## 5. A one-sided comparison in every rank

**Deduction.** In a fixed even degree, let \(r\) be the common rank
from (1), and write the sorted Smith exponent lists, including zeros,
as
\[
b_1\leq\cdots\leq b_r,\qquad
d_1\leq\cdots\leq d_r.
\]
Then
\[
\boxed{d_i\geq b_i\quad(1\leq i\leq r).}
\tag{7}
\]
Equivalently, the D square factors are termwise divisible by the
corresponding B factors.

### Proof

Fix a nonisotropic vector \(v\) in the odd-dimensional quadratic
space. Its orthogonal complement is nondegenerate of dimension
\(2n\), and gives a closed embedding
\[
i:X_D\hookrightarrow X_B,\qquad U\longmapsto U.
\]
It is the zero scheme of the section
\(U\mapsto(v,-)|_U\) of \(S_B^*\). It is smooth of the expected
codimension \(a\), and the tautological bundles agree under
restriction. The polynomial characteristic classes \(q_i\) agree
too, either by the squared elementary-series relation or by their
common polynomial in the quotient Chern classes.

Use the common integral source basis (1), and let \(M_B,M_D\) be its
two coefficient matrices in integral Schubert bases. Integral
restriction has an integer matrix \(T\), so
\[
M_D=T M_B.
\]
Consequently their coefficient-row lattices in \(\mathbb Z^r\)
satisfy
\[
R_D\subseteq R_B.
\tag{8}
\]
Both have full rank by Section 2. Transposition preserves nonzero
Smith factors, so the finite groups
\(\mathbb Z^r/R_B\) and \(\mathbb Z^r/R_D\) have exactly the Smith
factors under discussion. Equation (8) supplies a surjection from
the second group to the first. They are 2-primary: the B assertion is
in [TYPE_B_SMITH.md](TYPE_B_SMITH.md#3-only-the-prime-two-and-a-bound-on-each-exponent),
and the D assertion follows from the no-odd-factor argument in
[TYPE_D_SMITH.md](TYPE_D_SMITH.md#general-deduction-every-finite-smith-factor-is-a-power-of-2),
restricted to the square eigensector after inverting 2.

For every \(s\geq1\), the induced map on
\(2^{s-1}G/2^sG\) is surjective. Its dimension counts the exponents
at least \(s\). Thus D has at least as many such exponents as B at
every threshold. For sorted lists of equal length this is equivalent
to (7).

The unproved reverse inclusion in (8) is precisely what would turn
this bound into equality of the coefficient-row lattices. Equality
of the Smith multisets alone is weaker, but remains unproved in
general here as well.

## 6. A failed Gysin shortcut

The embedding above satisfies
\[
i_*i^*h=c_a(S_B^*)h.
\]
One might try to prove primitive restriction on the selected square
span by asserting that multiplication by \(c_a(S_B^*)\) preserves
its saturation. That assertion is false.

For \(a=k=2\), in polynomial degree \(4\), choose the source Schur
basis \(s_{(2)}(x_T^2),s_{(1,1)}(x_T^2)\). Its B coefficient matrix
has factors \((1,2)\), whereas the matrix after multiplying both
columns by \(x_1x_2=c_2(S_B^*)\), up to the immaterial Chern-root
sign, has factors \((2,2)\). Multiplication has introduced an
additional index of two in the selected rational span. Therefore
this purported primitive-multiplication lemma cannot be used to
prove square-sector equality.

The standalone experiment
[verify_type_d_square_comparison.py](experiments/verify_type_d_square_comparison.py)
checks this example using both the library's fundamental-weight BGG
expansions and an independent Cartesian type-C frontier with exact
B/C row rescaling. It also checks the asserted row inclusion,
rectangle ranks, primitive range, and finite B/D square comparisons
in its specified cases. Its
[report](experiments/type_d_square_comparison_verification.json)
records the precise scope. Finite equalities do not establish the
general reverse inclusion.

For the full D image, the separately proved Euler-sector comparison
and the integral eigenspace-gluing formula are still required; see
[TYPE_D_SMITH.md](TYPE_D_SMITH.md#exact-euler-sector-comparison-with-type-b)
and [ORTHOGONAL_SMITH.md](ORTHOGONAL_SMITH.md#type-d-an-exact-integral-gluing-formula).
Even a future square-sector equality theorem would not by itself
settle the full type-D partition formula.
