# Orthogonal Springer images: formulas, proofs, and the remaining question

**Status: audited repository deductions (2026-10-19), with a separate
type-D conjecture.** [PROOF_STATUS.md](PROOF_STATUS.md) identifies the
published inputs and the arguments checked here. This audit does not
establish novelty or constitute external mathematical review.

The notes derive several uniform statements for the normalized type B and D
tail generators. The [general type B proof](TYPE_B_GENERAL_SMITH.md)
now derives every normalized multiplicity in every rank, closing the
formerly missing integral comparison. A single
[closed partition formula](TYPE_B_CLOSED_CONJECTURE.md) describes all
layers. The type-D partition classification remains open.
The unit counting polynomial is already classical Birkhoff/Hall enumeration.
Its identification with the chosen Smith units is an additional argument
here, rather than a theorem attributed to Birkhoff.

## Conventions and results

Put $a=n-k$. The varieties are

\[
X_B=\operatorname{OG}(a,2n+1),\qquad
X_D=\operatorname{OG}(a,2n).
\]

Take $a,k\geq1$ in B and $a\geq1,k\geq2$ in D. The separate D note also
explains the $SO(2)$-tail boundary case.

The fixed normalized inputs are $q_i=e_i(x_{a+1}^2,\ldots,x_n^2)$.
In D include $p=x_{a+1}\cdots x_n$, with $p^2=q_k$. Thus the source
is $\mathbb Z[q_1,\ldots,q_k]$ in B and

\[
\mathbb Z[q_1,\ldots,q_{k-1},p]
\]

in D. Polynomial degree $d$ means cohomological degree $2d$;

\[
\deg q_i=2i,\qquad \deg p=k.
\]

All indices are taken inside the saturation of the rational image, which
can be smaller than the full cohomology group.

| Statement | Type B | Type D |
| --- | --- | --- |
| Possible normalized factors | $2^e$, $0\leq e\leq\min(a,k)$ | Powers of $2$ |
| Stable image, with $k$ fixed | Saturated in every degree | Saturated in every degree |
| Sufficient finite range for all factors $1$ | $d\leq\max(a,k)$ | $d\leq a$ |
| Complete closed families | Every pair $a,k$ | $a=1$ |
| Other uniform information | Block symmetry, rank formula, closed unit and largest-factor multiplicities, weighted-minor formula for every factor | Euler-sector comparison and the integral gluing formula below |
| Remaining classification | None for these normalized B generators | Closed partition multiplicities, including the integral gluing |

The detailed proofs are in [TYPE_B_SMITH.md](TYPE_B_SMITH.md) and
[TYPE_D_SMITH.md](TYPE_D_SMITH.md). The stable assertions strengthen
injectivity to primitivity by an explicit integral specialization into a
polynomial invariant ring. They do not follow from injectivity alone.

For the literal defining/vector polynomial-character forms specified in
[the normalization audit](INTEGRAL_LITERATURE.md), every degree-$d$
Smith factor is multiplied by $2^d$. This is a change of the chosen input
lattice, not a change in the cohomology basis.

## The type B formula at arbitrary finite rank

Write $d=2m$; odd-degree images are zero. Columns can be indexed by the
partitions $\lambda\subseteq a^k$ of size $m$. Let

\[
C_{w,\lambda}=
\bigl(\partial_w^C s_\lambda(x_{a+1}^2,\ldots,x_n^2)\bigr)(0),
\]

where rows run through degree-$d$ minimal parabolic representatives.
This is the integral type C matrix, whose column lattice is primitive.
If $\nu(w)$ is the number of negative entries of the signed permutation,
the type B matrix is exactly $B_{w,\lambda}=2^{\nu(w)}C_{w,\lambda}$.

Set $r=\#\{\lambda\subseteq a^k:|\lambda|=m\}$, $\delta_0=0$, and

\[
\delta_j=\min_{\substack{|I|=|J|=j\\\det C_{I,J}\ne0}}
\left(\sum_{w\in I}\nu(w)+v_2(\det C_{I,J})\right).
\]

Then every nonzero Smith factor, in divisibility order, is

\[
\boxed{s_j=2^{\delta_j-\delta_{j-1}},\qquad 1\leq j\leq r.}
\]

This is an exact all-rank determinantal formula with specified Schur and
Schubert indices. Its proof is row rescaling plus determinantal divisors.
It still requires a matrix calculation; it is not the desired closed
partition statistic for the multiplicities.

The number of factors equal to $1$ has the smaller description

\[
\operatorname{rank}_{\mathbb F_2}
\left(H^{2m}(\operatorname{Gr}(k,n);\mathbb F_2)
\xrightarrow{f\mapsto f^2}
H^{4m}(\operatorname{Gr}(k,n);\mathbb F_2)\right).
\]

It is verified independently using ordinary Schur multiplication in
[verify_type_b_frobenius.py](experiments/verify_type_b_frobenius.py), rather
than the library's BGG implementation.

### Closed type B partition formulas

The [two-row theorem](TYPE_B_TWO_ROW_SMITH.md) gives every factor when
$\min(a,k)=2$. For $A=\max(a,k)$, $u=\lceil A/2\rceil$,
$v=\lfloor A/2\rfloor$, attach to each partition $(r,s)\subseteq A^2$
the exponent

\[
e(r,s)=
\begin{cases}
0&r\leq u,\ s\leq v,\\
2&s>u,\\
1&\text{otherwise}.
\end{cases}
\]

In polynomial degree $2m$, take $2^{e(r,s)}$ for the pairs with $r+s=m$.
The [three-row theorem](TYPE_B_THREE_ROW_SMITH.md) gives all four
multiplicities when $\min(a,k)=3$, using closed Gaussian polynomials and
bounded three-part partition counts.

For every $l=\min(a,k)$, the [largest-factor theorem](TYPE_B_HIGHEST_SMITH.md)
gives the multiplicity of $2^l$ in every degree:

\[
\sum_m M_l(m)q^m
=q^{l(T+1)}{A-T-1+l\brack l}_q,\qquad
T=\left\lfloor\frac{A+l-1}{2}\right\rfloor.
\]

The [incidence and Schur-product filtration](TYPE_B_FILTRATION.md)
proves the common reduction and identifies the remaining question for
larger blocks.

The [unit-factor theorem](TYPE_B_UNIT_SMITH.md) supplies a positive
Gaussian-polynomial count in all ranks, strengthening the Frobenius-rank
description to a closed formula. The
[complete candidate](TYPE_B_CLOSED_CONJECTURE.md) organizes every factor
by one partition sum. The [general integral proof](TYPE_B_GENERAL_SMITH.md)
now establishes every layer, and [the high-layer note](TYPE_B_HIGH_LAYERS.md)
provides an independent elementary argument for the largest three.

## Type D: an exact integral gluing formula

In a fixed degree let $L_+$ be the lattice generated by the square-only
monomials, and $L_-$ the lattice generated by $p$ times square-only
monomials. The diagram involution fixes the first and negates the second.
Since $p^2=q_k$, the full image is $L=L_+\oplus L_-$.
The direct sum here is a sum of the two image lattices, not a claim that
their saturations exhaust the saturation of $L$.

The following lattice lemma gives every full Smith factor from the two
sector Smith lists and a finite-dimensional calculation over $\mathbb F_2$.
It applies in every rank. It is useful even when a closed partition formula
for the sector lists is unavailable.

The powers-of-two hypothesis holds for each D sector as well as for the
combined image: after inverting $2$, the combined image is primitive, and
the idempotents $(1\pm\iota)/2$ split it into the two sectors.

The Euler-sector list is already determined by type B:

\[
\operatorname{Smith}(L_-^{D}(a,k;d))
=\operatorname{Smith}(L^{B}(a-1,k;d-k)).
\]

The degree shift is the degree of $p$. For $a=1$ the variety on the right
is a point. The proof in [TYPE_D_SMITH.md](TYPE_D_SMITH.md) upgrades
[Buch--Kresch--Tamvakis, Proposition 2](https://math.umd.edu/~harryt/papers/oggiam.pdf)
using the explicit integral Schubert bases and their primitive difference
lattice. The analogous square-sector comparison with $B_n$ holds in the
recorded examples, but is not proved here.

### Lattice lemma

Let $H$ be a free abelian group with an involution $\iota$, and

\[
L_\pm\subset H^{\pm},\qquad
S_\pm=(L_\pm\otimes\mathbb Q)\cap H,\qquad
S=((L_+\oplus L_-)\otimes\mathbb Q)\cap H.
\]

Assume the sector Smith factors are powers of $2$, with exponent multisets

\[
E_+=(e_1^+,\ldots,e_{r_+}^+),\qquad
E_-=(e_1^-,\ldots,e_{r_-}^-).
\]

The quotient $W=S/(S_+\oplus S_-)$ is an $\mathbb F_2$-vector space.
It injects into each $S_\pm/2S_\pm$ by

\[
[s]\longmapsto s\pm\iota s\pmod{2S_\pm}.
\]

For $i\geq0$ define the intrinsic sector flags

\[
F_i^\pm=\operatorname{image}\left(
(2^{-i}L_\pm\cap S_\pm)\longrightarrow S_\pm/2S_\pm\right).
\]

Pull them back to flags $W_i^\pm$ on $W$, set negative-index flags
equal to zero, and put

\[
c_{ij}=\dim_{\mathbb F_2}(W_i^+\cap W_j^-),\qquad
m_{ij}=c_{ij}-c_{i-1,j}-c_{i,j-1}+c_{i-1,j-1}.
\]

Then $m_{ij}\geq0$, and the full Smith exponents are obtained as follows:

\[
\boxed{\text{For each of the }m_{ij}\text{ pairs, replace }(i,j)
\text{ in }E_+\sqcup E_-\text{ by }
\bigl(\min(i,j),\max(i,j)+1\bigr).}
\]

Leave the unpaired exponents unchanged and sort the resulting list.
In particular,

\[
\boxed{[S:L]=2^{\dim W}[S_+:L_+][S_-:L_-],\qquad
\dim W\leq\min(r_+,r_-).}
\]

### Proof

For $s\in S$, the elements $s+\iota s$ and $s-\iota s$ lie in
$S_+$ and $S_-$. Their sum is $2s$, so $2W=0$.
If either displayed map is zero, half of that component lies in the
appropriate $S_\pm$; subtracting it from $s$ gives an element of the
other sector. Thus both maps are injective. They identify $W$ with a
graph between two subspaces of the sector reductions modulo $2$.

Choose sector Smith bases. In such a basis $F_i^\pm$ is spanned by the
coordinates whose Smith exponent is at most $i$. The reduction modulo
$2$ of the group of integral automorphisms preserving both $S_\pm$
and $L_\pm$ is the full parabolic group preserving this flag. To see
surjectivity, use permutations and elementary matrices inside equal-exponent
blocks, and elementary matrices from a higher-exponent column to a
lower-exponent row. All have unimodular integral lifts preserving the
sector lattice, and their reductions generate that parabolic group.

Two flags on a vector space have a common adapted basis: choose successive
complements in their intersections. The number of basis vectors entering
the two flags at levels $i,j$ is precisely $m_{ij}$. Extend the images
of this common basis to adapted bases of $S_\pm/2S_\pm$, and lift the
changes of basis by the preceding parabolic-surjectivity observation.
We have therefore chosen genuine sector Smith bases in which the graph
for $W$ pairs individual coordinates with exponents $i,j$.

For one paired pair write its sector basis vectors as $u,v$. In $S$
the corresponding lattice is generated by $t=(u+v)/2$ and $v$.
The image-lattice generators $2^iu,2^jv$ have matrix

\[
\begin{pmatrix}2^{i+1}&0\\-2^i&2^j\end{pmatrix}
\]

in this basis. Its entry gcd is $2^{\min(i,j)}$, and its determinant
has absolute value $2^{i+j+1}$. Its Smith factors are therefore

\[
2^{\min(i,j)},\quad 2^{\max(i,j)+1}.
\]

The unpaired coordinates give the original sector factors. Taking the
direct sum of these blocks proves the formula. The quotient $W$ also
gives the stated index identity directly. This proof is integral; the
rational eigenspace decomposition alone would not give the answer.

### A concrete obstruction to taking a union

For $D_4$, $a=k=2$, in polynomial degree $4$, columns
$(q_2,q_1^2,pq_1)$ on four Schubert-support rows give

\[
M=\begin{pmatrix}
1&2&-1\\1&2&1\\2&2&0\\-2&-2&0
\end{pmatrix}.
\]

The square sector has factors $(1,2)$, the Euler sector has factor
$(1)$, and the combined factors are $(1,2,2)$. In these coordinates

\[
S=\{(b,c,z,-z)\},\quad
S_+=\{(x,x,z,-z)\},\quad
S_-=\{(y,-y,0,0)\}.
\]

The sum $S_++S_-$ imposes $b\equiv c\pmod2$, so $\dim W=1$.
The glued pair has exponents $(0,0)$, which the lemma replaces by
$(0,1)$; the remaining exponent $1$ gives the full list $(0,1,1)$.

For odd $k$, the two sectors occur in opposite polynomial-degree
parities. In every degree one sector is zero, so $W=0$: no gluing is
needed. For even $k$, this obstruction can occur.

## Reproducible checks and limits

The [finite search](experiments/orthogonal_smith_search.py) and
[its exact results](experiments/orthogonal_smith_search_results.json)
cover 24 parameter cases and 148 degree components: B ranks 2--6 and
D ranks 4--6, complete selected ranges through rank 4 and polynomial
degree at most 12 in ranks 5--6. All computed factors are powers of two,
and every one-head quadric prediction agrees.

The [independent Frobenius check](experiments/type_b_frobenius_verification.json)
covers all 78 recorded B components, including 35 stable-range checks
and 46 block-duality comparisons. It uses Schur Littlewood--Richardson
multiplication modulo two, not BGG coefficients.

The [sector experiment](experiments/orthogonal_sector_search.py) and
[its results](experiments/orthogonal_sector_search_results.json) cover 70 D
components in ranks 4--6, and record seven degrees with nontrivial gluing.
These are finite corroborations, not proofs of a partition classification.

The [gluing-formula check](experiments/orthogonal_gluing_flags.py) and
[its results](experiments/orthogonal_gluing_flags_results.json) recompute
all seven nontrivial examples, verify the 63 trivial cases from the sector
records, and compare the formula with exact Smith forms for 200 synthetic
overlattices. The synthetic models include positive sector exponents and
gluing dimension up to five. Four small flag-stabilizer patterns are also
checked by constructing all 224 unimodular lifts. All comparisons pass.

Reproduce the checks from `borel/` inside SageMath:

```bash
PYTHONPATH=src sage -python experiments/orthogonal_smith_search.py
PYTHONPATH=src sage -python experiments/verify_type_b_frobenius.py
PYTHONPATH=src sage -python experiments/orthogonal_sector_search.py
PYTHONPATH=src sage -python experiments/orthogonal_gluing_flags.py
```

The previously open type-B multiplicities are now covered by the
[general proof](TYPE_B_GENERAL_SMITH.md). External expert review and
priority assessment remain separate work.
For type D, the goal is to replace the
flag-intersection data by explicit statistics on partitions.
The established formulas constrain that search and rule out the naive union
of the two type D sector Smith lists.
