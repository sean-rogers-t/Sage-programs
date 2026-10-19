# A closed partition candidate for all type-B Smith multiplicities

**Status: general conjecture, with audited special-case deductions.**
See [PROOF_STATUS.md](PROOF_STATUS.md) for the published inputs, checked
arguments, and literature limits. Novelty has not been established.

This note proposes a positive partition formula for every normalized type-B
multiplicity. The endpoint cases and all three-row multiplicities have
derivations audited in the companion notes. The unrestricted assertion for
intermediate factors remains a conjecture; finite evidence does not prove it.
The polynomial $U_{n,l}$ is a classical Birkhoff/Hall counting polynomial.
The additional Smith interpretation requires the arguments cited below.

## The formula

Use the generators and integral lattice of [TYPE_B_SMITH.md](TYPE_B_SMITH.md).
Set

$$
n=a+k,\quad l=\min(a,k),\quad A=\max(a,k),\quad
r=\lfloor n/2\rfloor,\quad c=\lceil n/2\rceil.
$$

In polynomial degree $2m$, equivalently cohomological degree $4m$, let
$M_{n,l,e}(m)$ be the multiplicity of the nonzero Smith factor $2^e$.
Only $0\le e\le l$ can occur. Write
$E_{n,l,e}(q)=\sum_m M_{n,l,e}(m)q^m$.
Here $q$ is a formal counting variable.
The Gaussian coefficient ${v\brack u}_q$ counts partitions in $(v-u)^u$;
it is zero unless $0\le u\le v$. Define

$$
U_{n,l}(q)=\sum_{j=0}^{\lfloor l/2\rfloor}
q^{j(c-l+j)}{r\brack j}_q{c-j\brack l-2j}_q.
\tag{1}
$$

The proposed complete answer is

$$
\boxed{E_{n,l,e}(q)=q^{ec}{r\brack e}_q U_{n-2e,l-e}(q).}
\tag{2}
$$

The parameters inside $U$ are recomputed for its arguments. Equivalently,

$$
\boxed{
E_{n,l,e}(q)=
\sum_{j=0}^{\lfloor(l-e)/2\rfloor}
q^{ec+j(c-l+j)}
{r\brack e}_q{r-e\brack j}_q
{c-e-j\brack l-e-2j}_q.
}
\tag{3}
$$

Thus the candidate multiplicity counts tuples $(j,\lambda,\alpha,\beta)$
with

$$
\lambda\subseteq(r-e)^e,\qquad
\alpha\subseteq(r-e-j)^j,\qquad
\beta\subseteq(c-l+j)^{\,l-e-2j},
$$

and

$$
|\lambda|+|\alpha|+|\beta|+ec+j(c-l+j)=m.
$$

This requires partition counting, without Schubert matrices, Smith
reduction, or a rank calculation.

## Audited deductions and the remaining conjecture

* **$e=0$, all ranks:** (1) is the unit multiplicity polynomial.
  [TYPE_B_UNIT_SMITH.md](TYPE_B_UNIT_SMITH.md) proves this integrally
  through an affine Grassmannian model of an invariant-subspace variety.
* **$e=l$, all ranks:** (3) becomes $q^{lc}{r\brack l}_q$, the
  [largest-factor theorem](TYPE_B_HIGHEST_SMITH.md).
* **Every $e$ when $l\le3$:** the
  [two-row](TYPE_B_TWO_ROW_SMITH.md) and
  [three-row](TYPE_B_THREE_ROW_SMITH.md) notes prove these multiplicities
  by integral Schur straightening and primitive filtrations.
* **$l=4,e=3$:** the appendix of
  [TYPE_B_FILTRATION.md](TYPE_B_FILTRATION.md) gives the closed
  polynomial $q^{3c}{r\brack3}_q[c-3]$. All four-row prefixes are
  primitive; a closed first-prefix rank is still needed to separate
  the factor-two and factor-four multiplicities.
* **General integral reduction:** [TYPE_B_FILTRATION.md](TYPE_B_FILTRATION.md)
  identifies every weighted Schubert prefix with a specified
  Schur-product row lattice. Its rank increments give Smith
  multiplicities if all these row lattices are primitive.

The bullets describe deductions in these notes, not exact Smith theorems
located in the cited literature. The remaining claim in (2)--(3) is the
unrestricted interpretation of the intermediate layers as Smith
multiplicities. The geometry below proves
that these polynomials count natural strata. An integral comparison with
the selected Springer image still needs proof. No novelty claim is made.

## Why these polynomials arise geometrically

This is part of the existing subspace-profile counting framework:
for $N^2=0$, defect $e$ is profile $(l,e)$ in
[Ram, Definition 1.1 and Theorem 5.2](https://www.mat.univie.ac.at/~slc/wpapers/FPSAC2025/115.pdf).
That source gives a general finite-field count. The elementary Gaussian
specialization below and its compact-support interpretation are explained
here; neither is a proof of the additional integral Smith identification.

Let $V$ have complex dimension $n$ and choose a square-zero endomorphism
$N$ of maximal rank $r$. Put $I=\operatorname{im}N$ and $K=\ker N$,
of dimensions $r,c$. Choose $V=I'\oplus K$ with
$N:I'\longrightarrow I$ an isomorphism.
For $U\in\operatorname{Gr}(l,V)$ define its defect

$$
\delta(U)=\operatorname{rank}(N:U\longrightarrow V/U),
$$

and let $S_e$ be the locally closed stratum $\delta(U)=e$.
Independently of the Smith conjecture, its compactly supported even
Poincare polynomial is (3), and its integral compactly supported
cohomology is free with no odd part.

Indeed, put $D=N(U)$, $B=U\cap K$, and $j=\dim(D\cap B)$.
Then $\dim D=e+j$ and $\dim B=l-e-j$. Choose successively:

1. $D\in\operatorname{Gr}(e+j,r)$;
2. $C=D\cap B\in\operatorname{Gr}(j,D)$;
3. the image of $B/C$ in $K/D$, in
   $\operatorname{Gr}(l-e-2j,c-e-j)$.

The lifts of the third choice to $B$ form an affine space of dimension
$e(l-e-2j)$. With $D,B$ fixed, the lifts to $U$ form an affine space
of dimension $(e+j)(c-l+e+j)$. Their total dimension is
$ec+j(c-l+j)$.
These are locally trivial affine bundles over the indicated
Grassmannian bundles. Refining the bases by Schubert cells gives affine
pieces; filtering by intersection dimension supplies a closed filtration.
The contributions are therefore free and even. Their polynomial is

$$
\sum_j q^{ec+j(c-l+j)}
{r\brack e+j}_q{e+j\brack j}_q
{c-e-j\brack l-e-2j}_q.
$$

The flag identity ${r\brack e+j}_q{e+j\brack j}_q
={r\brack e}_q{r-e\brack j}_q$ proves (3).
For $e=0$, this is the projective variety of $N$-invariant
$l$-planes in the unit-factor proof.

In particular the geometric polynomials are nonnegative, and the closed
defect filtration of the ordinary Grassmannian gives

$$
\sum_{e=0}^l E^{\mathrm{geom}}_{n,l,e}(q)={n\brack l}_q.
\tag{4}
$$

This agrees with the known total selected-image rank.
At $q=1$ their total-degree counts have the simple identity

$$
\sum_{l,e}E^{\mathrm{geom}}_{n,l,e}(1)x^lt^e
=(1+x)^{n\bmod2}\bigl(1+(1+t)x+x^2\bigr)^r.
\tag{5}
$$

Here the geometric construction extends to $0\le l\le n$; terms with
impossible Grassmannian dimensions are zero. In that extension take
$0\le e\le\min(l,n-l,r)$ and sum over
$\max(0,l-c)\le j\le\min(r-e,\lfloor(l-e)/2\rfloor)$, so no negative
power of $q$ is used. For $n=8,l=4$ the candidate
counts over all degrees are $19,28,18,4,1$ for the factors
$1,2,4,8,16$. Their sum is ${8\choose4}=70$.

## The exact missing comparison

For the integral prefix row lattices $P_t$ in
[TYPE_B_FILTRATION.md](TYPE_B_FILTRATION.md), the desired bridge is

$$
P_t\text{ primitive},\qquad
\sum_m\operatorname{rank}(P_t)_m q^m
=\sum_{e=0}^t E^{\mathrm{geom}}_{n,l,e}(q).
\tag{6}
$$

Consider the closed locus $W_t=\{U:\delta(U)\le t\}$.
An incidence construction chooses an $N$-invariant $(l-t)$-plane inside
$U$, together with the balanced auxiliary subspace used to resolve that
invariant Grassmannian. Its iterated Grassmannian bundles reproduce
the Schur-product rows for $P_t$.

However, Chern classes generating the resolution cohomology do not alone
prove that the combined tautological Chern classes generate
$H^*(W_t;\mathbb Z)$. Primitivity of the resolution pullback also needs
an integral argument. For $t=0$ the affine Grassmannian proof supplies
both. For $t>0$, $U$ need not be $N$-invariant, so its inverse image is
not a lattice over $\mathbb C[[z]]$. The same proof cannot simply be
reused. A rational comparison would leave the Smith problem unresolved.

## Reproducible evidence

The [candidate verifier](experiments/verify_type_b_closed_candidate.py)
compares the closed polynomials with:

* **174** recorded exact Smith lists;
* **430** recorded prefix ranks, including **181** larger independent
  Littlewood--Richardson prefix computations;
* **108** parameter cases for polynomial consistency.

There are no counterexamples in
[the verification report](experiments/type_b_closed_candidate_verification.json).
The verifier imports no Borel or Weyl engine. The stored Smith data use
the library calculation and the independently checked
[Cartesian frontier](experiments/type_b_cartesian_frontier.py).
The [filtration verifier](experiments/verify_type_b_filtration.py)
separately checks **516 exact integer row-module equalities**.

From `borel/`:

```bash
PYTHONPATH=src sage -python experiments/verify_type_b_closed_candidate.py
PYTHONPATH=src sage -python experiments/verify_type_b_filtration.py
PYTHONPATH=src sage -python experiments/verify_type_b_partition_formulas.py
```

For the literal defining/vector Springer integral form described in
[INTEGRAL_LITERATURE.md](INTEGRAL_LITERATURE.md), a degree-$2m$ normalized
factor $2^e$ is multiplied by $2^{2m}$, becoming $2^{2m+e}$.
