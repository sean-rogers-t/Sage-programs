# Type B multiplicities and a Schur-product filtration

**Audited deduction (2026-10-19).** This is a repository derivation from
cited published inputs, not a published Smith classification. Novelty is
unestablished. See [PROOF_STATUS.md](PROOF_STATUS.md) for the audited scope
and dependencies.

The weighted-minor formula in [TYPE_B_SMITH.md](TYPE_B_SMITH.md) can be
reduced to a filtration of ordinary Grassmannian Schur spaces. This note
proves that reduction. It also identifies the remaining integral issue in
using its ranks as an unrestricted multiplicity formula.

## Conventions

By integral block duality put the smaller block at the head:

$$
l=\min(a,k),\quad A=\max(a,k),\quad n=A+l,\quad
X_C=\operatorname{IG}(l,2n).
$$

The source in polynomial degree $2m$ has the integral basis
$s_\lambda(y_1^2,\ldots,y_l^2)$, $\lambda\subseteq A^l$, $|\lambda|=m$.
It gives the same selected lattice as the tail generators, by inversion
of the two squared elementary series and Schur transposition.

Let $C$ be the primitive type C Schubert coefficient matrix in this basis,
and $B_{w,\lambda}=2^{\nu(w)}C_{w,\lambda}$ the type B matrix. Define

$$
P_t=\operatorname{span}_{\mathbb Z}
\{C_{w,\bullet}:\nu(w)\leq t\}\subseteq\mathbb Z^r,
\qquad P_{-1}=0.
$$

Here $r$ is the number of rectangle partitions of $m$. The full row
lattice is $P_l=\mathbb Z^r$.

Identify $\mathbb Z^r$ with the degree-$m$ part of the ordinary integral
Grassmannian ring

$$
R_{A,l}=\mathbb Z[e_1,\ldots,e_l]/(h_{A+1},\ldots,h_{A+l})
$$

by its Schur basis. The identification describes coefficient rows, not a
claim that the characteristic map itself changes cohomological degree.

The general Schur-product equality below also works with the actual head
size $l=H$ and complementary tail size $A=K$, without requiring $H\leq K$.
The smaller-head convention is chosen here using duality of the complete
Smith multiset. It does not assert that individual $\nu$-prefixes are
preserved by interchanging the blocks. The recorded-matrix verifier uses
the actual head size $a$ and conjugates the stored tail Schur column
partitions to obtain head Schur coordinates.

## A general description of every prefix

For $0\leq t\leq l$, put

$$
b=l-t,\quad b_0=\lceil b/2\rceil,\quad b_1=\lfloor b/2\rfloor,\qquad
A_0=\lceil n/2\rceil-b_0,\quad A_1=\lfloor n/2\rfloor-b_1.
$$

Then the following equality is integral:

$$
\boxed{P_t=
\operatorname{span}_{\mathbb Z}\{
s_\mu s_\alpha s_\beta\text{ in }R_{A,l}:
\mu\subseteq A^t,\quad
\alpha\subseteq A_0^{b_0},\
\beta\subseteq A_1^{b_1},\
|\mu|+|\alpha|+|\beta|=m\}.}
$$

An empty-height rectangle supplies only the empty partition.
Thus the general problem reduces to specified Schur-product spans in a
finite ordinary Grassmannian ring. This is not yet a rank-free partition
count for every $t$.

### Proof by incidence resolution

Fix a complex Lagrangian $L_0$. On a Schubert cell indexed by $w$,
the dimension of $\Sigma\cap L_0$ is $l-\nu(w)$: among the $l$ pivot
coordinates, precisely the negative signed entries lie beyond $L_0$.
Consequently the cells with $\nu(w)\leq t$ form the closed Schubert union

$$
Z_t=\{\Sigma:\dim(\Sigma\cap L_0)\geq b\}.
$$

Its incidence resolution chooses $S\subseteq\Sigma\cap L_0$ of dimension
$b$. It is a bundle over $\operatorname{Gr}(b,n)$ with fiber
$\operatorname{IG}(t,2(n-b))$. To give a section over a Schubert cell
explicitly, use its standard root-group coordinates
$z\longmapsto g(z)\Sigma_w$, with $g(z)$ in the Borel preserving $L_0$
and $\Sigma_w$ the coordinate isotropic plane. Choose a fixed
$b$-plane $S_w\subseteq\Sigma_w\cap L_0$, using any $b$ positive
coordinate pivots. Then $g(z)S_w\subseteq g(z)\Sigma_w\cap L_0$ varies
regularly and gives the required section.
The closure of its graph pushes a cycle forward with degree one onto
the corresponding Schubert closure. Homology pushforward is therefore
surjective. Since the Schubert union has free even homology, universal
coefficients show that pullback to the resolution is a primitive
injection on integral cohomology.

Choose a compatible quaternionic structure. The complex subspace $S$
spans a quaternionic rank-$b$ bundle $S\oplus JS$. Forgetting the complex
Lagrangian in the remaining quaternionic rank-$t$ tautological bundle
maps the resolution onto a quaternionic Grassmannian bundle $Q$ over
$\operatorname{Gr}(b,n)$. Its fiber is $\operatorname{LG}(t,2t)$.
The topological construction on each fiber is
[Pragacz--Ratajski, Appendix A, Section 10.5](https://www.impan.pl/~pragacz/download/comp3868.pdf);
applying it to the symplectic complement of $S\oplus JS$ gives this
bundle version. Their Section 10.3 gives the integral quaternionic Schur basis.
The Grassmannian bundle $Q$ is simply connected and has free even
cohomology, as do its base and quaternionic Grassmannian fiber.
Thus there is no local-coefficient issue for the Lagrangian fiber.
The integral even Serre spectral sequence shows that this pullback is
primitive as well. The incidence resolution itself has free even
cohomology by its isotropic-Grassmannian bundle over the ordinary base.

The classes $s_\mu(U)$, $\mu\subseteq A^t$, form a global integral
Leray--Hirsch basis for $Q$ over the complex Grassmannian: they restrict
to the integral Schur basis on every quaternionic Grassmannian fiber.
For the characteristic-class formulation and its relative version, see
[Panin--Walter, Theorems 11.2--11.3](https://arxiv.org/pdf/1011.0649).
If $x$ are the Chern roots of $S^*$, Schur branching gives the source
pullback

$$
s_\lambda(x^2+U)
=\sum_{\mu\subseteq\lambda}s_{\lambda/\mu}(x^2)s_\mu(U).
$$

All the relevant $\mu$ fit $A^t$; no fiber-basis reduction is needed.
The ordinary base has Schur basis indexed by
$\gamma\subseteq(n-b)^b$. The Littlewood two-quotient formula says
that the coefficient of $s_\gamma$ in $\psi_2(s_{\lambda/\mu})$
vanishes unless $\gamma$ has empty 2-core; in that case it is

$$
\epsilon_2(\gamma)\,
\operatorname{coeff}_{s_\lambda}(s_\mu s_{\gamma^{(0)}}s_{\gamma^{(1)}}).
$$

Padding $\gamma$ to $b$ rows, its beta numbers lie between $0$ and
$n-1$. The even and odd runners have respectively $b_0,b_1$ entries.
Their quotients therefore range independently over exactly
$A_0^{b_0}$ and $A_1^{b_1}$. Row signs do not change an integer span.

Finally, both primitive pullbacks above have integral linear left inverses.
Hence calculating the source coefficient rows after these pullbacks
preserves their row lattice, rather than just their rational rank.
This proves the boxed equality.

The two-quotient identity and runner argument are also detailed in
[TYPE_B_THREE_ROW_UNITS.md](TYPE_B_THREE_ROW_UNITS.md). The Schubert,
theta and B/C conventions are those of
[Buch--Kresch--Tamvakis](https://math.umd.edu/~harryt/papers/isogiamfinal.pdf).

## When the filtration ranks determine the factors

If every $P_t$ is primitive in $\mathbb Z^r$, then

$$
\boxed{\#\{\text{Smith factors }2^t\}
=\operatorname{rank}P_t-\operatorname{rank}P_{t-1}.}
$$

Indeed, a chain of primitive submodules of a free abelian group admits a
common adapted integral basis: each successive quotient is free. The
row lattice of $B$ is

$$
\sum_{t=0}^l2^t
\operatorname{span}\{C_{w,\bullet}:\nu(w)=t\}
=\sum_{t=0}^l2^tP_t.
$$

In the adapted basis, a coordinate first appearing in $P_t$ contributes
exactly the factor $2^t$. Transposition preserves Smith factors.

The primitivity assumption is essential for this argument. It is proved
for all prefixes when $l\leq3$ in the accompanying closed-form notes,
and when $l=4$ by the appendix below together with the all-rank endpoint
results. Closed multiplicity counts are proved here only for $l\leq3$;
primitivity alone does not supply missing closed rank counts. It is proved
for $P_0$ in every rank in
[TYPE_B_UNIT_SMITH.md](TYPE_B_UNIT_SMITH.md), and for $P_{l-1}$ in every rank in
[TYPE_B_HIGHEST_SMITH.md](TYPE_B_HIGHEST_SMITH.md). Its unrestricted
validity is still a question here. Finite checks alone do not establish it.
The proposed closed rank increments and their geometric motivation are
recorded in [TYPE_B_CLOSED_CONJECTURE.md](TYPE_B_CLOSED_CONJECTURE.md).

## Appendix: two integral prefix statements when the smaller block is four

Let $l=4$, $A\geq4$, and retain the identification of coefficient rows
with $R_{A,4}$. The complete-symmetric products

$$
h_\gamma=h_{\gamma_1}h_{\gamma_2}h_{\gamma_3}h_{\gamma_4},
\qquad \gamma\subseteq A^4,
$$

form an integral basis. Their expansions in the rectangle Schur basis
are Kostka unitriangular; terms of first part greater than $A$ vanish
in the Grassmannian quotient. Also $h_j=0$ for every $j>A$, by the
tautological quotient-bundle relation. We use this complete-symmetric
basis only in this appendix; it is a unimodular change from the head
Schur coordinates.

### The first prefix is primitive

For $t=1$, the Schur-product theorem gives rows

$$
h_r s_{(u,v)}h_w,
\quad 0\leq r\leq A,\quad
0\leq v\leq u\leq\lceil A/2\rceil,
\quad 0\leq w\leq\lfloor A/2\rfloor+1.
$$

Jacobi--Trudi rewrites each such row as

$$
h_rh_uh_vh_w-h_rh_{u+1}h_{v-1}h_w.
$$

For $v=0$, the second term vanishes and the row is a coordinate vector.
Otherwise the row is the difference of two coordinate vectors, after
sorting the four indices and discarding any vanished $h_j$ factor.
Sorting cannot create a coefficient of absolute value two: a collision
of the indices would cancel the row, rather than add its two terms.

Make a graph whose vertices are the surviving complete-symmetric basis
vectors, with an edge for each signed difference and an anchor for
each single coordinate row. On a connected component, its row lattice
is the full coordinate lattice if there is an anchor, and the
sum-zero coordinate lattice otherwise. A spanning tree proves both
claims integrally. In either case it is primitive. The direct sum over
components proves that $P_1$ is primitive, without making any
connectivity or closed-rank assumption.

### The second prefix has a closed coordinate description

Put

$$
U=\left\lceil\frac{A+2}{2}\right\rceil,
\qquad V=\left\lfloor\frac{A+2}{2}\right\rfloor.
$$

Then the stronger exact statement is

$$
\boxed{P_2=
\bigoplus_{\substack{\gamma\subseteq(A,A,U,V)\\|\gamma|=m}}
\mathbb Z h_\gamma.}
$$

Indeed, the Schur-product rows are

$$
s_{(r,s)}h_uh_v
=h_rh_sh_uh_v-h_{r+1}h_{s-1}h_uh_v,
\quad 0\leq s\leq r\leq A,\quad
0\leq u\leq U,\quad0\leq v\leq V.
$$

Every surviving monomial has sorted index in $(A,A,U,V)$: at most
two factors can exceed $U$, and at least one factor has index at most
$V$. The first term of a row reaches every basis vector in that shape
by choosing $(r,s)$ to be its two largest parts and $(u,v)$ its two
smallest parts.

Fix such a choice. Keep $u,v$ fixed and repeatedly transfer one unit
from $s$ to $r$. Each transfer is exactly one of the displayed row
differences. If $s$ reaches zero, the terminal coordinate is itself an
anchor row. If $r$ reaches $A$ first, the next Jacobi--Trudi term
contains $h_{A+1}=0$, so that terminal coordinate is again an anchor.
Back-substitution obtains the original coordinate integrally. Every
step stays in the same bounded four-part shape. This proves the boxed
equality and, in particular, primitivity of $P_2$.

Thus

$$
\operatorname{rank}P_2(m)
=\#\{\gamma\subseteq(A,A,U,V):|\gamma|=m\}
$$

is a closed partition count. Together with the all-rank unit-prefix
theorem, the highest-prefix theorem for $P_3$, and $P_4=\mathbb Z^r$,
the two arguments above prove primitivity of every prefix for $l=4$.
The rank-increment rule consequently applies to this family. A closed
count for $P_1$ is still needed here to turn that fact into a complete
closed multiplicity formula; no such count is asserted in this
appendix.

### The multiplicity of the factor eight is closed as well

The exact $P_2$ description determines one more layer even without a
closed $P_1$ count. Set

$$
c=\left\lceil\frac{A+4}{2}\right\rceil,
\qquad r=\left\lfloor\frac{A+4}{2}\right\rfloor.
$$

Only the primitivity of $P_2$ is needed to count exponents at least
three; no assertion about $P_0$ or $P_1$ is needed here. Let $M_m$ be
the degree-$m$ source coefficient lattice and $R\subseteq M_m$ the
weighted type-B row lattice. Choose
$M_m=P_2\oplus Q_2$, which is possible by the coordinate
description above. The rows of weights $\nu\leq2$ lie in $P_2$, and
their lattice contains $4P_2$. Every row of weight $\nu\geq3$ is
divisible by eight in $M_m$. Its $P_2$ component therefore
lies in $8P_2\subseteq R$ and can be subtracted integrally. Consequently
$R$ splits as a lattice in $P_2$, with factors dividing four, and a
lattice in $Q_2$, with all factors divisible by eight. The latter has
full rank in $Q_2$, since $16M_m\subseteq R$.

Thus the number of exponents at least three in each degree is the rank
deficit of $P_2$. Its generating polynomial counts partitions in $A^4$
outside $(A,A,U,V)$. We claim this is

$$
q^{3c}{r\brack3}_q[c-3]
 +q^{4c}{r\brack4}_q,
$$

where $[v]=1+q+\cdots+q^{v-1}$ is a q-integer. The second summand is
exactly the all-rank highest-factor formula for the factor sixteen.
Consequently the multiplicity of the factor eight has the closed
generating polynomial

$$
\boxed{\sum_m M_3(m)q^m=q^{3c}{r\brack3}_q[c-3].}
$$

Here is a direct partition proof of the claimed deficit. If $A$ is
even, $r=c$ and $U=V=c-1$. The excluded partitions have third part at
least $c$. Those with fourth part less than $c$ contribute
$q^{3c}{r-1\brack3}_q[c]$; those with fourth part at least $c$ contribute
$q^{4c}{r\brack4}_q$. The q-integer identity

$$
[r]{r-1\brack3}_q={r\brack3}_q[r-3]
$$

gives the result.

If $A$ is odd, $r=c-1$, $U=c-1$, and $V=c-2$. Divide the complement
into partitions with third part at least $c$ and fourth part at most
$c-2$, and partitions with fourth part at least $c-1$. Their generating
polynomials are respectively

$$
q^{3c}{r-1\brack3}_q[r],
\qquad q^{4c-4}{r+1\brack4}_q.
$$

Apply the preceding q-integer identity, the alternate Gaussian Pascal
identity

$$
{r+1\brack4}_q={r\brack3}_q+q^4{r\brack4}_q,
$$

and $[r-2]=[r-3]+q^{r-3}$. Their sum becomes exactly the displayed
rank deficit, since $r=c-1$.

Thus the compact general conjecture is proved for the factor eight
when the smaller block is four, using the audited $P_2$ and
largest-factor deductions independently of the all-rank unit argument.
This does not separate the factor-two
and factor-four multiplicities in that family; their closed count
still requires the missing $P_1$ rank calculation.

## What the new search establishes

The [Schur-column experiment](experiments/type_b_multiplicity_search.py)
stores the exact C and B matrices, row weights, integer prefix Smith
forms, timings and cross-checks in
[its results](experiments/type_b_multiplicity_search_results.json).
It certifies primitivity and the rank-increment formula for each recorded
component, without assuming their general validity.

Two tempting stronger assertions fail: C is neither always totally
unimodular nor always multiplicity free. For example, a prefix in the
$l=3,A=4,m=6$ calculation has a minor

$$
\begin{pmatrix}1&0&0\\0&1&-1\\0&1&1\end{pmatrix},
$$

of determinant $2$. In the $l=A=4,m=6$ calculation a coefficient itself
is $2$. The complete prefix lattices in those examples are nevertheless
primitive. A proof of primitivity cannot rely on either false assertion.

The [Cartesian frontier](experiments/type_b_cartesian_frontier.py) is an
experiment-only acceleration, independently compared with fundamental-weight
divided differences and the earlier elementary-generator matrices. It keeps
orthonormal monomials sparse; the terminal C divided difference removes an
odd terminal power. It changes the computation cost, not the input lattice
or the definition of the Smith factors.
