# A closed formula for type-B unit factors when one block has size three

This note determines the multiplicity of the Smith factor **1**, uniformly
when one block has size three. Together with two other integral prefix
calculations, it now gives the complete multiplicities of 1, 2, 4, and 8
in [TYPE_B_THREE_ROW_SMITH.md](TYPE_B_THREE_ROW_SMITH.md). The complete two-row answer is in
[TYPE_B_TWO_ROW_SMITH.md](TYPE_B_TWO_ROW_SMITH.md), and the general setup is
in [TYPE_B_SMITH.md](TYPE_B_SMITH.md). No novelty claim is made.

## Formula

By integral block duality, put the block sizes in the form \((A,3)\), with
\(A\geq1\). Use the normalized elementary squared-tail inputs. Write
\(d=2m\) for polynomial degree, so cohomological degree is \(4m\).
Let \(u_{A,m}\) be the number of nonzero Smith factors equal to 1.

Define \({N\brack r}_q\) to be the Gaussian binomial polynomial; its
coefficient of \(q^m\) counts partitions of \(m\) in the
\(r\)-by-\((N-r)\) rectangle. Then

\[
\boxed{
\sum_{m\geq0}u_{2b,m}q^m
={b+3\brack3}_q
+q^{b+1}\left({b+2\brack2}_q-q^{2b}\right)
-q^{2b+1}(1+q+\cdots+q^{b-1})
}\tag{1}
\]

for even \(A=2b\), and

\[
\boxed{
\sum_{m\geq0}u_{2b+1,m}q^m
={b+3\brack3}_q+q^{b+1}{b+2\brack2}_q
}\tag{2}
\]

for odd \(A=2b+1\). Formula (2) includes \(A=1\), where the answer is
\(1+q\). For example,

\[
\begin{aligned}
U_3(q)&=1+q+2q^2+2q^3+q^4,\\
U_4(q)&=1+q+2q^2+3q^3+3q^4+2q^5+q^6,\\
U_5(q)&=1+q+2q^2+3q^3+3q^4+3q^5+2q^6+q^7.
\end{aligned}
\]

All odd polynomial-degree images are zero. The rank of the selected image
in degree \(2m\), including its nonunit factors, is still
\([q^m]{A+3\brack3}_q\).

## 1. Frobenius and the two-quotient reduction

The proved Frobenius formula in [TYPE_B_SMITH.md](TYPE_B_SMITH.md) says that
\(u_{A,m}\) is the rank over \(\mathbb F_2\) of squaring

\[
H^{2m}(\operatorname{Gr}(3,A+3);\mathbb F_2)
\longrightarrow H^{4m}(\operatorname{Gr}(3,A+3);\mathbb F_2).
\]

On symmetric functions, squaring modulo 2 is the second Adams plethysm:
\(s_\lambda(y)^2=s_\lambda(y^2)\). The SXP rule gives its integral
Schur coefficient at \(\mu\) as zero unless \(\mu\) has empty 2-core;
otherwise that coefficient is

\[
\pm c^{\lambda}_{\alpha,\beta},
\]

where \((\alpha,\beta)\) is the 2-quotient of \(\mu\).
See [Wildon, *A generalized SXP rule proved by bijections and involutions*,
equation (1)](https://arxiv.org/abs/1508.07030). Row signs do not affect
the rank, including after reduction modulo 2.

Here the target partitions lie in \(A^3\). Pad \(\mu\) to three parts
and use its beta numbers \(\mu_i+3-i\), three distinct integers in
\([0,A+2]\). Empty 2-core means that two beta numbers are even and one
is odd, as for the empty partition's beta numbers \((2,1,0)\). Dividing
each runner by 2 shows that the allowable quotients are precisely

\[
\begin{array}{c|c|c}
A &\alpha &\beta\\ \hline
2b &\alpha\subseteq b^2 &\beta=(t),\quad0\leq t\leq b\\
2b+1 &\alpha\subseteq b^2 &\beta=(t),\quad0\leq t\leq b+1.
\end{array}
\]

Consequently the transpose of the nonzero Frobenius matrix is the Schur
coefficient matrix of

\[
s_{(r,s)}h_t,\qquad 0\leq s\leq r\leq b,
\quad0\leq t\leq c,\quad r+s+t=m,
\]

where \(c=b\) or \(b+1\). All Schur terms of these products have at most
three rows and first part at most \(b+c=A\), so the source rectangle
discards no terms. We can compute their span inside the full symmetric
function ring.

## 2. A graph computes the rank in every characteristic

The two-row Jacobi--Trudi identity is

\[
s_{(r,s)}h_t=h_rh_sh_t-h_{r+1}h_{s-1}h_t.
\tag{3}
\]

The products \(h_\gamma\), indexed by ordinary partitions \(\gamma\),
are an integral basis of the symmetric function ring. Regard each
\(h\)-monomial occurring in (3) as a graph vertex. A column with \(s>0\)
is the oriented incidence vector of an edge; when \(s=0\), it is a single
vertex vector, which we call an anchor. Zero parts are ignored, with
\(h_0=1\) and \(h_{-1}=0\).

A connected graph on \(v\) vertices has incidence rank \(v-1\) over every
field. If an anchor occurs, adjoining it gives rank \(v\). Both statements
also hold integrally with primitive image: a spanning tree supplies a basis
for the sum-zero lattice, and an anchor supplies the remaining coordinate.

In particular, this argument proves that the \(\nu=0\) prefix of the
integral type-C Schur matrix has primitive image for this three-row family.
The SXP correspondence uses only row signs and transpose; the Schur-to-\(h\)
basis change is integral unimodular. Embedding the rectangle Schur space in
the full symmetric-function group also has primitive image. Thus the graph
argument controls its nonzero integer Smith factors, as well as its rank
over \(\mathbb F_2\).

The internal vertices, with all three parts at most \(b\), form a connected
graph in each degree in which they exist. Indeed, transfers of one box from
a smaller positive part to a larger part connect every partition to the
lexicographically largest partition with that weight and cap \(b\).
Perform transfers to the first part until it reaches \(b\), and then to
the second part. Each transfer is an edge of (3). This gives the required
connectivity without a characteristic assumption.

## 3. Even width

For \(A=2b\), the full vertex set in degree \(m\) consists of:

1. partitions in \(b^3\);
2. partitions \((b+1,r,s)\) with \(0\leq s\leq r\leq b\) and
   \(s\leq b-1\).

There is no vertex with two parts \(b+1\). Every second-kind vertex is
joined to the internal vertex \((b,s+1,r)\) by an edge, so the graph is
connected in every degree \(0\leq m\leq3b\). The internal vertices are
nonempty throughout that interval. The anchors have weights \(r+t\),
\(0\leq r,t\leq b\); hence an anchor exists exactly when \(m\leq2b\).

The two vertex-count generating functions are

\[
{b+3\brack3}_q,\qquad
q^{b+1}\left({b+2\brack2}_q-q^{2b}\right).
\]

The removed monomial represents the excluded pair \((b,b)\). Subtract
one unanchored component in each degree \(2b+1,\ldots,3b\). This proves
(1).

## 4. Odd width

For \(A=2b+1\), the vertex set consists of:

1. partitions in \(b^3\);
2. partitions \((b+1,r,s)\), \(0\leq s\leq r\leq b\);
3. partitions \((b+1,b+1,j)\), \(0\leq j\leq b-1\).

Every second-kind vertex except \((b+1,b,b)\) attaches to an internal
vertex as above. Each third-kind vertex attaches to
\((b+1,b,j+1)\). At the final weight \(3b+1\), the only vertices are
\((b+1,b,b)\) and \((b+1,b+1,b-1)\), joined by an edge. Thus the graph
is connected in every degree with vertices. For \(b=0\), the two vertices
in different degrees are simply \(1\) and \(h_1\), both anchored.

The anchors now have all weights from 0 to \(2b+1\). There is therefore
one unanchored component in every degree \(2b+2,\ldots,3b+1\).
These degrees are exactly the weights of the third-kind vertices, one per
degree. Their vertex count cancels the component subtraction. The remaining
two vertex counts give (2).

## Verification and limits

The independent verification script
[verify_type_b_three_row_units.py](experiments/verify_type_b_three_row_units.py)
uses Schur Littlewood--Richardson multiplication over \(\mathbb F_2\),
and also checks the formulas against recorded type-B Smith calculations.
Its [report](experiments/type_b_three_row_unit_verification.json) records
the input hashes and exact scope. This theorem supplies the unit layer;
the nonunit layers are proved in
[TYPE_B_THREE_ROW_SMITH.md](TYPE_B_THREE_ROW_SMITH.md).
