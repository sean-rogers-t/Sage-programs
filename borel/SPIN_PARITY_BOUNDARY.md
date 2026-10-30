# Odd-spin kernel stability from adjacent-rank boundary lattices

Status, 22 October 2026: the boundary criterion and obstruction-rank formula
below are proved in every rank. They replace the degree-by-degree derivative
test by a comparison of two adjacent image lattices. They do **not** yet give
a closed numerical classification of the ranks n=0,3 modulo 4, or a closed
expression in n for the first failure degree. Those remain unresolved here.

## Setup

Use G=Spin(2n+1,C), P=B, R=Z_(2), and the full simply connected weight
lattice X. For S=Sym_R(X_R) let K be the actual integral Borel kernel.
The central-parity derivative on S/2S has D(z)=1, D(epsilon_i)=0, with
z=omega_n=(sum epsilon_i)/2. Stability means D(bar K) is contained in
bar K, where bar K=(K+2S)/2S. The quotient in question is the actual
image modulo twice itself; reduction in ambient cohomology is insufficient.
All degrees below are polynomial degrees.

The [symmetric-source reduction](SPIN_PARITY_RECURSION.md#1-exact-setup-and-classical-input)
identifies the decision, including the first failure degree, with that for
E_n=R[z,c2,...,cn], c1=2z, under z -> p1 and ci -> 2pi. Write

\[
 A_n=\bigoplus_{\lambda\subseteq\{1,\ldots,n\}}R P_\lambda,
 \quad M_n=\bigoplus_\lambda 2^{\ell(\lambda)}R P_\lambda,
 \quad N_n=M_n[p_1]\subset A_n.                                      \tag{1}
\]

A_n is the classical integral Schur-P/odd orthogonal Grassmannian ring.
Its presentation, Schubert basis and Pieri rule are published inputs; see
[Kresch--Tamvakis, introduction and Theorem 1 at q=0](https://math.umd.edu/~harryt/papers/qcogfinal.pdf).
The source under consideration here is the full Borel character source
after the stated free Artin-source reduction, rather than a newly chosen
representation or selected set of rational generators.

Let pi:A_(n+1) -> A_n kill strict partitions containing n+1. It restricts
to a surjection N_(n+1) -> N_n. Define the boundary lattice

\[
 L_{n,d}=\ker(N_{n+1,d}\longrightarrow N_{n,d}).                       \tag{2}
\]

On the ambient kernel of pi define a degree-minus-one linear map b_n by

\[
 b_n(P_{(n+1,\mu)})=
 \begin{cases}P_{(n,\mu)},&n\notin\mu,\\0,&n\in\mu.\end{cases}       \tag{3}
\]

Thus b lowers the boundary part once. It is not being asserted to be a
ring homomorphism or an operator preserving N_n.

## Boundary obstruction theorem

For every n and d, one has 2 b_n(L_(n,d)) contained in N_(n,d-1). The group

\[
 \mathcal O_{n,d}=
 \frac{b_n(L_{n,d})+N_{n,d-1}}{N_{n,d-1}}                            \tag{4}
\]

is therefore an elementary two-group. The central-parity derivative
preserves the integral full-flag kernel in degree d if and only if
O_(n,d)=0. In particular the entire kernel is preserved if and only if

\[
 b_n(L_{n,d})\subseteq N_{n,d-1}
 \quad\hbox{for every }n+1\leq d\leq T_n+1,
 \qquad T_n=n(n+1)/2.                                               \tag{5}
\]

Here the degree-d assertion does not assume stability in any earlier
degree. The exact first-failure formula is

\[
 f(n)=\min\{d:\mathcal O_{n,d}\ne0\},\qquad \min\varnothing=\infty.   \tag{6}
\]

Equation (6) is an effective boundary-lattice formula. It is **not** a
closed formula involving only elementary functions of n.

### Proof

Work first in the infinite Schur-P ring A_infinity, with its strict
partition basis and image subring N_infinity=R[p1,2p2,2p3,...]. The
integral derivation delta defined by

\[
 \delta(p_1)=1,\qquad \delta(p_i)=2p_{i-1}\quad(i>1)                  \tag{7}
\]

is well defined: for C(t)=1+2 sum pi*t^i it satisfies delta C(t)=2t C(t),
so it annihilates all defining coefficients of C(t)C(-t)-1. It preserves
N_infinity. On source generators delta(z)=1, delta(c2)=4z, and
delta(ci)=2c_(i-1) for i>2; consequently it reduces to D modulo twice
the actual image.

Its Schubert-basis formula lowers each part greater than one with
coefficient two, discarding repeated parts, and removes a part one with
coefficient one. For completeness, the Schur-Q generating expression is
the coefficient extraction from
product C(t_i) product_(i<j)(t_i-t_j)/(t_i+t_j), in the usual ordered
expansion. Differentiating multiplies this expression by 2 sum t_i.
Thus delta Q_lambda=2 sum Q_(lambda-e_i). Alternation kills a repeated
part. Passing from Q_lambda to 2^ell(lambda) P_lambda gives coefficient
two when length stays fixed and coefficient one when a final part one
is removed. This also follows directly from the two-row Pfaffian formula.

Let J_n be the span of partitions with largest part greater than n,
the kernel of A_infinity -> A_n. For g in N_infinity intersect J_n,
only parts exactly n+1 can survive after one application of delta and
truncation. Formula (7) in the Schubert basis gives

\[
             \pi_n\delta(g)=2 b_n(\pi_{n+1}g).                      \tag{8}
\]

Every element of L_n lifts to such a g: lift it from N_(n+1) to
N_infinity; its image in A_n is zero. Conversely every such g maps to
L_n. This explains why only adjacent ranks are needed. Since delta
preserves N_infinity, (8) also proves 2 b_n(L_n) is contained in N_n.

The map from the infinite polynomial source to N_infinity is surjective,
and its kernel is D-stable after reduction modulo two by (7). Every
finite-rank source relation lifts to a source polynomial whose infinite
image g is in J_n. Conversely every such g has a source representative,
and setting ci=0 for i>n gives a finite-rank source relation. This
specialization commutes with D modulo two: the infinite derivation sends
each killed ci to an even source element. Its derivative image is (8)
modulo 2N_n. Hence the degree-d derivative obstruction is exactly the
image of 2b_n(L_(n,d)) in N_(n,d-1)/2N_(n,d-1).

The map b+N -> 2b+2N identifies (4) with that obstruction. It is injective
because A_n is torsion-free: 2b in 2N implies b in N. This proves the
criterion in each degree. There are no boundary partitions in degrees
at most n, and the target degree d-1 vanishes for d>T_n+1. This proves
the finite range in (5), formula (6), and the full-flag assertions by
the unit Artin-source decomposition. QED.

## Explicit matrix and index formulas

The lattices in (2) need no polynomial ideal computation. In each degree
order the Schubert basis of A_(n+1) with partitions not containing n+1
first. Compute a lower triangular column Hermite basis from

\[
 H_d=\operatorname{Hermite}_{col}
 [\operatorname{diag}(2^{\ell(\lambda)})\mid U H_{d-1}],\quad H_0=(1).
                                                                        \tag{9}
\]

Write its block form as

\[
 H_d=\begin{pmatrix}F_d&0\\V_d&G_d\end{pmatrix}.                    \tag{10}
\]

Projection makes F_d a basis of N_(n,d); the last columns give the
entire boundary lattice L_(n,d), not just visible root generators in it.
Let W_d be the result of applying (3) to those last columns. Then

\[
 J_{n,d}=2F_{d-1}^{-1}W_d\ \text{is integral},\quad
 o(n,d)=\dim_{\mathbb F_2}\mathcal O_{n,d}
       =\operatorname{rank}_{\mathbb F_2}(J_{n,d}\bmod2).             \tag{11}
\]

The same number has a pure image-index expression. If H'_d is a lattice
basis of the columns [F_(d-1) | W_d], then

\[
 o(n,d)=v_2(\det F_{d-1})-v_2(\det H'_d).                           \tag{12}
\]

The index is a power of two because (4) has exponent two. These formulas
are uniform necessary-and-sufficient rank tests and give the first failure
by taking the smallest d with positive value. They do not require a
previously descended derivative matrix, unlike the earlier recurrence.

## A single-index criterion and a Hilbert formula

There is also a same-rank formulation that makes no reference to a
previously descended derivative. Set

\[
 C_n=M_n[p_1^2]\subseteq N_n=M_n[p_1],\qquad
 a_{n,d}=\#\{\lambda\subseteq\{1,\ldots,n\}:|\lambda|=d\},\qquad
 e_{n,d}=\log_2[N_{n,d}:C_{n,d}].                                   \tag{13}
\]

Both are actual embedded image subrings. Since 2p1 belongs to M_n,
2N_n is contained in C_n, so all these indices are powers of two and
the quotients have exponent two. With a_(n,d)=e_(n,d)=0 outside degrees
0,...,T_n, the exact obstruction formula is

\[
 \boxed{o(n,d)=a_{n,d-1}-e_{n,d}-e_{n,d-1}\ \geq0.}                  \tag{14}
\]

Consequently the following are equivalent:

1. D preserves the entire integral Borel kernel.
2. Every boundary obstruction in (4) vanishes.
3. The graded index polynomial is
   E_n(q)=sum e_(n,d) q^d = q product_(i=1)^n(1+q^i)/(1+q).
4. The single total image index satisfies
   log_2[N_n:C_n]=2^(n-1).

Here the total index refers to the underlying finite-rank R-lattices,
summed over degrees. Formula (14) gives

\[
 f(n)=\min\{d:e_{n,d}+e_{n,d-1}<a_{n,d-1}\},\qquad
 \sum_d o(n,d)=2^n-2\log_2[N_n:C_n].                                \tag{15}
\]

Equivalently, f(n) is the first degree d at which e_(n,d) differs from
the explicit partition count
# {lambda subset {2,...,n}: |lambda|=d-1}. At that first degree the actual
index exponent is strictly smaller than this count. The resulting
closed partition polynomial is the **stability target**, not an
unconditionally proved value of the image indices.

These are exact rank and first-failure criteria, rather than a guessed
congruence pattern. A closed evaluation of the actual indices in (13)
as a function of n is still required for a closed numerical answer.

**Proof.** In the F2-algebra A=N_n/2N_n, the subalgebra
B=C_n/2N_n has dimension b_d=a_(n,d)-e_(n,d) in degree d. Every source
polynomial modulo two splits into an even part and z times an even part.
The images of its even parts are exactly B: the extra root generator
2p1 in C_n maps to zero in A. Thus B_d direct-sum B_(d-1) maps
surjectively to A_d by (b,c) -> b+p1*c. Its kernel projects injectively
to the second summand, since B is embedded in A. That projection is
exactly the image of the derivative of source relations. Indeed the
derivative of an even part is zero and that of z times an even part
is the even part. Therefore
o(n,d)=b_d+b_(d-1)-a_d, which is (14), and this is the same actual
obstruction identified by the boundary theorem. Vanishing in all
degrees is equivalent to (1+q)E_n(q)=q product(1+q^i). Summing (14)
gives (15); its nonnegative summands show that the total-index equality
is also sufficient. QED.

The algebraic proof of (14) applies more generally to any finite free
graded image of a polynomial source R[z,x1,...,xs], with deg z=1 and
D=partial/partial z modulo two. Take M to be the image of
R[2z,x1,...,xs] and C=M[z^2]; replace a_(n,d) by the image rank in degree
d. The special odd-spin content is the Schur-P realization of those
image lattices and the explicit Pieri computation of their indices.
The dimension argument itself is elementary linear algebra, not a
special new property of spin groups.

The indices themselves have a two-step Pieri computation:

\[
 C_{n,d}=M_{n,d}+U^2C_{n,d-2},\qquad
 e_{n,d}=v_2(\det C_{n,d})-v_2(\det N_{n,d}).                       \tag{16}
\]

This provides a separate numerical check of (11), including after a
first failure. In particular the total dimension of all obstruction
groups must be even. Neither a rational Hilbert series nor the ambient
mod-two cohomology gives the indices in (13).

In this Grassmannian ring p1^2=p2, so C_n=M_n[p2]. For n>=2 this
strictly enlarges M_n already in degree two: M_(n,2)=2R P_(2), whereas
p2=P_(2) lies in C_n. Thus these factors compare the full spin image
with this specified intermediate image. They do not by themselves
give the Smith factors of N_n/M_n. In every derivative-stable rank,
N_n/C_n has only factors two, with multiplicity polynomial
q product_(i=2)^n(1+q^i). This applies in particular to every rank
n=1,2 modulo four by the earlier uniform theorem.

## Scope and next gap

The infinite derivation, Schur-Q identities, and Grassmannian presentation
are classical inputs. Equations (2)--(16) are deductions for the actual
central-parity kernel from those inputs. No literature-priority claim is
made. Results about integral Weyl invariants in H*(BSpin) concern a
different characteristic map; their invariant generators alone do not
settle preservation of this Borel-to-flag kernel. For example, the maps in
[Baek, Integral Weyl Invariants, introduction](https://arxiv.org/html/2607.18188v1)
and [the Spin(10) characteristic-image paper, introduction](https://arxiv.org/html/2607.23729v1)
start at the classifying space BG and restrict to BT. They are not the
map from the torus character algebra to the flag variety considered here.

The earlier proof already settles every n=1,2 modulo four. In the remaining
classes, exact finite calculations give stable ranks 4 and 11, and failures
at (n,d)=(3,4),(7,11),(8,16),(12,34),(15,37). The new rank-15 case means
Spin(31)/B, with first failure in cohomological degree 74. A closed
evaluation of (11)--(16) for arbitrary n=0,3 modulo four remains missing.
No proposed numerical
pattern for those ranks is promoted to a theorem or conjecture here.

### Independent checks

The [boundary computation](experiments/spin_parity_boundary.py) uses
FLINT integer Hermite bases and saves the entire adjacent-rank lattices,
not just their indices. It gives [ranks 3,4,7,8,11,12](experiments/spin-parity-boundary.json)
and the [rank-15 first failure](experiments/spin-parity-boundary-rank15.json),
with separate compressed certificate archives.

The [independent auditor](experiments/audit_spin_parity_boundary.py) uses
Python exact division and bit arithmetic, without a Hermite routine,
FLINT or Sage, to prove equality of the supplied lattices with the
generating lattices. Explicitly, the latter contain a diagonal
power-of-two sublattice. Exact inclusion followed by full row rank
modulo two of the inclusion coordinates therefore proves equality:
the possible lattice quotient is both two-primary and of odd order.
The auditor also replays the original derivative-factorization test on
the projected rank-n image bases and verifies the same obstruction rank
in every tested degree. This includes all degrees below 37 in rank 15.
It records concrete Schubert boundary witnesses for every first failure.
Reports: [initial cases](experiments/spin-parity-boundary-audit.json) and
[rank 15](experiments/spin-parity-boundary-rank15-audit.json).

Separately, 94 Schur-P lowering identities are checked by Pfaffian
expansion and monic elimination of the even special classes. This
corroborates the infinite generating-series proof without using its
coefficient extraction or the Pieri matrices. All checks share the
classical Grassmannian ring as a mathematical input; they are not
independent geometric constructions of that ring.

The next mathematical step is to describe the boundary lattice G_d in
(10) uniformly, including its off-diagonal dependence on V_d. Replacing
it by the diagonal root-character boundary misses cancellations with
z-multiples and cannot justify a numerical first-failure formula.

The [separate even-source calculation](experiments/check_spin_parity_even_source.py)
builds C and N using (16) in rank n, without reading the adjacent-rank
lattice certificates or applying the boundary map. It checks the same
first failures and every previously tested obstruction rank, and computes
the complete index polynomial for the six smaller ranks. It shares FLINT
arithmetic with the boundary computation; the Python lattice/derivative
audit described above supplies the separate arithmetic check. Its
[report](experiments/spin-parity-even-source-check.json) gives:

| n | log2[N_n:C_n] | Stability target 2^(n-1) | First failure |
| --- | --- | --- | --- |
| 3 | 2 | 4 | 4 |
| 4 | 8 | 8 | None |
| 7 | 48 | 64 | 11 |
| 8 | 112 | 128 | 16 |
| 11 | 1024 | 1024 | None |
| 12 | 1888 | 2048 | 34 |
| 15 | Not computed in all degrees | 16384 | 37 |

The rank-15 even-source check stops at degree 37 and does not assign a
total index. In that degree a_(15,36)=308, e_(15,36)=148 and
e_(15,37)=159, so (14) gives 308-148-159=1. Every earlier obstruction
vanishes. This supplies a separate exact check of the first-failure
degree, while the general-rank numerical problem remains unresolved.

## Reproduction

From borel, with python-flint and SymPy installed:

```sh
python experiments/spin_parity_boundary.py
python experiments/audit_spin_parity_boundary.py
python experiments/spin_parity_boundary.py --ranks 15 --stem spin-parity-boundary-rank15
python experiments/audit_spin_parity_boundary.py --stem spin-parity-boundary-rank15
python experiments/check_spin_parity_even_source.py
```
