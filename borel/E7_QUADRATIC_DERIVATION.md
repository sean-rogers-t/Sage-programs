# A structural route to the E7 112-case formula

Status, 21 October 2026: a general quadratic-derivation criterion, a
noncircular bounded E7 proof, and an arbitrary-rank type A test theorem.
The existing all-degree E7 multiplicity theorem remains proved by the
earlier case/transfer arguments. This note does **not** replace those
arguments by a uniform all-degree proof. Literature priority is unestablished.

## 1. The question in Springer and Borel notation

Let G be simply connected complex E7, X=X(T), Q the root lattice, and
R=Z_(2). Number nodes as in Bourbaki. Retained nodes I define W_I and
P_I. Write

\[
 M_I=\beta^{P_I}((\operatorname{Sym}_R Q_R)^{W_I}),\qquad
 N_I=\beta^{P_I}((\operatorname{Sym}_R X_R)^{W_I}).
\]

These are the full unscaled images, with polynomial degree d in
H^(2d). In the notation of the Kumar papers the map is
xi_rho^(P_I)=beta^(P_I) after (theta_rho^*)^(-1) on polynomial characters.
The previously proved literal-source formulas give
J_(56,I)^d=12^d M_I^d and J_(133,I)^d=36^d N_I^d. Thus the arithmetic
here compares the sources feeding those same Springer/Borel maps.
It does not introduce a new characteristic map or determine H/N_I.

The **112-case family** consists of I omitting at least one of
2,5,7. Equivalently X^(W_I) maps onto X/Q=Z/2. The already proved formula is

\[
 N_I^d/M_I^d=(\mathbb Z/2)^{m_{I,d}},\qquad
 \sum_d m_{I,d}q^d=\frac{qP_{W(E_7)}(q)}{(1+q)P_{W_I}(q)}.       \tag{1}
\]

The [D4 completion](E7_D4_REFINEMENT.md#2-common-quadratic-lifts-and-the-full-integral-d4-source)
already gives a uniform lattice split and quadratic lift, hence the
exponent two. The issue addressed here is the **multiplicities** in (1).

## 2. A general criterion that would force the multiplicities

Let A be a nonnegatively graded R-algebra, free of finite rank in each
degree, let B in A_2, and put

\[
 S=A[h]/(h^2-B),\quad \deg h=1.
\]

Suppose a graded map S -> E to a torsion-free R-algebra has image N
and kernel K. Put D=im(A -> N), M=D+2hD. Then N=D+hD, M is a subring,
and 2N is contained in M. These statements require no multiplicity formula.
Write bars for reduction of S, and define

\[
 \partial(a+hb)=b,\qquad \overline a,\overline b\in A/2A.
\]

This is a degree-minus-one derivation of S/2S: in characteristic two
it respects h^2-B. It squares to zero.

**Quadratic-derivation criterion.** If

\[
 \partial\big((K+2S)/2S\big)\subseteq(K+2S)/2S,                \tag{2}
\]

then it descends to N/2N. For
C=im(D -> N/2N)=M/2N there is an isomorphism of graded algebras

\[
 N/2N\simeq C[t]/(t^2-\overline B),\quad
 N/2N=C\oplus\overline h C,\quad N/M\simeq C[-1].             \tag{3}
\]

Here C[-1]_d=C_(d-1); the last isomorphism is an isomorphism of graded
C-modules. Consequently

\[
 \sum_d\dim_{\mathbb F_2}(N^d/M^d)q^d
       =\frac{q}{1+q}\sum_d\operatorname{rank}_R(N^d)q^d.     \tag{4}
\]

**Proof.** The source decomposition gives N/2N=C+bar(h)C. The descended
derivation kills C and sends bar(h)c to c. Thus this sum is direct and
ker(partial)=C. Multiplication satisfies bar(h)^2=bar(B), proving the
algebra assertion. Modulo 2N, M maps exactly to C; quotienting proves
the last assertion and the Hilbert identity. Conversely the displayed
quadratic presentation with the specified C and h defines partial and
implies (2). No invariant averaging or relative Smith count enters this
argument. The same argument works through a fixed degree if kernel
stability is known through that degree.

Condition (2) is essential. For example take the map from
R[a,h]/(h^2) onto R[a]/(a^2) sending h to a, with both variables of degree
one. Here h^2 is a root-image class and 2N is contained in M, but D=M=N.
The derivation kills a and sends h to 1, so it cannot respect h-a.
Thus a quadratic lift and the exponent bound alone do not imply (4).

## 3. What this says for all 112 E7 parabolics

Choose an omitted i in {2,5,7}; put K_i={1,...,7} minus {i},
eta=omega_i and V_i=sum_(j in K_i) R alpha_j. The odd determinants are
n_i=7,15,3, respectively. The earlier lattice calculation proves

\[
 X_R=V_i\oplus R\eta,\quad Q_R=V_i\oplus R(2\eta),\qquad
 B_i=-\frac{n_i}{2}\alpha_{K_i}^t C_{K_i}^{-1}\alpha_{K_i},
 \quad\beta(\eta)^2=\beta(B_i).                              \tag{5}
\]

B_i is integral and W_(K_i)-invariant. With
A=(Sym_R V_i)^(W_I), the actual full sources are A[eta] and A[2eta].
Reduction by the monic quadratic gives exactly A plus eta A, including
invariants; eta is fixed and V_i preserved. Thus the criterion applies
with h=beta(eta), provided (2) holds for this **invariant-source kernel**.

There is one canonical source derivation, independent of the choice of i:

\[
 \partial(\textstyle\sum a_j\omega_j)=a_2+a_5+a_7\pmod2.      \tag{6}
\]

It kills every root and is W-equivariant on Sym_(F2)(X/2X). In the
coordinates (5) it is d/deta. Equation (6), not a guessed modular
invariant ring, specifies the proposed mechanism.

**Bounded uniform theorem.** For every one of the 112 parabolics,
condition (2), the quadratic decomposition (3), and formula (1) hold
through polynomial degree 10 by the structural argument below,
independently of the previously computed relative multiplicities.

**Proof and exact coefficient inputs.** Work first in the full quadratic
source Sym_R(X)/(q_2), a free module on 1,omega7 over Sym_R(V_7).
Fresh minuscule weight moments supply primitive relations F6,F8,F10.
Their odd-omega7 coefficient valuations are respectively 2,1,2. Hence
partial(Fj)=0 modulo two. This assertion is coordinate independent:
q_2 is a primitive full-Weyl-invariant quadratic with a unit eta^2
coefficient in every split (5), and (6) is intrinsic.

In degrees 6,7,8,9,10 the columns of all Fj times quadratic-source
monomials have ranks 1,7,28,84,210 modulo two, with no saturation needed.
These equal the rational Borel-kernel dimensions. They are independent
over R and primitive, so they give the whole integral kernel in these
degrees. Below six the kernel is zero. This uses the established rational
Borel degrees 2,6,8,10,12,14,18 and torsion-free target, not (1).

Each Fj is a full-Weyl-invariant **class** in the quadratic source.
The expression of a kernel element as sum Fj*u_j is unique through
degree 10. If that element is W_I-invariant, uniqueness forces each
u_j to be W_I-invariant. This supplies the required integral invariant
lifting. It is legitimate for every I in the family: use its split (5)
to identify the invariant quadratic source with A plus eta A.

After reduction modulo two, each Fj lies in A/2A and is killed by
partial. Its invariant multiples therefore form a derivative-stable
kernel, proving (2). This argument also gives, in the indicated degrees,

\[
 C_I=(A/2A)/(\overline F_6,\overline F_8,\overline F_{10}),\qquad
 m_{I,d}=[q^{d-1}]\,(1-q^6)(1-q^8)(1-q^{10})H_A(q),          \tag{7}
\]

where H_A=1/((1-q)^(6-r) product_j(1-q^(d_j))) for rank-r W_I
with fundamental degrees d_j. The primitive independent columns justify
this truncated Hilbert computation. The omitted higher E7 factors start
at degree 12, so (7) agrees with (1) throughout the claimed range.

For the full flag alone, a second exact calculation proves (2) through
degree **13**, including the primitive half-sum corrections in degrees
11,12,13. The checks are:

| degree | full kernel rank | saturation index valuation | even projection rank | odd projection rank |
| --- | --- | --- | --- | --- |
| 10 | 210 | 0 | 148 | 62 |
| 11 | 462 | 1 | 314 | 148 |
| 12 | 924 | 8 | 610 | 314 |
| 13 | 1716 | 35 | 1106 | 610 |

The two projection ranks add to the full rank, so the reduced kernel
splits by eta parity. Derivatives of its columns lie in the preceding
degree's **primitive kernel**, checked separately. No relative image
table is read to prove these assertions.

The [primary calculation](experiments/e7-quadratic-derivation.json)
rebuilds the 56-weight orbit and saturates exact integer columns. The
[independent audit](experiments/e7-quadratic-derivation-audit.json)
instead constructs the weights from 27 E6 forms, replays the recorded
halving witnesses in integer congruences, and uses its own bit elimination
to check kernel rank, parity projections and derivative membership. It
also checks the three lattice splits and their common derivative (6).
The rational Borel degree list and the saturation witnesses are shared;
the weight construction and rank/membership arithmetic differ.

## 4. An arbitrary-rank type A model where the criterion is proved

This is a test of the proposed mechanism beyond E7, not a statement about
the all-even block family or about odd-primary factors.

**Type A theorem.** Let G=SL_(2r)(C), r odd, and let P have Levi block
sizes b_1,...,b_s summing to 2r, with at least one b_j odd. For the full
root and weight sources over R=Z_(2),

\[
 N_P=H^{2*}(G/P;R),\quad 2N_P\subseteq M_P,\qquad
 N_P^d/M_P^d\simeq(\mathbb Z/2)^{m_{P,d}},\qquad
 \boxed{\sum_d m_{P,d}q^d=
 \frac{q}{1+q}\frac{[2r]_q!}{\prod_j[b_j]_q!}.}              \tag{8}
\]

The quotient notation in (8) is 2-local; no claim that the global
integral quotient has only factors two is made.

**Proof.** Present X_R by x_1,...,x_n, sum x_j=0, n=2r. Each x_j has
the same nontrivial class in X_R/Q_R. Since n is even and some block is
odd, there is a block boundary k odd. Put eta=x_1+...+x_k. Define the
W_P-invariant R-linear functional ell by

\[
 \ell(x_j)=1/k\ (j\le k),\qquad
 \ell(x_j)=-1/(n-k)\ (j>k).
\]

It respects sum x_j=0, has ell(eta)=1, and all ell(x_j) are odd units.
For V=ker ell, X_R=V plus R eta and Q_R=V plus R(2eta).
Indeed Q_R is the kernel of coefficient-sum parity because n has
2-adic valuation one. The functional's reduction is that parity.
W_P fixes eta and preserves V.

Write v_j=x_j-ell(x_j)eta. The v_j sum to zero on each side of the
boundary. The ordinary quadratic Chern relation is

\[
 e_2(x)=e_2(v)-\frac{r}{k(n-k)}\eta^2,
 \quad\text{hence}\quad \eta^2=\frac{k(n-k)}r e_2(v).        \tag{9}
\]

All denominators are odd. The right side is a full root-source
W_P-invariant. Thus the quadratic setup, including 2N contained in M,
holds in every degree.

The full weight invariant source is the ring in block elementary
functions modulo global e1. To justify taking this quotient integrally,
choose an odd block of size b in R[x_1,...,x_n]. Subtract E/b from the
variables in that block, where E=sum x_j. The new variables sum to zero
and, together with E, give an equivariant polynomial splitting
R[x]=Sym_R(X_R)[E], with E fixed. Taking invariants coefficientwise and
then setting E=0 proves the assertion from ordinary symmetric polynomial
theory, without dividing by |W_P|.

The published Chern presentation now identifies its Borel quotient
with H^(2*)(G/P;R), with kernel generated by global e2,...,en. On
reduction modulo two our derivative is the simultaneous translation
derivative partial(x_j)=1. It satisfies

\[
 \partial(e_j(x))=(n-j+1)e_{j-1}(x).                         \tag{10}
\]

Consequently it preserves the **entire** Chern ideal. It is also well
defined modulo e1 since n is even. This proves (2) for the full invariant
source, in every degree. Apply the criterion and the classical flag
Poincare polynomial to obtain (8). No finite Smith data prove the theorem.

An [independent native Python calculation](experiments/type-a-quadratic-derivation-audit.json)
checks the complete 2-primary images on Gr(1,6), Gr(3,6), Gr(1,10),
Gr(3,10), Gr(5,10), of ranks 6,20,10,120,252. It generates the full
root image modulo two from centered elementary functions and computes
multiplication closure in a Schur basis using horizontal and vertical
Pieri rules. It compares with (8) only afterwards. This is independent
finite corroboration, not the arbitrary-rank proof. The audit inventories
28 odd-block parabolics for SL6 and 496 for SL10. The respective 3 and 15
proper all-even block cases are excluded.

## 5. What remains for the E7 structural project

The all-degree E7 formula already implies (3) by dimension counting, as
recorded in the D4 note. Defining partial from that known presentation is
a valid **consequence** of the existing theorem, not a new proof of it.
The contribution here is a precise kernel criterion, an independent
bounded verification, a uniform degree-10 invariant-lifting argument,
and a complete arbitrary-rank model illustrating the criterion.

The missing E7 step is to prove (2) for the full invariant-source kernels
in every degree, without importing (1). Starting in degree 11, rational
relation multiples require integral half-sums. Even full-flag kernel
stability does not by itself show those corrections lift through Levi
invariants. Neither reduction to the ambient mod-two cohomology ring
nor replacing A/2A by all modular invariants justifies that passage.
The source obstruction to a six-generator E6 polynomial invariant ring
and the nonfree D5-to-D4 transfer remain relevant.

A useful next target is an invariant-compatible generating or saturation
procedure for this kernel, proving derivative stability at every stage.
That would replace the earlier separate E6/D5/D4 image-index checks by
one mechanism. It is unresolved here. The other fifteen E7 parabolics
have no such fixed odd degree-one class and need a different mechanism.

Published foundations are rational Borel theory and type A symmetric,
Chern and Schur calculus. The ordinary integral partial-flag Chern
presentation used above is stated explicitly in
[Buch, section 2.3](https://sites.math.rutgers.edu/~asbuch/papers/qhpartial.pdf).
The criterion is elementary characteristic-two algebra; no claim is made
that it or (8) is absent from existing literature. This turn checked that
source for the Chern presentation, not the priority of the relative formula.

## Reproduction

From borel/ in Sage and native Python, respectively:

```sh
PYTHONPATH=src sage -python experiments/e7_quadratic_derivation.py
python experiments/audit_e7_quadratic_derivation.py
python experiments/audit_type_a_quadratic_derivation.py
```

The E7 audit uses NumPy and SymPy, no Sage or library imports. All three
reports bind their sources; the E7 archive preserves saturation witnesses.
