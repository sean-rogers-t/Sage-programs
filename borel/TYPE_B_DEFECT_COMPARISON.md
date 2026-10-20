# The integral comparison for the type-B defect filtration

**Status: repository deduction and a resolved comparison criterion,
19 October 2026.** The primitive incidence pullback below has an integral
proof. The subsequent [general argument](TYPE_B_GENERAL_SMITH.md)
uses the Chow envelope to construct integral triangular rows, resolves
the criterion, and deduces ambient Chern generation. This is unpublished
and makes no priority claim.

Use \(A\ge\ell\), \(n=A+\ell\), \(r=\lfloor n/2\rfloor\),
\(c=\lceil n/2\rceil\), \(\epsilon=c-r\le1\), and a maximal-rank
square-zero map \(N:V\to V\). Set \(I=\operatorname{im}N\), \(K=\ker N\).
For \(0\le t\le\ell\), let
\[
W_t=\{U\in\operatorname{Gr}(\ell,V):
\delta(U):=\operatorname{rank}(N:U\to V/U)\le t\}.
\]
The [candidate note](TYPE_B_CLOSED_CONJECTURE.md) counts its defect
strata integrally, giving free even cohomology and
\[
\sum_m\operatorname{rank}H^{2m}(W_t;\mathbb Z)q^m
=\sum_{e=0}^t E^{\rm geom}_{n,\ell,e}(q).
\tag{1}
\]
That counting statement alone does not identify Smith multiplicities.

## 1. The incidence variety has the required Schur-product rows

Put \(b=\ell-t\), \(s=\lfloor b/2\rfloor\), \(p=\lceil b/2\rceil\).
Consider the smooth projective iterated Grassmann bundle
\[
Y_t=\{(F,S,U):F\subset I,\ \dim F=s;\
F\subset S\subset N^{-1}F,\ \dim S=b;\
S\subset U,\ \dim U=\ell\}.
\tag{2}
\]
First choose \(F\in\operatorname{Gr}(s,I)\). Then choose
\(S/F\in\operatorname{Gr}(p,E)\), with
\[
E=N^{-1}F/F,\quad \operatorname{rank}E=c,\quad
0\to K/F\to E\to F\to0,\quad c(E)=1.
\]
Finally choose \(U/S\in\operatorname{Gr}(t,V/S)\), whose complementary
rank is \(A\). The integral Leray--Hirsch bases are therefore indexed by
\[
\beta\subseteq(r-s)^s,\quad
\alpha\subseteq(c-p)^p,\quad
\mu\subseteq A^t.
\]
Branching the Schur classes of \(U\) through
\(F\subset S\subset U\) shows that the rows of the pullback from the
ordinary Grassmannian to \(Y_t\) are precisely the three-factor
Littlewood--Richardson rows defining \(P_t\).
No reduction of an out-of-bounds fiber Schur class is needed:
the source partition already fits \(A^\ell\).

Since \(S\) is invariant, \(U\cap N^{-1}U\) contains \(S\), so the
forgetful map \(f:Y_t\to\operatorname{Gr}(\ell,V)\) factors through \(W_t\).

## 2. A generic-field lift proves integral primitivity

The following argument works over every extension field of \(\mathbb C\).
For \(U\in W_t\), put
\[
e=\delta(U),\quad d=\ell-e\ge b,\quad
C=NU\cap U,\quad j=\dim C,\quad J=U\cap I.
\]
The space \(T=U\cap N^{-1}U\) is invariant of dimension \(d\),
with image \(C\) and kernel \(U\cap K\) of dimension \(d-j\).
Thus \(j\le\lfloor d/2\rfloor\), \(C\subset J\), and
\[
\dim J\ge d-j-\epsilon.
\]
In particular \(\dim J\ge s\): if \(j\ge s\), use \(C\subset J\);
if \(j<s\), use
\(d-j-\epsilon\ge b-(s-1)-1=p\ge s\).

Choose an \(s\)-plane \(F\subset J\) whose intersection with \(C\)
has dimension \(\min(s,j)\). This means choosing \(F\subset C\) when
\(j\ge s\), or extending \(C\) inside \(J\) otherwise.
Because \(F\subset U\), \(F\cap NU=F\cap C\). It follows that
\[
\dim(U\cap N^{-1}F)
=d-j+\min(s,j)\ge b.
\]
For \(j\le s\) this is \(d\). For \(j>s\), it is at least
\(\lceil d/2\rceil+s\ge p+s=b\).
Also \(F\subset U\cap N^{-1}F\), since \(F\subset I\subset K\).
Extend \(F\) to a \(b\)-plane \(S\) in that space. This constructs a
point of (2) above \(U\), using only vector-space basis extension.

Apply this argument at the function field of any irreducible subvariety
\(Z\subset W_t\). It gives a rational section over a dense open subset
of \(Z\). The closure of its graph in the proper variety \(Y_t\)
maps birationally onto \(Z\). Hence every integral cycle of \(W_t\)
has a degree-one lift: \(f\) is a Chow envelope, and \(f_*\) is
surjective on integral Chow groups.

The defect stratification used for (1) also supplies the required
homology statement. Filter by \(e\), then within that stratum by
\(j=\operatorname{rank}(N|_U)-e\). Its pieces are iterated Grassmann
bundles and affine bundles: choose \(D=NU\), \(C=D\cap(U\cap K)\),
the image of \((U\cap K)/C\) in \(K/D\), and then the two graph lifts.
Refine the Grassmann bundles by Schubert filtrations; the bundles
trivialize over the resulting affine cells. This gives affine pieces
in a closed filtration. Localization and the compact-support exact
sequence show that the integral cycle classes generate the free
even homology. The same freeness holds on the Grassmann bundle \(Y_t\).
Therefore
\[
f_*:H_*(Y_t;\mathbb Z)\twoheadrightarrow H_*(W_t;\mathbb Z),
\qquad
f^*:H^*(W_t;\mathbb Z)\hookrightarrow H^*(Y_t;\mathbb Z)
\tag{3}
\]
is split injective degreewise. This proves primitivity of the incidence
pullback in every rank, independently of ambient Chern generation.

## 3. The comparison criterion and its subsequent resolution

For each \(t\), the following statements are equivalent:

1. Ambient restriction \(i^*:H^*(\operatorname{Gr}(\ell,n);\mathbb Z)
   \to H^*(W_t;\mathbb Z)\) is surjective.
2. Every \(P_t(m)\) is primitive and its rank equals the degree-\(2m\)
   Betti number of \(W_t\).

Indeed the composite \(f^*i^*\) has row lattice \(P_t\), by Section 1.
If \(i^*\) is onto, (3) makes its image primitive with the stated rank.
Conversely, primitivity and that rank make the image of \(f^*i^*\)
equal to \(f^*H^*(W_t;\mathbb Z)\): a primitive full-rank sublattice of
this free group is the group itself. Since \(f^*\) is injective,
\(i^*\) is onto. Equivalently, the tautological Chern classes must
generate \(H^*(W_t;\mathbb Z)\).

Proving this for all \(t\) would give the complete type-B formula through
(1) and the [primitive-prefix Smith lemma](TYPE_B_FILTRATION.md).
For \(t=0\), the [unit proof](TYPE_B_UNIT_SMITH.md) supplies it using
an affine-Grassmannian model. The [high-layer argument](TYPE_B_HIGH_LAYERS.md)
supplies the required prefixes at \(t=\ell-3,\ell-2,\ell-1,\ell\).
The [general proof](TYPE_B_GENERAL_SMITH.md) now covers every remaining
$t$: unbalanced incidence cycles lift through the Chow envelope above,
and their triangular pivots meet the Betti upper bound.
The equivalent criterion is a sufficient route to the complete Smith
formula; this note does not assert that it is necessary for every
possible proof of that formula.

## 4. Why a nearby published theorem does not close this step

Cerulli Irelli--Esposito--Fang--Fourier,
[*Specialization Map for Quiver Grassmannians*](https://publications.rwth-aachen.de/record/1036287/files/1036287.pdf),
2026, Corollary 1.2 (PDF p. 3), proves integral ambient cohomology
surjectivity for quiver Grassmannians of Dynkin type A.
The coefficient convention is integral cohomology throughout that paper.

The natural incidence space over \(W_t\) instead chooses an
\((\ell+t)\)-plane \(T'\) with \(U\subset T'\) and \(NU\subset T'\).
These are the two arrows \((\mathrm{id},N):V\rightrightarrows V\),
so it is a Kronecker-quiver Grassmannian, not a Dynkin-A one.
Its fiber over defect \(e\) is
\(\operatorname{Gr}(t-e,n-\ell-e)\).
The cited result therefore does not directly prove the required
restriction theorem. The [general proof](TYPE_B_GENERAL_SMITH.md)
instead obtains it from the integral triangular squeeze.

## 5. An algebraic formulation for further work

Homogenize the complete-symmetric coordinates: let
\(S=\mathbb Z[z_0,\ldots,z_A]\), all \(z_i\) of ordinary degree one
and weight \(i\), and set \(z_i=0\) outside \(0\le i\le A\).
The degree-\(\ell\), weight-\(m\) monomials identify with the
complete basis of \(R_{A,\ell}(m)\), by setting \(z_0\) to the
padding factor \(h_0=1\).
For \(d\ge u\), let \(J_{u,d}\) be the ideal of maximal minors of
\[
T_{u,d}=(z_{j-i})_{0\le i<u,\ 0\le j<d},\qquad J_{0,d}=S.
\]
Its minors are precisely the homogenized Schur functions
\(\lambda\subseteq(d-u)^u\), up to row/column signs.
By the spectator change (2) in the high-layer note,
\[
P_t(m)\cong
\bigl(J_{\lceil b/2\rceil,c}
      J_{\lfloor b/2\rfloor,r}\bigr)_{\ell,m},
\qquad b=\ell-t. \tag{4}
\]
In (4), juxtaposition means the **product of ideals**.
An integral standard-monomial or monic Gröbner theorem for these
particular products, with the required Hilbert series, is another
algebraic formulation of the integral prefix statement.

One cannot replace the product's initial module by the product of
individual initial modules without proof. For \((A,\ell)=(4,4)\),
the degree-four lattice from \(J_{2,4}^2\) has total rank \(19\);
the naive fourth power \((z_0,z_1,z_2)^4\) has rank \(15\) in that
degree. Generic Hankel-matrix theorems also require a separate
specialization argument when the negative-index variables are set to zero.
These cautions explain the next algebraic target, rather than supplying
a separate theorem.
