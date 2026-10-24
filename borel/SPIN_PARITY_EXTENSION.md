# Does the Spin7 obstruction extend to every odd spin group?

Status, 22 October 2026: **no**. The central-parity derivative preserves the
entire integral Borel kernel of Spin9/B, Spin11/B and Spin13/B. It fails on
Spin7/B, Spin15/B and Spin17/B, first in polynomial degrees 4, 11 and 16,
respectively. These are complete results for the specified ranks, with
independent certificates, not an arbitrary-rank classification.

There is also an arbitrary-rank statement: **failure in polynomial degree
four occurs if and only if n=3** for Spin(2n+1)/B. The general symmetric-source
reduction and the proof of this assertion are below. No literature priority
is claimed. The [earlier note](CENTRAL_PARITY_KERNEL.md) defines the question
and explains why ambient mod-two cohomology is insufficient.

Follow-up: the [arbitrary-rank recurrence](SPIN_PARITY_RECURSION.md) now
gives a recursive image-lattice presentation and exact derivative decision
criterion, proves all-degree stability for every n=1,2 modulo 4, and
stability through degree n in every rank. Independent checks add stable
rank 11 and a first failure in rank 12 at degree 34. A closed classification
in the remaining congruence classes is still unresolved here.

## 1. Lattice, map and known foundations

Let G=Spin(2n+1,C), n>=2, R=Z_(2), X=X(T), and Q the root lattice. Use
Cartesian root coordinates epsilon_i and the spin weight

\[
 z=\omega_n=\tfrac12(\epsilon_1+\cdots+\epsilon_n),\qquad
 X=\mathbb Z^n+\mathbb Zz,\quad Q=\mathbb Z^n.
\]

Here Z^n denotes the lattice on the epsilon_i. Put S=Sym_R(X_R),
K=ker(beta:S -> H^(2*)(G/B;R)), and N=beta(S). The nontrivial
character chi:X/Q -> F2 gives the derivation

\[
 D(z)=1,\qquad D(\epsilon_i)=0
 \quad\text{on }S/2S.
\]

We ask whether D(bar K) is contained in bar K, where bar K=(K+2S)/2S.
Equivalently, does D descend to **N/2N**? No representation or scaling of
its character source is chosen here; this concerns the full weight source.

The published inputs are the rational Borel presentation, the signed
permutation Weyl group, and the torsion-free Schubert basis of integral
flag cohomology. Consequently

\[
 K=S\cap (r_1,\ldots,r_n)S_{\mathbb Q},\qquad
 r_j=e_j(\epsilon_1^2,\ldots,\epsilon_n^2).
                                                        \tag{1}
\]

The torsion-free basis and classical presentation context are recalled by
[Duan--Zhao, introduction and Theorem 2.1](https://arxiv.org/html/0801.2444v12).
Equation (1) is an integral intersection, not a claim that the raw r_j
generate K over R.

## 2. A reduction valid for every rank

Write c_j=e_j(epsilon_1,...,epsilon_n), so c1=2z, and define the
**auxiliary symmetric subring** E=R[z,c2,...,cn] inside S. This argument
does not require identifying E with any full integral Levi invariant ring.
Give z degree 1 and c_j degree j.

The ordinary symmetric polynomial theorem and its integral Artin basis
give a free decomposition

\[
 S=\bigoplus_{b\in\mathcal B} E b,\qquad |\mathcal B|=n!,\quad1\in\mathcal B.
                                                        \tag{2}
\]

Indeed, start with R[epsilon_1,...,epsilon_n] free over R[c1,...,cn]
and base change by c1=2z. This is valid integrally even though 2 is not
invertible. In particular, E is a polynomial ring. Each b can be taken
to be a monomial in the epsilon_i and is killed by D modulo two.

Set K_E=E intersect (r_1,...,r_n)E_Q. The free decomposition implies

\[
 K=K_ES,\qquad K\cap E=K_E.                              \tag{3}
\]

To see the first equality, intersect the rational direct sum of copies
of (r_j)E_Q with the corresponding integral direct sum in (2).
The same unit basis vector gives contraction after reduction modulo two.
Therefore D preserves the entire full-flag kernel **if and only if**
partial/partial z preserves bar K_E. This reduces a full root-coordinate
calculation to a weighted symmetric-source calculation, in every rank.

The required rational relations in E have the simple formula

\[
 (-1)^j r_j=\sum_{k=0}^{2j}(-1)^k c_kc_{2j-k},\quad
 c_0=1,\ c_1=2z,\ c_i=0\ (i>n).                        \tag{4}
\]

The rational quotient E_Q/(r_j) has dimension 2^n and Hilbert polynomial
prod_(j=1)^n(1+q^j). This follows from the rational Borel presentation
and (2), or from the regular sequence of degrees 2,4,...,2n in a source
with generator degrees 1,2,...,n. These are known invariant-theoretic
inputs; (2)--(3) are the elementary reduction used here.

This also gives a terminating exact decision procedure in any specified
rank. Form the integer ideal (r_j,2t-1), eliminate t, reduce its generators
modulo two, and test their z-derivatives for membership in the reduced
ideal. Elimination gives the 2-saturation; localizing at all odd integers
identifies it with K_E over R. Reduction modulo two is unchanged by that
localization. Polynomial ideal algorithms over ZZ and F2 terminate, so
this decides the entire-kernel question, including the first failure
degree from homogeneous generators. This algorithmic deduction is valid
for every finite rank; it is not a closed rank classification or an
efficiency guarantee. The finite-flatness certificates below offer a
different way to audit its answers.

## 3. An arbitrary-rank theorem about the quartic obstruction

**Theorem.** For G=Spin(2n+1), n>=1, the central-parity derivative preserves
the reduced integral kernel through polynomial degree three. It fails
in degree four exactly when n=3.

**Proof.** For n>=2 the primitive quadratic is

\[
 q_2=2z^2-c_2=\tfrac12 r_1.
\]

It is monic, up to sign, in c2. Thus E/(q2) is a polynomial R-algebra,
and there is no further saturation in kernel degrees two or three.
Its reduction is c2 and D kills it.

The next rational relation satisfies

\[
 r_2=c_2^2-4zc_3+2c_4,\qquad
 2H_4=r_2-q_2^2+4z^2q_2,\qquad
 H_4=2z^4-2zc_3+c_4.                                    \tag{5}
\]

For n>=4, eliminating c2 with q2 and c4 with H4 shows that
E/(q2,H4) is a polynomial R-algebra. Hence this ideal is saturated and
its degree-four component is the entire integral kernel in that degree:
there are no other rational basic relations below degree six. Modulo two
these two generators are c2 and c4, both killed by D. Equation (3) gives
the full-flag assertion in every such rank.

For n=3, c4=0 and H4 has an additional factor two. Its primitive form is
F4=z^4-zc3. It belongs to K_E, whereas D(bar F4)=c3 is nonzero modulo
the entire degree-three kernel (q2). This recovers the Spin7 obstruction.
For n=2 the same primitive relation is z^4, whose derivative is zero;
the monic quotient by q2 and z^4 gives the entire quartic kernel.
For n=1, Spin3 is SL2, with source R[z] and kernel (z^2), which is stable
in every degree. This proves the theorem. The usual Spin5=Sp4
identification also gives all-degree stability in rank two. QED.

Thus adding a coordinate genuinely changes the saturation: the c4 term
prevents the extra division by two that created F4 in rank three.

## 4. A closed all-degree proof for Spin9/B

In rank four the entire symmetric-source kernel is generated over R by

\[
\begin{aligned}
 f_2&=2z^2-c_2,\\
 f_4&=z^2c_2-2zc_3+c_4,\\
 f_6&=2c_2c_4-c_3^2,\\
 f_8&=z^8+zc_3c_4-c_2^2c_4.
\end{aligned}                                           \tag{6}
\]

First, all four lie in the rational Borel ideal: 2f4=r2+c2f2,
f6=-r3, and modulo f2,f4,f6 one has f8=r4/4, with r4=c4^2.
These identities, including the last congruence, are also checked by
independent rational polynomial division in the audit.

Let L=(f2,f4,f6,f8). Its reduced quotient is exactly

\[
 E/(L+2E)=\mathbb F_2[z,c_3]/(z^8,c_3^2),                \tag{7}
\]

of dimension 16. Graded Nakayama implies E/L is a finite R-module:
every sufficiently high homogeneous component vanishes modulo two,
hence vanishes over R. Its rational rank is at least 16, since L lies
in K_E and the rational Borel quotient has rank 16. The structure
theorem for finite modules over R now forces rank exactly 16 and no
2-torsion, because its mod-two dimension is 16. Thus L=K_E.

The ideal (c2,c4,c3^2,z^8) is derivative-stable, so (3) proves stability
of the **entire** Spin9 full-flag kernel. This is not a bounded test.
It disproves uniform failure for all n>=3. It does not assert that N
is primitive or equal to the ambient integral cohomology.

## 5. Complete additional ranks and independent verification

The same finite-flatness argument supplies complete certificates for
the following ranks. First failure refers to the degree of the relation
being differentiated; cohomological degree is twice polynomial degree.

| n | Full flag | Entire kernel stable? | First failure degree |
| --- | --- | --- | --- |
| 1 | Spin3/B | Yes, classical SL2 presentation | None |
| 2 | Spin5/B | Yes, classical Sp4 presentation | None |
| 3 | Spin7/B | No | 4 |
| 4 | Spin9/B | Yes | None |
| 5 | Spin11/B | Yes | None |
| 6 | Spin13/B | Yes | None |
| 7 | Spin15/B | No | 11 |
| 8 | Spin17/B | No | 16 |

For example the reduced rank-seven kernel contains

\[
 g_{11}=z^8c_3+z^2c_3c_6+zc_3c_7+c_5c_6.
\]

Its derivative c3*c7 has nonzero normal form modulo the **whole**
reduced kernel. In rank eight the reduced relation

\[
 g_{16}=z^{16}+z^8c_8+z^2c_6c_8+c_3c_6c_7+zc_7c_8
\]

has derivative c7*c8, again a nonzero normal form. The
[certificate file](experiments/spin-parity-extension.json) supplies
actual integral lifts and full reduced ideals, not only modular
polynomials assumed to lift. These entire homogeneous generating ideals
also certify the stated first failure degrees: lower-degree generators
are derivative-stable, and higher-degree generators cannot contribute
in lower degrees.

The [primary calculation](experiments/spin_parity_extension.py) uses
integer elimination of (r_j,2t-1) and extracts short lists of integral
relations generating the resulting reduced ideals. The
[independent native audit](experiments/audit_spin_parity_extension.py)
does **not** repeat that saturation algorithm or use Sage/BGG/Smith data.
Using SymPy, it checks every supplied lift lies in the rational Borel
ideal, recomputes the modular Groebner bases and standard monomials, and
finds length 2^n and Hilbert polynomial prod(1+q^j) in each case. The
finite-flatness proof used for (6)--(7) then establishes that these
relations generate K_E in every degree. Finally it checks all generator
derivatives. See its [report](experiments/spin-parity-extension-audit.json).
Shared inputs are the root datum, rational Borel theorem and explicit
certificate polynomials. These are independent arithmetic/proof methods,
not independent external peer review.

## 6. What is known, what was deduced, and what remains

**Known inputs:** rational Borel theory, the Schubert basis, symmetric
polynomial freeness, and the signed-permutation root datum. Existing
integral spin invariant theory is substantial: Benson--Wood,
*Integral invariants and cohomology of BSpin(n)*, Topology 34 (1995),
13--28, is the invariant-ring input cited in
[Karpenko, Proposition 2.4](https://sites.ualberta.ca/~karpenko/publ/spin12-19.pdf).
Karpenko also distinguishes integral invariant reductions from full
modular invariants (Remark 3.3). His characteristic map is CH(BG) ->
CH(BT), whereas our kernel belongs to the map Sym(X) -> H*(G/B).
These maps and their kernels must not be identified.

The [2026 Baek abstract](https://arxiv.org/abs/2607.18188) concerns membership
of recursive invariants in the Chow characteristic image, also a different
question. Lyons, *Integral Weyl invariants and the Borel homomorphism for
Spin(n)*, Communications in Algebra 26 (1998), 2221--2239,
[DOI](https://doi.org/10.1080/00927879808826272), is particularly relevant
further reading; its full text was not retrieved in this check, so no
precise theorem from it is asserted here.

**Deductions recorded here:** the every-rank reduction (2)--(3), the
every-rank quartic theorem, the closed Spin9 proof, and complete exact
kernel certificates in ranks 3--8. These use known foundations; the
checked sources do not establish priority for their derivative formulation.

**Still unresolved here:** classify all ranks for entire-kernel stability
and determine the first failure degree when it fails. The observed mix
of stable ranks and later failures rules out both uniform failure for
n>=3 and uniform stability for n>=4. We do not propose a congruence or
power-of-two classification on this evidence alone. Parabolic descent
and a root-image Smith formula require their own hypotheses and are not
consequences of this full-flag calculation alone.

**Next step:** use the symmetric-source reduction to seek a recursive
integral presentation in arbitrary rank, with particular attention to
the c3*c7 and c7*c8 saturation obstructions. The goal is a rank criterion
and a first-obstruction formula, rather than additional E7 degree checks.

## Reproduction

From borel/:

```sh
PYTHONPATH=src sage -python experiments/spin_parity_extension.py
python experiments/audit_spin_parity_extension.py
```
