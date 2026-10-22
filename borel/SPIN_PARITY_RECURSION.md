# An arbitrary-rank recurrence for odd-spin Borel kernels

Status, 22 October 2026: a proved recursive presentation and exact decision
criterion in every finite rank, an infinite family with stable entire
kernels, and an arbitrary-rank lower bound for the first failure. A closed
classification of all ranks and a closed first-failure formula are **not**
proved. The remaining classes n=0,3 modulo 4 contain both behaviors.
This advances the [previous spin calculation](SPIN_PARITY_EXTENSION.md).

## 1. Exact setup and classical input

Use G=Spin(2n+1,C), B a Borel subgroup, R=Z_(2), the full simply connected
weight lattice X, and beta:Sym_R(X_R) -> H^(2*)(G/B;R). The central character
chi:X/Q -> F2 gives D(z)=1 and D(epsilon_i)=0 modulo two, where
z=omega_n=(sum epsilon_i)/2. Put c_i=e_i(epsilon_1,...,epsilon_n), c1=2z.
The question is preservation of bar K=(ker beta+2S)/2S, equivalently
descent to the actual image N/2N, rather than to ambient H/2H.

The earlier integral Artin-source decomposition reduces this exactly to
E_n=R[z,c2,...,cn]. Write K_n for its kernel. All degrees here are polynomial
degrees; cohomological degrees are twice these degrees.

Let

\[
 A_n=R[p_1,\ldots,p_n]/(\rho_1,\ldots,\rho_n),\qquad
 \rho_i=p_i^2+2\sum_{k=1}^{i-1}(-1)^kp_{i+k}p_{i-k}+(-1)^ip_{2i},       \tag{1}
\]

with p_j=0 for j>n. This is the **published classical** integral cohomology
presentation of OG(n,2n+1), whose special Schubert classes are p_i.
It is obtained at quantum parameter zero in
[Kresch--Tamvakis, introduction and Theorem 1](https://math.umd.edu/~harryt/papers/qcogfinal.pdf).
The same source identifies the odd maximal Grassmannian with the appropriate
even spinor variety. Its strict-partition Schubert basis is P_lambda,
lambda a subset of {1,...,n}, with degree |lambda|.

The characteristic-source substitution is

\[
 z\longmapsto p_1,\qquad c_i\longmapsto 2p_i.                         \tag{2}
\]

For our purposes it can also be justified without identifying a Levi
invariant ring: (1) implies the rational Borel relations after this
substitution, A_n is torsion-free, and the rational quotient has the same
dimension 2^n. Hence its source kernel is exactly K_n. Define

\[
 M_n=R[2p_1,\ldots,2p_n]\subset A_n,\qquad N_n=M_n[p_1].             \tag{3}
\]

Thus E_n/K_n=N_n. Classical Schur-Q straightening and the one-box Pieri rule
give, integrally,

\[
 M_n=\bigoplus_{\lambda\subseteq\{1,\ldots,n\}}
              R\,2^{\ell(\lambda)}P_\lambda,\qquad
 U P_\lambda=p_1P_\lambda=\sum_{\mu\gtrdot\lambda}P_\mu.              \tag{4}
\]

Here ell(lambda) counts parts; a cover adds 1 as a new part, or replaces
a part i by i+1, keeping distinct parts at most n. The identity
Q_lambda=2^ell(lambda)P_lambda and the Pieri rule underlying (4) are
classical symmetric-function inputs, not results newly claimed here.
For clarity about integrality, Q_lambda is an integer polynomial in the
c_i. Conversely products of the c_i straighten integrally into these
strict Q_lambda; terms with a part >n vanish. This proves both lattice
inclusions in (4), not just a rational spanning assertion.

## 2. The recursive presentation, in every rank

Let V_(n,d) be the free Schubert lattice on the strict partitions of d
with parts at most n. Put T_n=n(n+1)/2. Form the explicit diagonal matrix
M_(n,d)=diag(2^ell(lambda)) in this basis, and the 0/1 cover matrix U_d
from degree d to d+1. Then the image lattices satisfy

\[
 N_{n,0}=R,\qquad
 N_{n,d}=M_{n,d}+U_{d-1}N_{n,d-1},\qquad N_{n,d}=0\ (d>T_n).         \tag{5}
\]

This is an equality of embedded integral lattices, since N_n=sum_(k>=0)
U^k M_n. It is not a saturation guess. Integer Hermite bases compute it
without rational invariant ideals or multivariate Groebner elimination.

Choose a basis matrix H_d for (5), with H_0=(1). Let

\[
 A_d=[M_d\mid U_{d-1}H_{d-1}]=H_d T_d.                             \tag{6}
\]

The entries of T_d are integral: its columns express generating image
classes in their full image-lattice basis. They give a presentation of
N_(n,d) by root-character classes and z times previous basis classes.
Their integral relation modules, together with the classical source
relations for M_n, give a recursive module presentation of N_n and hence
of K_n. This is an explicit algorithmic presentation; it is not a closed
minimal list of polynomial ideal generators.

Suppose the derivative has descended through degree d-1, with matrix
delta_(d-1) in the image bases. Set delta_0=0 and

\[
 B_d=[0\mid H_{d-1}+U_{d-2}H_{d-2}\delta_{d-1}],\qquad
 C_d=H_{d-1}^{-1}B_d\pmod2.                                       \tag{7}
\]

Use zero for the second term when d=1. The inverse means integral
coordinates in the image lattice, not reduction in ambient cohomology.
Root-character classes have derivative zero; z*b has derivative
b+z*D(b), which proves (7). Choosing different integral representatives
of delta changes the answer by 2N and has no effect here.

**Recursive classification criterion.** The derivative extends through
degree d if and only if

\[
                 \ker(\overline T_d)\subseteq\ker(C_d).              \tag{8}
\]

If so, solve delta_d*bar(T_d)=C_d and continue. The first failed test
is exactly the first failure degree. If none fails through T_n+1,
the entire kernel is preserved in every degree.

**Proof.** T_d presents N_d/2N_d from the displayed generating classes.
Their proposed derivatives factor through that quotient exactly when
(8) holds. If a bit relation x has T_d*x=2y, the corresponding source
polynomial minus twice a representative of H_d*y is an actual integral
kernel element. Its derivative is C_d*x modulo 2N_(d-1). This proves
both necessity and sufficiency and produces an obstruction when the
test fails. Above T_n+1 both adjacent image groups vanish, so there is
nothing further to check. The unit Artin basis transfers this exact
decision, including the first failure degree, to Spin(2n+1)/B. QED.

Equations (4)--(8) are therefore a finite, explicit recursive answer for
each n, requiring only matrices indexed by strict partitions and integer
lattice arithmetic. A simplification of their output into a closed
description of all ranks remains separate work.

## 3. An infinite family: every n=1,2 modulo 4

**Theorem.** If n=1 or 2 modulo 4, the central-parity derivative preserves
the entire actual integral full-flag Borel kernel.

**Proof.** Work on the Schubert basis indexed by subsets of {1,...,n}.
Besides U, define a lowering operator L as follows. Removing a part 1
has coefficient T=T_n; replacing a part i+1 by i, when permitted, has
coefficient

\[
                b_i=(n-i)(n+i+1)=n(n+1)-i(i+1).                    \tag{9}
\]

Each b_i is even. A direct cover calculation gives, on degree d,

\[
                    [L,U]=(T-2d)\,\mathrm{id}.                    \tag{10}
\]

Here is an all-rank verification. Identify a subset with its ordered
exterior word, and let C_i,A_i be creation and annihilation. Adjacent
moves have positive signs, so U=C1+sum C_(i+1)A_i and
L=T*A1+sum b_i*C_i*A_(i+1). The standard exterior relations reduce
their commutator to
T+(b1-2T)N1+sum_(i=2)^(n-1)(b_i-b_(i-1))N_i-b_(n-1)N_n.
Substituting (9) gives T-2*sum i*N_i, which is (10).

Moreover L(M_n) is contained in 2M_n: same-length downward edges have
even coefficient, and removing 1 changes 2^ell into twice the generator
2^(ell-1). For homogeneous m in M_n of degree e, (10) implies

\[
 L(U^k m)=U^k L(m)+k(T-2e-k+1)U^{k-1}m.                           \tag{11}
\]

Thus L preserves N_n. If T is odd, the second coefficient is k modulo
two, and the first term vanishes modulo 2N. Consequently L induces
precisely D(z^k*m)=k*z^(k-1)*m on N_n/2N_n. This is the required
derivative, so the kernel is stable. Finally T_n is odd exactly when
n=1 or 2 modulo 4. QED.

This is an arithmetic deduction from the classical Schur-Q lattice and
Pieri operators. Raising/lowering constructions themselves are established
combinatorics; see, for context,
[Proctor, Representations of sl(2,C) on Posets and the Sperner Property](https://epubs.siam.org/doi/abs/10.1137/0603026).
The proof here explicitly checks the integral lattice condition needed
for our derivative, which a rational Lefschetz statement alone would not
establish. No novelty or necessity of this congruence condition is claimed.

## 4. A universal stable range in degree

**Theorem.** In every rank n, the entire kernel is derivative-stable
through polynomial degree n. In particular, any failure has degree >=n+1.

**Proof.** In the infinite classical Schur-P ring, let
delta(p1)=1 and delta(p_i)=2p_(i-1) for i>1. If
C(t)=1+2*sum_(i>=1)p_i*t^i, then delta C(t)=2t*C(t).
Thus delta(C(t)C(-t))=0, so delta annihilates every integral quadratic
relation rho_i, not merely the rational relations. It preserves
N_infinity=R[p1,2p2,2p3,...], and on its source its reduction is D:
delta(z)=1, delta(c2)=4z, delta(c_i)=2c_(i-1) for i>2.

Passing to A_n truncates special classes with indices >n. In degrees
at most n neither the source nor the strict-partition target changes,
so their kernels agree with the infinite ones. The infinite derivative
therefore descends in those degrees. The full-flag conclusion follows
from the same Artin-source decomposition. QED.

The [previous quartic theorem](SPIN_PARITY_EXTENSION.md#3-an-arbitrary-rank-theorem-about-the-quartic-obstruction)
is sharper at degree four: failure there occurs exactly at n=3.

## 5. Exact checks and limits of the classification

The recurrence gives the following entire-kernel decisions. Rank 11 is
a stable example outside the proved congruence family; rank 12 supplies
a later obstruction.

| Ranks n | Entire kernel stable? | First failure degree |
| --- | --- | --- |
| Every n=1,2 modulo 4 | Yes, uniform theorem | None |
| 3 | No | 4 |
| 4 | Yes, complete certificate and earlier closed presentation | None |
| 7 | No | 11 |
| 8 | No | 16 |
| 11 | Yes, complete certificate | None |
| 12 | No | 34 |

The [primary recurrence](experiments/spin_parity_recursion.py) records
all ranks 1--12 and supplies compressed image-basis and derivative
[certificates](experiments/spin-parity-recursion-certificates.json.gz).
The [native independent audit](experiments/audit_spin_parity_recursion.py)
recomputes Hermite indices with SymPy, verifies the embedded lattices by
exact triangular division, and checks each derivative factorization or
nonzero obstruction with bit arithmetic. It also checks the symbolic
coefficients in (10), all sparse commutators through rank 9, and nine
infinite quadratic derivation identities. Finite checks corroborate the
two written arbitrary-rank arguments; they are not their proofs.

A [separate quadratic calculation](experiments/check_spin_parity_quadratic.py)
uses the generators z and 2p_j and direct straightening of (1) into
squarefree special-class monomials. It does not use (4), Schur-Q diagonal
lattices, Pieri matrices or saved bases. This independently checks the
decisions in ranks 3,4,7,8,11,12. It shares Sage arithmetic with the primary
calculation; the native Hermite/bit audit does not. The classical
Grassmannian presentation is a shared mathematical foundation.

Reports: [recurrence](experiments/spin-parity-recursion.json),
[native audit](experiments/spin-parity-recursion-audit.json),
[quadratic check](experiments/spin-parity-quadratic-check.json).

**What remains unresolved here:** a closed classification within the
classes n=0,3 modulo 4 and a closed numerical formula for the first
failure. Neither class can be uniformly declared failing: n=4 and n=11
are stable. Nor does stability hold in all sufficiently small tested
degrees: n=12 first fails at degree 34. No untested rank is assigned a
status beyond the proved uniform family, and no congruence conjecture is
inferred from this table.

**Known versus deduced:** (1), Schur-P/Q theory, the one-box Pieri rule
and the general raising/lowering framework are established literature.
The explicit image-lattice recurrence and derivative factorization test,
the stated congruence-family deduction and the stable-range argument are
proved above. This check has not established priority for their precise
central-parity formulation. The full-flag results alone do not assert a
parabolic descent theorem or a root-image Smith formula.

**Next mathematical target:** simplify (8) for even T_n, using the
integral root-character lattice and the raising/lowering operators. An
extra integrality statement for a divided lowering operator is required;
dividing L by two and assuming it is the desired derivative is invalid.
In fact, when T is even, (11) shows L(N) is contained in 2N, but
(L/2)(z^2)=(T-1)z is nonzero modulo 2N, whereas D(z^2)=0. Thus this
particular divided operator cannot solve the remaining cases even when
it is integral. A different correction or a boundary-lattice argument
is required.
This isolates the missing step rather than declaring the classification
complete on finite evidence.

The [boundary follow-up](SPIN_PARITY_BOUNDARY.md) now proves that the exact
obstruction is obtained from ker(N_(n+1) -> N_n) by lowering the boundary
part once. It gives an adjacent-rank matrix-rank and image-index formula
for the obstruction in every degree, without propagating derivative
matrices. A closed numerical evaluation in n is still required.

## Reproduction

```sh
PYTHONPATH=src sage -python experiments/spin_parity_recursion.py
python experiments/audit_spin_parity_recursion.py
PYTHONPATH=src sage -python experiments/check_spin_parity_quadratic.py
```
