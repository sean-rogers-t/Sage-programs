# E8 adjoint fibers: certified Weyl classes and proof transfers

26 October 2026. **Seven genuinely distinct cubic Levi cases remain.**
The 93 cubic proper standard subsets form exactly ten Weyl-conjugacy
classes. Three existing E8 proofs cover three classes containing seven
standard subsets. Thus 169 of the 255 proper standard subsets are covered
by the currently recorded proofs; 86 subsets, forming seven classes,
remain unresolved here. These are counts for the adjoint generic
regular-source smoothness problem, not a classification of generic degrees.

## 1. Statement and normalization

Let G be the simple simply connected complex group of type E8, T a fixed
maximal torus, and alpha_1,...,alpha_8 its Bourbaki simple roots. For a
proper subset I of {1,...,8}, let K_I be the standard Levi whose roots are

    Phi_I = Phi(E8) intersect span_Z{alpha_i : i in I}.

The representation throughout is the full adjoint representation. Its
Springer map theta=theta_Ad is normalized by the Killing form B:

    B(theta(g),X) = tr(Ad(g) ad(X)).

Write z_I for the Lie algebra center of K_I, and

    Z_I^circ = {h in z_I : Z_G(h)=K_I}.

The assertion to transfer is that some dense open U_I in Z_I^circ has
smooth actual scheme fibers theta^(-1)(h) intersect G_reg for every
h in U_I. Empty fibers and smooth positive-dimensional fibers are
permitted. No reduction of the fiber scheme is implicit.

**Coverage theorem (GEO-86).** The three E8 seed results below imply this
assertion for every standard subset in their three certified Weyl classes.
Together with the existing quadratic and non-A obstruction theorems,
they cover 169 proper standard subsets. Exactly seven Weyl classes remain
outside this recorded coverage. This is a finite transfer theorem using
the existing seed proofs; it does not prove the seven remaining cases.

## 2. Why the transfer preserves the entire scheme

Suppose w(Phi_I)=Phi_J, and choose n in N_G(T) representing w. Then
n K_I n^(-1)=K_J and Ad(n)(z_I)=z_J. The Killing form and trace are
conjugation invariant. For every X in g,

    B(theta(n g n^(-1)),X)
      = tr(Ad(g) ad(Ad(n)^(-1) X))
      = B(Ad(n) theta(g),X).

Nondegeneracy of B proves equivariance as an identity of morphisms,

    theta composed with C_n = Ad(n) composed with theta,
    where C_n(g)=n g n^(-1).

Consequently C_n induces an isomorphism of the ACTUAL fiber schemes

    theta^(-1)(h) intersect G_reg
      -> theta^(-1)(Ad(n)h) intersect G_reg.

This is base change of the commuting diagram of isomorphisms, not just
a bijection on complex points. It preserves nilpotents, embedded
components, dimensions, smoothness and emptiness. Regularity of the
source is conjugation invariant, including at nonsemisimple points.
Also Z_G(Ad(n)h)=n Z_G(h) n^(-1), so Ad(n) maps Z_I^circ isomorphically
to Z_J^circ. The transported open is U_J=Ad(n)U_I. No new source
localization, quotient, saturation or finite-cover assumption is used.

Conjugacy of these standard Levis is exactly Weyl conjugacy of their root
subsystems. For the converse direction, if g K_I g^(-1)=K_J, then
g T g^(-1) and T are maximal tori of the connected reductive group K_J.
Conjugate the first to the second by an element of K_J; the resulting
representative lies in N_G(T) and sends Phi_I to Phi_J. Thus counting
root-subsystem orbits counts the genuinely distinct embedded Levi cases.

## 3. The complete cubic class table

Digits in a representative denote retained Bourbaki nodes. The orbit
size counts ALL embedded root subsystems in its Weyl orbit, not only
standard subsets. A matching abstract Dynkin type was not used as a
conjugacy test: full root-subsystem orbits certify the equivalences.
In this cubic census each displayed type happens to be one class.

| Derived type | Representative I | Standard subsets | Full Weyl orbit | Current status |
| --- | --- | ---: | ---: | --- |
| A1^3 | 125 | 21 | 37,800 | Unresolved |
| A2+A1^2 | 1235 | 28 | 302,400 | Unresolved |
| A1^4 | 1257 | 7 | 113,400 | Unresolved |
| A2^2+A1 | 12356 | 8 | 403,200 | Unresolved |
| A2+A1^3 | 12357 | 8 | 604,800 | Unresolved |
| A3+A1^2 | 12457 | 10 | 453,600 | Unresolved |
| A4+A1^2 | 123468 | 4 | 362,880 | Proved by transfer |
| A3+A2+A1 | 123567 | 4 | 604,800 | Unresolved |
| A2^2+A1^2 | 123568 | 2 | 604,800 | Proved by transfer |
| A4+A2+A1 | 1235678 | 1 | 241,920 | Previously proved |

The orbit sizes sum to 3,729,600. Both implementations enumerate all
these subsystems, including those not standard. The 93 standard targets
are partitioned without omissions or duplication. Each remaining
representative is therefore sufficient for its whole listed class.

The earlier count of 90 unclosed standard subsets subtracted only the
three individual seeds. Weyl transfer adds four subsets to the proved
coverage, leaving 86. This changes the proof workload to seven genuinely
distinct cases; it does not itself solve a new residual critical ideal.

## 4. Exactly which existing proofs transfer

1. **A4+A1^2:** seed 123468. The full relative critical surface is smooth
   and integral, including its w^2=1 boundary, by the argument at equation
   (46) in [the uniform proof note](UNIFORM_ADJOINT_PROOF_TARGET.md).
   Generic smoothness gives smooth-or-empty generic target fibers.
   Its certified standard class is
   **123468, 124568, 125678, 235678**.
2. **A2^2+A1^2:** seed 123568. The 25 October quasifiniteness theorem and
   six complete boundary exclusions in
   [the same note](UNIFORM_ADJOINT_PROOF_TARGET.md) give finite etale-or-zero
   generic critical quotients and hence the required regular-source
   smoothness. Its certified class is **123568, 123578**. The word
   **7,6,5,8,7,6**, applied in that order, transfers the first to the second.
   The proof does not establish generic nonemptiness or degree.
3. **A4+A2+A1:** seed 1235678 (node 4 omitted). This is
   [Theorem 3 of the coupled-multiplier note](MULTIPLIER_BRANCHING_FIBERS.md):
   six smooth regular K-orbits of dimension 28 over a generic nonzero
   central target. Its class has one standard representative.

The [machine-readable certificate](e8-weyl-transfers-verification.json)
contains a simple-reflection word from the indicated seed to EVERY one
of the 93 cubic standard subsets, grouped by class. For unresolved
classes these are equivalence certificates only. Listed words act in
order from left to right on roots; the corresponding group product is
s_last ... s_first. The four additional coverage transfers all start
from the two proved seeds above. E7 proofs are not transferred to E8
merely because a Levi has the same derived type: the ambient adjoint
characters differ.

## 5. Exact independent verification and limitations

The [Sage verifier](tools/verify_e8_weyl_transfers.py) constructs the 240
roots in simple-root coordinates, represents a subsystem by its 120-root
line indices, and closes its orbit under all eight simple reflections.
It independently counts all 255 proper subsets by restricted dominant
root degree and checks the 93 cubic subsets against the existing archive.
Breadth-first traversal retains explicit words to every standard member.
Exhaustive closure under the generators proves completeness, not merely
that the displayed pairs are conjugate.

The [independent audit](tools/audit_e8_weyl_transfers.py) does not import
Sage or the verifier. It builds E8 from the Euclidean roots
{+/-e_i +/-e_j} and half-sums with an even number of negative signs.
It uses all 240 signed roots, independently re-enumerates all ten orbits
in reverse generator order, and replays all 93 words. Each word is also
checked as an integral root-lattice isometry preserving the Killing
form and mapping the retained-root span to the target span; orthogonal
complements give the Lie center transfer. Counts and every orbit size
must agree. Its [report](e8-weyl-transfers-audit.json) records the result.
All arithmetic is exact; no finite-field sampling is used.

The geometry of the transfer is proved in Section 2. The three seed
proofs and their previous characteristic-zero audits are inputs; this
conjugacy audit does not rerun those polynomial calculations. In
particular it does not assert reducedness or smoothness for any of the
seven unresolved classes.

**Published inputs:** the equivariant trace-projection construction of
[Kostant and Michor](https://arxiv.org/abs/math/0109066), finite Weyl group
root theory, conjugacy of maximal tori in a connected reductive group,
and standard scheme-theoretic base change. The equivariance argument is
reproved above. **Repository deductions:** the three seed results and
their transferred coverage. **Computational certificates:** exhaustive
orbits and explicit words. No literature-priority claim is made for the
classification or coverage bookkeeping.

## 6. Next proof target

Work on the representative I=125, derived type A1^3. Its success would
settle all 21 standard subsets in that class. The mathematical goal is
to prove generic smoothness from the full E8 adjoint critical ideal,
including its central boundary, and to identify a criterion usable for
other classes. The existing E7 binary method is a starting technique,
not an E8 proof. A complete boundary argument or an exact obstruction
to that method is the next meaningful advance; further conjugate
examples are unnecessary.
