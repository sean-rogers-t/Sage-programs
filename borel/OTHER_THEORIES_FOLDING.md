# Fixed subgroups: compatible Springer maps, different integral sources

Exploration, 22 October 2026. This gives a controlled extension to other
groups and explains one limitation of diagram folding. Root-system folding
and fixed subgroups are classical. The trace-projection compatibility is
an elementary consequence of the defining equation. The source quotients
below are deductions checked from actual weights; their priority has not
been established.

## 1. A compatibility statement with exact hypotheses

Let G be a finite complex reductive group, H a reductive subgroup, and
rho a finite-dimensional representation whose trace form is nondegenerate
on both Lie algebras. Let pi_h:g -> h be orthogonal projection for that
trace form. Then

\[
                 \theta_{\rho|H}^{H}=
                       \pi_{\mathfrak h}\circ\theta_\rho^G|_H.    \tag{1}
\]

Indeed the equation defining theta_G pairs its projection error to zero
with every d rho(X), X in g, and in particular with every X in h. A second
orthogonal projection therefore satisfies exactly the defining equation
for theta_H. Nondegeneracy gives uniqueness. No extra Dynkin-index scalar
is inserted: these are the trace-projection normalizations with identity
differential.

Suppose now sigma is a finite-order pinned automorphism, H=(G^sigma)^0,
and rho composed with sigma is isomorphic to rho. The trace projection
is sigma-equivariant, so theta_G(h) belongs to g^sigma=h. Hence (1)
strengthens to

\[
                       \theta_G|_H=\theta_H.                       \tag{2}
\]

For compatible tori, the weight restriction res consequently gives
res(Lambda_(G,rho)) contained in Lambda_(H,rho|H). This follows by
restricting an actual integral virtual character, not just by comparing
rational root diagrams. Compatible flag pullbacks commute with the
ordinary, K-theoretic and specified oriented characteristic maps on
this restricted source. On partial flags impose the corresponding
compatible Levi/parabolic conditions; the full-flag inclusion has them
automatically for a stable Borel.

For root-system background see
[Stembridge, Folding by automorphisms](https://websites.umich.edu/~jrs/papers/folding.pdf).
Our matrices below are **character restriction** to the fixed torus.
One must distinguish this from a sum-of-orbit-roots convention, which
can interchange long and short roots in a folded root datum.

## 2. Compatibility does not make full integral sources equal

The three examples use simply connected G and its standard fixed group.
Representations are sums over the indicated diagram orbits, so the
sigma-stability hypothesis is satisfied.

| Folding and representation on G | Restricted representation | Restricted source | Full source on H | Quotient |
| --- | --- | --- | --- | --- |
| SL_(2r) -> Sp_(2r), V direct-sum V*, r>=2 | Two defining H modules | 2Q_C | 2X_C | Z/2 |
| E6 -> F4, 27 direct-sum 27* | Two 26 modules and two trivial modules | 12Q_F4 | 6Q_long,F4 | (Z/2)^2 |
| Spin8 -> G2, 8v direct-sum 8s direct-sum 8c | Three 7 modules and three trivial modules | 6Q_G2 | 6Z alpha_short + 2Z alpha_long | Z/3 |

Here X_C=Z epsilon_1+...+Z epsilon_r is the actual weight lattice
of Sp_(2r), and Q_C is its even-coordinate-sum sublattice. The quotient
column is full source divided by restricted source. Trivial summands
do not change a trace projection, and repeated summands have the same
normalized source as one copy.

For the first row, the weights of V and V* are disjoint for 2r>=4.
Their weighted coefficient lattice is still the full SL weight lattice,
while their trace matrix doubles. Thus their source is 2Q_A; root
restriction maps Q_A onto Q_C. The target defining source is 2X_C,
giving the index-two quotient in every rank. The exclusion r=1 matters:
there V and V* coincide and the disjoint-weight argument does not apply.

For E6 the two minuscule weight sets have different nonzero central
characters and are disjoint. Their common span is the full weight
lattice, their trace index is twelve, and hence the source is 12Q_E6.
Root restriction is onto Q_F4. In standard F4 coordinates Q_long is
the even-sum integer lattice D4 and
Q_F4=Z^4+Z*(1/2,1/2,1/2,1/2). It follows directly that
Q_long/(2Q_F4)=(Z/2)^2, giving the second quotient.

For triality the three disjoint eight-dimensional weight sets represent
the three nonzero central classes of Spin8. They span its weight lattice
and have total trace index six. Thus their source is 6Q_D4, which
restricts onto 6Q_G2. The already-established defining G2 source has
the stated unequal short/long factors, leaving a quotient Z/3.

These deductions use the existing moment-lattice source theorem. They
are not a classification of every representation under every folding.
In particular no claim is made that the quotient must always be nonzero
or determined just by the order of the diagram automorphism.

## 3. Exact checks and a useful research path

The [folding pilot](experiments/other_theories_folding_pilot.py) builds
the character-restriction matrix from the simple-root matrices. For
A3/C2, E6/F4 and D4/G2 it independently checks the full branching weight
dictionaries, every Laurent Springer coordinate in (2), and the actual
integer source quotient. It obtains (1,2), (1,1,2,2) and (1,3).
[Report](experiments/other-theories-folding-pilot.json).

This supplies a precise route to pullback comparisons in generalized
cohomology, with a measured integral defect rather than an assumed
identification. It could also inform a later twisted-affine program.
The [pointwise loop extension](OTHER_THEORIES_AFFINE.md#7-an-actual-pointwise-loop-map-already-exists)
does preserve the corresponding twisted loops. A central/rotation
extension, affine-flag source and completion remain separate questions;
finite compatibility alone does not provide their image theory.
It also does not supply the still-missing geometric explanation of the
earlier E7 even-block multiplicity formula.
