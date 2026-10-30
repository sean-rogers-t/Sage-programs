# Coupled character constraints and adjoint Springer fibers

Date: 2026-10-23. Characteristic-zero proofs with separate Sage and native
audits. The group application is E8, adjoint representation, Bourbaki node
4 omitted. General lemmas below concern specified potentials, not all
representations or all Levis. Novelty is unverified.

## 1. Exact simultaneous elimination

Work over R=C[z_1^(+/-1),...,z_s^(+/-1)]. Let a,tau be vectors of length
k and y further affine coordinates. Suppose the ACTUAL character is

    chi=tau^t(C(z)a-f(y,z))+G(a,y,z),                       (1)

where G is independent of tau and det C is a unit of R. For a central
target c, use all affine derivatives and D_j chi-c_j, where
D_j=z_j partial/partial z_j holds the original affine coordinates fixed.

**Theorem 1.** Put alpha=C^(-1)f. The complete critical scheme is
isomorphic, over the entire central torus and target space, to that of
G(alpha,y,z), with inverse formulas

    a=alpha, tau=-C^(-t)G_a(alpha,y,z).                    (2)

This is compatible with any finite central action preserving the original
character and the graph (2); the induced action on remaining coordinates
must be retained. It preserves nilpotents and commutes with base change.

Proof. The tau equations give a=alpha, and the a equations give tau in
(2). For a remaining derivation d, their substitution into d chi gives

    dG + tau^t((dC)alpha-df)
      = dG + G_a^t C^(-1)(df-(dC)alpha)
      = d(G(alpha,y,z)).

Thus both quotient rings, not just their closed points, are identical.
For d=D_j this retains the derivative of C and the target c_j. Graph
equivariance makes the isomorphism descend. If det C is not a Laurent
unit, this argument applies only off its zero locus; the critical scheme
on that locus requires a separate argument. In particular a matrix
invertible over C(z) does not certify global elimination.

### A branching criterion with an explicit certificate

In the central-charge expansion of a restricted character, choose tau
so each occurs only linearly and no monomial contains two tau variables.
Its coefficients must be affine linear in chosen paired variables a:
the a coefficients form C(z), with all remaining terms equal to -f.
Require det C to be a nonzero constant times a Laurent monomial. This
is a sufficient condition readable from the full branching coefficients.
Apply it in stages if needed, retaining every residual variable.

This condition is stronger than a matching in a monomial support graph:
coefficients can cancel in the determinant, and the remaining potential
still controls smoothness. The research helper rejects unsupported shapes
and nonunit determinants; it does not discard their boundary loci.

**Corollary 1.1 (unique-matching criterion).** Suppose every nonzero
entry of C is a nonzero scalar times a single central Laurent monomial.
Form its bipartite graph, with multiplier rows and paired-variable
columns. If that graph has a unique perfect matching, Theorem 1 applies.
Indeed the determinant expansion then contains precisely one nonzero
term, a Laurent unit. In particular, a forest with a perfect matching
suffices: two distinct perfect matchings would have an alternating cycle.
The condition also applies at each stage of a recursive elimination.
One must still check that the original multiplier coefficients have the
affine form in (1). A matching in the full cubic support graph alone
does not check that condition or determine the residual potential.

This is a root/charge criterion whenever the branching modules identify
these single-monomial coefficients. E7 has one such matched edge; E8
has two disjoint edges, with charges -2 and +2. It does not assume their
residual potentials are the same.

## 2. A general one-variable quotient theorem

Suppose exact eliminations leave

    psi=A(Y)+x B(Y)(z^m+z^(-m)), Y=x^2, m>=1.             (3)

Let A,B be complex polynomials, A nonconstant, B nonzero and squarefree,
B(0)!=0. Set C(Y)=B(Y)+2Y B'(Y), and assume gcd(A',C)=1. Take the finite
action z -> zeta z, x -> zeta^(-m)x, zeta^(2m)=1. It preserves (3).
There are no additional residual coordinates in this statement.

**Theorem 2.** For EVERY c!=0, the quotient of the critical scheme
(psi_x,D_z psi-c) by this action is

    Spec C[Y]/(F_c(Y)),
    F_c=4m^2 Y B(Y)^2 (Y A'(Y)^2-C(Y)^2)-c^2 C(Y)^2.     (4)

If a=deg A and b=deg B, its length is

    d=max(2a+2b,4b+1).                                   (5)

The family over c!=0 is finite flat of rank d and is generically etale.
Its exceptional targets are exactly the nonzero roots of disc_Y(F_c).
For the group consequence assume G is connected complex simply connected
semisimple, rho is an actual representation with nondegenerate trace
form, K is a Levi with one-dimensional center, and the nonzero central
target h has centralizer K. Require (3) to come from its restricted
character by the exact eliminations, with the displayed mu_(2m) as the
full residual central action and c=B_rho(h,H) for the chosen central
cocharacter. The general regular-quotient theorem then gives generically
d regular K-orbits, each smooth of dimension dim K-rank G. This is not
the degree of G -> Lie(G).

Proof. Write S=z^m+z^(-m), T=z^m-z^(-m). The equations are

    2x A'(Y)+C(Y)S=0,   m x B(Y)T=c,   S^2-T^2=4.       (6)

In the critical ring c!=0 makes x,B,T units. If C vanished, the first
equation would force A' to vanish. Bezout's identity for A',C therefore
makes C a unit in this ring. Substitute

    S=-2x A'/C, T=c/(m x B)

into the last equation to obtain (4). No radical was taken.

Conversely, (4) has no roots at Y=0, at a zero of B, or at a zero of C.
At Y=0 this follows from B(0)!=0 and c!=0. At a zero of B, C=2YB'
is nonzero by squarefreeness. At a zero of C, the other factors are
nonzero by the hypotheses. Consequently Y,B,C are units in its entire
Artinian quotient, including nonreduced fibers. Adjoin x with x^2=Y;
form S,T as above and k=(S+T)/2. Then k^(-1)=(S-T)/2 and adjoin z
with z^m=k. These two extensions form a finite etale torsor of degree
2m, with precisely the stated action. This reconstructs (6) and proves
the scheme quotient, including its nilpotents.

The two possible highest degrees of the numerator in (4) are 2a+2b
and 4b+1. They have different parity, so no cancellation is possible;
both exceed deg C^2=2b. Thus (5) holds at every c, with constant nonzero
leading coefficient. The rational function

    R(Y)=4m^2 Y B^2(Y A'^2-C^2)/C^2                       (7)

has coprime numerator and denominator by the same root exclusions and
has degree d. It is nonconstant and separable in characteristic zero.
The pullback of its generic fiber along c -> c^2 is reduced. This proves
generic etaleness and the discriminant statement.

This supplies a polynomial family of sufficient conditions. It does
not assert that arbitrary choices of A,B,m are representation characters.
For actual branching one must prove the reduction and identify the full
central action, not just a subgroup of it.

## 3. The E8 adjoint application

Let G be complex simply connected E8, rho=Ad, K the standard Levi omitting
Bourbaki node 4. Its derived group is SL3 times SL2 times SL5. Fix

    H=(10,15,20,30,24,18,12,6), h=tH,
    B_Ad(theta(g),X)=tr(Ad(g)ad(X)), c=B_Ad(h,H)=1800t.    (8)

Coordinates of H are in the simple-coroot basis. For t!=0, Z_G(h)=K.
Here dim K=36, rank G=8, and the multiplication cover has kernel mu30.
Write (a,b) for SL3 fundamental characters, x for SL2, and (p,q,r,s)
for SL5 in order along nodes 5,6,7,8. The fundamental coordinates in
retained-node order are (a,x,b,p,q,r,s).

The central-charge character is

    chi=ab+x^2+ps-2
       +z a x s+z^2 b r+z^3 x q+z^4 a p+z^5 x b+z^6 s
       +z^-1 x b p+z^-2 a q+z^-3 x r
       +z^-4 b s+z^-5 a x+z^-6 p.                       (9)

The positive-charge constituents, with their dimensions, are

| Charge | Character | Dimension |
| --- | --- | --- |
| 1 | a x s | 30 |
| 2 | b r | 30 |
| 3 | x q | 20 |
| 4 | a p | 15 |
| 5 | x b | 6 |
| 6 | s | 5 |

Negative charges are dual; the zero piece has dimension 36. Their sum
is 36+2*(30+30+20+15+6+5)=248. Root-weight reconstruction checks every
charge, not merely this dimension count.

### Two constraints and a paired elimination

The q,r derivatives are z^-2 a+z^3 x and z^2 b+z^-3 x.
Theorem 1 applies with C=diag(z^-2,z^2), determinant 1, giving

    a=-z^5 x, b=-z^-5 x,
    q=-z^2(b+z x s+z^4 p+z^-5 x),
    r=-z^-2(a+z^-1 x p+z^-4 s+z^5 x).                 (10)

Set

    U=z^6(x^2-1)+z^-9 x,
    V=z^-6(x^2-1)+z^9 x.

After (10), the potential is ps-pV-sU-2. Replacing p by p-U and
s by s-V is invertible over the entire torus. Its critical equations
eliminate the paired coordinates, leaving -UV-2, namely

    psi=-(x^2-1)(x^2+(z^15+z^-15)x)-3.                 (11)

Thus Theorem 2 applies with A=-Y^2+Y-3, B=1-Y, C=1-3Y, m=15.
All its coprimality hypotheses hold. There is no excluded source
determinant or unexamined boundary component on c!=0.

### Central descent and the complete nonzero-target answer

The central kernel acts on (a,x,b,p,q,r,s) with weights
(-10,-15,-20,-24,-18,-12,-6) modulo 30, and on z with weight 1.
These follow by canceling H(zeta) with the derived torus. The original
character and elimination graph are invariant. On the residual pair
the action is exactly (z,x)->(zeta z,zeta^-15 x), as in Theorem 2.

**Theorem 3.** For every t!=0 the critical quotient governing
theta_Ad^(-1)(tH) intersect G_reg is the length-six scheme

    Spec C[Y]/(900Y(1-Y)^2(4Y-1)(Y^2-3Y+1)
                     -c^2(1-3Y)^2), c=1800t.           (12)

Put

    Delta(k)=4096k^5+4485075k^4+1691730000k^3
             +292161500000k^2+23287500000000k
             +683437500000000.                          (13)

For c!=0 and Delta(c^2)!=0 there are exactly SIX regular K-orbits,
each smooth of dimension 28. Delta is squarefree and Delta(0)!=0,
so there are precisely ten exceptional nonzero complex targets. At
each of them (12) has one double point and four simple points. The
regular-source fiber has five underlying K-orbits there; along one
orbit its completed local ring is

    C[[epsilon,u_1,...,u_28]]/(epsilon^2),               (14)

and the other four orbits are smooth. In particular the generic fiber
is nonempty and smooth, but not every nonzero central fiber is reduced.

Proof. Formula (12) is (4), since
Y(1-2Y)^2-(1-3Y)^2=(4Y-1)(Y^2-3Y+1).
Its Y-discriminant is

    -30611001600000000 c^2 Delta(c^2).                   (15)

Exact polynomial arithmetic verifies (15) and gcd(Delta,Delta')=1.
At a simple zero of the discriminant a constant-degree polynomial has
exactly one double root and all other roots simple: two double roots or
a root of multiplicity at least three have discriminant vanishing order
at least two. This proves the assertion about (12).

The [regular-quotient theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md) identifies
the actual regular-source fiber as a smooth pullback of this critical
scheme, of relative dimension dim K-rank G=28. Every quotient point has
one regular K-conjugacy class of lifts. Equivariance puts Z_G(g) inside
K for a lift in the fiber, so regularity in K and in G agree. This proves
the orbit assertions and the local ring (14). It makes no claim about
the target zero, nonregular source points, or the generic degree of the
full map G -> Lie(G).

## 4. Structural scope, verification, and next gate

The earlier [scope theorem](ADJOINT_PENCIL_SCOPE.md) isolates E6/E7/E8
omitting node 4 as the only raw pair-pencil candidates under its stated
all-A and distinguished-A1 assumptions. It does NOT classify every
Levi admitting multiplier elimination. E6 already has the pair pencil;
E7 has one forced linear multiplier; E8 has the two simultaneous
constraints above. Thus this completes those three candidates by exact
methods, without turning a three-member list into an all-Levi theorem.

The genuinely general statements here are Theorems 1 and 2. The former
is elementary constrained critical-point elimination with explicit
global hypotheses. The latter gives an exact quotient and generic
orbit-count formula once a representation satisfies its reduction.
GEO-54 and GEO-59 remain the separate E6/E7 proofs; Theorem 2 is not
artificially imposed on their different residual forms.

The [Sage verifier](tools/verify_multiplier_branching_fibers.py) reconstructs
the actual character, checks the original affine and central derivatives
before substitution, verifies both directions of the residual quotient
over Q(c)[Y]/(F), checks the discriminant, and evaluates the original
adjoint-root differential at a split regular point over a finite field.
The recorded point is over F_32003, with z=343 and c=1389: the Cartan
rank is 8, the root-direction rank is 212, and the full kernel dimension
is 248-220=28. This only corroborates the characteristic-zero proof.
The [native audit](tools/audit_multiplier_branching_fibers.py) independently
rebuilds weights from Cartan reflections and subset models, checks the
elimination and target identities, central action, discriminant and root
differential. Its coordinate-ring checks retain the universal target,
including ramified fibers. Reports:
[verification](multiplier-branching-verification.json),
[audit](multiplier-branching-audit.json).

### Published input and literature status

- [Kostant--Michor, Proposition 2.4](https://www.mat.univie.ac.at/~michor/kostant.pdf)
  provides the trace-gradient identity and equivariance used here.
- [Kostant, Theorem 0.1](https://arxiv.org/pdf/0711.2809) gives the
  irreducibility of nonzero Levi central-root spaces, supporting the
  charge decomposition. Our weight checks identify the actual modules.
- [Steinberg, Regular elements of semi-simple algebraic groups](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf)
  supplies the regular conjugacy quotient input; the exact use is explained
  in the existing regular-quotient note.

Elimination, separability of rational functions in characteristic zero,
and discriminant criteria are standard algebra. The explicit E8 scheme
and its application are deductions written here. Targeted searches on
2026-10-23 for E8 adjoint Cayley fibers and the A4+A2+A1 Levi did not
identify this formula, but do not establish novelty or a complete survey.

**Next gate:** classify the central-charge constraint systems that are
unit-solvable, and the residual forms they leave, for the remaining
all-A exceptional Levis. First derive a diagram/charge criterion or an
explicit obstruction; then use a bounded classification to verify its
scope. A successful next result must explain a family or isolate the
precise residual obstruction. More unstructured point calculations are
not the next objective. The all-simple-adjoint generic-smoothness
question remains unresolved here.
