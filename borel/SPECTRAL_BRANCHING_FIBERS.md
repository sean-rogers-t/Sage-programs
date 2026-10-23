# A spectral criterion for arbitrarily many cubic couplings

Date: 2026-10-23. Written characteristic-zero deductions. Published
linear algebra and quotient geometry are inputs; novelty is unverified.
This supplies a sufficient branching criterion, not an exceptional-group
classification or a new assertion that every branching graph qualifies.

## 1. What the criterion reads from a character

Let G be complex simple simply connected, rho an actual nontrivial
representation, and normalize the Springer map by

    B_rho(theta_rho(g),X)=tr(rho(g) d rho(X)).

For a Levi K use the finite central cover [K,K] times T -> K, fundamental
characters of [K,K], and the full central Laurent coefficients in C[T].
For h in z(k), the target coordinates are c_j=B_rho(h,H_j), where H_j
are the chosen central cocharacters. We study theta_rho^(-1)(h) intersect
G_reg, for generic h with Z_G(h)=K, as an actual scheme.

Here is a directly checkable sufficient shape condition on chi_rho|K.
Choose one fundamental character x, typically the character of an A1
factor. Apart from a polynomial g(x,z), every occurring affine monomial
must contain exactly two distinct other variables, each to exponent one,
and at most one x. Variables absent from the polynomial are free w.
Require g to have degree two in x with invertible leading coefficient.
Draw an edge between the two other variables whenever such a monomial
occurs, retaining its entire Laurent coefficient and x exponent. Require
a bipartition of the active vertices into equally sized sets u and v.
These conditions are equivalent to a presentation

    chi=g(x,z)+u^t M(x,z)v,
    g=a(z)x^2+b(z)x+q(z),  M=A(z)+x B(z),                   (1)

with a invertible and u,v each of length m. All variables in (1) are
fundamental-character coordinates, or specified invertible changes of
them. The w coordinates are absent and remain free. No graph edge may
be inferred after setting z=1: central-charge cancellation loses data.

For simply laced adjoint maps, the previous
[norm theorem](BRANCHING_NORM_FIBERS.md#3-the-adjoint-cubic-obstruction-has-a-restricted-shape)
already bounds external terms by degree three. If there is exactly one
A1 factor, every cubic involves its character x. This explains why (1)
is a natural test. It does NOT force the bipartition, absence of linear
leaf terms, balanced sides, or invertibility of B. Odd cycles, other
affine squares, and unmatched linear terms are explicit failures of this
particular test. Failure does not imply singular Springer fibers.

The size m is unrestricted. The theorem below also handles repeated
semisimple eigenvalues, whose multiplicities determine excess dimension.

## 2. Spectral theorem and its exact hypotheses

Put X=V(partial_x chi, partial_u chi, partial_v chi) in
A^(1+2m) times T, and mu_j=D_j chi on X, D_j=z_j partial_z_j.
Free w will be reattached at the end. Let U be a nonempty open in T
where B is invertible and write

    L=-B^(-1)A,  M=B(xI-L).

Assume that after a finite etale cover of U the vector bundle splits
into eigensubbundles for L, with distinct regular eigenvalues lambda_j
and positive constant ranks d_j, sum d_j=m. Equivalently for the proof,
one may supply this splitting directly. Over the central function field,
a squarefree minimal polynomial ensures such a splitting after shrinking
U and passing to a finite etale cover. A squarefree CHARACTERISTIC
polynomial is sufficient but unnecessarily restrictive: d_j may exceed
one. A repeated eigenvalue with a Jordan block is not covered.

Set, on the spectral cover,

    x_0=-b/(2a), R_0=q-b^2/(4a),
    R_j=g(lambda_j,z), eta_j=g_x(lambda_j,z)
         =2a lambda_j+b.                                  (2)

Assume eta_j is not identically zero on any connected spectral component.
Finally impose the essential boundary hypothesis:

    closure(mu(X over T minus U)) is a proper subset of A^s. (3)

This refers to the ORIGINAL critical scheme and ORIGINAL target map.
Removing a bad central parameter divisor alone does not establish (3).

**Theorem 1 (spectral branching criterion).** Under these hypotheses
there is a dense open of targets c for which the entire critical scheme
mu^(-1)(c) is smooth, possibly empty. On the central cover its pieces are
described by the following scalar maps and affine quadrics:

| Branch | Scalar target equation | Affine part over each scalar solution | Dimension |
| --- | --- | --- | --- |
| Ordinary branch | D R_0=c | x=x_0, u=v=0 | 0 |
| Eigenvalue branch j | D R_j=c | x=lambda_j, u_j^t v_j=-eta_j; all other blocks zero | 2d_j-1 |

Only scalar maps dominant onto A^s contribute generically. Their generic
fibers are finite reduced. In particular the spectral rank d_j is NOT a
generic degree or an orbit count. Free coordinates add dim w throughout.
The description is etale-local in the central parameters; eigenspaces
need not have global frames and eigenvalue labels may have monodromy.

For an actual representation satisfying the hypotheses, the generic
regular Springer fiber is smooth, possibly empty, with contributing
component dimensions

    dim K-rank G+dim w                 (ordinary branch),
    dim K-rank G+dim w+2d_j-1          (eigenvalue branch).  (4)

It need not be pure-dimensional when several branches contribute.
Smoothness and these local dimensions descend through the finite central
cover. The table does not assert a global quadric decomposition downstairs.

### Proof, including scheme structure and collisions

Absorb B into u, use the eigensplitting, and choose local frames. The full
potential becomes exactly

    g(x,z)+sum_j (x-lambda_j(z)) u_j^t v_j.                 (5)

All changes are invertible on the stated chart. Derivatives in z of a
base-dependent affine change differ by linear combinations of the affine
gradient equations. Thus the induced mu on their quotient ring is
unchanged, provided the WHOLE potential is transformed.

The relative equations of (5) are

    (x-lambda_j)u_j=(x-lambda_j)v_j=0,
    g_x(x,z)+sum_j u_j^t v_j=0.                            (6)

Where all x-lambda_j are units, (6) sets u=v=0 and x=x_0
scheme-theoretically, and mu=D R_0.

Near x=lambda_j all other factors x-lambda_k are units. Eliminate their
variables. Where eta_j is invertible, the remaining equations force
x=lambda_j: at every underlying point u_j^t v_j=-eta_j is nonzero, so
some u or v coordinate is invertible locally, and its equation in (6)
sets x-lambda_j=0. These localizations cover the entire scheme there.
The remaining ideal is precisely

    (x-lambda_j, u_j^t v_j+eta_j),                         (7)

with the eliminated variables zero. This quadric is smooth of dimension
2d_j-1 because its level is nonzero. The original central derivatives are

    mu=D_z g|_(x=lambda_j) - (D lambda_j) u_j^t v_j
      =D_z g|_(x=lambda_j)+eta_j D lambda_j =D R_j.         (8)

The chain-rule term in (8) is necessary; differentiating g at fixed x
would give a different answer.

It remains to control eta_j=0. At such critical points x=lambda_j,
the remaining pairing is zero, and (8) still gives mu=D R_j. Therefore
the target image of this ENTIRE locus factors through the proper divisor
eta_j=0 in the s-dimensional spectral base. Its image has dimension at
most s-1, even when the quadric cone over it has large dimension. This
also covers intersections with the ordinary branch. Exclude its image
closure; no radical has been taken in the good-locus equations (6)-(7).

The remaining scalar maps have smooth s-dimensional sources, so
characteristic-zero generic smoothness gives finite reduced generic
fibers, or emptiness for nondominant components. Apply this to the finite
etale spectral cover, exclude (3), and descend. This proves the first
assertion, including all boundary support and possible nilpotents there.
The [regular-quotient presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md)
adds dim K-rank G and the free w, proving (4). Generic nonemptiness is
equivalent to at least one displayed scalar branch being dominant.

## 3. A sharp obstruction: the graph alone is insufficient

The persistent-collision hypothesis cannot simply be omitted. Compare
the two Laurent potentials, z in C*,

    chi_good=x^2+x+z+(x+1)uv,
    chi_bad =x^2+2x+z+(x+1)uv.                              (9)

They have exactly the same monomial support and the same coupling graph.
Both have B=1 and the simple eigenvalue lambda=-1. For the first,
eta=-1; the generic fibers consist of a point and the smooth curve uv=1.
For the second eta=0 identically. Put X=x+1; the critical ideal at target
c!=0 is

    (2X+uv, Xv, Xu, z-c),

and eliminating X gives

    (u^2 v, u v^2, z-c).                                  (10)

The nonzero class uv is nilpotent, since (uv)^2=0. Thus EVERY such target
fiber is nonreduced. This is an exact ideal calculation, not a sampled
singular point. It disproves a support-only criterion even for a
one-dimensional pencil with positive coefficients in (9).

There is also an honest representation realization, with a crucial
restriction on the group. Let H=SL2 times SL3 times C*, let V and W be
the defining modules of SL2 and SL3, and let C_1 be the one-dimensional
weight-one module of C*. Inflate all other summands trivially across
the central torus and set

    E=(V tensor V) direct_sum V^(direct_sum 2) direct_sum C_1
      direct_sum (W tensor W*) direct_sum (V tensor W tensor W*).
                                                               (11)

Its dimension is 4+4+1+9+18=36 and its character is exactly chi_bad,
where x=chi_V, u=chi_W, v=chi_(W*). The differential of this representation
is faithful and its trace form is nondegenerate. In simple-coroot
coordinates for SL2 and SL3, followed by the weight-one torus generator,
the trace Gram matrix has blocks 30, [[36,-18],[-18,36]], and 1.

**Corollary 2 (a reductive representation obstruction).** For this
representation of H, the regular part of the generalized Cayley fiber
over h=cZ is nonreduced for every c!=0, where Z generates Lie(C*).

Indeed the central derivative is z, so z=c. If omega is a primitive
cube root of unity, take the SL2 eigenvalues (omega,omega^(-1)) and
the SL3 eigenvalues (1,omega,omega^(-1)). This is regular semisimple,
with x=-1 and u=v=0. The character quotient is smooth there; its
relative critical ideal is precisely (10). The nonzero nilpotent uv
survives locally at the origin and under this smooth pullback. The
regular group fiber has local dimension 8+1=9 and tangent dimension
8+2=10 at this point. A direct weight-trace calculation confirms
theta_E(g)=(0,0,0,c) and differential rank two.

H is reductive but NOT simple, and E is not its adjoint representation.
We have NOT realized this as a Levi restriction of an ambient simple
group representation. Thus it does not refute the all-simple-adjoint
Springer conjecture. It does show that the coefficient relation is a
genuine generalized-Cayley obstruction, not merely an artificial
polynomial pathology. A root-theoretic extension must prevent it, or
prove its target image nondominant, in the intended simple-group scope.

Repeated semisimple eigenvalues are a different phenomenon. For every
m>=1 the model

    chi=x^2+z+(x+1) sum_(i=1)^m u_i v_i                    (12)

has, for every c!=0, a point and the smooth affine quadric
sum u_i v_i=2 of dimension 2m-1, both over z=c. Here L=-I_m,
eta=-2, so multiplicity creates smooth excess, not a collision. Formula
(12) is an arbitrary-m algebraic model, not a claimed family of actual
exceptional representations. Its proof is (6)-(7), with no rank sweep.

## 4. Recover the E6 suspension and check its boundary

For the adjoint E6 Levi with retained nodes {1,2,3,5,6}, use x for the
A1 character and put S=z^3+z^(-3). The charge-preserving formula is

    g=x^2+Sx-2,
    M=[[1, z^(-2)+zx], [z^2+z^(-1)x, 1]].                 (13)

The pair variables are (a,u) and (b,v). Thus

    det B=-1, det M=-x(x+S),
    eigenvalues of L: 0,-S,
    eta: S,-S, R_0=-2-S^2/4, R_1=R_2=-2.                 (14)

On U={S!=0} the spectral theorem applies. Both eigenvalue branches have
target mu=0, so neither contributes generically. Only the ordinary
branch survives, giving zero excess. This explains the previously
proved E6 dimension 14 through spectral critical values.

Check U's complement rather than discarding it. The existing invertible
completion of the pair variables gives

    chi=a'b'+P(1-uv)-2,  P=x(x+S).

On S=0, the relative equations include 2x(1-uv)=0 and the central
derivative is x (D S)(1-uv). Hence mu=0 in the boundary critical
quotient ring itself. This proves (3), even before passing to support.
No new E6 case or orbit count is asserted; this is a verified application
and explanation of the earlier result. The trace normalization remains
c=144t for h=t(2,3,4,6,4,2) in simple-coroot coordinates.

The other two E6 cubic refinements have additional A1 squares. They fail
the displayed bipartite-pair shape and retain their
[product/Fricke proofs](CUBIC_LEVI_REFINEMENTS.md). A correct recognizer
must report that limitation rather than force the spectral theorem.

## 5. Verification, published inputs, and the remaining question

The [Sage verifier](tools/verify_spectral_branching_fibers.py) reads the
full charge records and recognizes (13), rejects the two incompatible
shapes, checks spectral identities, checks exact localized critical
ideals, and preserves the nilpotent in (10). The
[native SymPy audit](tools/audit_spectral_branching_fibers.py) independently
checks the pencil, its eigenspaces, chain-rule target derivative, boundary
identity and the good/bad ideals. Both independently construct the 36
representation weights and check the trace Gram matrix, torus value,
and differential at the regular semisimple point of Corollary 2. Reports:
[verification](spectral-branching-verification.json),
[audit](spectral-branching-audit.json). These corroborate the written
arbitrary-size proof; they do not certify novelty or all embeddings.

Published inputs are the trace construction of
[Kostant--Michor](https://www.mat.univie.ac.at/~michor/kostant.pdf),
standard semisimple linear algebra, regular quotient geometry, and
characteristic-zero generic smoothness. Kostant's
[Levi-root-space theorem](https://arxiv.org/abs/0711.2809) supplies
branching structure, but does not by itself establish the spectral
noncollision and boundary hypotheses here. Targeted searches on
2026-10-23 found these sources, not a verified priority claim for this
criterion. The linear-algebra reductions are elementary deductions;
we have not resolved a named published open problem.

The next mathematical gate is now specific: for actual simply laced
adjoint restrictions having the bipartite shape, do root and central-charge
identities force spectral noncollision (eta_j not identically zero),
or at least nondominance of its target image? The same question is needed
on nonsemisimple/discriminant boundary strata. Proving that implication
would turn this conditional branching criterion into a group-theoretic
family theorem. The support-only counterexample explains exactly why
coefficients and embeddings, rather than more small examples, matter.
