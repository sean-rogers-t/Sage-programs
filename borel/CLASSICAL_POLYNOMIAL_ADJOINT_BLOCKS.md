# Classical block classification for polynomial-adjoint Springer maps

27 October 2026. GEO-97. Written characteristic-zero deduction, with
separate symbolic and full-root checks. This is an internal proof audit,
not external review or a literature-priority claim. No public API change.

## 1. What was already proved, and what is added

The all-rank ADJOINT nonemptiness and dimension classification was already
completed in [the classical theorem](CLASSICAL_ADJOINT_FIBERS.md), GEO-49.
The recent literature report incorrectly described this classification
itself as unfinished. The new question is its extension to arbitrary
actual characters polynomial in the adjoint character. The componentwise
[composition theorem](CHARACTER_COMPOSITION_FIBERS.md) supplies smoothness
and survival once constant critical values have been classified.

Let G be simple, complex, simply connected and classical, K a proper Levi,
chi=chi_Ad, and chi_sigma=f(chi) an ACTUAL representation character, with
f a nonconstant polynomial and f'(dim G) nonzero. Use trace normalization:

    theta_sigma = f'(chi)/f'(dim G) theta_Ad.

We consider the actual open fiber theta_sigma^(-1)(h) intersect G_reg,
for generic h in Z_K^o={h in z(k): Z_G(h)=K}. A generic assertion means
one dense target open for the fixed G,K,f. It does not include special
targets, nonregular source points, or full-factor zero targets.

### Theorem: complete generic nonemptiness and dimension

In A, with G=SL_N and block sizes b_i, the fiber is nonempty exactly when
all b_i<=2 and at least one block is a singleton. Its pure dimension is
sum b_i(b_i-1). This is the previously proved polynomial-character result.

In B/C/D write the defining-matrix Levi pattern as

    product GL_(b_i) times H_m,  sum b_i+m=n,

with H_m=SO_(2m+1), Sp_(2m), SO_(2m), respectively. In D absorb SO_2
into GL_1 blocks, so m=0 or m>=2, and require n>=4. In B,C require n>=2.
For Spin the actual Levi is the inverse image, not an identification
with this matrix product. The finite central isogeny is retained.

1. If any b_i>4, or m>2 in B, m>1 in C, m>3 in D, the regular fiber is
   empty for every centralizer-exact target, for every such f.
2. If every b_i=4, let q>=1 be their number and use the following SET V.
   The generic fiber is nonempty if and only if f'(v)!=0 for at least
   one v in V.

| Type | m | V: adjoint-character values on dominant components |
| --- | --- | --- |
| B | 0 | {-q/(2q-1)} |
| B | 1 | {-(3q-1)/(2q-1), (2q-1)/(4q)} |
| B | 2 | {-(2q+1)/(4q)} |
| C | 0 | {0} |
| C | 1 | {-1} |
| D | 0 | {0} |
| D | 2 | {-2, -1/q} |
| D | 3 | {-1} |

3. For the two additional patterns B_3: GL_1 times SO_5 and
   D_4: GL_1 times SO_6, generic nonemptiness is equivalent to f'(0)!=0.
4. Every other pattern satisfying the bounds has a dominant component
   on which chi is nonconstant. Its generic regular fiber is therefore
   nonempty for EVERY admissible f.

Whenever nonempty, the actual generic fiber is smooth and pure of dimension

    sum b_i(b_i-1)+2m^2          in B,C,
    sum b_i(b_i-1)+2m(m-1)       in D.

In particular the transverse excess is ZERO throughout this polynomial-
adjoint classical family. This conclusion uses GEO-49 and the dimension-
preservation part of GEO-95; it is not a new proof of their smoothness.

### A concrete consequence: pure adjoint tensor powers

For f(t)=t^k, k>=2, start with an adjoint-nonempty proper Levi. Its generic
regular fiber is erased exactly in these patterns:

* B_3: GL_1 times SO_5;
* C_(4q): GL_4^q with no residual symplectic factor;
* D_(4q): GL_4^q with no residual orthogonal factor;
* D_4: GL_1 times SO_6.

All other adjoint-nonempty classical patterns remain nonempty, including
all the type A patterns. No assertion about arbitrary irreducible
representations or minimum Cayley degree follows.

## 2. One reduction for every classical block pattern

Set epsilon=-1 in B,D and epsilon=+1 in C. Write T=tr_V(g), S2=tr_V(g^2).
The published tensor-square character identities give

    chi=(T^2+epsilon*S2)/2.

On a GL_b block use elementary symmetric coordinates e_1,...,e_(b-1)
and determinant d in C*. These are quotient coordinates; a finite central
cover gives the usual fixed-determinant coordinates, with no change to
the arguments below. Put s=d+d^-1, tau=-epsilon*T. The full relative
critical equations imply the following, with no radical taken to derive
the equations:

| b | Relative equations | Defining trace contribution | Square-trace contribution |
| --- | --- | --- | --- |
| 1 | none | s | s^2-2 |
| 2 | (1+d^-2)e_1=tau(1+d^-1) | tau F(d) | T^2 F(d)-2s |
| 3 | e_1=tau+d^-1, e_2/d=tau+d | 2tau+s | 2T^2-(s^2-2) |
| 4 | d=-1, e_1=e_3/d=tau; e_2 free | 2tau | 2T^2 |

The b=2 entries in the last two columns assume s!=0, and

    F(d)=(d+1)^2/(d^2+1)=1+2/s.

For example Newton identities give tr M^-2=(e_(b-1)/d)^2-2e_(b-2)/d.
Differentiating the FULL chi, including its linear terms, gives the table.
The b=4 equation is the derivative in e_2; its determinant equation has
a simple zero at d=-1. Thus e_2 is a genuine remaining coordinate, not
silently discarded by division. Blocks above four are excluded by GEO-49.

For the allowed residual factors, the zero-target residual equation is

    (g^2-1)(g^2-tau*g+1)=0.

The regular residual spectral branches give

    R=tr(residual)=-epsilon*c*T+b,   R2=c*T^2+d0,

with the following triples (c,b,d0):

| Type | m | Triples |
| --- | --- | --- |
| B | 0 | (0,1,1) |
| B | 1 | (0,-1,3), (1,1,-1) |
| B | 2 | (1,-1,1) |
| C | 0 | (0,0,0) |
| C | 1 | (1,0,-2) |
| D | 0 | (0,0,0) |
| D | 2 | (0,0,4), (1,2,0), (1,-2,0) |
| D | 3 | (1,0,2) |

These follow by using self-reciprocal eigenvalues +/-1 and at most one
reciprocal pair with sum tau. They agree with the independently proved
residual trace table in GEO-49. Work first away from T=+/-2, where the
residual quotient branches have locally smooth parametrizations by T
(possibly with finite choices). Section 5 explains why the excluded
collision locus cannot dominate. Central covers may split branches;
they do not change the constants or generic dimension.

There is also a direct quotient-coordinate check of these local scheme
assertions. Let L be the defining trace of all GL blocks and Phi_H its
adjoint tensor-square contribution. The full character is
Phi_H+L*R+chi_res. In residual fundamental characters one has

| Residual type | R | chi_res |
| --- | --- | --- |
| B1, coordinate u | u^2-1 | u^2-1 |
| B2, coordinates u,v (spin, vector) | v | u^2-v-1 |
| C1, coordinate u | u | u^2-1 |
| D2, coordinates u,v | uv | u^2+v^2-2 |
| D3=A3, coordinates u,v,w | v | uw-1 |

Differentiation gives B2: u=0,L=1; C1: u=-L/2; D3: u=w=0,L=0.
These are smooth equations, with the stated remaining T coordinate.
For B1 the equation is 2(L+1)u=0. Its branches u=0 and L=-1 meet
only at T=-2. For D2 the equations are 2u+Lv=2v+Lu=0. The branch
u=v=0 is regular for L!=+/-2; the other branches have L=+/-2 and
u=+/-v, with v nonzero away from T=+/-2. Thus each generic residual
branch really is smooth after the justified unit eliminations, not merely
a list of reduced spectral points. These standard low-rank character
identities also explain the two special constants directly: a singleton
H has exterior-square character 1; B2 forces L=1,chi_res=-v-1, giving
total 0, while D3 forces L=0,chi_res=-1, again giving total 0.

Let q be the number of four-blocks, p the number of three-blocks, and
k the number of two-blocks. For every non-four determinant put s_i=d_i+d_i^-1.
On the chart with every two-block s_i!=0 define

    A = 1+epsilon*(c+2q+2p+k)+2epsilon*sum_(b_i=2) 1/s_i,
    B = b+sum_(b_i=1 or 3) s_i,
    Q = d0+sum_(b_i=1)(s_i^2-2)
            -sum_(b_i=3)(s_i^2-2)-2sum_(b_i=2)s_i.

Then the single remaining defining-trace equation and the scalar value are

    A*T=B,              chi=A*T^2/2+epsilon*Q/2.             (1)

The only free derived coordinates in addition to T are the q variables
e_2 of four-blocks. Their determinants are fixed at -1. Thus this is a
uniform reduction, not an enumeration of ranks or Levis.

## 3. Ordinary chart A nonzero

Here T=B/A and the critical value is

    phi=B^2/(2A)+epsilon*Q/2.                                (2)

The critical locus is locally smooth on the generic residual branch:
the displayed GL equations eliminate by units, the four-block determinant
equations have nonzero differentials, and the trace equation eliminates T.
There are q free coordinates over the remaining determinant torus.
Consequently its dimension is the central rank.

When all blocks are four-blocks, (2) becomes

    b^2/[2(1+epsilon*(c+2q))]+epsilon*d0/2.                 (3)

Substitution of the residual triples gives exactly V in section 1.
Every such branch dominates: its determinant projection is a point, its
critical-value graph has finite horizontal stabilizer, and the component
is clean. Equivalently each e_2 supplies its independent central target
coordinate. The trace values from (1) avoid +/-2 for every listed q.

For other patterns, (2) is nonconstant in EACH remaining determinant
variable except for one allowed ordinary branch:

    D_5, GL_3 times SO_4, residual triple (0,0,4), phi=-3.

Here is a direct dependence check. For a one- or three-block variable s,
A is independent of s and B=s+B0. Its quadratic coefficient in phi is
1/(2A)+epsilon*sigma/2, with sigma=+1 for b=1 and -1 for b=3; its linear
coefficient is B0/A. Independence requires A=-epsilon*sigma and B0=0.
Any other one/three variable or any two-block prevents these identities.
With only that variable and four-blocks, the residual table forces
epsilon=-1, c=q=0, b=0. The only allowed-rank possibility is the D5
branch just stated. That branch has a constant graph over its entire
one-dimensional determinant torus and hence central image zero.
The two other SO4 residual branches give nonconstant values and dominate.

For a two-block variable s, write A=A0+2epsilon/s, with A0 and B
independent of s. The variable part is B^2*s/[2(A0*s+2epsilon)]-epsilon*s.
Its derivative cannot vanish identically unless A0=0 and B^2=4.
Another two-block prevents A0=0; a one/three-block prevents B^2=4.
With just one two-block and four-blocks the necessary residual conditions
would be epsilon=-1, c=q=0 and b=+/-2, absent from the table.

Finally, dependence in each variable proves dominance here, not just
nonconstancy: phi is invariant under independent inversions d_i->d_i^-1.
Its connected horizontal torus stabilizer is therefore invariant under
every coordinate sign change on Lie algebra. Such a subspace is spanned
by coordinate directions. Nonconstant dependence in each coordinate
excludes every one of those directions. The stabilizer is finite, so
the clean-component image theorem GEO-93 gives dominance. Four-block
coordinates cannot contribute a stabilizer because they are fixed.

## 4. Compatible singular chart A=B=0

If A and B are both nonconstant, they depend on disjoint sets of variables.
Their common locus has codimension two. Adding the free T and four-block
coordinates gives dimension one less than central rank. It cannot dominate.

If A is identically zero, the residual table forces epsilon=-1, c=1,
q=p=k=0: all GL blocks are singletons. The equation is

    sum s_i=-b,      chi=-[d0+sum(s_i^2-2)]/2.              (4)

For at least two singletons the value is nonconstant on every component.
Indeed the determinant-to-s map is finite flat, given by d_i^2-s_i*d_i+1,
and each component lies over the full hyperplane sum s_i=-b. The quadratic
sum is nonconstant on that hyperplane. For B with residual rank two and
D with residual rank three this is the only branch, so the established
adjoint nonemptiness theorem supplies a dominant component. In B with
residual rank one, the ordinary branch already supplies one.

For a single singleton, (4) gives dominant constant components:

    B2, GL1 SO3: value 1;
    B3, GL1 SO5: value 0;
    D4, GL1 SO6: value 0.

The fixed s is -1, 1, 0 respectively, so d!=+/-1 and dB has nonzero
derivative. The remaining T coordinate maps nontrivially to the one central
target; these are clean dominant components. In B2 there is also the
nonconstant ordinary branch. The B3/D4 patterns have no such alternative.
The putative SO4 cases have rank D3 and are outside the stipulated D range.

If B is identically zero and A is nonconstant, all blocks have sizes two
or four, and the residual branch has b=0. Write a0=1+epsilon*(c+2q+k).
The central equation and scalar value are

    a0+2epsilon*sum 1/s_i=0,
    chi=epsilon*d0/2-epsilon*sum s_i.                       (5)

For k>=2 the sum s_i is nonconstant on every component of this equation.
In v_i=1/s_i coordinates it is a reciprocal sum on an affine hyperplane.
Its only constant case for k>=2 is k=2 and sum v_i=0; this would require
a0=0, impossible for the allowed triples when k>=2. The finite flat
determinant cover does not create additional constant components.

For k=1, s=-2epsilon/a0 must be nonzero. The dominant constant components
are as follows; q counts four-blocks:

| Type and residual rank | q range | constant value |
| --- | --- | --- |
| C, m=0 | q>=0 | 1/(q+1) |
| C, m=1 | q>=0 | -(2q+1)/(2q+3) |
| D, m=0 | q>=1 | -1/q |
| D, m=2, triple (0,0,4) | q>=1 | -2-1/q |
| D, m=3 | q>=1 | -1-2/(2q+1) |

In these rows s!=+/-2, so dA is nonzero; at general T!=0 the actual
critical equations are smooth on this component. Its central projection
is finite and its critical value constant, giving a finite stabilizer
and dominance. Every listed pattern ALSO has a nonconstant ordinary
component, so these constants do not obstruct overall nonemptiness.

There is an essential excluded boundary: D5, GL2 SO6 (q=0,m=3) gives
s=-2, d=-1. Here A=-(d+1)^2/(d^2+1) has a double zero. Both dA and
dQ vanish, so its central image is zero, not dominant. It must not be
promoted to a dominant constant component by using only reduced support.
This boundary is checked directly in the original critical equations.

Sections 3-4 give all constant values on dominant classical critical
components, up to finite splitting of the central covers. They do not
count the components, their generic degrees, or individual group orbits.

## 5. Boundaries and actual schemes

No singular two-block denominator is deleted without checking its image.
If a set J of j>=1 two-blocks has s_i=0, then d_i=+/-i and their critical
equations force T=0. Their e_1 coordinates are free, subject to ONE
nontrivial linear trace equation (each coefficient 1+d_i^-1 is nonzero).
All other non-four derived coordinates are determined. Four-block e_2
coordinates remain free, and the residual quotient has finitely many
choices at T=0. If central rank is s, the support dimension is therefore

    (s-q-j)+(q+j-1)=s-1.

Its central image is nondense, including all intersections of these
denominator strata. This dimension statement concerns the image of its
underlying locus, not replacement of the actual scheme by its radical.

The residual collision values T=+/-2 also cannot dominate. At fixed
nonzero T, the quartic spectral lemma of GEO-49 forces every GL block
to have size four on a hypothetical dominant component. The trace equation
would be (1-2q-c)T=b in B,D, or (1+2q+c)T=0 in C. None of the residual
triples admits T=+/-2 with q>=1. Thus all collision strata are nondominant.

The map Y->target retains its ACTUAL structure throughout. We use reduced
components to classify image dimensions and constant functions, then use
the established actual adjoint smoothness theorem to remove the nonsmooth
image. GEO-95 proves that polynomial composition cannot turn any such
nondominant subvariety into a dominant one: (mu,chi) has the same image
dimension as mu by residue rigidity. It also proves f'(chi) is a unit in
every nonzero-target fiber, including its nilpotents. These inputs justify
the final scheme assertion; it is not inferred from reduced calculations.

## 6. Finish the classification

The bounds in section 1 are necessary by the old adjoint spectral theorem:
on a nonzero target, theta_Ad is a nonzero scalar multiple of that target,
so the same centralizer-exact bounds apply even when the scalar varies.

All-four patterns have exactly the dominant values (3). The only other
patterns with no dominant nonconstant component are the B3 and D4 cases
in (4). In every other allowed pattern, section 3 gives a clean dominant
nonconstant ordinary component, or (4) and old adjoint nonemptiness give
a dominant nonconstant component. The extra constants in (5) do not erase
these nonconstant alternatives. Apply the component-survival theorem:
a constant component is lost precisely at a root of f', whereas every
nonconstant dominant component survives with unchanged excess. This proves
the nonemptiness criterion and dimension assertion for all ranks.

The theorem is for actual polynomial-adjoint characters, not only
polynomials with nonnegative coefficients. If one restricts to nonnegative
coefficients, f' is positive at positive real values, so those components
cannot be erased by that narrower class of operations. The nonpositive
values can be erased by the explicit positive-coefficient construction
already proved in GEO-95. No claim that every finite subset of the table
can be selected independently by actual characters is needed here.

## 7. Independent checks, published inputs, and limits

The [native audit](tools/audit_classical_polynomial_blocks.py) derives the
block identities from companion-matrix/Newton data, verifies the symbolic
reduction and all-q constant formulas, and checks the double-zero boundary
in a separate D3=SL4 residual character presentation. The
[Sage check](tools/verify_classical_polynomial_blocks.py) independently
uses full ambient classical roots and exact number fields to check the new
constant branches and their central target equations. Neither imports the
other's formulas or reads its calculated answers as input.

These are independent arithmetic checks, not independent external refereeing
or formal verification of the general geometry. The all-rank proof is the
elimination and dependence argument above, not a finite sweep of examples.

Published inputs: the [Kostant--Michor trace construction and character
chain rule](https://www.mat.univie.ac.at/~michor/kostant.pdf), tensor-square
character/Newton identities, classical regular Jordan forms, and the
regular adjoint quotient. Repository inputs: GEO-49's full classical
classification, GEO-93's clean-component stabilizer theorem, and GEO-95's
actual-scheme polynomial-composition theorem. The table of constant
components and the assembled polynomial-adjoint classification are new
deductions within this repository. Searches on 27 October for adjoint
Springer/Cayley fibers and Levi classifications did not locate this exact
statement; that limited search does not establish publication priority.

Still open: arbitrary irreducible representations, exceptional polynomial-
adjoint constant components, special-target fibers, component/orbit counts,
and a uniform replacement of the exceptional smoothness proof. The focused
component-exhaustion audit below strengthens the justification of this
theorem; it is not another theorem about a larger representation family.

## 8. Focused component-exhaustion audit

27 October continuation of GEO-97. Outcome: the nonemptiness table and
dimension formula are unchanged. The following details were compressed in
the first version. They are needed to distinguish actual elimination from
a parametrization that merely supplies some solutions.

### 8.1 Equality of ideals and the local hypersurface

Adjoin a variable T together with T-tr_V(g)=0. This is an isomorphism of
coordinate rings, including nilpotents. In the differentiated full
character replace the OTHER trace by T minus the block trace. For every
GL block of size two, three or four, the resulting ideal equals the ideal
of the relative equations in section 2, with d invertible. This equality
holds before localizing at T, so the T=0 branches are retained. Newton's
identities give both inclusions directly: the three-block equations have
unit linear coefficients, and the four-block middle derivative is a unit
multiple of 1+d^-1. The latter first imposes d=-1, after which the other
two equations solve for e_1 and e_3; e_2 remains free.

The independent Sage check now verifies these six two-sided ideal
equalities (three block sizes, two signs), rather than only substituting
candidate solutions. For a two-block it also verifies exactly

    (full critical ideal, d^2+1) = (d*d_inverse-1, d^2+1, T).

The residual ideals are checked by equality with the INTERSECTION of the
branch ideals. In B1 these are (u,T-L+1) and (L+1,T-u^2+2). In D2 they are

    (u,v,T-L),
    (L-2,u+v,T-2+v^2),
    (L+2,u-v,T+2-v^2).

Here L is the GL trace. B2, C1 and D3 have the single ideals obtained by
the linear eliminations in section 2. These equalities are over Q in the
original fundamental-character coordinates, without taking radicals.
They also admit a direct hand check: in D2 use u+v,u-v to write the two
derivatives as (L+2)(u+v),(2-L)(u-v). The two exceptional values of L
cannot occur simultaneously. In B1 the single derivative factors into
the two distinct factors u and L+1. The T equation is then eliminated.

After removing T=+/-2 and the two-block loci s_i=0, the residual branches
are disjoint. Their remaining square roots have nonzero derivative.
Eliminating the displayed unit linear equations gives, locally in the
etale topology, precisely

    O(S)[T, w_1,...,w_q] / (A*T-B),                       (6)

where S is the remaining determinant torus with the stated denominators
inverted, and w_i are the free four-block e_2 coordinates. If a residual
branch requires a square root of a function of T, take its finite etale
pullback; this does not mean it has been chosen globally.

The equation A*T-B is not identically zero for any allowed proper Levi.
Indeed A identically zero forces only singletons and c=1, while B is then
nonconstant because there is at least one singleton. Thus (6) is a
hypersurface in a regular ring. It is Cohen--Macaulay and has no embedded
associated components. This is a justification for the height assertion,
not an assumption that the original derivatives form a regular sequence.
It does NOT assert reducedness: the D5 double-zero component persists.

### 8.2 Exhaustion of the support

Every point of (6) either has A nonzero, or has A=B=0. This dichotomy
exhausts the chart. If A and B are both nonconstant, their variables are
disjoint: A uses two-block determinants, B uses one/three-block determinants.
Each nonempty zero locus has codimension one in its own torus. Their
product therefore has codimension two, and adjoining T and the q free
coordinates gives dimension s-1, where s is the central rank. A constant
nonzero A or B makes the corresponding zero locus empty. The remaining
possibilities are exactly A identically zero or B identically zero, as
treated in section 4. Both cannot vanish identically.

For the nonconstancy arguments on these loci, the maps d_i to s_i are
finite free of degree two: d_i satisfies d_i^2-s_i*d_i+1, and is already
invertible in this algebra. The source over each integral hyperplane is
finite flat, so every minimal component dominates that hyperplane
(equivalently use going-down for flat extensions). Pullback cannot turn a
nonconstant function into a constant on one component. This argument uses
flatness, not merely finiteness, and allows ramification at s_i=+/-2.
For the reciprocal-sum locus use the same assertion after inverting s_i.

On A nonzero the graph (determinants,phi) is irreducible downstairs and
invariant under each determinant inversion. Its stabilizer Lie algebra
is invariant under all the resulting sign changes, hence is a coordinate
subspace. The dependence checks in section 3 exclude every unfixed
coordinate direction. Fixed four-block determinants exclude the others.
This proves dominance by GEO-93, rather than inferring it from scalar
nonconstancy alone. The constant D5 ordinary branch has its full remaining
torus as stabilizer and is correctly excluded.

Finally section 5 covers the complement of this chart: arbitrary
intersections of s_i=0, and T=+/-2. Those arguments bound the image of
every component, without a generic-position assumption on the other
parameters. Thus they also cover embedded support on that complement.
The actual nonsmooth locus is still removed using GEO-49 and GEO-95;
a dimension count on the support alone would not prove smoothness.

### 8.3 Central covers and regular group fibers

Two different cover operations must not be conflated. Taking an isogeny
of central tori to use fixed-determinant SL coordinates is finite etale
in characteristic zero. Residual Spin fundamental characters are already
used in section 2. We do NOT assert that every quotient map induced by
Spin -> SO is etale; it can be ramified. Where square roots occur in the
displayed residual branches, their etaleness follows from the explicit
nonzero derivative away from T=+/-2, not from the group isogeny alone.

On a splitting cover of the Levi the coordinate changes e_j=z^j a_j,
d=z^b are invertible with z a unit. Relative derivatives transform by an
invertible matrix. Changing the central derivative adds relative
derivatives, which vanish in the actual critical ring; the remaining
central differential is changed by the invertible isogeny matrix.
Consequently this operation preserves the ideal and dominance. Each
component of the finite etale pullback of the ordinary graph dominates
that graph. Inversions need not preserve an individual lifted component:
the stabilizer argument is made downstairs, and dominance and scalar
nonconstancy transfer by the cover. This avoids assuming an unproved
inversion symmetry on each Spin branch. Outer-automorphic type D Levi
embeddings have the same answer, since the adjoint trace construction is
preserved by the automorphism.

The return to group fibers also has a scheme-level justification. For
every representation Springer map the identity

    Ad(g) theta(g) = theta(g)

holds as an identity of regular functions. In the actual fiber ring at
h it implies Ad(g)h=h. Thus the entire fiber, not just its closed points,
factors through the scheme centralizer Z_G(h)=K. If g in K is K-regular
and theta(g)=h with h centralizer-exact, differentiation of equivariance
shows Lie Z_G(g) is contained in Lie Z_G(h)=k. It follows that
Lie Z_G(g)=Lie Z_K(g), of dimension rank G: the lift is G-regular.

Use the regular adjoint quotient on the splitting cover of K, whose
derived factors are simply connected. It is smooth and surjective with
relative dimension dim K-rank G. Pulling back the actual quotient fiber
and descending through the central cover gives the stated actual group
fiber and its dimension. This step uses the existing group-fiber
presentation, not an inference from eigenvalue multisets or reduced
conjugacy classes. In particular regular nonsemisimple lifts are included.

### 8.4 Audit verdict and remaining dependency

The original D5 ideal independently confirms (d+1)^2=0 with d+1 nonzero
after the indicated localization. Its central coordinate is likewise
nonzero nilpotent. Therefore its REDUCED target image is {0}; the actual
morphism is not the zero morphism. The classification uses only the former
statement for nondominance and keeps the latter in the verification.

No counterexample to the stated classification was found in this audit.
The added arguments close the previously compressed component-exhaustion
and cover steps. Six full block-ideal equalities, five residual-ideal
decompositions, six ambient-root checks, four tensor-weight convolutions,
and the original nonreduced D5 ideal check pass over characteristic zero.
These are arithmetic certificates and a written geometric audit, not
external mathematical review. The conclusion remains dependent on GEO-49,
GEO-93 and GEO-95; this audit does not independently reprove those theorems.

The next useful step is a review of that short dependency chain, especially
the stabilizer/dominance lemma and polynomial-composition image argument.
These are the general mechanisms that extend beyond classical blocks;
checking them has more value than adding another exceptional example.
