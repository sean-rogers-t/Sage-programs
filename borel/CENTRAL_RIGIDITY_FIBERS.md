# Character constancy and central rigidity of Springer fibers

24 October 2026. A uniform partial result in the attempt recorded in
[the active target](UNIFORM_ADJOINT_PROOF_TARGET.md). This is a deduction
from the trace identity and standard algebraic geometry. It does not
prove the all-simple-adjoint smoothness conjecture. No novelty claim.

## 1. The theorem, including its scheme limitation

Let G be a connected complex reductive algebraic group and rho a rational
representation whose trace form B_rho is nondegenerate. Write

    chi(g)=tr(rho(g)),
    B_rho(X,Y)=tr(d rho(X) d rho(Y)),
    B_rho(theta(g),X)=tr(rho(g) d rho(X)).

In particular, this includes every nontrivial representation of every
simple simply connected complex group, including the adjoint representation.
Fix a Levi K, with Lie algebra k=[k,k] direct-sum z, and put

    A_K=K/[K,K],    a:K -> A_K,
    Z_K^circ={h in z: Z_G(h)=K},    F_h=theta^(-1)(h).

Choose a Z-basis eta_1,...,eta_s of the character lattice of A_K, where
s=dim A_K.
View these as characters of K. Define target coordinates c_j(h) by

    B_rho(h,X)=sum_j c_j(h) d eta_j(e)(X),  X in k.       (1)

This is an invertible complex-linear coordinate system on z. In
particular c=0 if and only if h=0. It uses the actual abelianization
torus and its integral character lattice, not an unidentified central
cover. For a vector c define the saturated lattice and its subtorus

    Lambda_c={n in Z^s: sum_j n_j c_j=0},
    A_c subset A_K,  X_*(A_c)=Lambda_c.                   (2)

**Theorem.** For h in Z_K^circ and every irreducible component C of
the reduced fiber (F_h)_red:

1. chi is constant on C. Consequently chi takes only finitely many
   values on the complex points of F_h, and its scheme-theoretic image
   is finite, possibly nonreduced.
2. a(C) lies in a translate of A_c. Thus

       dim a(C) <= rank_Z Lambda_c.

3. For fixed rho and K there are finitely many nonzero integer vectors
   n^(1),...,n^(N) such that, outside the hyperplanes

       sum_j n_j^(i) c_j(h)=0,                            (3)

   the abelianization image of the complex points of F_h is finite,
   possibly empty. This is an actual
   Zariski-open conclusion, not only a very-general-target assertion.
4. At those targets the scheme-theoretic abelianization image is a
   finite scheme, possibly nonreduced. The ACTUAL fiber splits into
   finitely many open-and-closed subschemes, one over each supported
   abelianization point, retaining its possible infinitesimal directions.

If dim A_K=1, the image is finite for EVERY nonzero h in Z_K^circ.
For arbitrary dim A_K, Q-linear independence of c_1,...,c_s already
implies constancy of a on each C by part 2, without determining (3).

These assertions apply to the full semisimple-target fiber, including
nonregular source points. Constancy is asserted only on reductions;
part 4 retains the actual scheme. They do NOT assert smoothness, reducedness, finite
fibers, a finite morphism F_h -> A_K, or scheme-level constancy of a
or chi. Positive-dimensional components can remain inside a single
coset of [K,K]. The finite hyperplanes are not asserted to be the entire
singular-target discriminant or to be effective from this proof.

## 2. Proof from the trace identity

Equivariance gives scheme containment F_h subset K: the entries of
(Ad(g)-1)h equal the negatives of (Ad(g)-1)(theta(g)-h), so vanish in
the fiber ring. On K the defining trace formula gives

    d chi=B_rho(theta(g),g^(-1)dg).

Because h is central in k, B_rho(h,[k,k])=0. The restrictions of
B_rho to k and z are nondegenerate, and the characters of A_K span
the annihilator of [k,k]. Equation (1) is therefore well defined and
invertible. On the actual fiber we obtain

    d chi=sum_j c_j(h) d eta_j/eta_j.                    (4)

No choice of a logarithm is made. Pull (4) back to C, then to its
function field. The following elementary lemma applies.

### Exact versus logarithmic differentials

Let V be an irreducible complex variety, f a rational function and
z_1,...,z_s nonzero rational functions. Suppose

    df=sum_j c_j dz_j/z_j,   c_j in C.                   (5)

Then f is constant. Moreover, on a normal projective model of V,
every prime-divisor valuation vector

    v_D=(ord_D z_1,...,ord_D z_s)

belongs to Lambda_c.

Proof. A normal projective model exists by compactifying an affine
model and normalizing. At the generic point of a prime divisor it is
smooth. If f has a pole of order m>=1 there, then df has a nonzero
normal term of pole order m+1. Each dz_j/z_j has at most a simple
pole. Thus f has no divisorial poles. On a normal projective variety
such a rational function is regular and hence constant. The residue
of (5) at D is then sum_j c_j ord_D z_j=0, proving the assertion.
The same pole argument on normalized projective curves gives an
alternative proof of the constancy assertion.

For any m in the annihilator Lambda_c^perp subset Z^s, the rational
monomial z^m has order m.v_D=0 at every prime divisor. It and its
inverse are regular on the projective model, so z^m is constant.
Since Lambda_c is saturated, these monomials are exactly the characters
of the quotient torus A_K/A_c. The image of V therefore lies in one
coset of A_c. Applying the lemma to (4) proves parts 1 and 2. There
are only finitely many irreducible fiber components.

### From very general targets to finitely many hyperplanes

Use the entire incidence scheme

    Y={g in K: p_[k,k] theta(g)=0},
    nu:Y -> z,  nu(g)=p_z theta(g).

The restriction theta(K) subset k follows from equivariance. For
h in Z_K^circ, its fiber is F_h. Under (1) regard nu as a map to A^s.
The image I of (a,nu):Y_red -> A_K times A^s is constructible by
Chevalley's theorem. Hence

    E={c: dim I_c >=1}

is constructible: partition I into finitely many locally closed
subvarieties and use constructibility of fiber dimensions on each.
Here I_c is the ACTUAL image of the fiber, not a fiber of the closure
of I that might introduce additional points.

The preceding proof also applies to every fiber of nu, whether or not
its target has centralizer exactly K. If c belongs to E, at least
one fiber component has nonconstant central image. Part 2 implies
Lambda_c is nonzero. Consequently

    E subset union_(0!=n in Z^s) {c:n.c=0}.              (6)

Over the uncountable field C, an irreducible variety cannot be covered
by countably many proper closed subvarieties. Apply this standard
algebraic Baire fact to a finite irreducible locally closed partition
of E. Each member must be contained in one hyperplane in (6).
Finitely many such hyperplanes cover E, proving part 3. If E is empty,
take no hyperplanes. If s=0, the abelian quotient is a point and the
finite-image assertion is automatic.

Changing the integral character basis changes c by the dual integral
matrix. Thus rational hyperplanes, the lattice (2), and the conclusion
are intrinsic. Finite central covers change the lattice by finite index;
they do not change the rational nonresonance condition.

### Retaining the scheme structure in the finite central image

The fiber F_h is affine. Put R=O(F_h), including its nilpotents, and
let B be the image of O(A_K) in R. When the central point image is
finite, each eta_j has finitely many nonzero values lambda_(j,k) on
Spec R. Thus the monic polynomial product_k(T-lambda_(j,k)), evaluated
at eta_j, is nilpotent in R by the Nullstellensatz. A power annihilates
eta_j. The same holds for eta_j^(-1). Hence the finitely many generators
of B are algebraic over C, and B is a finite-dimensional C-algebra.
It defines precisely the scheme-theoretic abelianization image.

Decompose the Artinian algebra B as product_i B_i into its local factors.
Its idempotents lie in R, so R=product_i e_i R. This proves the
open-and-closed decomposition in part 4. It retains any nilpotent
central directions in B_i; setting eta_j equal to its reduced value
would be an additional quotient and is not justified.

The same argument applied to C[chi] proves the finite scheme-image
assertion in part 1, for every semisimple target. Empty fibers cause
no exception. This uses only elementary finite-type commutative algebra,
not reducedness or smoothness of the Springer fiber.

## 3. An independent geometric proof of generic central finiteness

On the incidence scheme Y, (4) holds with the variable functions nu_j:

    d chi=sum_j nu_j d log eta_j.

Taking its exterior differential and restricting to the smooth locus
of any reduced irreducible component gives

    sum_j d nu_j wedge d log eta_j=0.                    (7)

The right-hand side is the pullback of the standard symplectic form
on T^*A_K=A_K times A^s. In characteristic zero, the map to its image
is generically separable. Thus each irreducible image closure is
generically isotropic and has dimension at most s. Components dominating
the target have dimension exactly s, and have zero-dimensional generic
fibers over A^s. Remove the images of nondominant components and the
exceptional fiber-dimension loci. The remaining target open has finite
central image. This recovers generic central finiteness independently
of the residue/Baire argument, but does not identify rational hyperplanes.

This is the standard generating-family symplectic mechanism. We do not
assert that the singular incidence scheme is smooth, or that its whole
image is a smooth Lagrangian submanifold. Only its generic reduced image
and its dimension are used here.

## 4. What this does and does not resolve

For generic semisimple targets of any fixed Levi type, every reduced
positive-dimensional Springer-fiber component has fixed abelianization
coordinates. Any positive-dimensional freedom is within a derived-group
coset. The theorem therefore locates the possible excess geometrically,
uniformly across all groups and representations in its hypotheses.
The finite central algebra also gives a decomposition of the actual
fiber into pieces supported at individual central positions, without
discarding nilpotents.

It does not settle the active adjoint smoothness conjecture. The control

    phi=x^3+(z-1)y,  X=V(x^2,z-1),  mu=zy

has a double point for every target c, with z=1, y=c and phi=0. It
satisfies the rigidity conclusions while its singular incidence locus
still dominates the target. This is not an adjoint counterexample.

Even actual Springer maps require the reduction qualification. For the
standard representation of SL2 and regular semisimple target
h=diag(i,-i), the torus coordinate satisfies

    z^2-2i z-1=(z-i)^2=0.

Writing z=i+epsilon gives epsilon^2=0 and chi=z+z^(-1)=2epsilon.
Thus both z and chi can vary infinitesimally even when the reduced
fiber is a single point. This is the n=2 instance of the already
published standard SL_n fiber formula, not a new singularity result.

The remaining adjoint-specific problem is to exclude dominant singular
incidence components despite this now-proved central rigidity. Equivalently,
one must control the full transverse local algebras in these finitely
many central pieces, retaining nilpotents and excess directions. One
cannot simply set the central coordinates equal to reduced scalar values.
No new Levi cases or generic-smoothness coverage are claimed here.

## 5. Continuation: central nilpotents can persist at generic targets

A further attempted shortcut was to strengthen finite central image to
REDUCED finite central image after shrinking the target. The trace and
isotropy identities do not imply this. The following control disproves
that stronger inference even when the character value itself is exactly
constant on the actual fiber. It is NOT a representation character and
does not disprove the adjoint conjecture.

Put

    f=x^5+y^5+x^3 y^3,
    A=(C[x,y]/(f_x,f_y))[1/(25-9xy)].

This is the local critical algebra at the origin. Every other critical
point has x,y nonzero and xy=25/9, as follows by multiplying
5x^2+3y^3=0 and 5y^2+3x^3=0. The localization removes those points.
An exact presentation, obtained by saturating the Jacobian ideal, is

    A=C[x,y]/(5x^4+3x^2 y^3, 3x^3 y^2+5y^4, xy^4, y^6).

The displayed generators form a lexicographic Groebner basis. The
standard monomials number 16. In this algebra,

    f=-x^3 y^3/5=y^5/3 != 0,    f^2=0.                  (8)

These facts also have a direct linear-algebra check. On a vector space
with basis e_(i,j), 0<=i,j<=3, let X increment i when i<3 and let Y
increment j when j<3. At the boundary set

    X e_(3,0)=-(3/5)e_(2,3),   X e_(3,j)=0 for j>=1,
    Y e_(0,3)=-(3/5)e_(3,2),   Y e_(i,3)=0 for i>=1.

Then X,Y commute, satisfy f_x=f_y=0, and 25-9XY is invertible.
The vector e_(0,0) is cyclic. Moreover f(X,Y) has rank one, is nonzero,
and has square zero. This independently proves that the class in (8)
survives; it is not a numerical approximation or a radical computation.

Now take the smooth ambient space

    M=Spec C[x,y,1/(25-9xy)] times Gm_z,
    phi=z(1+f),    D=z partial_z.

Its relative critical incidence and central target function are

    Y=V(phi_x,phi_y)=Spec A times Gm,
    mu=D phi=z(1+f).

They satisfy d phi=mu d log z on the ACTUAL incidence. For every c!=0,
the actual target fiber is

    O(Y_c)=A,    z=c/(1+f)=c(1-f),    phi=c.             (9)

Thus the central point image is the single point z=c, but its
scheme-theoretic image is exactly

    Spec C[z,z^(-1)]/((z-c)^2).                          (10)

The whole fiber has length 16, whereas replacing z by its scalar value
c imposes f=0 and gives length 15. The lost direction is real scheme
structure. It persists for EVERY nonzero target, so no generic-target
argument removes it. Even a constant actual value of phi does not imply
constant actual z: the nonzero class f has df=0 in Omega_A/C because A
is a Jacobian algebra. One must not assume that the kernel of d on a
nonreduced algebra consists only of scalars.

This uses the familiar distinction between a local Jacobian (Milnor)
algebra and its quotient by the function (Tjurina algebra); no novelty
is asserted for that singularity-theoretic mechanism. All algebraic
claims needed here are proved by the displayed identities and the two
exact computations. This control is not evidence of a bad adjoint fiber.

### The obstruction survives a quartic degree bound

The preceding degree-six f is not essential. Take instead

    f_4=x^4+y^4+w^4+xyw,
    A_4=(C[x,y,w]/(partial f_4))[1/(1+64xyw)].

Multiplying the three derivative equations gives
(xyw)^2(1+64xyw)=0. A critical point with any zero coordinate is the
origin; every other critical point has xyw=-1/64 and is removed.
An exact lexicographic Groebner basis for the remaining local algebra is

    x^3+yw/4, xy+4w^3, xw+4y^3, y^4-w^4,
    y^2w, yw^2, w^5.

Consequently dim A_4=11 and

    f_4=xyw/4=-w^4 != 0,    f_4^2=0.

For an independent check use the cyclic basis

    1,x,x^2,x^3,y,y^2,y^3,w,w^2,w^3,q.

Multiplication obeys x^4=y^4=w^4=q, xy=-4w^3, xw=-4y^3,
yw=-4x^3; products of q with x,y,w and mixed products having one
exponent at least two vanish. These rules give commuting multiplication
matrices satisfying all three derivatives. The matrix of f_4 is
nonzero of rank one and square zero. This checks the surviving class
without using Groebner normal forms.

Applying the same construction phi=z(1+f_4) gives length-11 fibers,
central double-point images at every c!=0, and length 10 after the
invalid substitution z=c. In particular, the already-proved bound of
four affine variables and affine degree at most four is insufficient
BY ITSELF to establish generic smoothness. This polynomial is not
claimed to be a residual adjoint character. Additional constraints on
the actual root-derived coefficients are essential.

**Outcome of the uniform attempt.** The central-rigidity theorem remains
unchanged. There is no proof yet that adjoint branching excludes the
mechanism (8)-(10), or more generally excludes singular incidence over
generic central targets. The next proof step must work with the full
local critical algebra. Neither generic central finiteness nor removal
of exceptional target values justifies scalar specialization of the
central coordinates. No new Levi case or smoothness coverage is claimed.

## 6. Inputs and verification

- [Kostant--Michor](https://www.mat.univie.ac.at/~michor/kostant.pdf),
  sections 2 and 6.2: the trace construction and standard SL_n formula.
- Pole orders of exact and logarithmic rational differentials, normal
  projective compactification, and character/cocharacter duality are
  standard inputs; the proof above spells out their application.
- [Stacks, Lemma 37.30.3](https://stacks.math.columbia.edu/tag/05F9):
  constructibility of fiber-dimension level sets.
- [Katthaen--Michalek--Miller](https://arxiv.org/pdf/1706.03629),
  section 2 (Theorem 2.3 in arXiv v1): algebraic Baire theorem.
- [Tulczyjew--Urbanski](https://arxiv.org/pdf/1006.4183), sections 4 and 9:
  generating families and their symplectic interpretation. Their
  smooth-Lagrangian conclusions have hypotheses; (7) only uses the
  weaker generic isotropic-image argument proved here.

The [exact audit](tools/audit_central_rigidity.py) checks the universal
mixed-Hessian identity behind (7), divisorial residue/character-lattice
controls, a resonant subtorus family, excess and persistent-nilpotent
controls, the actual SL2 nonreduced fiber, and the generic nonreduced
central-image control (8)-(10), by saturation and a separate explicit
16-dimensional module, with a quartic three-variable variant checked
by a separate 11-dimensional module. Its
[report](central-rigidity-audit.json) is corroboration; it does not
replace the geometric proof or constitute independent human review.
The two written proof routes agree on generic central finiteness.
Targeted literature checks located the inputs above, not a priority
certificate for their application to Springer fibers. The result is
recorded as a deduction from standard results, not an original theorem claim.
