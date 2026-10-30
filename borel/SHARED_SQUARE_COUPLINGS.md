# Shared-variable quadratic couplings and their ramification boundary

26 October 2026. This continues the [uniform adjoint proof target](UNIFORM_ADJOINT_PROOF_TARGET.md).
The conclusion is a general algebraic family theorem, not a proof that all
remaining root-derived potentials belong to this family. No new exceptional
Levi coverage or literature priority is claimed.

## 1. Result and exact scope

Let T=(C*)^s with s>=1. Choose r>=1 blocks, each with m_i>=1 pairs of
affine variables u_ia,v_ia, and one shared affine variable x. Put

    w_i=1-sum_a u_ia v_ia, W=sum_i w_i,
    f=sum_i w_i*(x^2-delta_i(z))+Q(z),
    X=V(f_x, all f_uia, all f_via), mu_j=D_j f,
    D_j=z_j partial_zj.                                      (1)

Here Q and delta_i are Laurent polynomials. Assume each delta_i is
nonzero and its differential vanishes at every point of its zero locus:

    delta_i(z)=0 implies d delta_i(z)=0.                     (2)

In particular (2) holds if delta_i=S_i(z)^2 for a nonzero Laurent
polynomial S_i. More generally every irreducible factor may have
multiplicity at least two. Units satisfy (2) vacuously.

**Theorem (GEO-81).** The actual target fibers of mu are smooth over a
dense open subset of A^s. Empty fibers and positive-dimensional fibers
are allowed, as are different dimensions on disjoint components. The
relative critical ideal is NOT assumed radical, equidimensional, or a
complete intersection. The assertion includes every ramification point
and all coefficient-boundary intersections.

The shared x makes this a coupled family: this is not a disjoint product
of independently critical blocks. No equation f_i,x=0 for an individual
summand is imposed. The total equation is xW=0.

The same proof permits w_i=lambda_i-sum_a u_ia v_ia for constants
lambda_i with sum_i lambda_i!=0. For clarity the proof uses lambda_i=1.
Independent nonsingular quadratic blocks and free affine variables may
also be added. The central torus and derivatives must be the actual
ones; a z-dependent affine coordinate change is accompanied by the
chain-rule correction modulo the full derived gradient ideal.

For a Springer application, (1) must be derived from the FULL restricted
character by scheme-preserving changes and justified eliminations.
The one-block E6 A2+A1+A2 form already derived from ambient roots is an
example, with delta=(z^3+z^-3)^2/4 and Q=-2. There is no proof here that
the remaining E7/E8 couplings admit a shared-variable presentation (1).
In particular, (2) is a verified hypothesis for that existing E6 form,
not an assertion about every Schur discriminant in exceptional types.

## 2. Recover the quadratic dimension bound off the ramification locus

The exact derived ideal, ignoring nonzero scalar factors and signs, is

    J=(xW, (x^2-delta_i)u_ia, (x^2-delta_i)v_ia).             (3)

There are two kinds of charts. Their description uses the ORIGINAL
equations, not their radical.

**First, W!=0.** The equation xW=0 eliminates x exactly. The remaining
potential and target derivative are those of

    R=Q-sum_i w_i delta_i.                                  (4)

This is quadratic in all u,v, with arbitrary Laurent coefficients in z.
The [quadratic theorem](QUADRATIC_BRANCHING_FIBERS.md) gives smooth generic
target fibers. Restricting its fibers to the open W!=0 preserves this.

**Second, x!=0.** The equation xW=0 gives W=0. Since W=r-sum u_ia v_ia
and r is nonzero in C, at least one u_ja or v_ja is nonzero. Choose that
variable as a unit. Then w_j is a valid affine coordinate: solve for
its paired variable, retaining all other variables of that block as
free coordinates on this chart. Replace x by q=x^2. This change is
etale here, because 2x is a unit. Set

    a=q-delta_j, b=w_j+sum_(i!=j) w_i.

The full potential has the exact identity

    f=a*b+R_j,
    R_j=Q+sum_(i!=j) w_i*(delta_j-delta_i).                  (5)

The equations for a,b eliminate b,a respectively, without denominators
or regular-sequence assumptions. R_j is quadratic in all remaining
nonfree paired variables. The theorem again follows from the quadratic
case and etale invariance of smoothness. On the resulting critical
chart x^2=delta_j, so delta_j is a unit there; it is not inverted on
other charts or on the entire source. Free coordinates contribute
smooth factors. There are finitely many choices of j and of a unit
paired variable, so their target opens can be intersected.

For precision, (5)'s z-dependence changes D_z f by terms proportional
to f_a,f_b and the other derived gradients. These terms vanish in the
actual critical quotient. Thus this is an identification of target-fiber
schemes, not merely of reduced critical supports.

All points outside the closed incidence subset

    E=X intersect V(x,W)                                   (6)

are now covered. The two generic-smoothness arguments do not require
the separate summands of f to be critical.

## 3. The boundary image factors through a proper central divisor

At a point of E define an active block to be one containing a nonzero
u_ia or v_ia. Equation W=0 forces at least one active block. For an
active block (3), with x=0, forces delta_i=0. By (2), every D_j delta_i
also vanishes there. For an inactive block all pairs vanish and w_i=1.
Consequently the full target map on the underlying set of E is

    mu=D(Q-sum_i delta_i).                                 (7)

The central projection of E lies in

    B=union_i V(delta_i) subset T.

The product of the nonzero delta_i is nonzero, hence dim B<=s-1.
Equation (7) proves the explicit bound

    closure(mu(E)) subset closure(D(Q-sum_i delta_i)(B)),
    dim closure(mu(E)) <= s-1.                             (8)

Remove this proper closed target subset, then intersect the finitely
many target opens provided by the quadratic charts. This proves the
theorem. The same calculation gives D(Q-sum_i lambda_i delta_i) for
constant lambda_i with nonzero sum.

Only the IMAGE DIMENSION is computed on an underlying constructible
set. No critical ideal is radicalized. Every point of a retained actual
target fiber lies in a chart where its full scheme has been proved
smooth. Thus embedded and nonreduced structure cannot silently persist
over those targets. Conversely this proof does not supply an explicit
entire colon (J:sqrt(J)) or the sharp target annihilator for r>1.

This is the replacement for the quadratic argument in this family:
recover the quadratic bound in the charts, and use full coefficient
first derivatives to bound the image of their ramification boundary.
It is not a universal cubic-Hessian rank theorem.

## 4. Why aligned blocks cannot be glued without coefficient control

Here is a counterexample to dropping all control of the simultaneous
zero jets. It is an abstract potential, NOT an adjoint character:

    f=(1-u*v)(x^2-(z-1))+(1-p*q)(x^2+(z-1))+z.              (9)

Each individual summand, with half the constant term z, is in the
already proved one-block aligned family with nonzero delta. Nevertheless
their shared-variable sum has a persistent embedded nilpotent.

On the valid chart u,p!=0 set a=uv,b=pq. Over L=C(c), with c the
generic target coordinate, the full gradient and target ideal is EXACTLY

    (z-1, b-a-1+c, x^2, x*(1+c-2a)).                       (10)

The omitted coordinates u,p are invertible free variables. The class x
is nonzero and square-zero, with annihilator (x,1+c-2a) in the remaining
quotient. It is supported on an embedded component. Thus even generic
target localization does not remove it. The two delta functions vanish
at z=1 with opposite nonzero derivatives. Equation (7) fails precisely
because the target remembers the difference a-b, which can vary while
a+b=2 on E. Its image there is dominant.

This does not show (2) is necessary in every successful family. The
one-block theorem already works for arbitrary nonzero delta. It shows
that independent one-block smoothness is insufficient for gluing, and
identifies a concrete compatibility issue at shared boundaries.

The nonzero sum of the constants is also substantive. With one block
w=-uv and delta=1, f=-uv(x^2-1)+z, the chart near x=1 has coordinate
t=x^2-1 and critical ideal (uv,tu,tv). This is a reduced singular
three-axis crossing over EVERY generic target c=z. All coefficient
zero-jet conditions are vacuous, but the all-zero-pair locus is now
allowed in E and projects to the entire central torus. This explains
exactly where the proof uses sum lambda_i!=0.

## 5. Independent verification and existing controls

The [SymPy audit](shared-square-couplings-audit.json) checks (5) for
one through five blocks, (7) for every nonempty active subset in those
sizes, and the central chain rule. It checks (10) over QQ(c), including
the nonzero square-zero class, without numerical specialization.

The [independent Sage/Singular check](shared-square-couplings-verification.json)
imports no native audit. It starts from the ORIGINAL gradient ideal of

    f=(1-uv)(x^2-(z-1)^2)+(1-pq)(x^2-(z+1)^2)+z,

adjoins zr-1 to retain the torus, and imposes z*f_z=c over QQ(c).
Its exact ideal is the intersection of the following three pairwise
disjoint ideals (each also contains zr-1):

    (x,u,v,p,q, 4z^2-z+c),
    (p,q,uv-2,x^2-(z-1)^2,3z+c),
    (u,v,pq-2,x^2-(z+1)^2,5z-c).

Their dimensions are 0,1,1. Independent radical and Jacobian-ideal
computations prove each smooth over QQ(c); the entire original ideal
is reduced. This checks mixed dimensions and positive excess, not just
isolated critical points. The same verifier independently establishes
(10), its radical and the exact annihilator of x.

The existing E6 contact generic algebra L[z,w1,w2]/(z^2+1), with its
two-dimensional sheets, is retained. It is not claimed to have (1),
and no source localization removes it. The altered-Levi ideal
(u^2v,uv^2) remains nonreduced. The Sp4 special triple-point conclusion
is unchanged: generic smoothness never asserts that every special
target fiber is smooth. These small control checks are consistency
checks, not a fresh reconstruction of their root calculations.

## 6. An existing adjoint model obstructs a universal affine normal form

The theorem's hypotheses cannot be assumed for all residual root-derived
potentials. This can be tested against the EXISTING, already solved E8
model of derived Levi A2+A1+A2+A1; no new Levi census is involved.
In its archived residual variables x,u,v,y, over F=C(t,r), the top
homogeneous degree-four part is

    f_4=-(r*v*x+r*x*y+t*u*y)*(r*v*y+t*u*x+t*x*y)/(r*t).       (11)

The Hessian determinant of each quadratic factor is r^2*t^2!=0.
Each therefore has rank four and is irreducible even over an algebraic
closure of F: a reducible homogeneous quadratic is a product of linear
forms and has rank at most two. The two factors are not associates
(their monomial supports differ). Thus f_4 is squarefree.

In contrast, every potential (1), with any number of free variables or
independent degree-two summands, has top degree-four form

    -x^2*sum_(i,a) u_ia v_ia.                              (12)

It has a repeated linear factor. Invertible affine changes of the derived
variables over F transform the highest form by their invertible linear
parts and preserve factor multiplicities. Adding independent nonsingular
quadratic or free variables does not change this obstruction. Hence no
such affine change, even after those additions, puts the existing E8
residual potential into (1).

This is LIMIT-19: a counterexample to that UNIVERSAL AFFINE NORMAL-FORM
route, not to the shared-square theorem or the adjoint conjecture. It
does not exclude nonlinear triangular changes, further eliminations,
etale local normal forms, or other ways to bound the bad target image.
It explains why another successful star-shaped model cannot be presented
as a root-theoretic explanation for all coupled restrictions.

Both independent audits verify (11), the two Hessian determinants, and
the squarefree factorization over characteristic zero. They use the
previously audited root-derived potential as a shared INPUT; they do
not claim fresh independent reconstructions of its entire E8 branching.

## 7. Attribution, remaining root problem, and next gate

Published inputs are characteristic-zero generic smoothness and
[etale-local smoothness](https://stacks.math.columbia.edu/tag/03ZB).
The [quadratic proof](QUADRATIC_BRANCHING_FIBERS.md) is an earlier
repository deduction using the mixed-derivative rank bound and
[generic smoothness in the target](https://math.stanford.edu/~vakil/0506-216/216class5152.pdf).
The coupled chart argument and exact failure (10) are deductions here;
the finite tests corroborate, rather than prove, the arbitrary-size
statement. No exhaustive literature search or originality claim is made.

**Still unproved:** no nonlinear or local root-theoretic reduction puts
the remaining overlapping carriers into (1), independent paired blocks,
the same shared coordinate, or condition (2). The universal affine
version is now disproved by (11). Individual carrier formulas cannot be
summed with incorrect root multiplicities. The previous E8 determinant-only
obstruction also remains valid.

The next root-theoretic gate is a direct identity for the FULL coupling
that forces the singular-block boundary target image to be nondense,
without assuming a universal affine shared-variable form. A useful
first test is the simultaneous coefficient FIRST JETS, since
(9) shows exactly how they can otherwise produce dominant embedded
nilpotents. Another isolated successful specialization would not settle
that gate. The all-simple adjoint conjecture remains open.
