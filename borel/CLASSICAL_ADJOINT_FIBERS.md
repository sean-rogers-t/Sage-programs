# Classical adjoint Springer fibers: smoothness and dimension in every rank

Date: 2026-10-23. Written characteristic-zero deductions; finite symbolic
and original-matrix checks are corroboration, not the all-rank proof.
Literature priority is unverified. The public API is unchanged.

## 1. Statements and scope

Let G be simple complex simply connected of type A, B, C or D, let
rho=Ad, and use the trace-normalized Springer map

    B_Ad(theta(g),X)=tr(Ad(g) ad(X)),
    B_Ad(Y,X)=tr(ad(Y) ad(X)).                                  (1)

Thus G is SL_N (N>=2), Spin_(2n+1) (n>=2), Sp_(2n) (n>=2), or
Spin_(2n) (n>=4), respectively. Lower-rank isomorphisms are treated in
their A/B/C type, not as extra D cases. Fix a proper Levi K and put

    Z_K^circ={h in z(k): Z_G(h)=K},  r=rank G,
    F_h^reg=theta^(-1)(h) intersect G_reg.

All scheme assertions concern this open part of the fiber, not just its
reduction and not the whole fiber. The open set below can depend on K.

**Theorem A.** There is a dense open U in Z_K^circ such that, for h in U,
F_h^reg is smooth and either empty or of pure dimension

    dim K-r.                                                   (2)

Equivalently, the quotient critical scheme is finite and reduced over
generic central targets (possibly empty). Its regular group lift has
only the conjugacy-orbit dimensions; there is no additional quotient
parameter. This does not count those orbits or claim regularity preservation.

The type A nonemptiness criterion is already proved in
[the quadratic note](QUADRATIC_BRANCHING_FIBERS.md). For B/C/D, write the
image of a Levi in the defining matrix group as

    product_i GL_(b_i) times H_m,  sum b_i+m=n.                 (3)

Here H_m is SO_(2m+1), Sp_(2m), or SO_(2m). In D absorb SO_2 into
the GL_1 factors, so m=0 or m>=2. Finite central covers are retained.

**Theorem B (complete generic nonemptiness criterion).** The generic
regular fiber in B/C/D is nonempty if and only if every b_i<=4 and

    B: m<=2;    C: m<=1;    D: m<=3.                            (4)

If a bound fails, the regular fiber is empty for every h in Z_K^circ.
If the bounds hold, Theorem A gives its generic pure dimension explicitly:

    sum_i b_i(b_i-1) + 2m^2                 in B/C,
    sum_i b_i(b_i-1) + 2m(m-1)              in D.              (4a)

This is a statement for every block pattern in every permitted rank.
It is not a statement that the entire fiber is empty in the excluded
cases, or a classification of nongeneric target fibers.

## 2. A criterion for square substitutions

Work on A^ell_x times A^k_y times T_z, T=(C*)^s, and suppose

    chi=P(x,z)+sum_i A_i(x,z)y_i^2,                             (5)

where P is quadratic in x and each A_i is affine in x, with Laurent
coefficients in z. Let X=V(partial_x chi, partial_y chi) be the actual
relative critical scheme, and mu=(D_j chi):X->A^s, D_j=z_j partial_z_j.
Let E be the subset of X on which y_i=A_i=0 for at least one i.

**Square-substitution criterion.** If the closure of mu(E) is proper
in A^s, the generic fibers of mu are smooth, possibly empty.

Proof. At a point outside E, each vanishing y_i has A_i invertible.
The equation 2A_i y_i=0 eliminates that y_i scheme-theoretically.
For each remaining, nonzero y_i, replace it by q_i=y_i^2. This is an
etale substitution on that open set. The potential becomes

    P(x,z)+sum_(remaining i) A_i(x,z)q_i,

which is jointly quadratic in (x,q). Its exact critical equations are
the pullback equations: partial_y_i=2y_i partial_q_i, and the unit 2y_i
can be cancelled only on this open set. The quadratic theorem applies.
There are finitely many subsets of remaining indices. Intersect their
generic target opens and remove closure(mu(E)). This proves the claim.
The original incidence scheme has not been replaced by its radical.

The condition is necessary for this argument. For example

    chi=x^2/2+x y^2+z

has generic critical ideal (x+y^2,y^3,z-c); its y^3 structure persists
for c!=0. Here E dominates the target. This is a control potential, not
a representation character. Also replacing a square by y^d for d>2
does not preserve the argument: already y^3+z has a double critical
structure. Neither polynomial degree nor a ramified quotient change
alone is a smoothness theorem.

## 3. The adjoint matrix identity and a spectral lemma

For orthogonal groups put N=dim V, d=N-2, T=tr_V(g). The representation
Ad is exterior^2 V and B_Ad(X,Y)=d tr_V(XY). Then

    theta(g)=[T(g-g^(-1))-(g^2-g^(-2))]/(2d).                 (6)

For Sp_N put d=N+2. Ad is Sym^2 V and

    theta(g)=[T(g-g^(-1))+(g^2-g^(-2))]/(2d).                 (7)

These formulas are on the defining matrix image; for Spin they pull
back through its finite etale central isogeny. They follow by
differentiating chi_Ad=(T^2-tr(g^2))/2 or (T^2+tr(g^2))/2 and projecting
onto the form-skew matrices. Thus the constants, not only proportionality,
are fixed by (1).

Set tau=T in B/D and tau=-T in C. On a GL block with nonzero target
eigenvalue lambda, its matrix satisfies

    p_(tau,L)(g)=0,
    p_(tau,L)(u)=u^4-tau u^3+L u^2+tau u-1,                  (8)

where L=2d lambda in B/D and L=-2d lambda in C. The L_i are independent
coordinates on z(k). We use the dense open subset with L_i!=0 and
L_i!=+/-L_j, which lies in Z_K^circ. These are sufficient generic
conditions; a single zero GL_1 coordinate can still belong to Z_K^circ
in type D with m=0, so they are not claimed to be necessary there.

**Spectral lemma.** Fix tau!=0 over a characteristic-zero field and let
L be transcendental. For a subset of b roots of p_(tau,L), the sum

    S=sum_(chosen roots u) (u+u^(-1))

is algebraic over the constant field only when b=0 or b=4. In those
cases it is 0 or 2tau. Roots are simple for generic L.

Proof. The quartic has nonzero constant term. For b=1, a constant
u+u^(-1) makes u algebraic over constants, and (8) then makes L
algebraic, a contradiction. The case b=3 follows by taking the
complement, since the sum over all four roots is 2tau. For b=2 write
v for the chosen product and a for the chosen sum. The complementary
product is -1/v and its sum is tau-a. Comparing the cubic elementary
symmetric function gives

    a(v^2+1)=tau v(v+1),
    S=tau (v+1)^2/(v^2+1).                                   (9)

The denominator cannot vanish when tau!=0. The rational function on
the right is nonconstant in v; if S is algebraic over constants, so
are v,a and both chosen roots, again forcing L algebraic. Finally, a
multiple root satisfies the derivative of the nonconstant rational
function L=-u^2+tau u-tau/u+1/u^2. There are only finitely many such
L for fixed tau. This proves separability and the lemma.

When tau=0 a two-root set {u,-u} can have S=0. We do not apply the
lemma at tau=0. This exception is important for nonemptiness.

## 4. Close the two cubic families in every rank

Let H=sum_i(E_i+E_i^*) denote all the GL summands and write L_H=chi_H.
It is affine in the semisimple GL fundamental characters, including
the central contributions of GL_1 blocks. Write P_0=chi_(exterior^2 H),
which is quadratic in those characters with Laurent central coefficients.

For residual B_1=Spin_3=SL_2, let y be its fundamental character.
The residual vector and adjoint characters both equal y^2-1. Thus

    chi_Ad=P_0-L_H-1+(L_H+1)y^2.                             (10)

The bad subset in the square criterion has y=0, L_H=-1. Its defining
trace is T=L_H+(y^2-1)=-2, and its residual trace is -1.

For residual D_2=Spin_4=SL_2 times SL_2, let y_1,y_2 be the two
fundamental characters. The residual vector is y_1 y_2 and the adjoint
is y_1^2+y_2^2-2. Setting u=y_1+y_2, v=y_1-y_2 gives

    chi_Ad=P_0-2+(2+L_H)u^2/4+(2-L_H)v^2/4.                 (11)

The two square coefficients never vanish together. On the relative
critical scheme, either bad coefficient forces the other variable to
vanish too. Hence the bad subsets have y_1=y_2=0 and L_H=+/-2;
their trace is T=+/-2 and their residual trace is zero.

**Bad-image exclusion.** None of these bad subsets dominates A^s.
Suppose one did. Take independent generic target coordinates L_i and
a regular lift from the quotient to K. Such lifts exist by the regular
adjoint-quotient presentation; they are G-regular because equivariance
gives Z_G(g) subset K. Each GL block is cyclic, so (8) forces b_i<=4.
For these fixed nonzero T, its roots are simple for generic L_i.

The trace identity is

    T=sum_i S_i+R,                                           (12)

where R=-1 for the B bad locus and R=0 for D. Differentiate (12) in
each independent L_i, extending derivations to the finite root fields.
All the other summands are annihilated. Thus S_i is algebraic over
the constant field. The spectral lemma forces every b_i=4. If there
are t blocks, (12) now says

    B: -2=-4t-1;       D: T=2tT with T=+/-2.

These require t=1/4 and t=1/2, respectively, impossible for an integer.
This is an argument at the generic point of any hypothesized dominant
component, not a numerical check or a dimension guess.

The square-substitution criterion proves generic smoothness for the
entire cubic families. The previous quadratic theorem covers all other
classical adjoint Levis. It covers type A as well. Smoothness descends
through the finite central covers. This proves the smoothness part of A.

## 5. Residual bounds and the trace choices

On the zero-target residual block, (6)-(7) give

    (g^2-1)(g^2-tau g+1)=0.                                  (13)

A regular GL block has one Jordan block per eigenvalue, which proves
b_i<=4. For the residual block use the classical regular Jordan forms.
Away from tau=+/-2, the roots of (13) are simple. A non-self-reciprocal
pair has multiplicity at most one at a regular point. In SO, the +1
and -1 eigenspaces have dimension at most two, except that the odd
orthogonal +1 space has dimension one. In Sp no +1 or -1 eigenspace
is allowed at a regular semisimple point.

At tau=2 or -2 the corresponding root has multiplicity three. A
regular unipotent in SO_(2a+1), SO_(2a), or Sp_(2a) has partition
(2a+1), (2a-1,1), or (2a), respectively. Its largest block must have
size at most three. The other self-reciprocal eigenvalue still has
only simple blocks. These observations give (4), including at the
collision values. They also give the following complete residual
trace possibilities for regular residual solutions of (13):

| Family | m | Residual trace R |
| --- | --- | --- |
| B | 0 | 1 |
| B | 1 | -1 or T+1 |
| B | 2 | T-1 |
| C | 0 | 0 |
| C | 1 | -T |
| D | 0 | 0 |
| D | 2 | 0 or T+2 or T-2 |
| D | 3 | T |

At T=+/-2 choose the regular Jordan form just described, rather than
assuming the residual matrix remains semisimple. The table concerns
trace choices, not counts of Spin lifts or SO conjugacy classes.

There is also an independent character proof of the residual exclusions.
For B_m with m>=3 and D_m with m>=4, the residual adjoint character is
the second fundamental character, while its vector character is the
first. Its second-character derivative in the full restriction is 1.
For C_m with m>=2, the residual adjoint is x_1^2-x_2-1, and that
derivative is -1. The GL cross terms use only x_1. Thus there is no
relative critical point in these residual ranks, independently of any
Jordan-form computation.

## 6. No positive-dimensional generic quotient fiber

For fixed (T,L_1,...,L_s), (8) and (13) give finitely many eigenvalue
multisets and Jordan forms. Hence there are finitely many quotient
points; taking Spin lifts adds only a finite ambiguity. More formally,
all eigenvalues and their inverses are integral over C[T,L_1,...,L_s]
by the monic polynomials with unit constant terms. The reduced quotient
incidence is contained in a finite spectral cover with the trace equation.
Thus a positive-dimensional generic quotient fiber would force T to
be transcendental over C(L_1,...,L_s).

Assume this occurs on a component dominating the central target space.
Now work with T and the L_i all algebraically independent. Each S_i
is algebraic over C(T,L_i), and the residual trace is aT+b from the
table. Differentiate (12) in L_i while holding T and the other L_j
fixed. The same spectral lemma applies over C(T), since tau=+/-T
is nonzero there. Every block must have size four.

There is at least one GL block because K is proper. If their number
is t, in B/D their total trace is 2tT, while R=aT+b with a=0 or 1.
The identity (1-2t-a)T=b cannot hold for transcendental T and t>=1.
In C the total trace is -2tT and R=-aT with a=0 or 1; the identity
(1+2t+a)T=0 is equally impossible. This contradiction proves generic
quotient dimension zero. The already proved smoothness makes that
generic scheme finite and reduced. Its regular group lift is smooth
of pure dimension dim K-r. Together with the earlier type A result,
this proves Theorem A in full.

This argument does not assert that the quotient map is finite over
every special target or that a nonempty generic quotient exists.

## 7. Nonemptiness for every permitted block pattern

We prove dominance by deforming from target L_i=0, not by counting
parameters or extrapolating examples. At zero target the quartic roots
are +1, -1, alpha, alpha^(-1), with alpha+alpha^(-1)=tau. They are
distinct when T!=+/-2. Choose b_i of these four roots in each GL block.
Let epsilon_(i,+), epsilon_(i,-) be zero or one according as +1,-1
is selected, and let q_i in {0,1,2} count the alpha roots selected.
Set

    D=sum_i(epsilon_(i,+)-epsilon_(i,-)),  Q=sum_i q_i.

At L_i=0 the GL trace sum is 2D+QT in B/D and 2D-QT in C. Choose a
residual branch R=aT+b from the table. The trace equation becomes

    F(T,0)=(1-Q-a)T-(2D+b)             in B/D,
    F(T,0)=(1+Q+a_C)T-2D              in C,                  (14)

where a_C=0 for m=0 and a_C=1 for m=1. It suffices to choose roots so
that (14) has a solution T_0!=+/-2 with nonzero coefficient of T.
The roots then extend on an etale neighborhood in (T,L_1,...,L_s).
The nonzero T derivative implies that F=0 maps etale to the target
near this point. Hence the target map is dominant and its image
contains a dense open set. Remove the target hyperplanes L_i=0 and
L_i=+/-L_j. The resulting cyclic GL blocks and regular residual block
give K-regular solutions; equivariance makes them G-regular.

Here are choices proving the assertion for every block pattern. A
balanced choice means D=0 in each block: choose one alpha root for
b=1; either {+1,-1} or both alpha roots for b=2; {+1,-1,alpha} for
b=3; and all roots for b=4. If f is the number of odd blocks and t_4
the number of four-blocks, these choices allow

    Q=f+2t_4+2j,

where j is any number of two-blocks. In particular Q can be made
positive whenever there is a GL block. It can be made different from
one unless there is exactly one block and its size is one or three.

- C: use any balanced choice and T_0=0. The coefficient 1+Q+a_C
  is positive, so this always works for m=0,1.
- B, m=0 or 1: use R=1 or -1 and balanced Q!=1. Then
  T_0=+/-1/(1-Q), which is never +/-2. In the exceptional single-odd
  pattern, select {+1} for b=1 or {+1,alpha,alpha^(-1)} for b=3.
  These have D=1 and Q=0 or 2. They give T_0=3 or -3 for m=0,
  and T_0=1 or -1 for m=1; the T coefficient is 1 or -1.
- B, m=2: use R=T-1 and a balanced Q>0. Then T_0=1/Q, with
  nonzero coefficient -Q and T_0!=+/-2.
- D, m=0 or 2: use R=0 and balanced Q!=1, giving T_0=0. The
  remaining single-odd patterns for m=0 have n=1 or 3 and are outside
  D_n, n>=4. For m=2 the single b=1 pattern also has n=3. The one
  remaining pattern is m=2,b=3: use R=T-2 and select
  {+1,alpha,alpha^(-1)}. Then D=1,Q=2, F(T,0)=-2T, so T_0=0.
- D, m=3: use R=T and balanced Q>0. Then F(T,0)=-QT and T_0=0.

All residual choices in this construction are regular semisimple at
T_0, since T_0!=+/-2. The construction is in the matrix image, and
solutions lift through the finite etale Spin cover. At the starting
target zero, the matrix need only be K-regular, not G-regular; the
argument imposes Z_G(h)=K only after proving target dominance.
This distinction is essential. This proves the sufficiency in Theorem B.

For clarity, the D rank restriction is not disposable. Spin_6 is A_3;
its GL_3 Levi and its GL_1 times SO_4 Levi correspond to excluded
type A block patterns. The two failures of the balanced construction
above must not be silently advertised as D_3 nonemptiness cases.

### 7.1 A direct construction for four-blocks

Take t>=1 GL_4 blocks. For generic L_i use all four roots of (8) in
block i. In B/D their combined trace contribution is 2tT; in C it
is -2tT. Choose residual traces and T as follows:

| Family | m | Chosen R | T |
| --- | --- | --- | --- |
| B | 0 | 1 | 1/(1-2t) |
| B | 1 | -1 | 1/(2t-1) |
| B | 2 | T-1 | 1/(2t) |
| C | 0 or 1 | 0 or -T | 0 |
| D | 0, 2 or 3 | 0, 0 or T | 0 |

All these values avoid +/-2. The roots in each GL block are distinct
for generic L_i, and cannot be reciprocal partners since the reciprocal
polynomial corresponds to -L_i. Distinct target coordinates up to sign
exclude coincidences between blocks. The table supplies a regular
residual matrix. Embed the GL blocks with their duals; this enforces
the ambient form and determinant, and (12) ensures that the chosen T
is the actual defining trace. Equations (6)-(8) give exactly the target h.
Lift to Spin if applicable. Equivariance then gives G-regularity from
K-regularity. This independently constructs the four-block subfamily over
the algebraic closure of the generic target field. The torus case
instead follows from d theta_e=identity and torus dominance.

The next question is structural extension to representations other than
Ad. The proof used two distinct mechanisms: quadratic descent with
controlled square branching for smoothness, and a low-degree spectral
equation with a trace constraint for nonemptiness and dimension. An
extension must identify analogues of these mechanisms; the adjoint
classification alone is not evidence for an arbitrary-representation
theorem. Exceptional types require their own structural input.

## 8. Audit, sources and relation to the program

The [symbolic and matrix audit](tools/verify_classical_adjoint_fibers.py)
checks the Vieta identities, the square-substitution equations and a
nonreduced failure control. It verifies the residual B1/D2 characters
from explicit torus weights, checks the cubic bad-locus ideals in one
minimal model of each mechanism, verifies the symbolic deformation
coefficients (including the exceptional root choices), and computes the trace projection
independently in exterior/symmetric-square matrices. See its
[source-bound report](classical-adjoint-fibers-verification.json).
These bounded checks do not replace the variable-rank arguments above.

Published inputs include the [Kostant--Michor trace-projection construction](https://www.mat.univie.ac.at/~michor/kostant.pdf),
classical tensor-square character identities and regular Jordan forms,
and characteristic-zero generic smoothness. The existing
[critical presentation](GENERAL_CAYLEY_FIBER_CRITICAL.md) and
[quadratic theorem](QUADRATIC_BRANCHING_FIBERS.md) are repository
deductions used here. The square-substitution criterion, cubic exclusion,
dimension and deformation arguments are written deductions in this note. Their
originality is not established. Targeted searches on 2026-10-23 for
generalized Cayley adjoint fibers and quadratic regular Cayley fibers
did not identify these exact statements; that is not a complete priority audit.

The results concern the same Springer/Cayley map, not another cohomology
theory, integral Smith factors, or a new minimum-Cayley-degree computation.
Step 1 of the plan is closed by Theorem A's smoothness assertion. Step 2
has a complete generic dimension formula and an all-rank generic
nonemptiness criterion, including the mixed block patterns. Step 3
yields the reusable square-substitution criterion;
no exceptional-group conclusion or arbitrary-representation classification
is inferred from it. Whole fibers, special-target singularities and orbit
counts remain separate.
