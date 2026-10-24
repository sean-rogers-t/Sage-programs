# E7 refinements and a symmetric-pencil rank criterion

Characteristic zero; 2026-10-23. This note concerns the actual adjoint
Springer map of the simple complex simply connected group E7. Standard
inputs and written deductions are separated below; priority is unverified.

## 1. A criterion retaining all rank strata

Let T=(C*)^s, u in A^m, x in A^1, and

    phi=g(x,z)+(1/2)u^t Q(x,z)u,
    g=a(z)x^2+b(z)x+c(z), Q=Q_0(z)+x Q_1(z),                (1)

where Q is symmetric, all coefficients are Laurent polynomials, and a
is a Laurent unit. Set e=g_x and d=det Q. Keep the actual scheme

    X=V(Qu, e+(1/2)u^t Q_1 u), mu=D_z phi.

Dimensions below are dimensions of supports. The base for d and e is
A^1_x times T, of dimension s+1. Assume

    dim V(d,e) < s,
    dim V(d,d_x,e)+1 < s,
    dim {points of X with rank Q <= m-2} < s.               (2)

**Theorem 1.** Under (2), the generic fibers of mu are smooth and finite,
possibly empty. It is enough for each bad set in the proof to have
nondense target image; (2) gives readily checkable sufficient bounds.
The ordinary branch is

    u=0, x=-b/(2a), mu=D(c-b^2/(4a)), det Q!=0.             (3)

If its scalar gradient has an invertible logarithmic Hessian at one
point where det Q!=0, generic fibers are nonempty.

**Proof.** On det Q!=0, the equations give u=0 and e=0 exactly. The
ordinary branch is a smooth open of T. At a point of corank one, a
nondegenerate symmetric block can be eliminated by congruence, leaving

    phi=g(x,z)+(1/2)k(x,z)t^2

and zero for the other u coordinates. Here d is a unit times k. On k=0,
k_x is nonzero exactly when d_x is nonzero. The remaining equations are

    kt=0, e+(1/2)k_x t^2=0.                                (4)

At d_x!=0 and e!=0, t is nonzero, so kt=0 implies k=0 as a scheme
equation. The surface k=0 is etale over T; (4) then gives a nonzero square
root. This stratum of X is smooth of dimension s. The congruence may
depend on x,z: the extra derivatives vanish modulo the eliminated
relative equations, so these are the original equations and target map.

If d_x!=0 and e=0, the support has t=0 and maps to V(d,e), of dimension
less than s. If d_x=0, the relative equations require e=0; their support
has at most one free t and hence dimension at most dim V(d,d_x,e)+1.
The third condition handles every lower rank. Excluding the closures
of these bad target images retains the full generic fiber, including its
scheme structure. The remaining source is smooth of dimension s.
Characteristic-zero generic smoothness gives finite smooth fibers or
emptiness on a dense target open. A dominant ordinary branch proves
nonemptiness after shrinking that open. This proves the theorem.

In particular Q_1 need not be invertible and no global eigenbasis or
squarefree characteristic polynomial is assumed. This extends the useful
scope of the earlier [spectral criterion](SPECTRAL_BRANCHING_FIBERS.md),
but still imposes substantive bounds. It is not a theorem for every
quadratic pencil: persistent collisions can give nonreduced generic fibers.

## 2. The complete A1^4 E7 family

**Theorem 2.** Let K be an E7 Levi of the cubic A1^4 Weyl class, represented
by retained nodes 1257. For generic h in z(k) with Z_G(h)=K,

    theta_Ad^(-1)(h) intersect G_reg

is nonempty, smooth and pure of dimension 8. The two standard embeddings
are 1257 and 2357, related by the previously verified Weyl word (3,1).
This theorem is about regular-source fibers at generic central targets,
not global preservation of regularity or a full-map generic degree.

On the central cover let x,y,z,w be the A1 fundamental characters at
nodes 1,2,5,7. Write the original central torus parameters as p,q,r, for
omitted nodes 3,4,6, and make the unimodular change

    U=pq, V=pqr, W=pq^2r,
    A=U+U^(-1), B=V+V^(-1), C=W+W^(-1).

The FULL restricted adjoint character is

    phi=x^2+y^2+z^2+w^2
       +Axyz+Bxyw+Cxzw+BCyz+ACyw+ABzw
       +ABCx+A^2+B^2+C^2-7.                                (5)

Thus (1) holds with u=(y,z,w),

    g=x^2+ABCx+A^2+B^2+C^2-7,
    Q=[[2,BC+Ax,AC+Bx],[BC+Ax,2,AB+Cx],[AC+Bx,AB+Cx,2]].   (6)

It is useful to prove the dimension bounds over the independent affine
parameters A,B,C. The map (U,V,W)->(A,B,C) is finite, so its base change
preserves the dimension bounds. Its ramification does not invalidate
the theorem: on the good stratum d_x!=0 the spectral equation remains
etale over the ORIGINAL torus, and all original target derivatives are
retained. No derivative in A is silently substituted for D_U.

### 2.1 Collisions and lower ranks

Here e=2x+ABC. Let

    P(A,B,C)=d(-ABC/2,A,B,C),
    H(A,B,C)=d_x(-ABC/2,A,B,C).

Exact polynomial arithmetic gives P!=0 and gcd(P,H)=1. In particular
P(2,2,2)=-216. For reproducibility,

    H=(ABC/2)(3A^2B^2C^2-4A^2B^2-4A^2C^2-4B^2C^2
               +8A^2+8B^2+8C^2-24).                       (7)

The complete P is defined by the displayed 3 by 3 determinant (6), so
no omitted table is needed to recover it. Coprimality in Q[A,B,C] gives
dim V(d,e)=2 and dim V(d,d_x,e)<=1, also over C. These are the first
two bounds in (2) for s=3.

Q has nonzero diagonal, so rank zero is impossible. At rank one its
off-diagonal entries are 2 epsilon_1, 2 epsilon_2,
2 epsilon_1 epsilon_2. Sign changes in A,B,C and in the leaf coordinates
reduce this to all off-diagonal entries equal to 2; they preserve g.
The rank-one base is therefore determined by

    Ax+BC=Bx+AC=Cx+AB=2.                                   (8)

Subtracting equations yields (A-B)(x-C)=0 and its permutations.
Consequently (8) is the union of four curves: either A=B=C=t and
x=2/t-t, or two of A,B,C equal x=t and the remaining one equals 2/t-t.
Here t!=0; there are no missing t=0 solutions. Thus its dimension is one.

On the normalized rank-one locus, ker Q is y+z+w=0. The remaining
quadratic term in the x-equation restricts to

    Ayz+Byw+Czw=-B y^2-C z^2+(A-B-C)yz.                    (9)

It cannot vanish identically: that would force A=B=C=0, contradicting
(8). The x-equation therefore cuts a curve (or the empty scheme) in
this two-dimensional kernel for every base point. The entire low-rank
critical locus has dimension at most 1+1=2<3. This proves the last
condition in (2), including intersections of the four rank-one curves.

### 2.2 Nonemptiness, dimensions and normalization

The ordinary potential is

    R_0=A^2+B^2+C^2-7-A^2B^2C^2/4.                         (10)

At U=V=W=1, A=B=C=2, x=-4, y=z=w=0. The relative Hessian is
invertible: det Q=-216, so its determinant is -432. The logarithmic
Hessian of R_0 in U,V,W is -24 times the identity. This proves dominance
of the ordinary branch. The witness is used for quotient dominance;
regular lifts at generic centralizer-K targets follow from the
[regular-quotient theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md).
The smooth orbit factor has dimension dim K-7=4*2=8, proving Theorem 2.

For h=sum h_i H_i in the original primitive central cocharacters, the
target is c_i=B_Ad(h,H_i), with trace matrix computed from the 126 roots
in the verification report. In the unimodular torus change above,

    (c_p,c_q,c_r)^t = [[1,1,1],[1,1,2],[0,1,1]]
                     (D_U R_0,D_V R_0,D_W R_0)^t           (11)

on the ordinary branch. The other branches retain the derivatives of
the full character (5). Finite central descent and Weyl conjugacy
preserve the smoothness and dimension conclusions.

## 3. Status of the companion A2+A1^3 refinement

For retained nodes 12357, eliminate the globally paired A2 endpoint
characters as in the earlier quartic reduction. With x,y,w the remaining
A1 characters, put U=q, V=p^3q and

    A=U+U^(-1), B=V+V^(-1), C=UV+(UV)^(-1).

The exact residual is

    psi=x^2+y^2+w^2-x^2y^2-x^2w^2-y^2w^2
       -xyw(Ax+By+Cw)-ABxy-ACxw-BCyw-ABC-3,                (12)

with A^2+B^2+C^2-ABC=4. Unlike (5), these three trace coefficients are
not independent. This difference prevents copying the matrix-pencil
proof verbatim. The map from the two-dimensional torus to this Cayley
cubic surface is the quotient by simultaneous inversion. Its four
ramification points have U,V in {+1,-1}; all three trace differentials
vanish there, so their entire critical fibers have zero central target.
This torus-quotient description is standard: see
[Rebelo--Roeder, section 4](https://arxiv.org/pdf/2307.10962), allowing
their coordinate signs and replacing one torus variable by its inverse.
The surface itself is not a new construction in this work.

The remaining task is to certify a dimension bound for the nonsmooth
locus of the relative-critical scheme of (12) on the smooth part of
the cubic surface. A computed dimension alone is insufficient for the
independent proof audit. This section records the exact reduction, not
a completed generic-smoothness theorem for this companion family.

For a concrete continuation, form the ideal generated by psi_x,psi_y,
psi_w,h and the four by four minors of their Jacobian in x,y,w,A,B,C,
where h=A^2+B^2+C^2-ABC-4. The
[reproducible computation](tools/explore_e7_cayley_quartic.py) gives
dimension one over Q using Singular's modular standard-basis algorithm
with its exactness option, with a separate prime-101 corroboration.
The [exploration record](e7-quartic-boundary-exploration.json) deliberately
does not promote this to an independently audited companion theorem.
Six selected polynomials in the computed basis have leading monomials

    x^5, y^5, w^2 B^4, w^3 C^3, A^6, B^5 C.

Their membership in the ORIGINAL Jacobian-minor ideal is the precise
remaining certificate. Once those six memberships are independently
verified, no independent Groebner-basis computation is needed: the
initial ideal contains these monomials, whose support ideal is
(x,y,A,wB,wC,BC), of dimension one. The candidates and original
generators are saved in a hash-bound archive. Reduction of the original
generators by the computed basis verifies the opposite inclusion and
would not, by itself, close this audit gap.

That bound would exclude the bad locus from generic two-dimensional
targets. On its complement, the four equations in six variables have
full Jacobian rank, giving a smooth two-dimensional incidence scheme.
The ordinary critical branch x=y=w=0 has torus potential -ABC-3; at
U=V=1 its relative Hessian has determinant -216 and its logarithmic
torus Hessian is [[-16,-8],[-8,-16]], of determinant 192. Thus dominance
would give nonempty generic fibers and the orbit factor would have
dimension 6+2+2+2=12. These are conditional consequences of the pending
certificate, not additions to the established coverage count.

## 4. Attribution, checks and next step

The trace-gradient construction is the published generalized Cayley map
of Kostant--Michor. Regular quotient geometry, symmetric congruence,
dimension of images and generic smoothness are standard inputs. The
rank-stratum criterion and its E7 application are written deductions;
no literature-wide novelty assertion is made.

The [Sage verification](tools/verify_e7_symmetric_pencil.py) and
[native root/SymPy audit](tools/audit_e7_symmetric_pencil.py) independently
check the full character, determinant coprimality, rank-one identities,
dominance witness, original torus transformation and trace normalization.
The original target equations are retained throughout.
The records are [verification](e7-symmetric-pencil-verification.json)
and [independent audit](e7-symmetric-pencil-audit.json).

The next gate is the independent dimension certificate for (12). Its
purpose is to establish that the bad critical locus has nondense target
image, the same geometric mechanism used in Theorem 1. The A1^4 result
adds two standard embeddings, bringing E7 coverage to 101 of 127;
26 cubic embeddings in four Weyl families remain unresolved here.
