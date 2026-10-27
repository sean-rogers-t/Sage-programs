# All proper Levi strata for the F4 adjoint Springer map

26 October 2026. Actual adjoint character of complex simply connected F4;
Killing trace normalization. Case proof with independent exact arithmetic;
no literature-priority claim.

**Theorem.** For every proper Levi K of F4 there is a dense open subset
U_K of z(k), consisting of h with Z_G(h)=K, such that the actual scheme
theta_Ad^(-1)(h) intersect G_reg is smooth for every h in U_K. Empty fibers
are allowed. This theorem does not classify which strata are nonempty.

## 1. Exhaustive reduction

Use Bourbaki nodes 1,2 long and 3,4 short. The full character is expressed
in derived fundamental characters and the primitive central cocharacters
on the finite central cover. Root--coroot pairings are used directly,
avoiding a transpose error in nonsimply-laced Cartan conventions.

There are fifteen proper standard subsets. The ten subsets

    empty, 1, 2, 3, 4, 12, 23, 34, 123, 234

have characters of degree at most two in the derived fundamental
characters, with all central Laurent coefficients retained. They are
covered by the [quadratic theorem](QUADRATIC_BRANCHING_FIBERS.md), which
allows rank drops, nilpotents in incidence and smooth excess fibers.
The remaining subsets are 13,14,24,124,134. The formulas below handle all
five. Character identities are checked against all 48 ambient roots and
four zero weights, independently by Sage roots and native Cartan
reflections with exterior/spin models for the Levi fundamental characters.

## 2. The three orthogonal long/short A1 pairs

Let x be the long-root and y the short-root fundamental character. The
full potential in each of 13,14,24 is

    f=x^2+x*[A*y^2+B*y+A*(B^2-3)]
       +(B^2-1)*y^2+A*B*y+A^2-B^2.                     (1)

Here A=u+u^(-1), B=v+v^(-1). For the primitive central coordinates p,q
in the certificate, (u,v) are respectively

    (p,pq), (pq^2,pq), (pq^2,q).

Each change is a torus isomorphism. Complete the unit square in x. If
L=A*y^2+B*y+A*(B^2-3), the exact residual is obtained by x=-L/2.
Twice its y-derivative is

    J=-2*A^2*y^3-3*A*B*y^2
      +(-2*A^2*B^2+6*A^2+3*B^2-4)*y+A*B*(5-B^2).       (2)

It is primitive as a polynomial in y over C[u^+/-1,v^+/-1]: any common
factor of its coefficients must divide A, from the leading coefficient;
on either divisor u=+i or -i, its y coefficient is 3*B^2-4, not identically
zero. At u=2,v=3 (A=5/2,B=10/3), it specializes, without losing degree, to

    -25*y^3/2-25*y^2-1297*y/18-1375/27,

which is coprime to its derivative. Hence J is squarefree over C(u,v),
and primitivity makes it squarefree in the full Laurent polynomial ring.
The ACTUAL incidence is therefore a reduced hypersurface of pure dimension
two. Its nonsmooth locus has dimension at most one, so cannot dominate the
two-dimensional target. No trace-cover ramification divisor is removed.
Generic smoothness on the remaining smooth source proves the assertion
for these three cases, allowing components with nondominant target image.

## 3. Long A2 plus short A1: retained nodes 124

Write a,b for the A2 fundamental characters, x for A1, and t for the
primitive central coordinate. The full character is

    f=ab+a*(x^2/t^2+xt+t^4-t^(-2))
         +b*(x^2*t^2+x/t-t^2+t^(-4))
         +x^2+x*(t^3+t^(-3))-1.                        (3)

Eliminate the unit-paired a,b globally. With S=t^3+t^(-3), minus the
remaining derivative is

    M=4*x^3+3*S*x^2+2*(S^2-4)*x-S.                     (4)

Its leading coefficient is the unit 4. At t=2 it is the squarefree cubic
4*x^3+195*x^2/8+3969*x/32-65/8. Thus its actual Laurent hypersurface is
a reduced curve. A reduced complex curve has only finitely many nonsmooth
points; their target images can be excluded. Generic smoothness on its
smooth part gives the required generic smooth fibers, including possible
emptiness. There is no eliminated source divisor in this argument.

## 4. Long A1 plus short A2: retained nodes 134

Let x belong to A1 and y,z to A2, with primitive central coordinate t.
The full character has a quadratic block in v=(y,z):

    f=Q+(1/2)v^t H v+ell^t v,
    Q=x^2+(t^3+t^(-3))*x-1,
    H=[[2*(x/t+t^2),1],[1,2*(xt+t^(-2))]],
    ell=(-xt-t^(-2),-x/t-t^2)^t,
    D=det H=4*x^2+4*(t^3+t^(-3))*x+3.                  (5)

On D!=0 eliminate v exactly. Set R=Q-ell^t adj(H)ell/(2D). Then

    N=t^9*D^2*R_x

is a polynomial of x-degree five, with unit Laurent leading coefficient
32*t^9. At t=2 its coefficients, in descending x-degree, are

    16384, 316160, 1917376, 4606160, 1063168, 1061060.

This quintic is coprime to its derivative. Therefore the actual incidence
on D!=0 is a reduced curve, handled as in the previous case.

It remains to retain D=0. The matrix has rank one everywhere on this
divisor, since its off-diagonal entry is 1. Consistency requires
adj(H)ell=0. For k=t^3*D and c1 the numerator of its first component,
the resultant in x is, with the numerator convention used in the audit,

    -4*t^6*(t^6-1)*(2*t^4+t^2+2)
        *(4*t^8-2*t^6-3*t^4-2*t^2+4).                  (6)

It is nonzero. Consequently the consistent base (x,t) is finite. Where
D_x!=0, eliminating the nonzero block leaves an x-equation quadratic in
the kernel variable, with nonzero quadratic coefficient. The full critical
boundary there is finite, including nonreduced points.

There are NO consistent points where D_x=0. Here x=-(t^3+t^(-3))/2 and
t^12-t^6+1=0. Write H=[[a,1],[1,b]]; then

    a=(t^6-1)/t^4, b=-(t^6-1)/t^2,
    ab=1, ell=(-b/2,-a/2).

Consistency implies a=b^2 and b=a^2, hence a^3=b^3=1. Since b/a=-t^2,
this forces t^6=-1, contradicting t^12-t^6+1=0 (its value would be 3).
This covers every singular-block point, without an unjustified unit.
Its finite target image can be excluded, completing the fifth cubic case.

## 5. Scheme transfer, checks and scope

Every elimination above is an isomorphism of the actual critical schemes
on its stated chart. The original logarithmic central derivatives become
those of the residual by the chain rule; the eliminated affine derivatives
vanish. For nodes 134 the entire missing chart D=0 was analyzed separately.
After excluding the indicated proper closed target sets, generic smoothness
of the smooth source gives a dense good target open. Persistent nilpotents
are not discarded: the surviving hypersurfaces were proved reduced, and
every excluded scheme component is supported over the excluded target set.

The [regular-quotient theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md), finite
central descent and conjugacy transfer this statement to actual regular
group fibers. Intersect with the nonempty centralizer-exact target open.
Together with the ten quadratic cases this proves the theorem for all
proper F4 Levis. Genuinely special central targets may still be singular.

The [Sage verification](tools/verify_f4_adjoint_fibers.py) and
[native independent audit](tools/audit_f4_adjoint_fibers.py) check full
characters, the cubic models, squarefreeness and boundary identities in
characteristic zero. Root characters, classical fundamental-character
models, the quadratic theorem and generic smoothness are existing inputs.
The finite case completion is a repository deduction, not a claim that
the statement is new in the literature. E7 is completed separately in
[its final-family note](E7_FINAL_FAMILIES.md); E8 remains outstanding.
