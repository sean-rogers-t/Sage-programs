# All proper Levi strata for the G2 adjoint Springer map

26 October 2026. Complex simply connected G2, actual adjoint representation,
Killing trace normalization. This completes the proper-Levi coverage audit;
it does not claim literature priority or a new full-map generic degree.

**Theorem.** For every proper Levi K of G2 there is a dense open subset of
z(k), consisting of targets with centralizer exactly K, on which the actual
scheme theta_Ad^(-1)(h) intersect G_reg is nonempty and smooth. Its dimension
is zero for K=T and two for either rank-one Levi.

## Full restricted characters

Use a short simple root alpha and a long simple root beta. Their positive
roots in this basis are (1,0),(0,1),(1,1),(2,1),(3,1),(3,2). The root--coroot
pairing matrix, with roots indexing rows, is [[2,-1],[-3,2]]. This orientation
matters: it is the transpose of Sage's Cartan-matrix convention here.

For the short-root Levi take central cocharacter alpha^vee+2 beta^vee;
for the long-root Levi take 2 alpha^vee+3 beta^vee. Both cover maps from
SL2 times C* have finite kernel of order two. Let x be the fundamental SL2
character, z the central coordinate, and S=z+z^(-1). The ENTIRE restricted
adjoint characters are, respectively,

    f_short=x^2+S*(x^3-2*x)+S^2-2,
    f_long =x^2+(S^3-2*S)*x+S^2-2.                       (1)

These identities follow by summing the twelve root weights and two zero
weights. In particular no central-only or linear terms are omitted. The
target coordinate is c=B_Ad(h,H)=z*f_z. The primitive central trace norms
are B_Ad(H,H)=16 for the short-root Levi and 48 for the long-root Levi.
The quotient incidence is the actual equation f_x=0 in A^1 times C*.

## Short-root Levi: the incidence is already smooth

Write F=3S*x^2+2*x-2S. A singular point of F=0 would require

    F_x=6S*x+2=0,
    z*F_z=(z-z^(-1))*(3*x^2-2)=0.                    (2)

If 3*x^2-2=0, F=0 implies 2*x=0, a contradiction. Otherwise (2) forces
z=+1 or -1 and S=+2 or -2. But the simultaneous equations F=F_x=0
force S^2=-1/6, again a contradiction. Thus the entire incidence is a
smooth curve, including all points where f_xx vanishes.

At (x,z)=(-1,1), F=0, f_xx=-10, the mixed x/log(z) Hessian entry vanishes,
and the logarithmic central Hessian equals 10. Hence its central target
map is dominant on a component. A nonconstant map from a smooth complex
curve to A^1 has smooth finite fibers on a dense target open; any constant
components contribute only finitely many targets, which can be excluded.
This is a statement about the actual incidence, not its reduced support.

## Long-root Levi and torus

Here f_x=2*x+S^3-2S, so the incidence is globally the graph
x=-(S^3-2S)/2 over C*. The target is the logarithmic derivative of

    R=S^2-2-(S^3-2S)^2/4.

At z=1 its logarithmic second derivative is -32, proving dominance.
The same curve argument gives nonempty smooth finite generic quotient
fibers. For K=T there are no relative critical equations: the incidence
is T itself. The logarithmic Hessian of the full adjoint character at the
identity is the nondegenerate Killing matrix, so the target map is dominant
and generically etale. Generic fibers are nonempty, finite and smooth.

For both rank-one Levis, exclude c=0 to enforce centralizer-exactness.
For the torus exclude the root hyperplanes. The established
[regular-quotient theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md) then adds the
smooth regular orbit factor, of dimension two for either A1 Levi. It retains
regular nonsemisimple points and the actual fiber scheme. The finite central
covers are etale, so smoothness descends. Every proper Levi is conjugate to
one of these three standard Levi subgroups, proving the theorem.

## Checks and attribution

The reduction to character gradients and regular quotient geometry are the
standard inputs already recorded in the program. The elementary character
and smooth-curve argument above completes this group's coverage. The
[verifier](tools/verify_g2_adjoint_fibers.py) compares Sage ambient roots with
an independent enumeration from the displayed simple-root coefficients,
checks both Laurent characters, trace normalization, singular-incidence
ideals over QQ and dominance Hessians. This is an exact finite check
supporting the written argument, not an assertion of originality.
