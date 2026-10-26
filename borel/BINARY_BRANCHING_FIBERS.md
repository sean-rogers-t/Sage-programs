# A binary critical-equation criterion and two E7 Weyl families

Date: 2026-10-24. Characteristic zero. The general statement is a
sufficient criterion for specified characters, not an all-simple-group
classification. Generic smoothness, hypersurface geometry and polynomial
separability are standard inputs. Literature priority is unverified.

## 1. A common criterion

Let Z be a product of affine spaces and algebraic tori of dimension s>=2,
with coordinate ring R over C. Let

    phi(x,y,z)=a(y,z)x^2+b(y,z)x+c(y,z),   a,b,c in R[y],
    X=V(2ax+b, a_y x^2+b_y x+c_y) in A^2 x Z.

Put

    P=4a^2 c_y-2ab b_y+b^2 a_y,
    E=V(a,b,a_y,b_y,c_y) in A^1_y x Z.                    (1)

Assume:

1. a is nonzero and gcd(a,b)=1 in the UFD R[y].
2. P is nonzero and squarefree in R[y,1/a].
3. dim E<=s-2 (empty sets satisfy this bound).

**Theorem 1.** Outside a closed subset of X of dimension at most s-1,
X is smooth of pure dimension s. Consequently any morphism
mu:X -> A^s has smooth finite generic fibers, possibly empty.
Dominance is an additional condition for nonemptiness.

Proof. On a!=0 eliminate x=-b/(2a) exactly. The remaining potential is
R_0=c-b^2/(4a), and its derivative is P/(4a^2). Thus this open part of
X is the hypersurface V(P) in Spec R[y,1/a], with its actual scheme
structure. Squarefreeness makes it reduced, of pure dimension s. Its
nonsmooth locus has dimension at most s-1 in characteristic zero.

On a=0 the first equation also gives b=0. Coprimality implies
dim V(a,b)<=s-1. Except over E, the remaining equation in x is a
nonzero polynomial of degree at most two, so its fibers are finite or
empty. Over E at most one free x remains. Hence

    dim(X intersect {a=0}) <= max(s-1,dim E+1) <= s-1.      (2)

Remove the closures of the target images of this boundary and the
nonsmooth part of the hypersurface. All remaining source points lie
in a smooth s-dimensional scheme. Characteristic-zero generic
smoothness, with nondominant components excluded, proves the assertion.
This removes supports of bad loci, not their radicals from the original
equations; no nilpotent boundary structure is used as a good component.

The criterion also applies on a finite cover of Z that is etale away
from a locus whose full target image is nondense. The dimension bounds
survive finite pullback. A ramified divisor cannot simply be deleted:
its target-image assertion is a separate required hypothesis.

### Checking squarefreeness without a large critical-locus basis

It suffices that P be primitive in R[y] and squarefree over Frac(R).
A specialization that preserves degree in y and gives a squarefree
polynomial proves the latter: a repeated positive-degree factor would
retain its degree and repetition at that specialization. Primitivity
excludes repeated factors depending only on the parameters. A successful
specialization alone, without this content check, is insufficient.

This uses the geometry of a reduced hypersurface, not interpolation of
Smith factors or point counts. It is valid in any polynomial degree in y.
The two applications below have degree five after elimination.

## 2. Torus traces and the actual Springer targets

Work with the adjoint Springer map of simple complex simply connected E7,

    B_Ad(theta(g),H)=tr(Ad(g) ad(H)).

Here B_Ad(Y,H)=tr(ad(Y)ad(H)). Write G_reg for the elements whose
group centralizer has dimension rank G=7. In Theorems 2 and 3, a
"generic regular-source fiber" means the actual scheme

    theta^(-1)(h) intersect G_reg,
    h in a nonempty Zariski-open subset of {h in z(k): Z_G(h)=K}.

For a proper Levi K use a finite central cover [K,K] x T -> K. Its
fundamental-character and central-torus coordinates describe K//K
etale-locally. Target coordinates are c_j=B_Ad(h,H_j), using the
cocharacters of the chosen cover. They are independent linear coordinates
on z(k); changing a finite-index cocharacter basis changes them invertibly
over C. All original logarithmic derivatives are retained.

When a commuting A1 factor is replaced by its torus, its character is
w=t+t^(-1) and

    D_t phi=(t-t^(-1)) phi_w.                              (3)

Thus every relative-critical point on t=+/-1 has c_t=0, including
nonreduced or positive-dimensional branches. For generic c_t these
ramification points are absent. Away from them t -> w is finite etale.
The same statement applies to each of several independent A1 traces.
This is the elementary Weyl-reflection target constraint; it does not
assert smoothness transfers automatically under Levi refinement.

After Theorem 1, the regular-quotient presentation adds dim K-7 smooth
directions. Smoothness and local dimensions descend through the finite
central cover. Neither the degree-five critical equation nor its number
of roots at fixed central parameters is a full Springer-map degree or
a count of regular K-orbits at fixed target.

For completeness, this passage uses the trace-critical presentation in
the [general fiber theorem](GENERAL_CAYLEY_FIBER_CRITICAL.md). Equivariance
gives g in Z_G(h)=K whenever theta(g)=h. On K_reg the quotient K -> K//K
is smooth and surjective, of relative dimension dim K-7. Vanishing of
the derived-coordinate character derivatives and the prescribed central
logarithmic derivatives is exactly theta(g)=h there. Finally
Z_G(g) is contained in Z_G(theta(g))=K, so such a K-regular solution is
G-regular. Conversely a G-regular solution in K is K-regular. Thus this
presentation captures precisely the regular-source fiber, and quotient
nonemptiness supplies regular solutions. It makes no assertion about
additional nonregular solutions.

## 3. The cubic A1^3 class

Take K represented by retained nodes 125, in Bourbaki numbering. There
are ten standard embeddings in its cubic Weyl class. Use a cover with
seven mutually orthogonal coroot directions, three derived and four
central. Put A=U+U^(-1), B=V+V^(-1), C=W+W^(-1), w=t+t^(-1).
In derived A1 characters X,x,y, the full adjoint character is

    Phi=X^2+x^2+y^2+w^2+A Xxy+B Xxw+C Xyw
        +BCxy+ACxw+AByw+ABCX+A^2+B^2+C^2-7.               (4)

It follows by restricting the earlier A1^4 character; both implementations
also check it against the original E7 roots. The seven coroot columns,
in the order X,x,y,w,U,V,W, are

    e1, e2, e5, e7,
    (1,1,2,2,1,0,0), (1,1,2,2,2,2,1), (1,2,2,4,3,2,1).

Their Gram matrix is 2I and their lattice index is 8. They specify the
actual cover, not a replacement of E7 by its adjoint isogeny form.

Eliminate the globally nondegenerate square in X. Set

    L=Ay+Bw, M=Cyw+ABC,
    a=1-L^2/4,
    b=BCy+ACw-LM/2,
    c=y^2+w^2+AByw+A^2+B^2+C^2-7-M^2/4.                 (5)

The residual is a x^2+b x+c and X=-(Lx+M)/2. The coordinate change is
invertible over the entire parameter space, so all target derivatives
agree modulo the eliminated equation, by the chain rule.

Apply Theorem 1 over Z=A^4_(A,B,C,w), s=4. The following give small
exact certificates for all its hypotheses:

- a is primitive in y. At A=B=C=1,w=0, a=1-y^2/4 and b=y/2 are
  coprime, and the degree of a stays two. Thus gcd(a,b)=1 globally.
- The coefficient of y^5 in P is A^4/2. At A=B=C=w=0, P=8y.
  Any common parameter factor of the coefficients would divide A,
  which the second identity excludes. Hence P is primitive.
- At A=B=C=1,w=0,

      P=y(y^2-3)(y^2-5)/2.

  Its degree is five and its roots are distinct, proving separability
  over the parameter function field and therefore squarefreeness.
- On a=a_y=0, necessarily A=0 and Bw=+/-2. Each such parameter locus
  (with y included) has dimension three. There

      b_y=-CB(w^2-2)/2

  is not identically zero on either irreducible component. Thus E has
  dimension at most two, as required. No division by A is made.

Pull back to the four central torus coordinates and use (3) to exclude
the four ramification loci for generic targets. Theorem 1 gives finite
smooth generic quotient fibers.

For nonemptiness take X=x=y=-4 and A=B=C=w=2 (U=V=W=t=1).
The relative Hessian of (4) has determinant -216; in (5), a=-3.
At this point (Phi_A,Phi_B,Phi_C,Phi_w)=(-108,36,36,36).
All first derivatives of the trace coordinates in logarithmic torus
coordinates vanish. Hence the logarithmic Hessian on the local relative
critical branch is diagonal (-216,72,72,72), and is invertible.
This proves dominance. Although the witness lies over target zero and
on the trace ramification locus, dominance is an open condition on this
smooth local branch. It therefore also holds away from those divisors.
Only after establishing dominance do we restrict to targets with Z_G(h)=K.

**Theorem 2.** For this cubic A1^3 Weyl class, generic central targets
with centralizer K have nonempty smooth regular-source adjoint Springer
fibers of pure dimension 6. All ten standard embeddings are covered by
the verified Weyl words. Other embeddings of the same abstract Dynkin
type are not identified with this class merely by their type.

## 4. The cubic A2+A1^2 class

Take K represented by retained nodes 1235. This class has twelve standard
embeddings. Restrict the exact A2+A1^3 character at nodes 12357 by
replacing the node-7 character w by t+t^(-1). Globally eliminate the
paired A2 endpoint characters. In the remaining A1 characters x,y,
the full residual potential is

    psi=x^2+y^2+w^2-x^2 y^2-x^2 w^2-y^2 w^2
        -xyw(Ax+By+Cw)-ABxy-ACxw-BCyw-ABC-3,              (6)
    A=U+U^(-1), B=V+V^(-1), C=UV+(UV)^(-1).

Use Z=(C*)^2_(U,V) x A^1_w, s=3. The torus variables are independent;
A,B,C are not independent. In the original central variables p,q,
U=q, V=p^3q, a finite etale isogeny. Consequently

    D_p=3 D_V,  D_q=D_U+D_V.                              (7)

The additional central direction is e7. The p,q cocharacters are
(4,6,8,12,9,6,3) and (2,3,4,6,5,4,2), respectively. Together with
e1,e2,e3,e5,e7 these specify the cover. The two eliminated A2 coordinates
are recovered on this cover; no unsupported descent of their graph is used.

Writing (6) as a x^2+b x+c gives

    a=1-y^2-w^2-Awy,
    b=-By^2w-Cyw^2-ABy-ACw,
    c=y^2+w^2-y^2w^2-BCyw-ABC-3.                          (8)

Again all three hypotheses of Theorem 1 have short certificates:

- a has leading y coefficient -1. At U=V=1,w=0, a=1-y^2 and
  b=-4y, so a,b have no common factor.
- The y^5 coefficient of P is

      2((V-V^(-1))^2 w^2+4).

  It is independent of U. The y^0 coefficient of P, as a Laurent
  polynomial in U, has highest U power U^5 with coefficient -V^2 w^3.
  Any nonunit content must therefore divide both the displayed leading
  coefficient and w^3. Their gcd is one in C[V^(+/-1),w]. Thus P is
  primitive in C[U^(+/-1),V^(+/-1),w][y].
- At U=V=1,w=0,

      P=8y(y^2-3)(y^2+1),

  a squarefree degree-five polynomial. Hence P is squarefree globally.
- On a=a_y=0 one has y=-Aw/2 and

      (U^2-1)^2 w^2+4U^2=0.

  Thus U^2!=1 and w!=0. Over C this locus has two irreducible surfaces,
  w=+/-2iU/(U^2-1), with free U,V and y determined. Modulo its equation,

      b=-w(U^2+1)(U^4+1)(V^2-1)/(2U^2 V(U^2-1)).          (9)

  This is nonzero on each surface. The further condition b=0 cuts their
  dimension to at most one. In particular dim E<=1=s-2.

Use (3) for the single w-trace ramification locus. No passage to the
singular Cayley cubic surface is needed: U,V stay torus coordinates.
This argument does NOT use the unfinished A2+A1^3 membership certificate.

At x=y=-4 and U=V=t=1, w=2, the relative Hessian has determinant -108
and a=-3. In trace parameters (psi_A,psi_B,psi_C,psi_w)=(108,108,-36,36).
The logarithmic Hessian on the local relative critical branch is

    [[144,-72,0],[-72,144,0],[0,0,72]], determinant 1119744.

It proves dominance as in section 3, and hence generic nonemptiness.

**Theorem 3.** For this cubic A2+A1^2 Weyl class, generic central targets
with centralizer K have nonempty smooth regular-source adjoint Springer
fibers of pure dimension 10. Weyl transfer covers all twelve embeddings.

## 5. Audit scope and next step

The [Sage verifier](tools/verify_binary_branching_fibers.py) reconstructs
the characters, checks exact eliminations, content/separability and
boundary identities, and verifies the dominance Hessians and Weyl words.
The [native audit](tools/audit_binary_branching_fibers.py) uses independent
Cartan reflections and explicit A1/A2 weights, plus SymPy polynomial
arithmetic. Their reports are [verification](binary-branching-verification.json)
and [audit](binary-branching-audit.json). Finite checks corroborate the
written theorem; neither specialized Hessians nor root counts replace it.

Known inputs are the [Kostant--Michor trace map](https://www.mat.univie.ac.at/~michor/kostant.pdf),
[Steinberg regular quotients](https://www.numdam.org/item/PMIHES_1965__25__49_0.pdf),
and standard reduced-hypersurface/generic-smoothness arguments as recorded
in the [review](CAYLEY_FIBER_REVIEW.md). The conditional binary criterion
and its applications are deductions written here; no novelty claim is made.

The 22 newly covered standard embeddings raise E7 coverage from 101 to
123 of 127. Two cubic Weyl families remain: A2^2+A1 (three embeddings)
and A2+A1^3 (one). The latter's separate three-variable certificate is
still pending. The all-simple-adjoint conjecture is not proved.

The next useful test is whether the A2^2+A1 residual admits this same
binary reduction with a controlled boundary. A proof would extend one
mechanism to a third family; failure should identify a structural
obstruction, rather than trigger an unrestricted new computation.
