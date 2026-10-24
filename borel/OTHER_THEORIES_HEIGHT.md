# A formal-group-height test for the symplectic Springer source

Exploration, 22 October 2026. This is a complete elementary deduction for
one arbitrary-rank family and the specified trace-coordinate lift. Formal
group height and the formulas for complex conjugation/formal inverse in
Morava K-theory are established literature, not new results claimed here.
The point here is their consequence for the **same literal Springer
character source** used in the finite program.

## 1. The source and the space

Let G=Sp_(2r,C), rho its defining representation, and P=P1. Then
G/P=CP^(2r-1). Write epsilon_i for the standard torus weights. With the
normalization in Part I,

\[
 \theta^*(\epsilon_i)=(t_i-t_i^{-1})/2,\qquad
 \Lambda_\rho=2\mathbb Z\epsilon_1+\cdots+2\mathbb Z\epsilon_r.
\]

The full integral Levi source is

\[
 A=\mathbb Z[z,q_1,\ldots,q_{r-1}],\qquad
 z=2\epsilon_1,\quad
 q_j=e_j((2\epsilon_2)^2,\ldots,(2\epsilon_r)^2).
\]

Its ordinary image is Z[2h] inside Z[h]/h^(2r), with relations
q_j=(-z^2)^j and z^(2r)=0. These follow by setting to zero the positive
global signed-symmetric invariants in the squares of 2epsilon_i.
The literal torus/Levi sources and Cayley coordinates are published in
[Kumar, Part I, Sections 6--7](https://kumar.math.unc.edu/papers/kumar76.pdf).
This source is not the unscaled normalized generator ring of Part II.

For a complex orientation with formal group law F and inverse iota_F,
put x=c1^E(O(1)). The projective bundle theorem gives
E^*(CP^(2r-1))=R[x]/x^(2r). The
[trace-coordinate construction](OTHER_THEORIES_FINITE.md#4-a-uniform-oriented-lift-with-explicit-hypotheses)
sends

\[
                            z\longmapsto b=x-\iota_F(x).            \tag{1}
\]

In complex cobordism the global invariant relations above vanish by the
integral torsion-free argument in that note. Naturality of the complex
orientation then preserves those identities on specialization to the
characteristic-two theories below. Consequently the full Levi image is
R[b]; the higher q_j add no extra independent classes.

## 2. Height and the complete image

Assume R is a characteristic-two coefficient domain and

\[
                 [2]_F(x)=u x^{2^s}+\text{higher terms},\quad u\in R^\times.
                                                                    \tag{2}
\]

This is the usual finite-height situation, with s>=1; for Morava K(s)
one can use R=F2[v_s,v_s^(-1)]. The even projective-space groups are
commutative in the usual complex-oriented calculation. No identification
of arbitrary mod-two invariant rings is used.

**Theorem.** Put m=2r and

\[
                         L=\left\lceil\frac{m}{2^s}\right\rceil.
\]

For the base change of the literal source above, its entire oriented
image and complete source kernel are

\[
 \boxed{\operatorname{im}\xi_{\rho,E}^{P1}\cong R[z]/(z^L),\qquad
 \ker\xi_{\rho,E}^{P1}=(q_j-(-z^2)^j\ (1\le j<r),\ z^L).}          \tag{3}
\]

In particular its R-rank is L, although the ambient projective-space
rank is m. The linear Springer generator is nonzero exactly when m>2^s.
Its first nonzero support-filtration term is u*x^(2^s).

**Proof.** Since F(x,iota_F(x))=0, formal division in the second variable
gives

\[
 F(x,y)=(y-\iota_F(x))U(x,y),\qquad U(0,0)=1.
\]

Setting y=x shows [2]_F(x)=(x-iota_F(x))*U(x,x). The factor U is a
unit, so (1) has order 2^s with unit leading coefficient. In R[x]/x^m,
the powers 1,b,...,b^(L-1) have distinct leading orders and are R-linearly
independent; b^L=0. The already-vanishing q_j relations reduce every
source element to a polynomial in z. This proves both assertions in
(3), including completeness of the kernel. QED.

At an odd prime, the first term of b is 2x with unit coefficient, so the
image is the whole projective-space theory by a finite inverse-coordinate
change. For the additive theory in characteristic two the inverse is x,
so b=0 and only constants remain. Thus (3) expresses a genuine dependence
on the formal group law, not just a change of Schubert-basis notation.

For example on Sp4/P1=CP3, height one gives rank two (1,b), while height
two gives only constants. On Sp6/P1=CP5, height two gives rank two. These
are examples of the proved formula, not its evidence in place of proof.

## 3. What is known and what is being deduced

The height/formal-inverse mechanism is classical. In particular the
formal inverse in Morava K-theory and its interpretation as complex
conjugation are explicitly studied in
[Kriz--Laures, Real structures and Morava K-theories, Section 2](https://www.ruhr-uni-bochum.de/topologie/KL.pdf)
and [Kitchloo--Wilson, The Morava K-theory of spaces related to BO, Section 3](https://math.jhu.edu/~nitu/papers/submit3.pdf).
The source identification is Part I's Springer construction. Formula
(3) applies those established facts to the full specified Springer image;
its literature priority has not been established. It is not a new
computation of the ambient Morava K-ring of projective space.

The [Sage series check](experiments/other_theories_height_pilot.py) verifies
320 height/rank cases using 2-typical logarithms at heights one through
five, to precision 130. It separately computes the formal inverse of the
supersingular elliptic curve y^2+y=x^3 over F2, and verifies the height-two
image ranks in another 64 cases. The comparison uses the common leading
order, rather than asserting that the formal group coordinates are
identical. [Report](experiments/other-theories-height-pilot.json).

## 4. A concrete next generalization

Investigate Levi sources whose generators have higher signed/symmetric
degrees and whose images are not generated by one class. The exact
questions are the image rank over the coefficient field, the new modular
source kernel and the first nonzero support filtration. The universal
torsion-free lift supplies a starting presentation, while reduction can
introduce additional relations. It would be unjustified to infer a
general symplectic or exceptional height formula from (3).

## 5. An arbitrary-rank integral K-theory consequence

For this same source on CP^(2r-1), ordinary integral K-theory gives

\[
 K=\mathbb Z[v]/v^{2r},\quad L=1+v,\qquad
 J_K=\mathbb Z[b],\quad b=L-L^{-1}.
\]

Its abstract source image is still Z[z]/z^(2r), but its ambient quotient
has different factors from the cohomological list 1,2,4,...,2^(2r-1).
The following is a complete recursive presentation of that quotient,
plus a closed formula for its exponent. It does not give a closed list
of all Smith multiplicities. That remaining even-dimensional-lattice
case is now completed by the [closed ordinary Hankel formula](SYMPLECTIC_K_ORDINARY_SMITH.md).
The older argument below also treats odd lattice lengths m, which are
outside that symplectic projective theorem.

**Proposition.** For m>=2 put c_0=0, c_1=1 and

\[
 c_{j+1}(b)=(b-2)c_j(b)+b c_{j-1}(b).
\]

Then for K_m=Z[v]/v^m and b=(1+v)-(1+v)^(-1),

\[
 K_m/\mathbb Z[b]\cong
           \mathbb Z[b]/(c_m(b),\ b c_{m-1}(b))                    \tag{4}
\]

as modules over Z[b]. Set j=floor((m-1)/2). The largest nonunit Smith
factor is 2^e, where e=1 if m=2 and otherwise

\[
                         \boxed{e=4j-s_2(j).}                     \tag{5}
\]

Here s_2(j) counts binary ones. That largest factor occurs twice if
m is even and m>=4, and once if m is odd or m=2. The number of unit
factors is ceil(m/2), and the total index is 2^(m(m-1)/2).

For our symplectic flag m=2r, (5) is 4(r-1)-s_2(r-1) for r>=2.

### Proof

The coordinate relation is v^2=(b-2)v+b. Write
v^k=a_k(b)+c_k(b)v with integral polynomials; then a_k=b*c_(k-1)
and the displayed recurrence follows. Regard the quotient by v^m
as a module with generators 1,v over Z[b]. Its relations come from
v^m and v^(m+1). Quotienting further by the source copy Z[b] kills
the generator 1, leaving relations c_m and a_m+(b-2)c_m. This is
exactly (4).

Over Q the inverse coordinate is

\[
 v(b)=\frac b2+\sqrt{1+b^2/4}-1
     =\frac b2+\sum_{k\ge1}\binom{1/2}{k}\frac{b^{2k}}{4^k}.
\]

The k-th even coefficient has 2-denominator exponent
4k-s_2(k), using the Catalan identity for binomial(1/2,k) and
Legendre's factorial formula. These exponents strictly increase with k.
The least power of two making the inverse inclusion matrix integral
is consequently (5): every v^k=a_k+c_k*v has no worse denominator
than v itself, and the last even coefficient of v attains the bound.
This least denominator is exactly the exponent of the finite cokernel.

To count the largest factors, multiply the inverse inclusion matrix by
2^e and reduce modulo two. Only the last even coefficient of v survives.
For odd m it is in degree m-1, so the matrix has rank one. For even
m>=4 it is in degree m-2; multiplication by the integral c_k can only
also reach degree m-1. The columns v and v^2 have c_1=1 and c_2=b-2,
so those last two rows have rank two. This rank counts exactly the
Smith factors equal to 2^e. For m=2 the assertion is immediate.

Modulo two, b=v^2/(1+v), so its powers have distinct first degrees
0,2,4,...<m. Their rank is ceil(m/2), proving the unit count. Finally
the original inclusion matrix is triangular with diagonal 2^k,
0<=k<m, giving the stated total index. QED.

For example r=2 gives Smith factors (1,1,8,8), and r=3 gives
(1,1,1,2,128,128). These are complete calculations for those ranks;
the proposition establishes the exponent and largest multiplicity for
every r. The [independent native pilot](experiments/other_theories_symplectic_k_pilot.py)
checks complete Smith forms through m=64 using truncated v-powers and
separate projective Euler characteristics, checks (4) through m=32,
and checks 1,000 exact inverse-coordinate denominators.
[Report](experiments/other-theories-symplectic-k-pilot.json).

The general quotient presentation and exponent are deductions here;
their literature priority has not been established. The
[subsequent ordinary Hankel theorem](SYMPLECTIC_K_ORDINARY_SMITH.md)
now completes the full multiplicity list for even m=2r, the symplectic
projective family. The analogous full list for odd truncation length m
is not asserted in that theorem.
