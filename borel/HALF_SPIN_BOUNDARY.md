# Half-spin Springer maps: the complete underlying zero fiber

Date: 22 October 2026. Status: uniform written deductions, with independent
exact matrix and character corroboration. This note classifies the underlying
zero set and its dimension. It does not establish the boundary scheme structure,
all irreducible components, generic degree, or publication priority.

## 1. Goal, conventions and result

The goal is to explain exceptional fibers of representation-dependent
generalized Cayley maps by a mechanism that works across ranks. The
[preceding half-spin note](HALF_SPIN_CAYLEY.md) classified the zero fiber
inside a Cayley chart but left the boundary unresolved. In particular its
33-dimensional Spin_10 family was only a lower bound. The present argument
settles that boundary and proves the bound is sharp in every rank.

Work over C with G=Spin_(2r), r>=3, the standard cover
pi:G -> SO_(2r), and either half-spin representation S_sigma, sigma=+/-1.
Its weights are epsilon/2 with product epsilon_i=sigma. Use the actual
trace normalization

\[
 B_\sigma(\theta_\sigma(g),H)
 =\operatorname{tr}(S_\sigma(g)dS_\sigma(H)).
\]

Let c=chi_(S_+ direct-sum S_-)/2^r. On c!=0 use X=(pi(g)-I)(pi(g)+I)^(-1)
and the normalized Pfaffian p(X)=product x_i on the split Cartan. The
preceding note proves theta_sigma=2c(X+sigma J(X)), J=pX^(-1) extended
polynomially. Define the closed rank locus in SO_(2r)

\[
 \mathcal Z_r=\{a:\dim\ker(a-I)\ge4,\quad
                        \dim\ker(a+I)\ge4\}.                 \tag{1}
\]

These are ordinary eigenspaces, not generalized eigenspaces: equivalently
there are at least four Jordan blocks at each of the eigenvalues +1 and -1.

**Theorem.** As underlying closed sets,

\[
 \boxed{\theta_\sigma^{-1}(0)=Z(G)\ \sqcup\
       \bigsqcup_{\substack{u^{r-2}=1\\u\ne1}}
                  (\mathcal O_{\sigma,u,+}\sqcup\mathcal O_{\sigma,u,-})
       \ \sqcup\ \pi^{-1}(\mathcal Z_r).}                   \tag{2}
\]

Here Z(G) consists of four points. The chart orbits are precisely those
in the preceding note:

\[
 X^2=uI,\qquad p(X)=-\sigma u,\qquad c^2=(1-u)^{-r},
\]

with the two signs distinguishing the two nonzero values of c. Each orbit
has dimension r(r-1). The last term lies on c=0 and is the same for both
half-spin representations. In particular

\[
 \boxed{\dim\theta_\sigma^{-1}(0)=r(2r-1)-12\quad(r\ge4).} \tag{3}
\]

For r=3 the zero fiber consists of the four central points and has dimension
zero. The chart fiber is known to be smooth, but no reducedness assertion
is made here for the boundary fiber scheme. The rank equations in (1)
describe its underlying set; they are not asserted to generate the
scheme-theoretic fiber ideal.

## 2. Which semisimple boundary elements can be zeros?

Write torus parameters z_i=a_i^2, and set
A=product(a_i+a_i^(-1)), B=product(a_i-a_i^(-1)). Then
chi_sigma=(A+sigma B)/2. Let k be the number of z_i=-1 and l the number
of z_i=1. The boundary c=0 means k>=1.

If k=1, the derivative of A in that coordinate is nonzero, while the
derivative of B there is zero. Thus the Springer image cannot vanish.
If k>=2, every first derivative of A is zero. For l=1, the derivative
of B in that unique +1 coordinate is nonzero. For l=0, the derivative
of B in any coordinate with z_i!=-1 is nonzero, so vanishing forces
k=r. Finally, if k,l>=2 all first derivatives vanish.

Consequently the semisimple boundary zeros are exactly:

- the two central lifts of -I;
- the elements with -1 and +1 eigenspaces each of dimension at least four.

This also rules out extra semisimple cancellation cases at the boundary.

## 3. Reduction to the unipotent factors

Write the multiplicative Jordan decomposition g=su, with u unipotent
and commuting with s. A cocharacter in the reductive centralizer of s
contracts u to 1. By equivariance and continuity, theta_sigma(g)=0
therefore implies theta_sigma(s)=0. Also c(g)=c(s), since characters
take the same value on an element and its semisimple part.

When s is central, S_sigma(s) is a nonzero scalar, so
theta_sigma(su)=scalar*theta_sigma(u). A unipotent u has c(u)=1 and
lies in the preceding Cayley chart. Its zero-fiber classification
forces u=1. Thus central semisimple parts yield only central points.

It remains to treat s whose orthogonal image has -1 and +1 eigenspaces
V_- and V_+ of dimensions 2k and 2l, with k,l>=2. The orthogonal
complement V_R contains no +/-1 eigenvalues. Its centralizer factors
are GL_(m_j) for the remaining reciprocal eigenvalue pairs. Write
u_-,u_+,u_R for the corresponding restrictions of pi(u).

For an even-dimensional orthogonal block W, let A_W be its full-spin
character and B_W its difference of half-spin characters. Restriction
to the orthogonal block subgroup gives

\[
 \chi_\sigma=\tfrac12(A_-A_+A_R+\sigma B_-B_+B_R).        \tag{4}
\]

One may work on the product of the block Spin groups; its finite central
kernel does not affect Lie derivatives. Orientations fix the sign in (4).
When V_R=0 set A_R=B_R=1.

At a unipotent element in Spin_(2h), A_W=2^h and B_W=0. A lift of -I
acts as a nonzero scalar times the chirality operator, interchanging
A_W and B_W. Therefore at su we have

\[
 A_-=0,\quad B_-\ne0,\qquad A_+\ne0,\quad B_+=0.
\]

Both A_R and B_R are nonzero: their squared values are nonzero constant
multiples of det(I+pi(s)|_(V_R)) and det(I-pi(s)|_(V_R)), respectively.
These character values are unaffected by u_R.

Differentiate (4) in the centralizer Lie algebra. On the minus block
only the derivative of A_- survives, and it is a nonzero multiple of
dB_(V_-)|_(u_-). On the plus block only dB_(V_+)|_(u_+) survives.
Every derivative in a remaining GL block is zero because its coefficient
contains A_-=0 or B_+=0. There are no further conditions from directions
outside the centralizer: tr(S_sigma(su)dS_sigma(H)) is invariant under
Ad(s), hence vanishes on its eigenspaces with eigenvalue different from 1.
Thus

\[
 \theta_\sigma(su)=0\quad\Longleftrightarrow\quad
 dB_{V_-}|_{u_-}=0\ \hbox{and}\ dB_{V_+}|_{u_+}=0.        \tag{5}
\]

This reduction uses the actual chiral characters and is independent of
sigma. In particular the two block contributions cannot cancel each
other: they lie in different direct summands of the centralizer.

## 4. The four-eigenvector condition

For a unipotent v in Spin_(2h), h>=2, put
Y=(pi(v)-I)(pi(v)+I)^(-1). Character differentiation on the torus and
equivariant density on this Cayley chart show that the gradient of B_W,
under the defining vector trace pairing, is a nonzero constant
multiple of J_h(Y), the Pfaffian cofactor. The full-spin scalar equals
one at v. This assertion remains valid for h=2: it uses the difference
of the two characters, not a nondegenerate individual half-spin form
on so_4.

Because Y is nilpotent, its Pfaffian is zero. The cofactor vanishes
exactly when rank Y<=2h-4: a skew matrix of rank 2h-2 has a nonzero
rank-two Pfaffian cofactor, and all cofactors vanish at lower rank.
Since pi(v)+I is invertible,

\[
 dB_W|_v=0\quad\Longleftrightarrow\quad
 \operatorname{rank}(\pi(v)-I)\le2h-4
 \quad\Longleftrightarrow\quad\dim\ker(\pi(v)-I)\ge4.   \tag{6}
\]

On V_- the kernel in (6) is the -1 eigenspace of pi(su), and on V_+
it is its +1 eigenspace. Equations (5)--(6) prove exactly (1). Together
with the previous chart theorem and the central cases, this proves (2).

## 5. Dimension: why the earlier family is maximal

Stratify the boundary by the Jordan partitions lambda_- of 2k and
lambda_+ of 2l, and the partitions nu_j of the GL_(m_j) blocks.
Both lambda partitions must have at least four parts. Orthogonal
partitions have each even part occurring with even multiplicity.
Very-even splitting can produce distinct orbits of the same dimension;
it does not change this dimension argument.

For an orthogonal partition lambda of 2h, write lambda' for its
transpose and o(lambda) for the number of odd parts, with multiplicity.
The standard centralizer dimension is

\[
 z_{\rm O}(\lambda)
 =\tfrac12\left(\sum_j(\lambda'_j)^2-o(\lambda)\right).  \tag{7}
\]

For clarity, this formula also follows by writing the natural module
as the sum of sl_2 modules of dimensions lambda_i and using
so(V)=wedge^2 V. The number of irreducible summands, hence the dimension
of the nilpotent centralizer, is
sum_i floor(lambda_i/2)+sum_(i<j) min(lambda_i,lambda_j), which is (7).
Thus the dimension formula is a standard input, not a new orbit theorem.

If b=length(lambda)>=4, then lambda'_1=b and o(lambda)<=b. Since the
remaining column lengths sum to 2h-b and squares dominate integers,

\[
 z_{\rm O}(\lambda)\ge h+\frac{b(b-2)}2\ge h+4.          \tag{8}
\]

Equality is attained by (2h-3,1,1,1). For a GL_m block the centralizer
dimension is sum_j(nu'_j)^2>=m. Write m=sum m_j=r-k-l and let t be the
number of distinct remaining reciprocal eigenvalue pairs, so t<=m.
The stratum dimension is the conjugacy-orbit dimension plus its t
eigenvalue parameters. With d=dim G=r(2r-1), it is at most

\[
 d-(k+4)-(l+4)-m+t\le d-k-l-8\le d-12.                 \tag{9}
\]

Equality requires k=l=2, m=t=r-4 and trivial unipotents in the two
four-dimensional blocks. These are precisely the previously constructed
semisimple families with four +1 and four -1 eigenvectors and all
remaining reciprocal pairs distinct. They attain d-12. Chart zero
orbits have dimension r(r-1)<d-12 for r>=4; the center is finite.
This proves (3) and identifies the generic Jordan data of every stratum
attaining that dimension. It does not prove that every smaller stratum
lies in the closure of this family, or that its Spin inverse image is
irreducible.

## 6. Spin_10 explicitly

For either 16-dimensional half-spin of Spin_10 the whole underlying
zero fiber consists of the four central points, the four chart classes
of dimension 20, and the inverse image of Z_5. The possible boundary
Jordan data and their dimensions are:

| -1 block partition | +1 block partition | Free reciprocal eigenvalue pairs | Orbit dimension | Stratum dimension |
| --- | --- | --- | --- | --- |
| (1,1,1,1) | (1,1,1,1) | 1 | 32 | 33 |
| (1,1,1,1) | (1,1,1,1,1,1) | 0 | 24 | 24 |
| (1,1,1,1) | (2,2,1,1) | 0 | 30 | 30 |
| (1,1,1,1) | (3,1,1,1) | 0 | 32 | 32 |
| (1,1,1,1,1,1) | (1,1,1,1) | 0 | 24 | 24 |
| (2,2,1,1) | (1,1,1,1) | 0 | 30 | 30 |
| (3,1,1,1) | (1,1,1,1) | 0 | 32 | 32 |

These are Jordan strata in SO_10, whose finite Spin inverse images have
the same dimensions. The table does not assert a count of Spin orbits
or connected components over each row. In particular the global zero
fiber has dimension exactly 33, and nonsemisimple points do occur.

## 7. Consequence for localization, and the next goal

There is also a direct necessary condition for harmonic separation by
pullback along this map. Suppose a nonzero invariant f is nonzero on
one of the positive-dimensional zero conjugacy classes. Every function
in the image of

\[
 \mathbb C[G]^G[f^{-1}]\otimes\theta_\sigma^*\operatorname{Harm}(\mathfrak g)
 \longrightarrow\mathbb C[G][f^{-1}]
\]

is constant on that class, whereas regular functions on the affine
open G_f separate its distinct points. Multiplication is therefore not
surjective. For r>=4 the chart zero classes have both c!=0 and
chi_sigma!=0. Consequently localizing at c alone, or at the individual
half-spin character alone, cannot give this harmonic separation.
This is an elementary consequence of the collapsed classes, not a
classification of all sufficient localizations.

The bounded boundary investigation has reached its intended goal:
an all-rank criterion, the full underlying zero set and a sharp dimension
formula. The next priority should be the generic fiber for Spin_10,
starting from the coupled quadratic equations in the preceding note.
That addresses the number of inverse branches directly and can suggest
a degree law across ranks. Scheme multiplicities and incidence between
zero-fiber components remain legitimate questions, but they are not
prerequisites for that degree calculation and are not the next priority.

## 8. Verification and source distinctions

The [verifier](tools/verify_half_spin_boundary.py) uses actual chiral
Clifford matrices at nonsemisimple points, including negative controls.
It constructs orthogonal Jordan representatives over Q(i), exponentiates
their spin Lie matrices and verifies the covering action on every Clifford
generator before testing the trace projection. It checks all admissible
partition pairs for the Spin_10 boundary semisimple types, with both
chiralities, plus Spin_8 controls and higher-rank examples with GL
unipotents, conjugated points and a very-even partition.

Separate checks compare centralizer dimensions from partitions, sl_2
weights on wedge^2 V and direct commutator kernels. They audit the
semisimple classification from original weight sums and the dimension
bounds over finite ranges. The [report](half-spin-boundary-verification.json)
records these ranges, results and source hashes. The uniform proofs do
not rely on extrapolating these finite calculations.

Published inputs include the trace construction and Clifford/Cayley
chart in [Kostant--Michor, Sections 7.7 and 7.14--7.15](https://arxiv.org/pdf/math/0109066),
Jordan decomposition and sl_2 theory, the half-spin tensor restriction,
and the standard orthogonal nilpotent classification. For a source
recalling the partition classification and very-even qualification, see
[Goodwin--Pengelly, Section 2.2](https://link.springer.com/article/10.1007/s00031-022-09722-y).
Formula (7) is rederived above to make the dimension argument inspectable.

The boundary rank criterion, its application to the full half-spin
zero fiber, and the sharp dimension statement are written deductions
here. The limited source search does not establish that they are new
to the literature. There is no claim of external review, a boundary
scheme theorem, generic degree, or a new runtime API in this step.
