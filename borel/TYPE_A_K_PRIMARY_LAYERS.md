# Prime-primary lower bounds from a uniform coefficient family

27 October 2026. IMG-38. Continuation of the
[faithful-specialization theorem](TYPE_A_K_SPECIALIZATION.md), with the
same full literal defining SL_m/P1 source. The earlier result and its
dependencies were pushed in commit 6b1279e. This note does not assert a
closed formula for all factors of the full defect Delta_m.

## 1. The advance and its exact scope

**Theorem.** For every prime power m=p^s>=3 there is a surjection

    Delta_m -> direct-sum_(e=1)^s (Z/p^e)^((p-1)p^(e-1)).

For an odd prime degree this is Delta_p -> (Z/p)^(p-1). In general the
full p-primary defect has at least m-1 cyclic factors, exponent at least
m, and order at least p^(s*m-(m-1)/(p-1)). These are lower bounds, not a
direct-summand assertion or a formula for the full defect.

The proof constructs one uniform coefficient family with exactly this
kernel, proves its entire image-module presentation, and proves that its
defect is a quotient of the original defect. It is not a rank-by-rank
enumeration. For m=5 the full kernel has twenty factors 5 while this
quotient has four; the lower bound need not be sharp. For m=4 it matches
the previously certified full kernel Z/2 plus (Z/4)^2.

**Still conjectural:** m Delta_m=0 for every m>=2. Complete full-image
calculations through m=6 satisfy this. The preceding theorem does not
prove an upper bound of that strength. Section 5 explains exactly why the
previous determinant mechanism alone cannot do so.

## 2. Base change gives quotient defects, not upper bounds

Let S,A,J and their integer augmentation be as in the earlier note.
For any augmented coefficient change S -> B -> Z_(p), let J_B be the
IMAGE of J tensor_S B in A tensor_S B. Do not identify this image with
J tensor_S B unless injectivity has been proved.

There is a surjection J tensor_S B -> J_B, hence after augmentation a
surjection from J tensor_S Z_(p) to J_B tensor_B Z_(p). Both map onto the
same ordinary image J_0 tensor Z_(p). Therefore their specialization
kernels fit into a surjection

    (Delta_m)_(p) -> Delta_B.

Indeed lift a class of Delta_B through the first surjection. Its ordinary
image is zero, so the lift already belongs to the original kernel. No
flatness of B over S and no commuting-of-intersections assertion is used.
This allows smaller coefficient families to prove lower bounds and
disprove proposed annihilators, but not to prove an annihilator of Delta_m.

## 3. A uniform prime-power trace family

Fix m=p^s>=3 and work over B=Z_(p)[a]. Specialize the symmetric
coefficients by

    c_j=binomial(m,j)+binomial(m-1,j-1)a  (1<=j<m), c_m=1.

Augmentation is a=0. Set v=L-1 and epsilon=(-1)^m. The ambient
projective polynomial becomes

    F(v)=v^m-a v^(m-1)-epsilon*a.

Comparing coefficients in F(v+u), the full source generators become

    f_j=binomial(m,j)v^j-binomial(m-1,j-1)a v^(j-1)
       =binomial(m-1,j-1) h_j,
    h_j=(m/j)v^j-a v^(j-1),           1<=j<m.

Every binomial(m-1,j-1) is a p-local unit: the base-p digits of m-1
are all p-1. Also m/j lies in Z_(p), with valuation s-v_p(j)>=1.
Thus the full image is B[h_1,...,h_(m-1)] in A_B=B[v]/F(v). Put r=m-1.

**Complete image-module theorem.** This image is precisely the B-span of

    1, h_1,...,h_r, a^2 v,...,a^2 v^r.                 (1)

Here is a proof of both inclusions and multiplication closure. Define
h_k by the same rational polynomial formula for every k needed below.
For 1<=i,j<m and i+j<=2m-2, the identity is

    h_i h_j = [m(i+j)/(ij)] h_(i+j)+a^2 v^(i+j-2).    (2)

All coefficients used here are p-integral. Indeed
v_p(m(i+j)/(ij))>=s-max(v_p(i),v_p(j))>=1, and
m/(i+j) is p-integral in this range, including i+j=m.
In A_B, h_m=epsilon*a. For 1<=k<=m-2, use (2) with i=1,j=k+1
to express a^2 v^k using input generators and the scalar h_m.
The remaining power follows by reducing h_(m+1):

    h_2 h_(m-1) - [epsilon*m/(2(m-1))] a h_1
       = [(m-2)/(2(m-1))] a^2 v^(m-1).

The coefficient on the right is a p-local unit: for odd p, m-2 and 2
are units; for p=2, m>=4 and both m-2 and 2 have valuation one.
The other coefficient is p-integral. Consequently a^2 A_B is contained
in the image. This is the step that fails for m=2.

Modulo a^2 A_B, for 1<=k<=m-2 the defining relation gives

    h_(m+k) = epsilon*[m/(m+k)] a v^k
            = epsilon*[k/(m+k)] a h_k.

The ratio k/(m+k) is a p-local unit, since v_p(k)<s. Thus every h_(i+j)
appearing in (2), after reduction, belongs to the proposed span.
Equation (2) proves closure under products of inputs. Products involving
a^2 A_B stay in that ideal. This proves the FULL module presentation (1).

## 4. Exact defect of the family

Let C_B=A_B/J_B. Remove the constant coordinate and denote the classes
of v,...,v^r by e_1,...,e_r. Put q_j=m/j and w_j=v_p(q_j)=s-v_p(j).
The presentation (1) gives ALL relations

    q_1 e_1=0,
    q_j e_j=a e_(j-1)  (2<=j<=r),
    a^2 e_j=0         (1<=j<=r).                       (3)

Because A_B is free and B/(a)=Z_(p), the specialization defect is

    Delta_B=Tor_1^B(C_B,B/(a))=ker(a:C_B->C_B).

Put z_j=a e_j. As a Z_(p)-module, relations (3) say

    q_j z_j=0, q_1 e_1=0, q_j e_j=z_(j-1) (2<=j<=r).

Eliminating z_1,...,z_(r-1) therefore gives

    C_B = Z_(p)/(q_1) e_1 direct-sum
          direct-sum_(j=2)^r Z_(p)/(q_(j-1)q_j) e_j
          direct-sum Z_(p)/(q_r) z_r.

These relations are complete: polynomial coefficients reduce to their
constant and a terms using a^2 C_B=0. Multiplication by a sends e_j to
q_(j+1)e_(j+1) for j<r, sends e_r to z_r, and kills all z_j.
Its kernel is the independent direct sum generated by

    q_2 e_2,...,q_r e_r,z_r,

with orders p^w_1,...,p^w_r, respectively. Hence

    Delta_B = direct-sum_(j=1)^(m-1) Z/p^(s-v_p(j)).

There are (p-1)p^(e-1) indices j with s-v_p(j)=e. This proves the closed
formula in section 1; summing the exponents gives
s*m-(m-1)/(p-1). The quotient lemma transfers it to the full Delta_m.
For m=2 the image-module claim was not established and the original
specialization defect is zero; it is explicitly excluded.

## 5. Why the stronger annihilator does not follow formally

There is a uniform counterexample to an OVERGENERAL module argument.
It is not a counterexample to m Delta_m=0 for Springer images.

Fix a prime p and any r>=2. Over B=Z_(p)[a], let A=B^r and let J be
generated by the r columns of M=pI+aN, where N has ones immediately
below its diagonal, together with the additional column a e_1.

At a=0 the ordinary image is p Z_(p)^r, so every ordinary Smith factor
is p. At the generic characteristic-p point the image is all of A:
the shift columns and a e_1 span. Nevertheless

    ker(J/aJ -> p Z_(p)^r) = Z/p^r.

To prove this, the full relation module among the r+1 displayed columns
is free on the vector whose first r entries are

    -a p^(r-1), a^2 p^(r-2), ..., (-1)^r a^r,

and whose last entry is p^r. The recurrence from the bottom row, using
that a is prime to p, forces the last coefficient of any relation to
be divisible by p^r. This proves completeness. At a=0 the relation
has only its last entry p^r nonzero. The claim follows.

Thus generic characteristic-p image recovery plus ordinary exponent p
does not bound the specialization exponent by p. A proof of the proposed
Springer bound must use additional multiplicative/source relations.
The input J in this control is a module, not a Springer image algebra.

## 6. Checks, computational limits and next proof target

The native audit verifies the source shift, product identities, membership
and multiplication closure independently of the standard-basis engine.
It also constructs the finite abelian presentation of C_B and the map a
and checks its kernel size and factors for several prime powers. The
Sage calculation starts from the original specialized f_j, closes their
integral image module under multiplication, computes all syzygies, and
recovers the predicted kernels for m=3,4,5,7,8,9 in two monomial orders.

Two additional c_1-only coefficient slices in m=7,8 have kernels
(Z/7)^6 and Z/2 plus (Z/4)^5 plus (Z/8)^2 respectively. They are only
quotients of the full kernels and are not promoted to full-rank results.
The attempted full m=7 standard-basis run was stopped near 10 GB of
memory before completion; it yielded no certified full answer. This
motivates the structural coefficient-family proof rather than a larger
uncontrolled matrix computation.

The algebra tools here (Tor, base change, finite presentations and Smith
theory) are standard. The specific coefficient-family presentation and
its resulting lower bound are repository deductions, not externally
reviewed results or established novelty claims. No new literature
classification of these defects is asserted.

Next: prove or disprove m Delta_m=0 from the full multiplicative source
relations. The precise missing step is a mechanism controlling the full
augmentation kernel, not merely its quotients under coefficient changes.
The prime-power family calculation gives a rigorous lower bound and a compact
model to test such a mechanism, while the module counterexample rules
out relying only on the previous determinant criterion.
