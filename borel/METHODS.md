# The mathematical meaning of the methods

This guide connects the public API to the maps and algebras it computes.
Read [MATHEMATICS.md](MATHEMATICS.md) for the Borel/Springer framework and
Schubert conventions, and [PRESENTATIONS.md](PRESENTATIONS.md) for the full
coinvariant and completeness proofs. The executable companion is
[examples/mathematical_methods.py](examples/mathematical_methods.py).

## The map being studied

Fix a flag variety `G/P`. Let `W_P` be the Weyl group of its Levi subgroup, and
put

$$
S=\mathbb Q[\omega_1,\ldots,\omega_r],\qquad
\beta^P:S^{W_P}\longrightarrow H^*(G/P;\mathbb Q).
$$

Choose homogeneous Levi-invariant polynomials `f_1,...,f_m` of positive
polynomial degrees `d_1,...,d_m`. `GeneratorAlgebra` constructs the map

$$
\phi:R=\mathbb Q[y_1,\ldots,y_m]\longrightarrow H^*(G/P;\mathbb Q),
\qquad y_j\longmapsto\beta^P(f_j),\qquad \deg y_j=d_j.
$$

Its answers concern the selected image

$$
\mathcal A=\mathbb Q[\beta^P(f_1),\ldots,\beta^P(f_m)]
\cong R/\ker\phi.
$$

The full Borel map is surjective over `QQ`; a selected image can be a proper
subalgebra. The chosen inputs may be algebraically dependent, and a nonzero
input polynomial may have zero cohomology image.

| Object in code | Mathematical object |
| --- | --- |
| `b.ring`, `b.weights` | `S` and its fundamental-weight variables |
| `A.generator_polynomials` | The actual Cartan polynomials `f_j` |
| `A.domain` | The separate formal polynomial ring `R` |
| `A.generator_degrees` | The positive weights `d_j` of the formal variables |
| `A.evaluate(p)` | The Schubert expansion of `phi(p)` |
| `A.presentation().kernel_ideal()` | The entire ideal `ker(phi)` in `R` |

Polynomial degree `d` corresponds to cohomological degree `2*d`. Formal
variables use their input degrees: a squared tail elementary generator of
index `i` has polynomial degree `2*i` and cohomological degree `4*i`.

`A.domain` uses ordinary degrevlex for its **term order**. The mathematical
weights are recorded separately in `A.generator_degrees`; calling Sage's
ordinary `total_degree()` on a formal expression does not compute its weighted
degree. For example, with weights `(1,2)`, `c**4+c**2*e+e**2` has weighted
degree four, despite its terms having different ordinary total degrees.

## Borel and Schubert methods

`BorelMap(type, omitted_nodes=...)` specifies the root datum and the standard
parabolic. Retained nodes generate `W_P`; omitted nodes give the available
fundamental-weight divisors. `None` selects the full flag, an integer a maximal
parabolic, an iterable a general partial flag, and `()` the point `G/G`.

The normalization is `beta(omega_i)=sigma_(s_i)` on the full flag. For `G/P`
this divisor statement applies to omitted nodes. Schubert classes are indexed
by minimal **right** coset representatives in `W/W_P`.

| Method | Mathematics and result |
| --- | --- |
| `b.expand(f)` | Computes `beta^P(f)=sum_w a_w sigma_w`; checks Levi invariance |
| `b.coefficient(f, word)` | Computes the single scalar `a_w` by divided differences |
| `b.reflect(f,i)` | Applies the Weyl reflection `s_i` to a Cartan polynomial |
| `b.divided_difference(f,i)` | Computes `(f-s_i(f))/alpha_i` by exact division |
| `b.chevalley(i,word)` | Computes the product of `sigma_w` with the divisor `beta(omega_i)` |
| `b.divisor_product(indices)` | Repeated Chevalley products, equal to `beta(prod omega_i)` |
| `b.multiply_divisor(expansion,i)` | Extends divisor multiplication by linearity |
| `b.basis(d)` | Lists minimal right coset representatives of length `d` |
| `b.poincare_dual(word)` | Returns the index `w_0*w*w_(0,P)` of the dual Schubert class |
| `b.longest_representative()` | Returns the top index `w_0*w_(0,P)` |

For a reduced word `w=s_i1*...*s_id`, the coefficient is

$$
a_w=(\partial_{i_1}\cdots\partial_{i_d}f)(0),\qquad
\partial_i f=\frac{f-s_i(f)}{\alpha_i}.
$$

The rightmost operator acts first. `poincare_dual` returns an index; it does
not itself return an intersection number. General multiplication of arbitrary
Schubert expansions is outside this API. Polynomial representatives can be
multiplied before `expand`, or formal representatives within a selected image
can be multiplied and then reduced as below.

The classical helpers construct specified symmetric polynomials: `tail_polynomial`
and `head_polynomial` return literal inputs in `S`; their `...elementary`
counterparts return Schubert expansions after the characteristic map.
`euler_polynomial` returns the type D coordinate product, and `euler` its image.
The [exceptional catalogue](EXCEPTIONAL.md) similarly separates named input
polynomials, normalization/source metadata and computed images.

## Degreewise images and relations

For a nonnegative integer `d`, the formal monomial basis of `R_d` consists of
`y^alpha` with `sum(alpha_j*d_j)=d`. The degreewise map is the exact linear map

$$
\phi_d:R_d\longrightarrow H^{2d}(G/P;\mathbb Q).
$$

`A.degree(d)` builds its coefficient matrix `C_d`. Columns are formal
monomials, and rows are Schubert classes appearing in their images. A relation
`sum_alpha v_alpha*y^alpha` is therefore represented by a vector in the
**right kernel** of `C_d`.

| Method or property | Mathematical interpretation |
| --- | --- |
| `A.monomials(d)` | The weighted monomial basis of `R_d` |
| `D=A.degree(d)` | The exact homogeneous map `phi_d` |
| `D.matrix`, `D.monomials`, `D.schubert_indices` | The matrix and its column/row labels |
| `D.rank` or `D.image_dimension` | `dim_Q(im(phi_d))` |
| `D.relations` | A vector-space basis of `(ker(phi))_d` |
| `D.image_basis_monomials` | Formal monomials with independent image columns |
| `A.substitute(p)` | The literal polynomial `p(f_1,...,f_m)` in `S`, before applying `beta` |
| `A.evaluate(p)` | `phi(p)` in Schubert coordinates |
| `A.through_degree(D)` | The homogeneous reports in degrees `0,...,D` |
| `A.hilbert_function(D)` | The prefix `(dim A_0,...,dim A_D)` |
| `A.first_relation_degree(D)` | First positive kernel degree through `D`, or `None` if none is found there |

For example:

```python
from sage_borel import BorelMap, GeneratorAlgebra

b = BorelMap("A2", omitted_nodes=1)  # G/P is the projective plane
A = GeneratorAlgebra(b, {"h": b.weights[0]})
h, = A.domain.gens()

D = A.degree(2)
assert D.rank == 1 and D.kernel_dimension == 0
assert A.degree(3).relations == (h**3,)
assert A.hilbert_function(4) == (1, 1, 1, 0, 0)
```

Only nonzero Schubert rows are stored; adding the omitted zero rows would not
change rank or kernel. `include_ambient=True` additionally computes the full
ambient degree dimension. A few homogeneous kernel spaces alone do not
establish generators of the entire ideal. Likewise, `None` from the bounded
first-relation search is a statement only about its requested range.

## Complete presentations and normal forms

`K=A.presentation()` computes a complete defining ideal rather than a degree
prefix. Write `J=K.kernel_ideal()` for the ideal it returns. The first
isomorphism theorem gives

$$
R/J\cong\mathcal A.
$$

| Method or property | Mathematical interpretation |
| --- | --- |
| `K.relations` | The reduced monic Groebner basis of `ker(phi)` in the formal domain's term order |
| `K.kernel_ideal()` | A fresh Sage ideal with those generators |
| `K.standard_monomials` | Formal monomials outside the initial ideal; their images form a basis of the selected image |
| `K.dimension` | The total vector-space dimension of the selected image |
| `K.hilbert_coefficients` | Its complete Hilbert function, ending at the last nonzero degree |
| `K.hilbert_series` | `sum_d dim(A_d)*t^d`, represented by a finite polynomial in `QQ[t]` |
| `K.reduce(p)` | The unique remainder of a **formal** polynomial modulo the returned Groebner basis |
| `K.backend` | The coordinate system used for exact image comparisons |
| `K.as_dict()` | Exact source data, grading, relations, basis and certificate in JSON-safe form |

The Groebner basis need not be a minimal set of ideal generators. Changing the
formal variable order can change that basis and the standard monomials while
leaving the abstract image algebra unchanged. Its weighted Hilbert series is
invariant under a change of quotient basis. Here “Hilbert series” is not the
eventual Hilbert polynomial of a graded module.

Continuing the projective-plane example:

```python
K = A.presentation()
assert K.relations == (h**3,)
assert K.standard_monomials == (1, h, h**2)
assert K.dimension == 3
assert K.kernel_ideal() == A.domain.ideal(h**3)
assert K.reduce(h**4 + 2*h**2) == 2*h**2
assert K.reduce((1+h)*(1+h+h**2)) == 1 + 2*h + 2*h**2
assert K.verify().certified
```

Thus equality and multiplication inside the selected image can be calculated
in formal coordinates: `phi(p)=phi(q)` exactly when `K.reduce(p-q)==0`, and
`K.reduce(p*q)` represents their product. An arbitrary Schubert class outside
the selected image has no formal representative in this quotient.

## The two backends compute the same kernel

The BGG backend compares Schubert-coordinate vectors. The coinvariant backend
compares normal forms in

$$
S/I_W,\qquad I_W=(S^W_+)S.
$$

The pullback `H*(G/P;QQ) -> H*(G/B;QQ)` is injective, and its image is the
`W_P`-invariant subalgebra. Consequently the two backends detect precisely the
same relations on the selected Levi-invariant inputs. Choosing a backend
changes coordinates used for the calculation, not `phi` or its kernel.

| Choice | Behavior |
| --- | --- |
| `A.presentation()` or `backend="auto"` | Coinvariant normal forms for standard-numbered A/B/C/D; BGG otherwise |
| `backend="coinvariant"` | Explicitly requests that classical quotient; unsupported types raise `NotImplementedError` |
| `backend="bgg"` | Explicitly requests Schubert-coordinate calculation |
| `K.verify()` | Rechecks with `K.backend` |
| `K.verify(backend="bgg")` | Performs a separate Schubert-coordinate check |

For small examples, the latter is a useful independent comparison:

```python
assert K.backend == "coinvariant"
assert K.verify(backend="bgg").certified
assert A.presentation(backend="bgg").relations == K.relations
```

There are three different bases involved:

1. The **ambient coordinate monomials** of `S/I_W` are internal comparison
   coordinates for the coinvariant backend.
2. `K.standard_monomials` are **formal** monomials of `R/ker(phi)`.
3. The **Schubert basis** is the basis used by `expand`, `evaluate` and `degree`.

`K.reduce` uses the second basis. It neither returns an internal Cartan normal
form nor converts an expression to Schubert coordinates. Use `A.evaluate` for
that conversion, when needed.

## Why verification proves the full kernel

Let `B` be the retained formal monomials and `G` the returned relations. The
certificate checks that every relation has zero image; the images of `B` are
independent; `1` is retained; and every variable multiple of an element of `B`
is either retained or belongs to the monomial ideal generated by the leading
terms of `G`. It also checks reduced-basis and output-grading consistency.

The border condition proves that `B` spans the proposed quotient. If a
polynomial in the true kernel had a nonzero remainder in that basis, the
remainder would be a nontrivial linear combination of independent images with
zero image. Hence every true relation belongs to the returned ideal. This is
why the certificate establishes completeness rather than just vanishing of
some candidate relations. [PRESENTATIONS.md](PRESENTATIONS.md#why-the-certificate-proves-completeness)
gives the detailed proof and the ordering and termination arguments.

`verify()` repeats exact checks using the selected evaluator and the recorded
inputs, with the geometric degree bound recomputed. It returns a
`KernelCertificate`; inspect `.certified` on that object. A single same-backend
check does not claim that a second algorithm or a proof assistant was run.

## The rank-nine calculation, derived

For type `C9` with omitted node `8`, the flag variety is the isotropic
Grassmannian `IG(8,18)`. The tail factor is `Sp(2)`. Set `q=x_9^2`, which has
polynomial degree two and cohomological degree four.

In the full symplectic coinvariant ring, a Groebner basis of the invariant
ideal is

$$
g_i=h_i(x_i^2,\ldots,x_9^2),\quad 1\le i\le9,
\qquad \operatorname{LM}(g_i)=x_i^{2i}.
$$

In particular, `g_9=x_9^18=q^9`, so `q^9=0`. Every `q^j=x_9^(2*j)` with
`0 <= j <= 8` is a standard ambient monomial, and therefore is nonzero. These
nine classes have distinct degrees and are independent. It follows that

$$
\ker(\mathbb Q[q]\longrightarrow H^*(IG(8,18);\mathbb Q))=(q^9),
\qquad \mathcal A\cong\mathbb Q[q]/(q^9).
$$

Its dimension is nine and its Hilbert series is
`1+t^2+t^4+...+t^16`. The ambient variety has complex dimension 52; the
relation already occurs in polynomial degree 18 because of the invariant
ideal. In general, the same argument gives `(q^n)` for the one-tail-generator
family in type `C_n`, with omitted node `n-1` and `n >= 2`.

```python
b9 = BorelMap("C9", omitted_nodes=8)
A9 = GeneratorAlgebra(b9, {"q": b9.tail_polynomial(1, 1)})
q, = A9.domain.gens()
K9 = A9.presentation()

assert K9.relations == (q**9,)
assert K9.standard_monomials == tuple(q**j for j in range(9))
assert K9.dimension == 9
assert K9.hilbert_coefficients == (1, 0)*8 + (1,)
assert K9.reduce(q**8) == q**8 and K9.reduce(q**9) == 0
assert K9.verify().certified
```

This calculation only needs reduced coordinate products. The
[benchmark](benchmarks/coinvariant-results.json) records zero Schubert images
and no ambient Schubert basis enumeration.

## Weighted mixed relations and the type D Euler generator

For the smallest two-step symplectic example, write `c=c_1` and `e=e_1` with
weights `(1,2)`. Its selected image has presentation

$$
\mathbb Q[c,e]/(c^4+c^2e+e^2,\ e^3).
$$

The standard monomials are `c^i*e^j`, with `0 <= i < 4` and `0 <= j < 3`.
Their images give dimension 12 and Hilbert series

$$
(1+t+t^2+t^3)(1+t^2+t^4)
=1+t+2t^2+2t^3+2t^4+2t^5+t^6+t^7.
$$

The ambient flag cohomology has total dimension 24. The mixed relation has
weighted degree four and describes how middle and tail classes interact.
The [two-step derivation](experiments/TWO_STEP_DERIVATION.md) proves the
general family presentation from classical Borel theory.

In type D, square-symmetric inputs alone do not contain the full orthogonal
invariant information. For a two-coordinate tail `x_3,x_4`, write

$$
q=x_3^2+x_4^2,\qquad u=x_3x_4,\qquad
u^2=e_2(x_3^2,x_4^2).
$$

The last equality is an identity of actual input polynomials. If the chosen
formal ring has only variables `q,u`, it does not introduce a third independent
formal variable `e_2`. The executable example retains `u`, verifies its
independence from `q` in polynomial degree two, and recovers a nine-dimensional
selected D4 image. The full type D coinvariant ideal likewise includes the
Euler product of **all** coordinates, as explained in
[MATHEMATICS.md](MATHEMATICS.md).

## Scope and reproducibility

All calculations here are rational. Integral image lattices, equivariant
parameters and quantum products are outside the implemented API. The
Springer connection identifies which invariant polynomials to feed into the
characteristic map; the current library does not automatically invert an
arbitrary representation character into such a polynomial.

Caches affect runtime, not the map or its answers. `A.clear_cache()` releases
the cached degreewise reports and presentations; a presentation is cached
separately for each resolved backend. Returning all relations can still be
expensive when the selected image or its coordinate support is large.

Run the mathematical companion from the `borel/` directory:

```bash
PYTHONPATH=src sage -python examples/mathematical_methods.py
```

The verification records linked from [PRESENTATIONS.md](PRESENTATIONS.md)
distinguish independent invariant-ideal checks, BGG cross-checks and complete
regression runs. The [documentation record](documentation-verification.json)
checks the Python blocks in these guides and their local links. Reproduce it
with `PYTHONPATH=src sage -python tools/verify_documentation.py`.
These methods implement established Borel and invariant
theory; the documented computational results do not assert novelty of those
underlying theorems.
