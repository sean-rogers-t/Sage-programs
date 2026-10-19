# Complete kernels of chosen Borel generators

`GeneratorAlgebra.presentation()` computes the entire rational kernel ideal of
the chosen generator map, a monomial basis of its image, and its finite Hilbert
series. Its completeness certificate uses exact image coordinates and a finite
monomial border. Classical types normally use coinvariant quotient normal
forms; other supported types use Schubert images. No user-supplied degree
cutoff is needed.

For the mathematical meaning of each method and returned field, see
[METHODS.md](METHODS.md). The geometric Borel map, Springer transport,
Schubert conventions and classical invariant ideals are developed in
[MATHEMATICS.md](MATHEMATICS.md).

All degrees below are **polynomial degrees**, half the cohomological degrees.
The formal variables have the positive degrees of their input Cartan
polynomials. The computation applies to any `GeneratorAlgebra`, including
classical and exceptional types and arbitrary parabolics. Its cost depends on
the selected image, which can still be large.

## A first example

```python
from sage_borel import BorelMap, GeneratorAlgebra

b = BorelMap("A2", omitted_nodes=1)  # projective plane
A = GeneratorAlgebra(b, {"h": b.weights[0]})
h, = A.domain.gens()
P = A.presentation()

assert P.domain is A.domain
assert P.generator_degrees == (1,)
assert P.relations == (h**3,)
assert P.standard_monomials == (1, h, h**2)
assert P.dimension == 3
assert P.hilbert_coefficients == (1, 1, 1)
print(P.hilbert_series)              # t^2 + t + 1
assert P.reduce(h**4 + 2*h**2) == 2*h**2
assert P.kernel_ideal() == A.domain.ideal(h**3)
assert P.backend == "coinvariant"
assert P.certificate.certified
assert P.verify().certified
```

The quotient `A.domain / P.kernel_ideal()` is the **selected image algebra**.
It need not be all of `H*(G/P; QQ)`. The basis consists of formal monomials whose
cohomology images are independent. The coinvariant backend constructs no
Schubert expansions; neither backend requires the complete ambient Schubert
basis to be enumerated.

## API

`A.presentation(backend="auto")` returns a frozen `KernelPresentation`.
The available backends are:

* `"auto"` uses `"coinvariant"` for standard-numbered types A/B/C/D and `"bgg"`
  for other supported Cartan types, including relabelled classical types.
* `"coinvariant"` computes exact normal forms in the full flag coinvariant
  quotient. Explicitly requesting it for an unsupported Cartan type raises
  `NotImplementedError`.
* `"bgg"` uses exact Schubert expansions from divided differences, preserving
  the original method for independent comparisons and exceptional types.

Presentations are cached separately by their resolved backend. Thus the default
call and an explicit `backend="coinvariant"` call share a cached presentation
when auto selects that backend. `A.clear_cache()` clears these presentations
together with the algebra's other computed data.

| API | Result |
| --- | --- |
| `P.domain` | The existing formal polynomial ring `A.domain` |
| `P.generator_degrees` | Positive polynomial degrees of the formal variables |
| `P.backend` | The actual backend used: `"coinvariant"` or `"bgg"` |
| `P.relations` | Tuple containing a reduced monic Groebner basis of the complete kernel |
| `P.standard_monomials` | Tuple of monomials giving a rational basis of the quotient and image |
| `P.dimension` | Total rational dimension of the selected image |
| `P.hilbert_coefficients` | Image dimensions in polynomial degrees zero through the last nonzero degree |
| `P.hilbert_series` | The Hilbert series, represented by a finite polynomial in `QQ[t]` |
| `P.certificate` | The exact completeness checks recorded for this presentation |
| `P.kernel_ideal()` | The defining ideal in `P.domain` |
| `P.reduce(p)` | The unique normal form of a formal polynomial in the standard monomial basis |
| `P.verify(backend=None)` | Recheck the presentation with its original backend, or explicitly select another backend |
| `P.verify().certified` | Whether the completeness and output-consistency checks pass |
| `P.as_dict()` | JSON-safe exact source, relation, basis and certificate data, using schema version 2 and recording the backend |

The Groebner basis uses the domain's existing ordinary degree reverse
lexicographic term order. The mathematical grading is supplied separately by
`generator_degrees`: ordinary total degree and the grading degree need not
agree. Inputs can have algebraic dependencies or zero cohomology image; neither
prevents the construction.

For example, a small classical presentation can be checked using both exact
coordinate systems:

```python
P = A.presentation(backend="coinvariant")
assert P.verify().certified                 # coinvariant normal forms
assert P.verify(backend="bgg").certified    # independent Schubert coordinates
Q = A.presentation(backend="bgg")
assert P.kernel_ideal() == Q.kernel_ideal()
```

Cross-backend verification may cost substantially more than the original
coinvariant calculation. The bounded `A.degree(d)` API and public Schubert
expansion operations retain their existing behavior.

This complements the bounded methods in [IMAGE_RELATIONS.md](IMAGE_RELATIONS.md).
`A.degree(d).relations` is a vector-space basis of one homogeneous kernel
component. `P.relations` instead generates the **entire kernel ideal**. The
number of Groebner basis elements need not be the minimum number of ideal
generators.

## Two-step symplectic example

```python
from sage_borel import two_step_symplectic

F = two_step_symplectic(n=3, b=1, k=1)
A = F.joint_algebra()
c1, e1 = A.domain.gens()
P = A.presentation()

assert P.generator_degrees == (1, 2)
assert set(P.relations) == {c1**4 + c1**2*e1 + e1**2, e1**3}
assert P.dimension == 12
assert P.hilbert_coefficients == (1, 1, 2, 2, 2, 2, 1, 1)
assert P.reduce(c1**4) == -c1**2*e1 - e1**2
assert P.reduce(e1**3) == 0
assert P.verify().certified
```

Here the selected image has presentation

$$
\mathbb Q[c_1,e_1]/(c_1^4+c_1^2e_1+e_1^2,\ e_1^3),
\qquad \deg(c_1)=1,\quad \deg(e_1)=2.
$$

Its standard monomials are `c1**i * e1**j` for `0 <= i < 4` and
`0 <= j < 3`, arranged in discovery order. Thus its Hilbert series is

$$
(1+t+t^2+t^3)(1+t^2+t^4)
=1+t+2t^2+2t^3+2t^4+2t^5+t^6+t^7.
$$

The ambient two-step flag has cohomology dimension 24. This selected
middle-plus-tail image has dimension 12. The independent
[two-step derivation](experiments/TWO_STEP_DERIVATION.md) supplies a separate
mathematical presentation for the family, rather than being an assumption of
the generic algorithm.

## The exact staircase construction

Write

$$
R=\mathbb Q[y_1,\ldots,y_r],\qquad \deg(y_i)=d_i>0,
\qquad \phi:R\longrightarrow H^*(G/P;\mathbb Q).
$$

The map substitutes the chosen Cartan polynomials and applies the Borel map.
Its kernel is an ideal homogeneous for the positive weights `d_i`.

Start with the candidate monomial `1`. Process candidates in increasing
weighted degree, breaking ties by the domain's ordinary degree reverse
lexicographic order. Retain a monomial when its exact image is
linearly independent of the retained images in the same weighted degree.
Only retained monomials enqueue their products by the formal variables.

If a candidate image is dependent, exact rational linear algebra expresses it
in the earlier retained images. Record the monic relation

$$
m-\sum_{b\in B}q_b b,
\qquad \phi(m)=\sum_{b\in B}q_b\phi(b),
$$

where the sum uses retained monomials of the same weighted degree. Its leading
monomial is `m`. Candidates divisible by an already discovered leading
monomial are discarded. Continue until no candidate remains. This constructs
the standard monomial staircase and the relations at its boundary.

The construction is FGLM-style: it discovers an ideal basis from exact linear
dependencies of monomial images in a finite-dimensional algebra. It is not
literally the original change-of-order FGLM algorithm, whose input includes an
existing Groebner basis.

### Coinvariant image coordinates

Let `S` be the rational polynomial ring on the Cartan subalgebra, `W` the full
Weyl group and `W_P` the Levi Weyl group. The Borel presentation identifies

$$
H^*(G/B;\mathbb Q)=S/(S^W_+S),\qquad
H^*(G/P;\mathbb Q)=S^{W_P}/(S^W_+S^{W_P}).
$$

The inclusion of invariants induces an injection from the second quotient to
the first. Indeed, averaging over `W_P` is an `S^W`-linear retraction from `S`
to `S^{W_P}`. Applying this retraction to an expression in the invariant ideal
proves the contraction identity

$$
(S^W_+S)\cap S^{W_P}=S^W_+S^{W_P}.
$$

All chosen input generators are already checked to be Levi invariant.
Consequently, computing their relations and image ranks in the full flag
coinvariant quotient gives exactly the same answers as computing in the
original partial flag cohomology ring. The completeness proof below applies
to either coordinate system.

The backend converts fundamental-weight polynomials to standard coordinate
polynomials before multiplying them. Type `A_r` uses `r+1` coordinates and the
trace-zero fundamental weights; the relation `x1+...+x_(r+1)=0` is included.
Types B/C/D use `r` coordinates, with the half-integral spin-weight factors
included in B and D. All coefficients remain rational.

For `N` coordinates, the type A invariant ideal has the triangular basis

$$
g_i=h_i(x_i,\ldots,x_N),\qquad 1\leq i\leq N.
$$

For types B/C, use `N=r` and replace every `x_j` by `x_j^2`. To verify equality
with the usual elementary-symmetric invariant ideal, write `y_j=x_j` in A and
`y_j=x_j^2` in B/C. The generating-function identity

$$
\left(\prod_{j=1}^N(1-y_jt)\right)
\left(\prod_{j=i}^N(1-y_jt)^{-1}\right)
=\prod_{j=1}^{i-1}(1-y_jt)
$$

has zero coefficient of `t^i`. Hence

$$
g_i=\sum_{j=1}^{i}(-1)^{j+1}e_j(y_1,\ldots,y_N)
             h_{i-j}(y_i,\ldots,y_N).
$$

This is an invertible triangular change of ideal generators: the coefficient
of `e_i` is `(-1)^(i+1)`. In ordinary degree reverse lexicographic order with
`x1 > ... > xN`, the leading monomials of the `g_i` are `x_i^i` in A and
`x_i^(2*i)` in B/C. They are pairwise coprime, so Buchberger's product criterion
proves that these polynomials already form a Groebner basis. No general
Groebner basis calculation is needed in these three types.

Type D uses the invariant generators

$$
e_1(x^2),\ldots,e_{r-1}(x^2),\qquad x_1\cdots x_r,
$$

including the Euler product. The backend computes a Groebner basis from these
generators using Sage. Omitting the Euler product would give the wrong
coinvariant quotient.

Generator images and intermediate products are reduced and cached. A newly
requested monomial is obtained from a cached parent by multiplying by one
reduced generator and reducing again; the algorithm avoids expanding the raw
high powers in fundamental-weight variables. The quotient's complete ambient
monomial basis is never enumerated.

### Why the existing term order is preserved

The discovery order first compares positive weighted degree and then compares
ordinary degree reverse lexicographic order. It is a global monomial order:
it is well ordered and multiplication preserves comparisons. Inside one
weighted homogeneous component, it agrees with the domain's term order.

Every recorded relation is weighted homogeneous, so it has the same leading
monomial in both orders. More generally, let `f` belong to the kernel and take
the weighted homogeneous component containing its ordinary-degrevlex leading
monomial. That component also belongs to the kernel, and its leading monomial
in the discovery order is the ordinary-degrevlex leading monomial of `f`.
Consequently, a Groebner basis for the discovery order is also a Groebner basis
for the existing domain order. No new formal polynomial ring or conversion of
the public variables is required.

The tails of the recorded relations use only standard monomials. Their
leading monomials do not divide one another, and every leading coefficient is
one. The resulting basis is therefore reduced and monic.

`P.standard_monomials` belongs to the formal generator ring, and its images
form a basis of the selected image. This is different from both the ambient
coordinate monomials used internally by the coinvariant backend and the
Schubert basis used by BGG. In particular, `P.reduce(p)` accepts a formal
polynomial and returns its representative in this quotient basis. Use
`A.evaluate(p)` to obtain its Schubert expansion.

## Why the certificate proves completeness

Let `B` be the reported standard monomials, `G` the reported relations, and
`L` the monomial ideal generated by the leading monomials of `G`. The
`KernelCertificate` records five conditions:

| Condition | Meaning |
| --- | --- |
| `relation_images_zero` | Every relation in `G` has zero exact image in the selected backend |
| `basis_images_independent` | The images of `B` are linearly independent in every weighted degree |
| `standard_monomials_reduced` | No monomial in `B` is divisible by a leading monomial in `G` |
| `border_closed` | For every `b` in `B` and variable `y_i`, the monomial `b*y_i` lies in `B` or in `L` |
| `contains_unit` | The monomial `1` belongs to `B` |

The implementation additionally checks `groebner_basis_reduced` (monic
homogeneous relations, minimal leading terms and tails in the standard basis)
and `summary_data_correct` (the formal domain, weights and Hilbert data match
the reported basis). `P.verify()` recomputes the geometric top degree rather
than trusting the recorded value. `certified` is true when all these checks
hold. The five mathematical conditions above imply
completeness, rather than merely checking that some proposed relations vanish.

First, `B` is precisely the set of monomials outside `L`. To see this, suppose
there were a monomial outside both `B` and `L`, and choose one of smallest
ordinary total degree. It is not the unit. Remove one variable factor. The
parent is outside `L`, because monomial ideals are closed under multiplication;
by minimality the parent belongs to `B`. Border closure now contradicts the
choice of the original monomial. The reverse inclusion follows from
`standard_monomials_reduced`.

Division by the monic relations terminates in a linear combination of `B`.
Let `J` be the ideal generated by `G`. The zero-image condition gives
`J <= ker(phi)`. For any `f` in `ker(phi)`, its remainder is a linear
combination of `B` with zero image. Independence forces that remainder to be
zero, so `f` belongs to `J`. Therefore

$$
J=\ker(\phi),\qquad R/J\cong\operatorname{im}(\phi).
$$

The same argument shows that `G` is a Groebner basis: a nonzero kernel
polynomial with a standard leading monomial would retain that term after
division and could not have zero remainder. This proof also gives the unique
normal form returned by `P.reduce(p)`.

The certificate uses the selected exact evaluator and records its backend.
`P.verify()` defaults to `P.backend`; passing a different backend recomputes
the image checks in that coordinate system. The test suite also compares
complete ideals against independent elimination calculations. A single
certificate is an exact completeness check, rather than a claim that two
independent evaluators were run.

## Termination and cost

Let `D` be the dimension of the selected image. Each retained monomial adds an
independent image, so there are exactly `D` retained monomials. Only these
monomials enqueue variable multiples. At most `1 + r*D` distinct candidates
are therefore created, including the initial unit; pruning can reduce the
number of actual image evaluations.

There is also a geometric degree bound. Put

$$
N=\dim_{\mathbb C}(G/P)=\ell(w_0w_{0,P}),
$$

the highest nonzero polynomial degree of the cohomology ring. A retained
monomial has nonzero image and hence weighted degree at most `N`. Every
candidate is the unit or a variable multiple of a retained monomial, so no
candidate has degree greater than `N + max(generator_degrees)`. Images above
`N` are zero by the grading. These bounds prove finite termination; completion
is detected by the exhausted border, not by assuming a guessed degree is
sufficient.

This is output-sensitive in the number of selected image basis monomials,
not a promise of fast calculations in every rank or exceptional type.
Coinvariant normal forms or Schubert supports, rational coefficients and input polynomials can grow, and
the selected image itself can approach the full cohomology ring. In
particular, supporting type E8 as an input does not make every E8 presentation
practical. The degreewise API remains useful when only a bounded portion is
needed.

In version 0.4.0, an exploratory `C9/P8` calculation with the single tail
generator `x9**2` exceeded a 45-second process limit while other validation was
running. That historical observation was a timeout, not a completed timing or
a controlled speed comparison.

Version 0.5.0 computes this example using coinvariant normal forms. The
expected complete answer is

$$
\mathbb Q[e]/(e^9),\qquad \deg(e)=2,\qquad
\operatorname{Hilb}=1+t^2+\cdots+t^{16},\qquad \dim=9.
$$

The [benchmark record](benchmarks/coinvariant-results.json) records completed
calculations and verification timings with their exact outputs and runtime
context. Timings there distinguish calculation inside an initialized Sage
process from any process-startup cost. The optimized rank-9 regression test
also rejects any attempted call to the Schubert expansion evaluator, checking
that the new backend actually avoids the previous bottleneck.

## Verification records

The version 0.4.0 [independent suite](presentation-verification.json) checked 13 complete
elimination ideals in types A/C/D/G2, three exceptional ample-divisor kernels,
65 normal forms and 32 vanishing relations. It also rejected forged degree
bounds, incomplete bases and inconsistent Hilbert data. The complete 39-test
regression suite passed for that release under SageMath 9.5; its package and
historical-file checks are recorded in
[presentation-validation.json](presentation-validation.json). These retained
records describe the original BGG implementation.

The version 0.5.0 [backend verification record](coinvariant-verification.json)
adds comparisons between BGG and coinvariant kernels, normalization and
fallback checks, and the rank-9 test that forbids Schubert expansion calls.
All 49 regression methods pass; the
[release validation record](coinvariant-validation.json) also records the
isolated wheel import and installed-package examples.

The [runnable example](examples/kernel_presentations.py) writes
[current sample presentations](coinvariant-examples.json) for projective space,
Gr(2,4), the two-step C3 image, a D4 Euler example, the full G2 flag ring and an
F4 divisor image. Run it with `--output PATH` to reproduce the JSON data.
The [version 0.4.0 samples](presentation-examples.json) remain available.

## Mathematical context

The rational Borel presentation and finite-dimensional flag cohomology are
classical. This feature computes a certified presentation of a user-selected
subalgebra of that known cohomology ring. It makes no novelty claim for the
underlying Borel theorem or for the staircase linear-algebra method.

* A. Borel, *Sur la cohomologie des espaces fibres principaux et des espaces
  homogenes de groupes de Lie compacts*, Annals of Mathematics 57 (1953),
  115-207, especially Theorem 26.1. [Original paper](https://www.sas.rochester.edu/mth/sites/doug-ravenel/otherpapers/Borel-Sur.pdf).
* J.-C. Faugere, P. Gianni, D. Lazard and T. Mora, *Efficient Computation of
  Zero-dimensional Groebner Bases by Change of Ordering*, Journal of Symbolic
  Computation 16 (1993), 329-344. [Publisher and DOI](https://www.sciencedirect.com/science/article/pii/S0747717183710515).
