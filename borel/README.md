# sage-borel

A small SageMath library for exact Borel characteristic maps in a Schubert basis.
The historical programs in the repository root are preserved. This package is a
separate implementation incorporating the corrections found during their audit.

The [exceptional generator catalogue](EXCEPTIONAL.md) adds three-parameter
calculations such as `exceptional_image("F4", 1, "e2")`, with named polynomial
inputs and computed Schubert images. It includes the corrected Part III inputs,
the missing E6/P6 Euler generator, and the displayed E8 dissertation cases.

It computes

\[
\beta:\operatorname{Sym}(\mathfrak t^*)^{W_P}\longrightarrow H^{2*}(G/P;\mathbb Q),
\qquad \beta(\omega_i)=\sigma_{s_i}.
\]

Here `i` is an omitted node; on the full flag variety every node is allowed.

This is the characteristic-map stage of the research calculations. For a
representation-ring element written as the Springer pullback `q(Theta(t))`, pass
the invariant polynomial `q(omega)` in fundamental-weight coordinates to `expand`.
The independent exceptional Springer calculations remain in
[`../articles/calculations`](../articles/calculations/).

## Run it

Requires SageMath and Python 3.9 or later; verified with SageMath 9.5 / Python 3.10.
Run from this `borel/` directory without installing anything:

```bash
PYTHONPATH=src sage -python examples/basic.py
```

The exceptional catalogue has its own example:

```bash
PYTHONPATH=src sage -python examples/exceptional.py
```

Or install into your Sage environment and import from any directory:

```bash
sage -pip install .
sage -python examples/basic.py
```

On this Windows machine these commands run inside the Ubuntu WSL installation.
Use `**` for exponentiation in Python files. In a Sage notebook, `^` is also
preparsed by Sage. This package requires the Sage runtime; ordinary Python alone
does not provide its root systems or exact polynomial arithmetic.

```python
from sage_borel import BorelMap, classical_borel

# Your old type C parameters: rank n=12, tail size k=10.
b = classical_borel("C", 12, 10)
answer = b.tail_elementary(10, 5)  # beta(e_5(x_3^2,...,x_12^2))
print(answer)
print(answer.to_terms())          # JSON-safe words and exact coefficients
print(b.last_stats)               # frontier sizes and operation counts

# The corrected type D example from the audit:
d = classical_borel("D", 4, 3)
assert d.tail_elementary(3, 1).to_terms() == [
    {"word": [2, 1], "coefficient": "-1"}
]
```

`classical_borel("A", n, k)` uses **GL(n) dimension**, matching the old type A
programs. `BorelMap("A3")` instead uses **Dynkin rank** and means SL(4).
The factory omits node `n-k`. Its ranges are A: `1 <= k < n`; B/C:
`0 <= k < n`; D: `2 <= k < n`. Use `BorelMap` directly for spin or other
parabolics.

## General inputs and parabolics

```python
b = BorelMap("A3", omitted_nodes=(1, 3))
w1, w2, w3 = b.weights
answer = b.expand((w1 + 2*w3)**3 + 7)
print(b.coefficient((w1 + 2*w3)**3, [1, 2, 3]))

f = BorelMap("F4")
print(f.expand(f.weights[3]**3))

e = BorelMap("E7", omitted_nodes=7)
print(e.divisor_product([7, 7, 7]))
```

Supported types are finite irreducible crystallographic A, B, C, D, E, F, G.
`omitted_nodes=None` gives the full flag variety; an integer gives a maximal
parabolic; an iterable gives an arbitrary parabolic; `()` gives G/G. Retained
nodes generate the Levi Weyl group. Inputs to a parabolic map must be invariant
under those nodes; the library checks this and rejects invalid inputs.

Weights use Sage's Dynkin numbering. Polynomial degree `d` means cohomological
degree `2*d`. A squared tail elementary polynomial of index `i` has polynomial
degree `2*i`, hence cohomological degree `4*i`. Type A coordinate weights have
sum zero; the GL(n) central character has zero characteristic-map image.

## API

| Operation | Meaning |
| --- | --- |
| `b.ring`, `b.weights`, `b.W` | Exact QQ polynomial ring, fundamental weights, Weyl group |
| `b.expand(f)` | Full Schubert expansion, including zero, constants and mixed degrees |
| `b.coefficient(f, word)` | One coefficient without generating other classes |
| `b.reflect(f, i)`, `b.divided_difference(f, i)` | Simple reflection and exact BGG divided difference |
| `b.chevalley(i, word)` | Multiply a Schubert class by the divisor `omega_i` |
| `b.divisor_product(indices)` | Repeated divisor product, with empty product equal to 1 |
| `b.multiply_divisor(expansion, i)` | Multiply an existing sparse expansion by a divisor |
| `b.classical_coordinates()` | Standard coordinate weights for A/B/C/D |
| `b.tail_polynomial(k, i)` | Literal tail elementary input polynomial |
| `b.tail_elementary(k, i)` | Its Borel image, using the complementary shortcut when cheaper |
| `b.head_polynomial(r, i)`, `b.head_elementary(r, i)` | Unsquared head elementary polynomial and image |
| `b.euler_polynomial(k)`, `b.euler(k)` | Type D tail Euler product and image |
| `b.basis(d)` | Minimal right-coset representatives in degree `d`, generated on demand |
| `b.poincare_dual(word)` | Dual index `w0 * w * w0P` |
| `b.longest_representative()` | Top Schubert index, without enumerating W |
| `b.clear_cache()` | Release generated bases, words and per-class products |
| `to_type_a(w)` | Structural signed-permutation embedding for B/C/D |
| `exceptional_image(type, parabolic, generator)` | Compute a named exceptional generator image |
| `exceptional_generators(type, parabolic)` | List names, degrees and source metadata |
| `exceptional_polynomial(type, parabolic, generator)` | Inspect its fundamental-weight polynomial |
| `exceptional_borel(type, parabolic)` | Reusable map with named generator methods |

`tail_polynomial` and `tail_elementary` default to unsquared coordinates in A,
squared coordinates in B/C/D; `squared=True/False` overrides this. Classical
coordinate helpers and the signed embedding require standard Dynkin numbering;
generic Borel operations also accept integer relabelings of Dynkin nodes.

An `Expansion` is an immutable mapping from Weyl elements to rational
coefficients. `answer[word]` returns zero for an absent valid word;
`answer.as_dict()` returns a copy; `answer.to_terms()` serializes reduced words
and rational numbers as strings. Addition and scalar multiplication are provided.
General multiplication of two Schubert expansions is outside this small API;
use `expand(f*g)` when polynomial representatives are available. Weyl elements
must have the map's parent; passing a list of Dynkin indices converts a word
explicitly. Parabolic divisor operations require omitted nodes.

## Why this is faster

The original routines repeatedly enumerate every Weyl element, reconstruct
roots, expand all subsets and accumulate duplicate term lists. At type C rank 8,
W already has 10,321,920 elements.

The replacement follows nonzero divided differences only. Braid-equivalent
paths merge once, rather than being counted as distinct contributions. Divisor
products use cached positive-root covers. Elementary polynomials use dynamic
recurrences; tail images can use
`e_i(tail^2) = (-1)^i h_i(head^2)` in the coinvariant quotient. This often makes
the input much smaller. [The mathematical conventions and algorithm](MATHEMATICS.md)
explain why these optimizations preserve the result.

Measured results, including setup times and complete-degree coefficient checks,
are in [`benchmarks/results.json`](benchmarks/results.json). The benchmark script
records actual timings on the local machine; these depend on hardware, degree,
input size and output size. Large dense expansions can still be expensive.

```bash
PYTHONPATH=src sage -python benchmarks/run.py
```

## Verification and migration

The independent test oracle uses explicit coordinate swaps/sign changes and
Schubert representatives derived from the positive-root discriminant. It does
not reuse the library's Chevalley root/coroot product formula. The exhaustive
suite checks all **1,088** divisor products in B2, C2, B3, C3 and D4, plus
classical generators, mixed degrees, parabolics, duality/intersection pairings,
historical counterexamples, exceptional-type sanity checks, and rank-10 word
handling. See [`verification.json`](verification.json) for counts and environment.

```bash
PYTHONPATH=src sage -python -m unittest discover -s tests -v
PYTHONPATH=src BOREL_EXHAUSTIVE=1 BOREL_TEST_REPORT=verification.json \
  sage -python -m unittest discover -s tests -v
```

[`MIGRATION.md`](MIGRATION.md) maps your old functions to the new API and records
the corrected examples. [`legacy_manifest.json`](legacy_manifest.json) records
canonical Git and local checkout SHA-256 hashes of the 35 unchanged historical
program files. `python verify_legacy.py` verifies repository contents;
`python verify_legacy.py --checkout-bytes` also verifies this machine's recorded
checkout bytes. These are recorded separately because Windows Git converts line
endings on checkout.
