# Moving calculations from the historical scripts

All historical `.sage`, `.py` and `.txt` program files are kept unchanged. Import
the new package inside SageMath; no historical loader is needed.

```python
from sage_borel import BorelMap, classical_borel, to_type_a
b = classical_borel("C", n, k)
```

| Old function | New calculation |
| --- | --- |
| `cohom_0(n,k,h)` | `classical_borel("C",n,k).tail_elementary(k,h)` |
| `cohom_0_b(n,k,h)` | `classical_borel("B",n,k).tail_elementary(k,h)` |
| `cohom_0_d(n,k,h)` | `classical_borel("D",n,k).tail_elementary(k,h)` |
| `cohom_0_spec(n,k)` | `classical_borel("D",n,k).euler(k)` |
| `chevalley(n,a,word)` | `BorelMap(["C",n]).chevalley(a,word)` |
| `chevalley_b/d` | `BorelMap` with the corresponding Dynkin type |
| `finchev(n,indices)` | `BorelMap(["C",n]).divisor_product(indices)` |
| `finchev_b/d` | `divisor_product` on the corresponding full flag map |
| `mincp(n,k,m)` | `classical_borel("C",n,k).basis(m)` |
| `cp_long(n,k)` | `classical_borel("C",n,k).longest_representative()` |
| `PD(n,k,w)` and variants | `b.poincare_dual(w)` |
| `sr/sr_b` | `to_type_a(w)`; also supports the corrected D fork |
| `redec(w)` | `list(w.reduced_word())` |

The new result is an `Expansion`, rather than a mutable list of `[coefficient,w]`
pairs. For an old-style consumer, use `[[c,w] for w,c in answer.items()]`.
Inputs and outputs retain their Cartan and parabolic context. If a historical
Weyl element has a different parent, convert it using
`b.W.from_reduced_word(list(old_w.reduced_word()))`.

The historical type A files supply parabolic and duality helpers rather than
a complete elementary-image pipeline. The new package adds that pipeline:
`classical_borel("A",n,k).tail_elementary(k,h)` computes the unsquared tail image.

For the analogous exceptional workflow, use
`exceptional_image("F4", parabolic=1, generator="e2")`. The three arguments
select the group, omitted Dynkin node and named Levi-invariant polynomial.
See the [exceptional guide](EXCEPTIONAL.md) for available generators and
normalization conventions.

## Corrections included

The type D elementary construction now includes every tail coordinate. The old
`range(1,k-1)` omitted one coordinate and duplicated another; the corrected
historical range would be `range(2,k)`. The new coordinate implementation avoids
that indexing pattern entirely. In D4, `k=3,h=1`, the correct image is
`-sigma_(s2*s1)`.

Repeated type D divisor multiplication stays in D throughout. The historical
`finchev_d` called the B routine after the first two factors. For example, the
correct D4 product of three `omega_1` factors is
`sigma_(s3*s2*s1) + sigma_(s4*s2*s1)`. Both fork classes are retained. The Euler
input is unchanged mathematically; its image is now evaluated in the correct
type. For D4, `k=3`, it is `-sigma_(s3*s2*s1) + sigma_(s4*s2*s1)`.

The A longest-class helper now stays in A, and duality uses the correct order
`w0*w*w0P` for every parabolic, including D spin parabolics. For A2 with omitted
node 2, the dual of the identity is `s1*s2`.

Reduced words use Sage's structural API, so rank 10 and custom generator names
are handled correctly. The D signed embedding uses
`s_n -> r_n*r_(n+1)*r_(n-1)*r_n`. Zero, empty and single-factor products have
explicit behavior; negative or fractional indices are rejected. Rational
coefficients remain exact.

The old `open.sage` uses Python 2 syntax, and the old mathematical `.py` copies
contain Sage exponent syntax. `Schubert3.py` is a saved HTML page containing a
separate intersection-theory toolkit; it is preserved as historical material and
is not loaded by this package.

## Preservation record

[`legacy_manifest.json`](legacy_manifest.json) records each historical program's
canonical Git SHA-256 hash, local checkout hash and the Git baseline. Windows
checkout line endings are recorded separately from repository bytes. Run
`python verify_legacy.py` to confirm unchanged Git contents, or add
`--checkout-bytes` to check this machine's recorded files exactly. The root README
is updated only to link here.
