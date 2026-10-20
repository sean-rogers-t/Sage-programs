# Two-step symplectic experiment

The completed grid uses `b,k` in `{1,2}`, `a=1,...,5`, `n=a+b+k`,
and polynomial degrees `0,...,8` (cohomological degrees `0,...,16`).
The selected image uses the joint `GL(b) x Sp(2k)` generators.

All **20 cases and 180 homogeneous ranks** agree with the independent
formal presentation, and all **23 presentation relations** within the
degree bound have zero BGG image. The separate test suite also compares
whole kernel spaces with coordinate-coinvariant Groebner reductions.

| n | b | k | a | First relation found through degree 8 | Degree-8 domain | Degree-8 image |
| ---: | ---: | ---: | ---: | --- | ---: | ---: |
| 3 | 1 | 1 | 1 | 4 | 5 | 0 |
| 4 | 1 | 1 | 2 | 6 | 5 | 2 |
| 5 | 1 | 1 | 3 | 8 | 5 | 4 |
| 6 | 1 | 1 | 4 | None through 8 | 5 | 5 |
| 7 | 1 | 1 | 5 | None through 8 | 5 | 5 |
| 4 | 1 | 2 | 1 | 4 | 9 | 2 |
| 5 | 1 | 2 | 2 | 6 | 9 | 6 |
| 6 | 1 | 2 | 3 | 8 | 9 | 8 |
| 7 | 1 | 2 | 4 | None through 8 | 9 | 9 |
| 8 | 1 | 2 | 5 | None through 8 | 9 | 9 |
| 4 | 2 | 1 | 1 | 4 | 15 | 5 |
| 5 | 2 | 1 | 2 | 6 | 15 | 11 |
| 6 | 2 | 1 | 3 | 8 | 15 | 14 |
| 7 | 2 | 1 | 4 | None through 8 | 15 | 15 |
| 8 | 2 | 1 | 5 | None through 8 | 15 | 15 |
| 5 | 2 | 2 | 1 | 4 | 22 | 11 |
| 6 | 2 | 2 | 2 | 6 | 22 | 18 |
| 7 | 2 | 2 | 3 | 8 | 22 | 21 |
| 8 | 2 | 2 | 4 | None through 8 | 22 | 22 |
| 9 | 2 | 2 | 5 | None through 8 | 22 | 22 |

The first computed relation degrees are 4, 6 and 8 for a=1,2,3.
For a=4,5 there is no relation through degree 8. Finite computations
alone establish only these bounded facts. The
[invariant-theory derivation](TWO_STEP_DERIVATION.md) separately proves the
full ideal `(h_(a+1)(Z),...,h_n(Z))` over QQ, and hence the sharp first
relation degree `2(a+1)`. No novelty claim is made.

For n=3,b=k=1, a mixed relation is `c1^4+c1^2*e1+e1^2=0`, and
`e1^3=0`. The image Hilbert polynomial has coefficients
`(1,1,2,2,2,2,1,1)`.

The [complete exact report](two_step_symplectic_results.json) contains input
polynomials, generator degrees, every computed kernel basis, structured
rational coefficients and timing information. No complete ambient Weyl
group or ambient Schubert basis was enumerated for this grid.

Run from `borel/`:

```bash
PYTHONPATH=src sage -python experiments/two_step_symplectic.py
python tools/summarize_two_step.py
```

The recorded SageMath 9.5 run took 106.97 seconds; timing is machine dependent.

The script accepts degree/rank-grid bounds; edit or use its `--help` options
for a new experiment. This summary requires the recorded degree-8 grid.
