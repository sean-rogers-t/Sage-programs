# Recorded performance

Run from `borel/`:

```bash
PYTHONPATH=src sage -python benchmarks/run.py
```

The exact recorded results are in [`results.json`](results.json), including
setup time, all output coefficients, Weyl-group orders and frontier statistics.
These are single-run local measurements with SageMath 9.5 on Ubuntu under WSL2.

| Type | Tail size | Elementary index | Calculation | Setup | Output terms |
| --- | ---: | ---: | ---: | ---: | ---: |
| C8 | 6 | 2 | 0.015 s | 0.173 s | 3 |
| C12 | 10 | 5 | 0.166 s | 0.091 s | 6 |
| C20 | 18 | 5 | 0.529 s | 0.349 s | 6 |
| D12 | 10 | 5 | 0.188 s | 0.089 s | 7 |
| A15 | 13 | 5 | 0.038 s | 0.245 s | 1 |

The C20 calculation starts from just 11 fundamental-weight monomials and never
holds more than six frontier classes. Its Weyl group has about `2.55 * 10^24`
elements; the algorithm does not list them.

For the bounded identical-input comparisons, the historical C3 `k=2,i=1`
calculation took 0.306 seconds; the new calculation took 0.00186 seconds.
The C4 `k=3,i=1` calculation took 2.330 seconds historically and 0.00161 seconds
with the library. Exact results agree. New-map setup is separated in the JSON.
Both improvements combine the sparse engine with the complementary-function
shortcut. They are examples, not uniform speedup promises for every polynomial.

The benchmark also independently computes the literal tail polynomial in
orthonormal coordinates and checks every coefficient in the complete degree
basis, including zeros. Those checks happen outside the calculation timer.
Historical comparisons are limited to ranks 3 and 4 to keep the script bounded.

Higher degree and larger intermediate polynomial or Schubert support can still
increase cost substantially. Use `coefficient(f, word)` when a full expansion
is unnecessary; inspect `last_stats` to see the work performed; reuse a map for
repeated divisor products so its cover cache is retained.
