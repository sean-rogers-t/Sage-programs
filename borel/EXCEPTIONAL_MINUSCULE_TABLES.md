# Complete E6/P1 and E7/P7 integral image tables

Generated from [the exact dataset](experiments/exceptional-minuscule-results.json).
See [the mathematical guide](EXCEPTIONAL_MINUSCULE.md) for the source forms,
published presentations, integer generator identities and independent checks.

Degree d means H^(2d). Each tuple includes all Smith factors, including 1.
Every image has full rational rank. Original catalogue images are already
integral; the second column of factors divides each generator image by its
own Schubert content before forming products. Higher target degrees vanish.

## E6/P1

Complex dimension: 16. Total Schubert classes: 27.

| d | Original catalogue | Individually primitive images |
| --- | --- | --- |
| 0 | (1) | (1) |
| 1 | (3) | (1) |
| 2 | (3) | (1) |
| 3 | (9) | (1) |
| 4 | (3, 288) | (1, 16) |
| 5 | (3, 432) | (1, 16) |
| 6 | (3, 144) | (1, 16) |
| 7 | (9, 144) | (1, 16) |
| 8 | (3, 144, 6912) | (1, 16, 256) |
| 9 | (9, 432) | (1, 16) |
| 10 | (9, 144) | (1, 16) |
| 11 | (9, 432) | (1, 16) |
| 12 | (9, 432) | (1, 16) |
| 13 | (54) | (2) |
| 14 | (18) | (2) |
| 15 | (54) | (2) |
| 16 | (54) | (2) |

## E7/P7

Complex dimension: 27. Total Schubert classes: 56.

| d | Original catalogue | Individually primitive images |
| --- | --- | --- |
| 0 | (1) | (1) |
| 1 | (1) | (1) |
| 2 | (1) | (1) |
| 3 | (1) | (1) |
| 4 | (1) | (1) |
| 5 | (1, 933120) | (1, 81) |
| 6 | (1, 933120) | (1, 81) |
| 7 | (1, 933120) | (1, 81) |
| 8 | (1, 933120) | (1, 81) |
| 9 | (1, 186624, 76187381760) | (1, 81, 19683) |
| 10 | (1, 186624, 10883911680) | (1, 81, 6561) |
| 11 | (1, 186624, 10883911680) | (1, 81, 6561) |
| 12 | (1, 186624, 10883911680) | (1, 81, 6561) |
| 13 | (1, 186624, 10883911680) | (1, 81, 6561) |
| 14 | (2, 186624, 10883911680) | (1, 81, 6561) |
| 15 | (2, 186624, 10883911680) | (1, 81, 6561) |
| 16 | (2, 186624, 10883911680) | (1, 81, 6561) |
| 17 | (2, 186624, 10883911680) | (1, 81, 6561) |
| 18 | (2, 186624, 2176782336) | (1, 81, 6561) |
| 19 | (6, 373248) | (3, 81) |
| 20 | (6, 373248) | (3, 81) |
| 21 | (6, 373248) | (3, 81) |
| 22 | (6, 373248) | (3, 81) |
| 23 | (6) | (3) |
| 24 | (6) | (3) |
| 25 | (6) | (3) |
| 26 | (6) | (3) |
| 27 | (6) | (3) |

Regenerate with `python tools/render_exceptional_minuscule_tables.py`.
