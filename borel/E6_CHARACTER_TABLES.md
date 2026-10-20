# E6/P1: complete integral character-image comparison

See the [proof and methods](E6_CHARACTER_FORM.md) and
[exact dataset](experiments/e6-character-form-results.json).

Polynomial degree d means cohomological degree 2d. All image lattices
have full rank in their target degree. Factors 1 are retained.
The literal character image is contained in each of the two prior images;
containments are checked in common Schubert coordinates. Degrees above 16 vanish.

| d | Literal Levi characters | Original catalogue | Individually primitive images |
| --- | --- | --- | --- |
| 0 | (1) | (1) | (1) |
| 1 | (18) | (3) | (1) |
| 2 | (108) | (3) | (1) |
| 3 | (1944) | (9) | (1) |
| 4 | (3888, 23328) | (3, 288) | (1, 16) |
| 5 | (23328, 419904) | (3, 432) | (1, 16) |
| 6 | (139968, 839808) | (3, 144) | (1, 16) |
| 7 | (2519424, 5038848) | (9, 144) | (1, 16) |
| 8 | (5038848, 15116544, 90699264) | (3, 144, 6912) | (1, 16, 256) |
| 9 | (90699264, 544195584) | (9, 432) | (1, 16) |
| 10 | (544195584, 1088391168) | (9, 144) | (1, 16) |
| 11 | (3265173504, 19591041024) | (9, 432) | (1, 16) |
| 12 | (19591041024, 117546246144) | (9, 432) | (1, 16) |
| 13 | (705277476864) | (54) | (2) |
| 14 | (1410554953728) | (18) | (2) |
| 15 | (25389989167104) | (54) | (2) |
| 16 | (76169967501312) | (54) | (2) |

## Quotients of the previous images by the character image

These are quotients between embedded lattices, not ratios of individual
ambient Smith factors. They use the integer containment matrices.

| d | Catalogue / character | Primitive / character |
| --- | --- | --- |
| 0 | (1) | (1) |
| 1 | (6) | (18) |
| 2 | (36) | (108) |
| 3 | (216) | (1944) |
| 4 | (81, 1296) | (486, 11664) |
| 5 | (972, 7776) | (2916, 209952) |
| 6 | (5832, 46656) | (17496, 419904) |
| 7 | (34992, 279936) | (314928, 2519424) |
| 8 | (6561, 209952, 1679616) | (19683, 1889568, 45349632) |
| 9 | (1259712, 10077696) | (11337408, 272097792) |
| 10 | (7558272, 60466176) | (68024448, 544195584) |
| 11 | (45349632, 362797056) | (408146688, 9795520512) |
| 12 | (272097792, 2176782336) | (2448880128, 58773123072) |
| 13 | (13060694016) | (352638738432) |
| 14 | (78364164096) | (705277476864) |
| 15 | (470184984576) | (12694994583552) |
| 16 | (1410554953728) | (38084983750656) |

## Removing the uniform trace-index scaling

The character form in degree d is 6^d times the invariant root-lattice
form Sym_Z(Q)^(W_L). The following factors are therefore the character
Smith factors divided by 6^d, not an alternative character normalization.

| d | Root-lattice image factors |
| --- | --- |
| 0 | (1) |
| 1 | (3) |
| 2 | (3) |
| 3 | (9) |
| 4 | (3, 18) |
| 5 | (3, 54) |
| 6 | (3, 18) |
| 7 | (9, 18) |
| 8 | (3, 9, 54) |
| 9 | (9, 54) |
| 10 | (9, 18) |
| 11 | (9, 54) |
| 12 | (9, 54) |
| 13 | (54) |
| 14 | (18) |
| 15 | (54) |
| 16 | (27) |

Regenerate with `python tools/render_e6_character_tables.py`.
