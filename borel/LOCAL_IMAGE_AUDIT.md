# Audit of the first local-image and equivariant work package

Date: 20 October 2026. This is a stepwise internal proof/source audit by
the current research session, with a separate SymPy computation importing
neither Sage nor the library. It is not external expert review, formal
verification or a literature-priority certification.

The audited notes are the [local bound and projective family](LOCAL_IMAGE_BOUND.md)
and the [equivariant specialization counterexample](EQUIVARIANT_IMAGE_PILOT.md).
No gap was found in the arguments under their stated hypotheses.

## Claims and dependencies

| Claim | Integral argument checked | Scope and verdict |
| --- | --- | --- |
| Lambda_rho is contained in X(T) | Differentiate the integral Laurent pullback at 1 and use d theta_1=id | Valid for the normalized construction; no irreducibility assumption |
| The full literal Levi source has integral Borel image | Degree-one source theorem, integral full-flag characteristic map, primitive Schubert pullback | Valid in the simply connected simple characteristic-zero scope of the source theorem; arbitrary selected rational generators are not included |
| A single t_G supplies every parabolic Schubert lift | F_w=partial_(w^-1 w0)F is Levi invariant by left descents; nil-Coxeter multiplication gives exactly t_G sigma_w | Valid; no Reynolds average or extra Levi-prime exclusion |
| Every degree-d cokernel is killed by t_G e_rho^d | e_rho X lies in Lambda_rho; apply scalar clearing to the homogeneous invariant lifts | Valid for the full literal source; not necessarily sharp |
| Good primes give the whole integral cohomology after localization | Invert t_G e_rho in the preceding bound | Standard characteristic-map consequence applied to these sources; no novelty claim |
| Good primes also recover the full equivariant image | Integral equivariant Schubert pullback, localization of invariant kernels and triangular induction on Schubert length | Valid using the standard free equivariant Schubert basis; no universal bad-prime faithfulness assertion |
| Fixed-degree specialization kernels are finite with primes dividing t_G e_rho | Rational and good-prime equivariant surjectivity plus finite graded integer modules | Valid; the sharp bad-prime exponent is not determined |
| Defining SL_n has Lambda=Q | Actual Laurent coefficients c_i-average(c) and the derivative integrality lemma | Valid with X=Z^n/Z(1,...,1); source polynomial generation uses the existing theorem |
| Projective Levi invariants are the elementary functions of differences | The differences form a Z-basis of Q and the Levi permutes that basis | Full invariant ring, not only rational generators |
| Ordinary factors are n/gcd(n,d) | Shift identity gives binomial generator coefficients; lower and attainable upper p-valuations prove the monomial gcd | Valid for 0<=d<n on SL_n/P_1; higher target degrees vanish |
| SL4 equivariant source formulas have the stated normalization | Projective bundle relation and all fixed-point restrictions | Valid with the coefficient ring/sign convention explicitly fixed |
| The first equivariant specialization defect is Z/2 in degree 5 | Complete low-degree product reductions, a primitive odd coefficient in a^2, and J_5=M_5+Z nu | Valid; a full description of higher-degree kernels is not asserted |

The elementary valuation assertion can be checked without a separate
Lucas-theorem hypothesis: when 1<=r<p^a and p^a divides n, every factor
(n-i)/i for 1<=i<=r-1 has valuation zero. This proves
v_p binomial(n-1,r-1)=0. The same argument applies at r=p^a.

## Independent calculation coverage

The [Sage report](experiments/local-image-equivariant-pilot.json) records:

- Top characteristic-index polynomials for A1, A2, A3, B2, B3, C3 and G2,
  including 427 parabolic Schubert lifts.
- BGG checks of all elementary-difference generators for SL_n/P_1,
  2<=n<=8, and source-lattice indices/exponents in those cases.
- The gcd formula in 3160 positive-degree rows for 2<=n<=80.
- Twelve fixed-point checks and the displayed SL4 multiplication identities.
- Exact integer image and parameter-times-image lattices in degrees 0--5.
  Degree 5 has ambient rank 52 and a single nonunit quotient factor 2.

The [independent report](experiments/local-image-independent.json) uses
SymPy 1.12 and Python stdlib. It independently constructs the quotient
polynomials and the two rank-52 degree-five lattices, obtains their
Hermite bases, and Smith-reduces the inclusion. Its result is 51 unit
factors and one factor 2. It checks the witness membership, failure of
membership in parameter times image, and membership after multiplication
by two. It also rechecks all 3160 ordinary-family coefficient gcd rows.

Both implementations use the mathematical projective-bundle presentation;
their agreement independently checks its arithmetic, not the published
geometric theorem behind the presentation. The written arguments identify
the source and prove the general statements independently of the finite grid.

Reproduce from borel/:

```sh
PYTHONPATH=src sage -python experiments/local_image_equivariant_pilot.py
python experiments/verify_local_image_independent.py
```

The second command requires SymPy. No general equivariant runtime API,
new wheel or runtime version change is supplied by this work package.

## Literature attribution and search limits

The normalized weight formula is Rogers's published input, and the
polynomial-character/Borel construction is Kumar's. The degree-one
source theorem is the earlier repository argument.

The torsion-index polynomial and divided-difference characteristic-map
machinery are classical. See
[Demazure, Inventiones 21 (1973), 287--301](https://doi.org/10.1007/BF01418790)
and [Calmès–Petrov–Zainoulline, §5.2 and Theorem 6.4](https://bcalmes.perso.math.cnrs.fr/pro/recherche/publis/cobordismflag/cobordismflag.pdf).
The all-Levi bound is best described as an explicit application of that
machinery and a lattice-exponent inclusion, not a newly discovered
good-prime surjectivity theorem.

The projective and symmetric-polynomial presentations are standard.
The root-characteristic interpretation for the adjoint root datum is
an important literature comparison for the n/gcd(n,d) image formula;
we have not completed that comparison. It should not be advertised as
an original type A theorem on the strength of a targeted search.

Equivariant characteristic maps and fixed-point methods are established;
see [Calmès–Zainoulline–Zhong](https://arxiv.org/abs/1409.7111).
Integral invariants versus divided-difference conditions are also
treated by [Holm–Sjamaar](https://arxiv.org/abs/math/0607069).
The Tor/base-change mechanism is general commutative algebra.
The counterexample is the explicitly identified literal SL4 source
inside that framework, not a new general theory of equivariant cohomology.

Targeted searches on 20 October 2026 concerned torsion-index annihilators,
root-lattice/projective characteristic images and equivariant
characteristic-subring specialization. They did not identify a precise
prior statement of this projective image formula or this exact SL4
specialization witness. The search was limited and does not establish
their absence from the literature or mathematical priority.

## Remaining program gaps

The conservative ordinary/equivariant good-prime benchmark and first E1
pilot are complete in the forms above. OPEN-02/09 still ask for sharper
bad-prime structure and uniform invariant/image descriptions. OPEN-14
now asks for conditions eliminating or controlling the bad-prime
specialization kernel, since unconditional faithfulness has a
counterexample. Shared-prime representation comparisons outside the
projective family remain OPEN-05.
