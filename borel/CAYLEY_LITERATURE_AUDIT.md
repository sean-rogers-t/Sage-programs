# Literature audit of the Springer and Cayley research brief

Audit date: 22 October 2026. This is a claim-specific audit of all eight
bibliographic entries in the supplied Springer_Cayley_Research_Brief.docx,
its linked notebook, and a targeted search of related and citing work.
It checks the proposed research boundary. It is not a line-by-line
verification of every proof or printed exceptional table in those papers,
an external review, or a certificate of originality.

The earlier Smith project is paused. No PDF, LaTeX article or source paper
is changed by this audit. The original brief is preserved.

## Conclusions

The identification of the two trace-projection constructions is correct.
The proposed invariant pullbacks, differential, field-degree comparison
and harmonic-localization criterion have established foundations. The
brief generally distinguishes them correctly from proposed research.
Its starting calculations overlap the repository and, in rank one, the
source literature itself. A useful next target is a representation-level
classification of regularity preservation, with its priority still
unverified. A fresh mathematical proof of such a classification is not
being asserted here.

Subsequent work on the same date completes the first pilot family: the
[exterior-power note](EXTERIOR_REGULARITY.md) proves that the maps for
exterior^k(C^n) preserve torus regularity, equivalently global regularity,
exactly for k=1,n-1. This is a separate written deduction with a uniform
proof, not a novelty certification supplied by this literature audit.
The follow-on [type A proof](TYPE_A_REGULARITY.md) extends the deduction
to every nontrivial irreducible SL_n highest weight, with exactly the same
two positive cases. The source and originality qualifications still apply.

Two additions materially improve the bibliography: the original
Lemire-Popov-Reichstein *Cayley groups* paper, and Agricola-Goodman's
invariant-vector-field and harmonic-separation paper. In particular,
the latter already supplies a general basis of equivariant regular maps
and an explicit differential-operator harmonic complement for SL2.

There are also concrete cautions: Part I has an apparent scalar
normalization error in its adjoint SL2 example, and the current notebook
is a mutable, partially configured working file. These observations are
separate from the validity of the main published theorems.

## 1. Coverage of each reference in the brief

PDF page numbers below refer to the named version, not necessarily the
printed journal page. The complete source files were obtained and their
relevant constructions, theorem statements, proofs and references inspected.
Unrelated proofs elsewhere in a long paper were not independently re-proved.

| Brief entry | Version and locators checked | What is established and what does not follow |
| --- | --- | --- |
| KM | [Kostant-Michor, arXiv v3](https://arxiv.org/pdf/math/0109066), 33 pages; Sections 2.1, 2.4-2.5, 3.2-3.3, 4.5, 4.8, 4.11, 4.13, 5, 6 and 7 | Trace projection; identity differential; centralizer inclusion and equality when the differential is invertible; localized harmonic separation and its regular-orbit criterion; quotient-degree equality; Jordan behavior; explicit examples. These are published inputs, not new consequences of our code. |
| R | [Rogers, arXiv v1](https://arxiv.org/pdf/1701.01538v1), 10 pages; Section 2.1, Theorems 1-2 and their Sections 3-4 proofs | Explicit weight-moment coordinate formula and symmetrized Cartan Gram matrix. This is the preprint, not a newly obtained final journal copy. Multiplicities are essential. The word invariant in some prose is read as equivariant, consistently with the construction. |
| K | [Kumar, published Part I](https://kumar.math.unc.edu/papers/kumar76.pdf), 21 pages; Definition 1, Lemma 2, Example 3, Definition 4, Theorem 5 and Sections 5-10 | Definition of polynomial characters, inverse-on-image construction and complex finite-stage surjectivity, with classical source descriptions and maps. Neither a uniform integral image theorem nor general regularity preservation is asserted by Theorem 5. See the Example 3 normalization issue below. |
| KR | [Kumar-Rogers, published Part II](https://kumar.math.unc.edu/papers/kumar81.pdf), 23 pages; Definitions 1 and 3, Theorem 4, Propositions 7, 11, 20 and 24, Theorems 16 and 29, Remarks 17 and 30 | Explicit defining C/B maps, selected integral forms, and injectivity from the fixed classical factor into stable cohomology. The finite map from the full Levi source is surjective; the restricted stable map is not. This distinction is correctly made in the brief. |
| KX | [Kumar-Xie, published Part III](https://kumar.math.unc.edu/papers/kumar86.pdf), 30 pages; introduction, Theorems 3.6 and 3.8, Remark 3.9, Sections 4-8, reference 16 | Type D stable injectivity with its Euler sector, and minimal-representation exceptional maximal-parabolic computations. E8 is explicitly omitted because the expressions are too long; the method is said to apply in principle. E8 implementation alone would extend a computation, without automatically supplying a new theorem or method. |
| LPR | [Lemire-Popov-Reichstein, degree paper](https://arxiv.org/pdf/math/0608473), 8 pages; all six sections | Minimum Cayley degree, isogeny/product bounds, degree-preserving reduction to W-equivariant torus maps, the SLn upper bound and Cay(G2)=2. Section 6 also discusses the character-difference construction giving a degree-at-most-six G2 map. These minimum-degree statements must not be confused with our specified trace maps. |
| P | [Popov, 2007 problem list](https://arxiv.org/pdf/math/0702533), 3 pages; Section 2 and Problems 2.1-2.4, with the other sections screened for scope | Exact historical source for the SL5, degree-greater-than-two and unboundedness questions. It is not evidence by itself that these remain open in 2026. |
| LY | [Lu-Ye, v2 of 18 July 2025](https://arxiv.org/pdf/2411.02071v2), 22 pages; introduction, Definition 4.1, Theorems 4.9 and 5.6, Section 6 and Theorems 6.13/6.15 | Applicability of the rational matrix transform to a represented group/image, characterized through odd powers and weight configurations. This is a different property from regularity preservation of trace projection. Theorem 6.13 is stated for classical complex simple groups; exceptional exclusion is discussed at the start of Section 6. Finite kernels and the represented image must be retained. |

The journal versions of K, KR and KX were used, rather than transferring
their claims solely from abstracts or preprints. R's inaccessible final
journal copy remains outside this audit; the user's earlier preference
not to pursue that access issue is retained.

## 2. The common map and the exact literature boundary

Assume G is a complex connected semisimple group and rho has faithful
differential. Put

\[
 B_\rho(X,Y)=\operatorname{tr}(d\rho(X)d\rho(Y)),\qquad
 B_\rho(\theta_\rho(g),X)=\operatorname{tr}(\rho(g)d\rho(X)).
\]

This is the construction in KM Section 2.1 and K Definition 1, recalled
in KR Definition 1 and R Section 2.1. For a rational representation of
a complex reductive group the needed trace nondegeneracy is established
in KM Section 2.2; arbitrary real representations require the additional
hypotheses recorded there. No general real statement is silently inferred.

With fundamental-weight coordinates for weights and simple-coroot
coordinates for the output, R's formula is

\[
 S=\sum_\mu m_\mu\mu\mu^t,\quad
 \theta_\rho(t)=S^{-1}\sum_\mu m_\mu\mu e^\mu(t),\quad
 J_T(t)=S^{-1}\sum_\mu m_\mu\mu\mu^t e^\mu(t).
\]

The last formula is elementary differentiation, also compatible with
KM Proposition 2.4(7). It is not a new differential theorem. The full
root-direction determinant in the [repository geometry note](SPRINGER_GEOMETRY.md)
follows from equivariance; a source critical divisor is still distinct
from its target branch locus and from the nonproperness locus.

K Definition 4 and Theorem 5 give

\[
 R_\rho(L)=\theta_T^*\mathbb C[\mathfrak t]^{W_L},\qquad
 \xi_\rho^P=\beta^P\circ(\theta_T^*)^{-1}.
\]

The inverse is an isomorphism from the image algebra, not a geometric
inverse to a generally many-to-one map. For every Levi, the brief's
fraction-field index is deg(theta_T): cancel the two faithful W_L
extensions of degree |W_L|. This is an elementary finite-group field
argument of the kind in LPR Section 3, not a new Levi degree theorem.
Generic finite field degree does not imply a finite or free extension
of coordinate rings.

The familiar projective and Grassmannian results concern beta after
taking invariants. The geometry of theta before beta contains information
that beta kills, so the integral Borel-image and generic-degree problems
remain different questions.

## 3. Important related work missing from the brief

### Original Cayley classification

[Lemire-Popov-Reichstein, Cayley groups](https://arxiv.org/abs/math/0409004),
JAMS 19 (2006), 921-967, is the foundational classification cited by the
brief's shorter degree paper. Its torus/birational reduction, character
lattice obstruction and Section 10 regular-map construction were inspected.
Its Example 10.5 explicitly includes the trace-projection construction.
One must distinguish the character lattice used there to study a W-torus
from the representation-dependent literal source lattice in our integral
image calculations. Having lattices in both settings does not identify
the two classification problems.

### Invariant vector fields and an explicit harmonic complement

[Agricola-Goodman](https://arxiv.org/pdf/math/0207161v2), Michigan Math. J.
51 (2003), 607-630, is directly relevant and should be added to the program.
Theorem 4.1 says that, for simply connected semisimple G, the gradients of
fundamental characters form a free basis over C[G]^G of the module of
conjugation-invariant algebraic vector fields. Under left trivialization,
these are exactly equivariant regular maps G -> g. Therefore merely
proving such a module basis or character-gradient description would
repeat a known theorem.

For simple G, fix a common invariant form B and write B_rho=I_rho B.
If chi_rho=F(chi_1,...,chi_r), ordinary differentiation gives

\[
 I_\rho\theta_\rho
 =\sum_i I_i\frac{\partial F}{\partial\chi_i}\theta_i.
\]

This is the standard chain rule applied to known gradients. The new
question would concern regularity, critical factors or degrees for the
maps it expresses, rather than the existence of the expression.

AG Theorem 5.3 constructs an explicit invariant differential operator
and harmonic separation for SL2; Theorem 5.2 handles its Weyl-torus
action. AG's introduction expressly says its proof does not extend
obviously to arbitrary conjugation actions. The known existence of some
harmonic complement, KM's particular pullback harmonic space, and an
explicit differential-operator kernel are three different assertions.

### Later birational work screened

[Borovoi-Kunyavskii-Lemire-Reichstein](https://arxiv.org/abs/1207.1329v5),
IMRN 2014, 5340-5397, Theorems 1.3-1.4, gives the quasi-permutation
character-lattice criterion and classification for stable Cayleyness in
characteristic zero. [Borovoi, with Dolgachev's appendix](https://arxiv.org/abs/1212.1065v4),
J. Algebra 436 (2015), 35-60, Theorems 1.2-1.3, concerns real groups of
absolute rank at most two and stabilization. These address birational
or stable questions over specified fields, not the trace-map regularity
classification proposed here.

The citing literature also led to
[Hassett-Tschinkel's low-dimensional quadrics paper](https://arxiv.org/abs/2411.00226)
(author manuscript dated 23 May 2025), which studies stable linearization
of finite-group actions, and
[Popov's Jordan-groups survey](https://arxiv.org/abs/1307.5522).
Their introductions and relevant Cayley passages were screened; no
resolution of the three minimum-degree questions or of trace-map
regularity was identified there. No claim is made to have audited all
proofs in these peripheral papers.

## 4. Two normalization and reproducibility findings

### Part I Example 3

On printed page 398 (PDF page 4), the published adjoint SL2 example
has diagonal coefficient (z^2-z^(-2))/2. With the paper's trace-projection
definition and simple coroot h=diag(1,-1), the weights are 2,0,-2, so

\[
 S=8,\qquad
 \theta_{2\omega_1}(\operatorname{diag}(z,z^{-1}))
 =\frac{z^2-z^{-2}}4\operatorname{diag}(1,-1).
\]

A direct matrix projection with d(rho)(h)=diag(2,0,-2) gives the same
coefficient. Its logarithmic derivative at 1 is one; the printed
coefficient gives two. KM Section 6.1's general SL2 formula also gives
the denominator four for this representation. The PDF page was visually
inspected. This is a proposed scalar correction, not an author-approved
erratum. Over C, nonzero scalar changes do not change the generated
polynomial-character subalgebra or generic degree; an integral source
calculation must retain the correct normalization.

Part III's E6 scalar issue was already recorded in the
[earlier article review](../articles/README.md#verification-note-for-shrawan);
it is not a new discovery of this audit. No exceptional table was
silently corrected or recertified here.

### The public Xie notebook

The [notebook in KX reference 16](https://colab.research.google.com/drive/1WdrZseCHoRlWADOng7SwfefXqgbmz3CL)
was retrieved through its public download endpoint and all 30 cells
inspected. No notebook cell, Drive mount or embedded command was executed.
The exact downloaded-file SHA256 and literal matrix fixtures are in the
[source record](experiments/cayley-literature-source-record.json).

Using zero-based cell numbers:

- Cell 3 sets an E8 Cartan matrix A and an SInv matrix. SInv's last row
  duplicates its second row. It has rank seven, and is neither symmetric
  nor the inverse of A, whose determinant is one. Even treating a missing
  scalar separately cannot make a singular matrix an inverse. The correct
  last row of A^(-1) is (2,3,4,6,5,4,3,2).
- Cell 8 explicitly defers multiplication by the trace factor 2/x.
  Consequently the displayed SInv is intended as the Cartan inverse up
  to the separately tracked scale. That convention does not resolve the
  duplicate-row defect.
- Cell 16 requires a separate mounted file E8MCR.txt. The notebook itself
  does not contain that file.
- Cell 20 assigns seven-column E7-style vectors while A remains 8-by-8
  in a linear read of the saved file. The matrix product has incompatible
  dimensions without manual reconfiguration.
- All saved execution counts are null, despite retained outputs. This
  is not evidence of a fresh linear execution of the saved settings.

These are reproducibility defects of the retrieved working file. They
do not identify which historical settings produced a published table,
prove an error in every table, or establish that an E8 computation was
completed. The independently rebuilt repository calculations remain the
appropriate reference for our checked formulas.

## 5. What the regularity proposal can legitimately claim

KM Theorem 3.2 already proves

\[
 \mathbb C[G]=\mathbb C[G]^G\otimes\theta_\rho^*\operatorname{Harm}(\mathfrak g)
 \quad\Longleftrightarrow\quad
 \theta_\rho(G_{\rm reg})\subseteq\mathfrak g_{\rm reg}.
\]

A classification of the rho satisfying the right side was not located
in the inspected sources and searches. This is a candidate research
target with unverified priority, rather than a documented open problem
explicitly posed in KM. KM's original question is an explicit harmonic
construction, and its criterion is already a theorem.

For regular semisimple t, an elementary root test is

\[
 e^\beta(t)\ne1\ (\beta\in\Phi),\qquad
 \theta_\rho(t)\text{ regular}\iff
 \alpha(\theta_\rho(t))\ne0\ (\alpha\in\Phi).
\]

This gives counterexamples and a first necessary condition. A positive
classification still has to handle regular nonsemisimple elements;
KM Section 4's centralizer and Jordan results are relevant. An invertible
torus derivative alone is not sufficient.

The brief's Sp4 witness t=(2,-1/2) passes this audit: its positive root
characters are -4,-1,4,1/4, while both trace coordinates are 3/4. The
torus derivative has determinant -25/16, but the image is nonregular.

There is an equally small way to separate LY's classification from this
question. For SL4 acting on exterior^2(C^4), LY Theorem 6.13(iii) says
the matrix Cayley transform is applicable to the represented image.
The torus element diag(2,-1/2,1,-1) is regular semisimple, yet direct
trace projection gives diag(0,0,3/4,-3/4), which is nonregular. The
representation has central kernel {1,-1}. Our
[independent exact check](tools/verify_cayley_literature.py) uses the
six representation eigenvalues and the simple-coroot trace matrix.
This counterexample is an elementary deduction illustrating the distinction,
not a claimed new classification theorem.

For a first uniform family, exterior powers of SLn remain a reasonable
pilot. The standard and dual endpoint representations preserve centralizer
dimensions; interior exterior powers have additional elementary-symmetric
root factors. The subsequent [uniform proof](EXTERIOR_REGULARITY.md)
settles the complete exterior-power statement, including the endpoint
nonsemisimple cases. Its literature priority remains unverified. The later
[positive-character argument](TYPE_A_REGULARITY.md) settles arbitrary
irreducible SL_n highest weights. Other groups and reducible
representations remain follow-up questions.

## 6. Present status of the historical degree questions

| Question | Verified historical statement | Status after this audit |
| --- | --- | --- |
| Cay(SL5) | Popov Problem 2.2; LPR proves 2 <= Cay(SL5) <= 3 | No primary resolution located. Current openness remains unverified. |
| Some G, or some simple G, with Cay(G)>2 | Popov Problem 2.3 | No primary resolution located. The degree of a particular representation does not answer it. |
| Unbounded minimum Cayley degrees | Popov Problem 2.4, both general and simple versions | No primary resolution located. An unbounded family of trace-map degrees would not prove it. |

Searches included the exact problem phrases, SL5/SL_5/SL(5) variants,
regular-orbit and representation formulations, the authors' publication
lists, OpenAlex citing-work discovery for the eight entries and the
foundational Cayley-groups paper, and a Semantic Scholar check of the
degree preprint. The citing lists are incomplete: LY cites KM and Cayley
groups but was absent from both returned OpenAlex lists. A zero count or
an unsuccessful phrase search therefore certifies neither completeness
nor openness. Graph-embedding papers using the same term Cayley map were
excluded as unrelated. Exact source versions, hashes and discovery limits
are retained in the source record.

## 7. Updated recommendation

Use the published gradient basis and KM regularity criterion as starting
theorems. Seek a uniform classification for a specified representation
family, beginning with exterior powers of SLn, and retain the existing
G2/E6/E7 geometry as benchmarks. Do not restart the already published
all-SL2 coordinate calculation, the known SLn defining degree, a general
equivariant-gradient basis theorem, or an SL2 harmonic-complement existence
project under a new name. Minimum Cayley degrees and E8 table completion
remain separate projects, with different novelty and difficulty questions.

The [audit check report](cayley-literature-verification.json) verifies the
source-record consistency, notebook matrix defects, the Part I scalar
discrepancy and the two regularity witnesses. It does not mechanically
prove a literature-status claim. No new runtime API is added in this audit.
