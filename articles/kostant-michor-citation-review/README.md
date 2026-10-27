# Kostant–Michor citation review

Prepared 27 October 2026 for Sean Rogers.

**Roadmap correction and subsequent progress:** the ordinary classical adjoint nonemptiness and dimension classification was already completed at GEO-49. The unfinished extension concerned polynomial-adjoint representations. The subsequent [GEO-97 theorem](../../borel/CLASSICAL_POLYNOMIAL_ADJOINT_BLOCKS.md) now supplies that extension, with its own independent checks and explicit review limitations. The report describes the preceding literature-review checkpoint; its roadmap wording has been corrected accordingly.

- [Report PDF](output/pdf/Kostant_Michor_Citation_Review.pdf)
- [LaTeX main document](main.tex), with [framework](framework.tex), [paper profiles](profiles.tex), [research assessment](research.tex), and [references](references.tex)
- [Source manifest and research-note snapshot](source-manifest.json)
- [Exact comparison checks](check_comparisons.py) and [results](comparison-checks.json)
- [PDF verification](report-verification.json)

The report covers twelve identified citing works, the original Kostant–Michor paper, and two graph-theory follow-ups. It distinguishes substantive map uses, local-coordinate uses, inequality uses, and a bibliography-only citation. Coverage is complete for this identified corpus, not asserted exhaustive worldwide.

It compares the published results with the current repository proofs, explains what does and does not follow, and recommends a structural classical nonemptiness/excess classification after external review of the adjoint fiber theorem. It also records two precisely located dissertation corrections. Historical source files and existing proof certificates are unchanged.

Compile `main.tex` in the Codex LaTeX editor or with Tectonic. The delivered PDF uses US Letter pages and SinglePage viewer layout. The exported file is named `output/pdf/Kostant_Michor_Citation_Review.pdf`.

Run `python check_comparisons.py` with SymPy to reproduce the exact matrix comparisons. These checks do not rerun the E8 certificate or establish publication priority. Public source PDFs are stored outside this report directory and are identified by URLs and hashes in the manifest.

Current mathematical status is described from [PROOF_STATUS.md](../../borel/PROOF_STATUS.md), [the adjoint theorem](../../borel/ADJOINT_GENERIC_FIBER_THEOREM.md), [critical values](../../borel/ADJOINT_NONEMPTINESS_CRITICAL_VALUES.md), [character compositions](../../borel/CHARACTER_COMPOSITION_FIBERS.md), [the regularity classification](../../borel/REGULARITY_CLASSIFICATION.md), and [harmonic localization](../../borel/SYMPLECTIC_HARMONIC_LOCALIZATION.md). These are internally audited repository arguments, not externally refereed results. The manifest binds the exact note contents read, including uncommitted work.
