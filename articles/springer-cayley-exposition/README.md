# Springer-Cayley maps: exposition and unresolved questions

An eleven-page research exposition dated 24 October 2026. Each PDF page
contains one full-size manuscript page. The note distinguishes published
inputs, internally audited repository deductions, exact computational
corroboration, and conjectures. It makes no literature-priority claim and
does not constitute a fresh audit of every historical proof.

Subsequent direction change: the user stopped automatic Levi-case completion.
The PDF's final recommendation is now superseded by the
[direct uniform proof target](../../borel/UNIFORM_ADJOINT_PROOF_TARGET.md).
The mathematical results and open-case counts have not changed.

- [Read the PDF](output/pdf/Springer_Cayley_maps_progress_and_open_questions.pdf)
- [Editable LaTeX](main.tex)
- [Exact source-version manifest](source-manifest.json)
- [Artifact and status verification](verification.json)

The exposition defines the trace-normalized map, explains Rogers's torus
formula and the Kumar polynomial-character/Borel composite, and separates
regularity, harmonic separation, generic degrees, subgroup coordinates,
and fiber geometry. It gives the character-critical presentation, the
classical adjoint theorem, exceptional progress, and a proof of the common
binary critical-equation criterion.

The immediate E7 gaps are the cubic A2^2+A1 family (three standard subsets,
seed 12356) and A2+A1^3 (one, seed 12357). Coverage is 123/127. The latter's
six ideal memberships remain unverified. The general simple-adjoint
generic-regular-fiber smoothness conjecture remains open, as do the
separate nonregular-fiber, exceptional-target, degree-boundary, and minimal
harmonic-localization questions described in the note.

The [proof-status guide](../../borel/PROOF_STATUS.md) and
[research ledger](../../borel/RESEARCH_LEDGER.md) remain authoritative for
the evolving program. This exposition is a dated snapshot.

## Rebuild and verification

Compile `main.tex` with Tectonic to `tmp/pdfs/main.pdf`, preserving its log.
Render each page with Poppler at 85 dpi into `tmp/pdfs/page-*.png`, then run
`python verify_note.py` with pypdf, pdfplumber, and Pillow installed. The
checker copies the checked PDF to `output/pdf`, refreshes the source hashes,
checks current E7 status against both independent reports, and creates
temporary contact sheets for visual review. It does not rerun the historical
mathematical computations. The source also opens in the desktop LaTeX editor.
