"""Editorial/provenance checks; this does not independently prove the mathematics."""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import fitz

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
BASELINE = '9f6c000a4a97a5f7540424ba1f8f2881cc55a2ae'

CLAIMS = [
    ('Published maps', 'published', 'MATHEMATICS.md', 'Rogers formula and Kumar polynomial-character construction; original papers cited separately.'),
    ('Full torus and Levi source', 'repository-proof', 'TRACE_SOURCE_THEOREM.md', 'Simple simply connected, faithful differential; torus generation is not Levi linear generation.'),
    ('Irreducible source classification', 'repository-proof', 'IRREDUCIBLE_CHARACTER_LATTICES.md', 'All finite simple types; direct sums excluded from the closed formula.'),
    ('Full-source exponent bound', 'repository-deduction-from-published-inputs', 'LOCAL_IMAGE_BOUND.md', 'Literal full Levi source; torsion-index machinery is published.'),
    ('Normalized B factors', 'repository-proof', 'TYPE_B_GENERAL_SMITH.md', 'All allowed single-tail blocks; index in selected saturation.'),
    ('D square row equality', 'repository-proof', 'TYPE_D_SQUARE_THEOREM.md', 'Normalized square sector, all allowed blocks.'),
    ('Full D factors', 'repository-proof', 'TYPE_D_EVEN_GLUING_THEOREM.md', 'Euler gluing included; literal source requires separate scaling.'),
    ('Adapted B/D classes', 'repository-proof', 'CLASSICAL_ADAPTED_SMITH_BASIS.md', 'Constructive deduction from primitive prefixes, with finite independent controls.'),
    ('Full D basis', 'repository-proof', 'TYPE_D_EXPLICIT_INTEGRAL_BASIS.md', 'Unit gluing and multiplication retained.'),
    ('Denominator formula', 'repository-proof-and-conjecture', 'CLASSICAL_SATURATION_GENERATORS.md', 'Formula, inclusion and top layer proved; intermediate higher-block equality conjectural.'),
    ('Two-variable generators', 'repository-proof', 'CLASSICAL_TWO_HEAD_GENERATORS.md', 'All larger blocks; multiplication rules are not a complete defining ideal.'),
    ('C joint images', 'repository-deduction-from-published-inputs', 'TYPE_C_JOINT_PRIMITIVITY.md', 'Specified normalized two-step generators.'),
    ('B/D joint images', 'repository-proof-with-open-ambient-problem', 'ORTHOGONAL_JOINT_IMAGES.md', 'Relative sum primitive; full ambient factors open.'),
    ('E6 relative images', 'computer-assisted-repository-proof', 'E6_REMAINING_CASES.md', 'Together with E6_ROOT_WEIGHT_COMPARISON: all 63 proper standard parabolics, not ambient quotients.'),
    ('E7 relative images', 'computer-assisted-repository-proof', 'E7_ALL_PARABOLICS.md', 'All 127 proper standard parabolics; not ambient quotients.'),
    ('G2 full flag', 'exact-finite-calculation-and-kernel-proof', 'G2_FULL_FLAG.md', 'Absolute literal minimal/adjoint full-flag image.'),
    ('Filtered K comparison', 'repository-proof', 'OTHER_THEORIES_FINITE.md', 'Full literal source; support filtration, not gamma filtration.'),
    ('Symplectic ordinary K factors', 'repository-deduction-from-published-inputs', 'SYMPLECTIC_K_ORDINARY_SMITH.md', 'Projective defining family; nonequivariant K-theory.'),
    ('Specialization exponent', 'repository-proof', 'TYPE_A_K_INTEGRAL_RECURRENCE.md', 'All m; full multiplicative jet source, not just its linear span.'),
    ('Stable partition factors', 'repository-proof-in-bounded-degree-range', 'TYPE_A_K_STABLE_MULTIPLICITIES.md', 'd<=p, also p+1 for odd p; transported grading.'),
    ('First correction', 'repository-proof-in-bounded-degree-range', 'TYPE_A_K_FIRST_CORRECTION.md', 'p+2 for odd p; root boundary distinct.'),
    ('Next elementary degree', 'repository-proof-in-bounded-degree-range', 'TYPE_A_K_NEXT_ELEMENTARY.md', 'p+3 for p>=5; full first Frobenius band remains open.'),
    ('Quantum images', 'repository-deduction-from-published-inputs', 'OTHER_THEORIES_QUANTUM.md', 'Specified projective lifts, not all partial flags.'),
    ('Affine images', 'repository-deduction-from-published-inputs', 'OTHER_THEORIES_AFFINE.md', 'Evaluation-and-level cup-cohomology construction, distinct from Peterson transport.'),
    ('Exceptional degrees', 'computer-assisted-repository-proof', 'EXCEPTIONAL_DEGREES.md', 'Specific normalized maps and correct torus lattices; not minimum Cayley degree.'),
    ('Regularity classification', 'repository-proof', 'REGULARITY_CLASSIFICATION.md', 'Actual representations of simple groups; external review still needed.'),
    ('Symplectic harmonic localization', 'repository-deduction-from-published-inputs', 'SYMPLECTIC_HARMONIC_LOCALIZATION.md', 'Fixed pulled harmonic space and invariant principal localizations.'),
    ('Adjoint generic regular fibers', 'computer-assisted-repository-proof', 'ADJOINT_GENERIC_FIBER_THEOREM.md', 'Actual schemes, generic centralizer-exact targets; emptiness and positive excess allowed.'),
    ('E8 certificate audit', 'internal-audit-with-explicit-engine-dependency', 'E8_DEGREE_TRANSFER_AUDIT.md', 'Native checks do not independently reproduce every membership-engine result.'),
    ('Classical polynomial-adjoint blocks', 'repository-proof', 'CLASSICAL_POLYNOMIAL_ADJOINT_BLOCKS.md', 'Actual polynomial characters, generic proper-Levi targets; not arbitrary representations.'),
]

def digest(path):
    data = path.read_bytes()
    if path.suffix not in {'.pdf', '.gz'}:
        data = data.replace(b'\r\n', b'\n')
    return hashlib.sha256(data).hexdigest()

def main():
    sources = sorted(HERE.glob('*.tex'))
    text = '\n'.join(p.read_text(encoding='utf-8') for p in sources)
    refs = sorted(set(re.findall(r'\\record\{([^}]+)\}', text)))
    for ref in refs:
        assert (REPO / 'borel' / ref).is_file(), ref
    for inc in re.findall(r'\\input\{([^}]+)\}', text):
        assert (HERE / (inc + '.tex')).is_file(), inc
    cite_keys = set()
    for group in re.findall(r'\\cite(?:\[[^]]*\])?\{([^}]+)\}', text):
        cite_keys.update(group.split(','))
    bib_keys = set(re.findall(r'\\bibitem\{([^}]+)\}', text))
    assert cite_keys <= bib_keys, cite_keys - bib_keys
    assert BASELINE in text
    log = (HERE / 'output/pdf/main.log').read_text(encoding='utf-8', errors='replace')
    for bad in ('Overfull \\hbox', 'Undefined control sequence', 'There were undefined references', 'multiply defined'):
        assert bad not in log, bad
    pdf = HERE / 'output/pdf/Springer_Borel_Cayley_Research_Statement.pdf'
    visual_path = HERE / 'visual-review.json'
    if visual_path.exists():
        visual = json.loads(visual_path.read_text(encoding='utf-8'))
        assert visual['pdf_sha256'] == digest(pdf), 'Refresh render review for the changed PDF'
    document = fitz.open(pdf)
    assert 25 <= len(document) <= 35, len(document)
    all_text = '\n'.join(p.get_text() for p in document)
    assert '\ufffd' not in all_text
    page_words = [len(p.get_text().split()) for p in document]
    assert min(page_words) >= 75, page_words
    inventory = []
    for title, status, source, scope in CLAIMS:
        path = REPO / 'borel' / source
        assert path.is_file(), source
        inventory.append(dict(claim=title, status=status, scope=scope,
                              source='borel/' + source, sha256=digest(path)))
    (HERE / 'claim-inventory.json').write_text(json.dumps(dict(
        date='2026-10-29', baseline=BASELINE, claims=inventory,
        caveat='Status inventory and source cross-check, not external review or fresh independent proof.'
    ), indent=2) + '\n', encoding='utf-8')
    proof_hashes = {ref: digest(REPO / 'borel' / ref) for ref in refs}
    for ref in refs:
        committed = subprocess.check_output(['git', 'show', BASELINE + ':borel/' + ref], cwd=REPO)
        if Path(ref).suffix not in {'.gz', '.pdf'}:
            committed = committed.replace(b'\r\n', b'\n')
        assert hashlib.sha256(committed).hexdigest() == proof_hashes[ref], ref
    result = dict(
        result='passed', date='2026-10-29', baseline=BASELINE,
        pages=len(document), extracted_words=len(all_text.split()),
        page_word_counts=page_words, citation_keys=len(cite_keys),
        local_proof_links=len(refs), claim_inventory_entries=len(inventory),
        pdf_sha256=digest(pdf), source_sha256={p.name:digest(p) for p in sources},
        verifier_sha256=digest(Path(__file__)),
        claim_inventory_sha256=digest(HERE / 'claim-inventory.json'),
        proof_sha256=proof_hashes,
        compiler='Tectonic 0.17.0; desktop compiler also checked',
        compiler_errors=0, overfull_boxes=0,
        underfull_boxes=len(re.findall(r'Underfull \\hbox', log)),
        viewer_layout='SinglePage; one letter page per PDF page',
        visual_review='Render review recorded separately in visual-review.json',
        scope='Editorial, citation-key, pinned-source, and artifact checks. Not a mathematical proof audit or exhaustive novelty search.'
    )
    (HERE / 'verification.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','proof_sha256','page_word_counts')}, indent=2))

if __name__ == '__main__':
    main()
