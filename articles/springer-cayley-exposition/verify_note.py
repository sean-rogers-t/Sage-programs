"""Artifact QA and source/status provenance; not a new mathematical proof audit."""
from hashlib import sha256
from pathlib import Path
import json
import re
import subprocess

import pdfplumber
from PIL import Image, ImageDraw
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
BOREL = REPO / 'borel'
PDF = ROOT / 'tmp/pdfs/main.pdf'


def digest(path):
    data = path.read_bytes()
    return sha256(data if path.suffix == '.pdf' else data.replace(b'\r\n', b'\n')).hexdigest()


def main():
    source = ROOT / 'main.tex'
    tex = source.read_text(encoding='utf-8')
    log = (ROOT / 'tmp/pdfs/main.log').read_text(encoding='utf-8', errors='replace')
    assert not re.search(r'Overfull \\[hv]box|undefined references|Citation .* undefined', log)
    citations = {key for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', tex)
                 for key in group.split(',')}
    assert citations <= set(re.findall(r'\\bibitem\{([^}]+)\}', tex))
    reader = PdfReader(PDF)
    assert str(reader.trailer['/Root'].get('/PageLayout')) == '/SinglePage'
    words = []
    with pdfplumber.open(PDF) as pdf:
        for page in pdf.pages:
            assert (page.width, page.height) == (612, 792)
            text = page.extract_text() or ''
            assert '\ufffd' not in text and '??' not in text
            words.append(len(text.split()))
            for char in page.chars:
                assert -0.5 <= char['x0'] <= char['x1'] <= page.width + 0.5
                assert -0.5 <= char['top'] <= char['bottom'] <= page.height + 0.5
    assert len(words) == 11 and min(words) > 100, words
    images = sorted((ROOT / 'tmp/pdfs').glob('page-*.png'))
    assert len(images) == len(reader.pages)
    for start in range(0, len(images), 4):
        sheet = Image.new('RGB', (1280, 1720), '#e8edf0')
        draw = ImageDraw.Draw(sheet)
        for offset, filename in enumerate(images[start:start+4]):
            page = Image.open(filename).convert('RGB')
            page.thumbnail((620, 820))
            x, y = 640 * (offset % 2), 860 * (offset // 2)
            draw.text((x+10, y+5), f'Page {start+offset+1}', fill='black')
            sheet.paste(page, (x+10, y+27))
        sheet.save(ROOT / f'tmp/pdfs/contact-{start//4+1}.png')

    records = sorted(set(re.findall(r'\\record\{([^}]+)\}', tex)) - {'borel'})
    for name in records:
        assert (BOREL / name).is_file(), name
    current = json.loads((BOREL / 'binary-branching-verification.json').read_text())
    audit = json.loads((BOREL / 'binary-branching-audit.json').read_text())
    assert current['result'] == audit['result'] == 'passed'
    assert current['resolved_proper_standard_Levis'] == audit['resolved_proper_standard_Levis'] == 123
    assert current['remaining_embeddings'] == audit['remaining_embeddings'] == 4
    assert current['remaining_Weyl_classes'] == 2
    assert current['regular_fiber_dimensions'] == {'A1^3 cubic class': 6, 'A2+A1^2 cubic class': 10}
    assert not audit['companion_memberships_verified']
    assert not audit['all_simple_adjoint_conjecture_proved']
    for record in (current, audit):
        for name, expected in record['source_sha256'].items():
            assert digest(BOREL / name) == expected, name

    files = [BOREL / name for name in records]
    files += [source, Path(__file__), REPO / 'articles/PROVENANCE.json']
    manifest = dict(date='2026-10-24',
        base_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=REPO, text=True).strip(),
        includes_uncommitted_sources=True,
        scope='Dated synthesis of repository proofs and status; no fresh all-proof or exhaustive literature audit.',
        source_sha256={p.relative_to(REPO).as_posix(): digest(p) for p in sorted(set(files))},
        primary_sources_consulted=[
            'https://www.mat.univie.ac.at/~michor/kostant.pdf',
            'https://arxiv.org/abs/1701.01538',
            'https://www.numdam.org/item/PMIHES_1965__25__49_0/'])
    (ROOT / 'source-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    final = ROOT / 'output/pdf/Springer_Cayley_maps_progress_and_open_questions.pdf'
    final.write_bytes(PDF.read_bytes())
    report = dict(result='passed',date='2026-10-24',pages=len(words),words_by_page=words,
        total_words=sum(words),single_manuscript_page_per_pdf_page=True,
        page_layout='SinglePage',unresolved_references=False,overfull_boxes=False,
        source_records_checked=len(records),E7_resolved=123,E7_unresolved=4,
        E7_unresolved_Weyl_families=2,companion_memberships_verified=False,
        all_simple_adjoint_conjecture_proved=False,
        scope='PDF/text/layout and citation checks, source provenance, and current E7 status consistency; not a new proof verification.',
        source_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in (
            source,Path(__file__),ROOT/'source-manifest.json',final)})
    (ROOT / 'verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))


if __name__ == '__main__':
    main()
