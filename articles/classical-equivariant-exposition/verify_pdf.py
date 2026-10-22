"""Check this exposition's compiled PDF and prepare page contact sheets.

Run with Python containing pypdf, pdfplumber and Pillow after compiling
to tmp/pdfs/main.pdf and rendering page-*.png with pdftoppm.
This is artifact QA, not an external mathematical proof audit.
"""
from hashlib import sha256
import json
from pathlib import Path
import re

import pdfplumber
from PIL import Image, ImageDraw
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "main.tex"
PDF = ROOT / "tmp/pdfs/main.pdf"


def main():
    tex = SOURCE.read_text(encoding="utf-8")
    reader = PdfReader(PDF)
    assert str(reader.trailer["/Root"].get("/PageLayout")) == "/SinglePage"
    log = (ROOT / "tmp/pdfs/main.log").read_text(encoding="utf-8", errors="replace")
    assert not re.search(r"Overfull \\[hv]box|undefined references|Citation .* undefined", log)
    tags = [int(n) for n in re.findall(r"\\tag\{(\d+)\}", tex)]
    assert tags == list(range(1, 44)), tags
    keys = set(re.findall(r"\\bibitem\{([^}]+)\}", tex))
    cites = {
        k for group in re.findall(r"\\cite(?:\[[^\]]*\])?\{([^}]+)\}", tex)
        for k in group.split(",")
    }
    assert cites <= keys, cites - keys
    notes = re.findall(r"\\record\{([^}]+)\}", tex)
    for relative in notes:
        assert (ROOT.parents[1] / "borel" / relative).is_file(), relative
    words = []
    outside = []
    with pdfplumber.open(PDF) as pdf:
        for index, page in enumerate(pdf.pages, 1):
            text = page.extract_text() or ""
            assert "\ufffd" not in text
            words.append(len(text.split()))
            for char in page.chars:
                if char["x0"] < -0.5 or char["x1"] > page.width + 0.5:
                    outside.append(index)
                if char["top"] < -0.5 or char["bottom"] > page.height + 0.5:
                    outside.append(index)
    assert not outside, outside
    pages = sorted((ROOT / "tmp/pdfs").glob("page-*.png"))
    assert len(pages) == len(reader.pages)
    for start in range(0, len(pages), 4):
        sheet = Image.new("RGB", (1120, 1500), "#e8edf0")
        pen = ImageDraw.Draw(sheet)
        for index, path in enumerate(pages[start:start + 4]):
            x, y = (index % 2) * 560, (index // 2) * 750
            pen.text((x + 8, y + 3), "Page " + str(start + index + 1), fill="black")
            sheet.paste(Image.open(path).convert("RGB"), (x, y + 23))
        sheet.save(ROOT / "tmp/pdfs" / ("contact-" + str(start // 4 + 1) + ".png"))
    report = {
        "result": "passed",
        "pages": len(reader.pages),
        "page_layout": "SinglePage",
        "manuscript_pages_per_pdf_page": 1,
        "words_extracted": sum(words),
        "words_per_page": words,
        "numbered_equations": len(tags),
        "bibliography_items": len(keys),
        "resolved_citation_keys": sorted(cites),
        "local_companion_targets_checked": len(notes),
        "characters_outside_page": 0,
        "overfull_boxes_or_undefined_references": 0,
        "source_sha256": sha256(SOURCE.read_bytes().replace(b"\r\n", b"\n")).hexdigest(),
        "pdf_sha256": sha256(PDF.read_bytes()).hexdigest(),
        "source_baseline_commit": "2bc20702f142d6124a22ab97aaa511caf227f916",
        "scope": "Compilation, bibliography, local companion links and PDF geometry; rendered pages require separate visual inspection. This is not an external proof audit or priority claim.",
    }
    (ROOT / "verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
