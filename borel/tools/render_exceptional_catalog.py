"""Render the discoverable exceptional generator index without SageMath.

Run: python tools/render_exceptional_catalog.py
"""

import json
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    data = json.loads((root / "src/sage_borel/data/exceptional_generators.json").read_text(encoding="utf-8"))
    lines = [
        "# Exceptional generator index", "",
        "Generated from the packaged input catalogue. Degrees below are polynomial",
        "degrees; cohomological degrees are twice these values. See the",
        "[exceptional guide](EXCEPTIONAL.md) for normalization and coverage.", "",
        "Call `exceptional_image(type, parabolic, name)` using a row below. Names",
        "with `squared` distinguish elementary polynomials in squared coordinates.",
        "The E8 rows are the displayed dissertation subset, not complete generator lists.", "",
        "| Type | Parabolic | Generator name | Degree | Source |",
        "| --- | --- | --- | ---: | --- |",
    ]
    type_order = {kind: i for i, kind in enumerate(("G2", "F4", "E6", "E7", "E8"))}
    ordered = sorted(enumerate(data["generators"]),
                     key=lambda pair: (type_order[pair[1]["kind"]], pair[1]["parabolic"], pair[0]))
    for _, entry in ordered:
        source = entry["source"]
        if source["key"] == "part_iii":
            citation = "Part III " + source.get("theorem", "supplement")
            if "omitted_euler_product" in source.get("row_id", ""):
                citation += " (supplement)"
        else:
            citation = "Xie dissertation p." + str(source["printed_page"])
        if source.get("pdf_page"):
            citation += "; PDF p." + str(source["pdf_page"])
        lines.append("| {} | P{} | `{}` | {} | {} |".format(
            entry["kind"], entry["parabolic"], entry["name"],
            entry["polynomial_degree"], citation))
    lines.extend(["", "This file can be regenerated with `python tools/render_exceptional_catalog.py`.", ""])
    (root / "EXCEPTIONAL_GENERATORS.md").write_text("\n".join(lines), encoding="utf-8")
    print("Rendered {} named generators.".format(len(data["generators"])))


if __name__ == "__main__":
    main()
