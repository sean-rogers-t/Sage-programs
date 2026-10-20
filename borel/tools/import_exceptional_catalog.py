"""Import the frozen exceptional-generator audit into separate input/test data.

This script needs ordinary Python, not SageMath. It never evaluates a Borel map.
Run from any directory, passing the audit directory from the research review:

    python tools/import_exceptional_catalog.py /path/to/papers/table_audit

The installed package receives polynomial specifications and provenance only.
The complete corrected reference expansions are written under tests/data.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path


INPUT_FILES = (
    "part_iii_corrected_rows.json",
    "e6_omitted_generator.json",
    "e8_dissertation_rows.json",
)

SOURCES = {
    "part_iii": {
        "title": "Representation ring of Levi subgroups versus cohomology ring of flag varieties III",
        "authors": ["Shrawan Kumar", "Jiale Xie"],
        "publication": "Journal of Pure and Applied Algebra 228 (2024), 107539",
        "url": "https://kumar.math.unc.edu/papers/kumar86.pdf",
        "doi": "10.1016/j.jpaa.2023.107539",
        "pdf_sha256": "973adeabbff3fef2be81daeb870d6ea33ba3aefa4187bd9ace91b786884645c6",
        "version": "Published journal PDF; not the arXiv coefficient/scaling convention",
    },
    "xie_dissertation": {
        "title": "Polynomial Representation Ring versus Cohomology Ring of Flag Varieties for Classical and Exceptional Groups",
        "author": "Jiale Xie",
        "publication": "University of North Carolina at Chapel Hill dissertation, 2022",
        "url": "https://drive.google.com/file/d/1XcFdCYXQjDD_rKGnQZejXrkOk8l5zEjx/view",
        "pdf_sha256": "95550100dda01bf99b5b10f49ae877701459d34dee1f0ad7c1fde81c0158963e",
        "version": "109-page dissertation PDF; displayed E8 expansions only",
    },
}

CONVENTIONS = {
    "coefficient_ring": "QQ",
    "cartan_numbering": "SageMath finite Cartan types (Bourbaki exceptional numbering)",
    "coordinates": "Each coordinate vector lists fundamental-weight coefficients omega_1,...,omega_rank.",
    "springer_substitution": "Normalized Springer coordinates Theta_i are replaced by omega_i before applying the Borel map.",
    "cohomological_degree": "Twice polynomial_degree",
    "schubert_index": "Minimal right coset representative w in W/W_L; beta(omega_i)=epsilon_(s_i).",
    "e6_normalization": {
        "location": "Part III, Lemma 7.1, PDF pages 19-20",
        "printed_outer_denominator": 6,
        "corrected_outer_denominator": 18,
        "effect": "The printed Springer coordinates are three times the normalized coordinates. Retaining that scaling multiplies a degree-d image by 3^d.",
    },
    "psi_e7": {
        "x_simple_root_coefficients": [
            [5, 0, 4, 3, 2, 1, 0],
            [-1, 0, 4, 3, 2, 1, 0],
            [-1, 0, -2, 3, 2, 1, 0],
            [-1, 0, -2, -3, 2, 1, 0],
            [-1, 0, -2, -3, -4, 1, 0],
            [-1, 0, -2, -3, -4, -5, 0],
        ],
        "z_simple_root_coefficients": [-3, -6, -6, -9, -6, -3, 0],
        "formula": "sum_i (x_i+z)^d + sum_i (x_i-z)^d + sum_(i<j) (-x_i-x_j)^d",
        "scaling": "Unscaled journal psi_d, not the divided arXiv inputs",
        "source": "Part III, Theorem 8.2(g), PDF pages 29-30",
    },
}


def journal_page(row):
    """One-based PDF page of the displayed image (verified against the PDF)."""
    kind, p = row["kind"], row["parabolic"]
    if kind == "G2":
        return 17
    n = int(row["id"].rsplit("row", 1)[1])
    if kind == "F4":
        return 18 if p <= 2 else 19
    if kind == "E6":
        if p == 1:
            return 20 if n <= 3 else 21
        if p == 2:
            return 21
        if p == 3:
            return 21 if n <= 4 else 22
        if p == 4:
            return 22
        if p == 5:
            return 22 if n == 1 else 23
        return 23
    if p == 1:
        return 25 if n <= 5 else 26
    if p == 2:
        return 26
    if p == 3:
        return 27
    if p == 4:
        return 27 if n <= 6 else 28
    return {5: 28, 6: 29, 7: 30}[p]


def is_euler(row, siblings):
    spec = row["input"]
    if spec["operation"] != "elementary" or spec.get("squared"):
        return False
    if spec["degree"] != len(spec["variables"]):
        return False
    return any(
        other["input"]["operation"] == "elementary"
        and other["input"].get("squared")
        and other["input"]["variables"] == spec["variables"]
        for other in siblings
    )


def name_and_aliases(row, siblings):
    spec = row["input"]
    operation = spec["operation"]
    if operation == "omega":
        return "omega" + str(spec["index"]), []
    if operation == "linear":
        return "y" + str(spec["variable"]), []
    if operation == "psi":
        return "psi" + str(spec["degree"]), []
    if operation in ("g2quadratic", "g2product"):
        return "quadratic", [operation]
    if operation != "elementary":
        raise ValueError("Unsupported audited operation: " + operation)
    short = "e" + str(spec["degree"])
    alphabet = "_".join("y" + str(v) for v in spec["variables"])
    long = short + "_" + ("squared_" if spec.get("squared") else "") + alphabet
    if is_euler(row, siblings):
        return "euler", [long, "euler_" + alphabet]
    collisions = [
        other for other in siblings
        if other["input"]["operation"] == "elementary"
        and other["input"]["degree"] == spec["degree"]
        and not is_euler(other, siblings)
    ]
    if len(collisions) == 1:
        return short, [long]
    return long, []


def source_location(row):
    if row["kind"] == "E8":
        return {
            "key": "xie_dissertation",
            "row_id": row["id"],
            "pdf_page": row["pdf_page"],
            "printed_page": row["printed_page"],
            "scope": "One of the nine displayed E8 expansions; not a complete generator catalogue",
        }
    thm = {"G2": "5.2", "F4": "6.2", "E6": "7.2", "E7": "8.2"}[row["kind"]]
    result = {
        "key": "part_iii",
        "row_id": row["id"],
        "theorem": thm + "(" + chr(96 + row["parabolic"]) + ")",
    }
    if "omitted_euler" in row["id"]:
        result.update({
            "pdf_page": 23,
            "scope": "Supplement to the printed list: the Euler generator is described but its image was omitted",
        })
    else:
        result["pdf_page"] = journal_page(row)
    return result


def normalization(row):
    result = {"convention": "Normalized Springer coordinates; Theta_i -> omega_i"}
    if row["kind"] == "E6":
        result.update({
            "springer_coordinate_denominator": 18,
            "printed_springer_coordinate_denominator": 6,
            "note": "Corrected Lemma 7.1 normalization; all polynomial coefficients below use the normalized coordinates.",
        })
    if row["input"]["operation"] == "psi":
        result.update({
            "psi_scaling": "unscaled journal convention",
            "arxiv_input_divisor": {2: 144, 5: 5760, 6: 3456, 8: 9216, 9: 387072, 12: 663552}[row["input"]["degree"]],
        })
    return result


def correction_metadata(row):
    result = []
    for item in row.get("corrections", []):
        result.append({"category": "schubert_coefficient", **deepcopy(item)})
    if "normalization_correction" in row:
        result.append({
            "category": "springer_normalization",
            "location": "Lemma 7.1, PDF pages 19-20",
            "description": row["normalization_correction"],
        })
    if "typesetting_correction" in row:
        result.append({"category": "missing_sign", "description": row["typesetting_correction"]})
    if "omitted_euler" in row["id"]:
        result.append({
            "category": "omitted_generator_image",
            "description": "The P6 Euler product is a generator in the description of the algebra, but its image is omitted from the journal generator list.",
        })
    if "printed_label" in row:
        result.append({
            "category": "generator_label",
            "printed": "e3",
            "corrected": "e4",
            "description": "The last displayed E8/P2 input is e4; its printed repeated e3 label is a typo. The numerical coefficients are unchanged.",
        })
    return result


def polynomial_degree(spec):
    if spec["operation"] in ("omega", "linear"):
        return 1
    if spec["operation"] in ("g2quadratic", "g2product"):
        return 2
    return spec["degree"] * (2 if spec.get("squared") else 1)


def build(audit_directory):
    inputs = {
        name: json.loads((audit_directory / name).read_text(encoding="utf-8"))
        for name in INPUT_FILES
    }
    rows = inputs[INPUT_FILES[0]] + [inputs[INPUT_FILES[1]]] + inputs[INPUT_FILES[2]]
    assert len(rows) == 128
    siblings = defaultdict(list)
    for row in rows:
        siblings[(row["kind"], row["parabolic"])].append(row)
    data = []
    expected = []
    for row in rows:
        name, aliases = name_and_aliases(row, siblings[(row["kind"], row["parabolic"])])
        record = {
            "kind": row["kind"],
            "parabolic": row["parabolic"],
            "name": name,
            "aliases": aliases,
            "polynomial_degree": polynomial_degree(row["input"]),
            "spec": deepcopy(row["input"]),
            "coordinates": deepcopy(row.get("y_definitions", {})),
            "source": source_location(row),
            "normalization": normalization(row),
            "corrections": correction_metadata(row),
        }
        data.append(record)
        # Frozen independently of the runtime data. Tests reconstruct the input
        # from these audited source fields and compare every Schubert coefficient.
        expected.append({
            "kind": row["kind"],
            "parabolic": row["parabolic"],
            "name": name,
            "source_id": row["id"],
            "input": deepcopy(row["input"]),
            "y_definitions": deepcopy(row.get("y_definitions", {})),
            "terms": deepcopy(row["published_terms"]),
        })
    for key, group in siblings.items():
        names = [n for r in data if (r["kind"], r["parabolic"]) == key for n in [r["name"]] + r["aliases"]]
        if len(names) != len(set(names)):
            raise ValueError("Ambiguous generator name or alias: " + repr(key))
    counts = Counter(r["kind"] for r in rows)
    coverage = {
        kind: {
            "maximal_parabolics": sorted(p for t, p in siblings if t == kind),
            "complete_generator_lists": kind != "E8",
            "catalogued_rows": counts[kind],
            "description": (
                "All maximal-parabolic journal generator lists, with the omitted E6/P6 Euler image restored"
                if kind == "E6" else
                "All maximal-parabolic journal generator lists"
                if kind != "E8" else
                "Only the displayed P1/P2 dissertation inputs: P1 has y1 and e1-e4; P2 has e1-e4. Further E8 generators and parabolics are not catalogued."
            ),
        }
        for kind in ("G2", "F4", "E6", "E7", "E8")
    }
    provenance = {
        "importer": "tools/import_exceptional_catalog.py",
        "input_files": {name: sha256((audit_directory / name).read_bytes()).hexdigest() for name in INPUT_FILES},
        "verification": "Complete corrected coefficient vectors from the separate table audit; exact BGG and positive-root Chevalley crosschecks.",
    }
    catalogue = {
        "schema_version": 1,
        "conventions": CONVENTIONS,
        "sources": SOURCES,
        "provenance": provenance,
        "coverage": coverage,
        "generators": data,
    }
    fixture = {"schema_version": 1, "provenance": provenance, "rows": expected}
    return catalogue, fixture


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit_directory", type=Path)
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    catalogue, fixture = build(args.audit_directory)
    destinations = (
        (args.output_root / "src/sage_borel/data/exceptional_generators.json", catalogue),
        (args.output_root / "tests/data/exceptional_expected.json", fixture),
    )
    for path, content in destinations:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(content, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
        print(path, len(content.get("generators", content.get("rows", []))), "rows")


if __name__ == "__main__":
    main()
