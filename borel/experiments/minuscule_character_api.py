"""Full E6/E7 API replay and new E6 adjoint image comparison.

Fresh API source lattices are compared with the previously audited minimal
tables. Every rational Borel column is checked with the separate published
ordinary-ring/Chevalley implementation. No saved images feed the API.
"""
from pathlib import Path
from hashlib import sha256
from time import perf_counter
import json
from sage.all import QQ, ZZ, matrix, prod, CartanMatrix
from sage_borel import e6_character_source, e7_character_source
from e6_character_form import certify_preimage
from exceptional_published_rings import PublishedMinusculeRing

ROOT = Path(__file__).resolve().parents[1]


def rows(M):
    return [[str(c) for c in r] for r in M]


def same_lattice(A, B):
    T = A.change_ring(QQ).inverse()*B
    assert T.denominator() == 1 and abs(T.det()) == 1


def coefficient_matrix(source, monomials):
    fs = [prod(g**int(e) for g, e in zip(source._inputs, m.exponents()[0]))
          for m in monomials]
    dictionaries = [f.dict() for f in fs]
    support = sorted(set().union(*(set(D) for D in dictionaries)))
    return matrix(QQ, sorted(set(tuple(D.get(e, 0) for D in dictionaries)
                                for e in support)))


def main():
    started = perf_counter()
    output = ROOT/"experiments/minuscule-character-api-results.json"
    report = dict(result="running", cases=[])
    baselines = {}
    files = [Path(__file__).resolve()]
    for kind in ("e6", "e7"):
        p = ROOT/("experiments/"+kind+"-character-form-results.json")
        audit = ROOT/("experiments/"+kind+"-character-form-audit.json")
        data, check = json.loads(p.read_text()), json.loads(audit.read_text())
        assert data["result"] == check["result"] == "passed"
        assert sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest() == check["source_sha256"][p.relative_to(ROOT).as_posix()]
        baselines[kind.upper()] = data
        files += [p, audit]
    sources = [("E6", "minimal", e6_character_source()),
               ("E6", "adjoint", e6_character_source(representation="adjoint")),
               ("E7", "minimal", e7_character_source())]
    libraries = {kind: PublishedMinusculeRing(kind) for kind in ("E6", "E7")}
    comparisons = []
    for kind, representation, source in sources:
        algebra = source.rational_algebra
        case = dict(kind=kind, representation=representation, degrees=[])
        pub = libraries[kind]
        for old in baselines[kind]["degrees"]:
            d = old["degree"]
            row = source.degree(d, max_matrix_entries=10000000)
            im = row.image()
            powers = [list(map(int, m.exponents()[0])) for m in row.monomials]
            assert powers == old["monomial_exponents"]
            rational = algebra.degree(d)
            words = [list(w.reduced_word()) for w in rational.schubert_indices]
            columns = []
            for p in powers:
                polynomial = prod(pub.generator_expressions[n]**e
                                  for n, e in zip(algebra.generator_names, p))
                check = {tuple(w.reduced_word()): c for w, c in pub.expand(polynomial, d).items()}
                # Convert Weyl elements, since reduced words need not be unique.
                converted = {algebra.borel_map.W.from_reduced_word(w): c for w, c in check.items()}
                columns.append([converted.get(w, 0) for w in rational.schubert_indices])
            assert rational.matrix == matrix(QQ, columns).transpose()
            entry = dict(degree=d, monomial_exponents=powers, source_rank=len(powers),
                         source_basis=rows(row.basis_matrix), schubert_words=words,
                         rational_image_matrix=rows(rational.matrix),
                         image_matrix=rows(rational.matrix*row.basis_matrix),
                         image_basis=rows(im.basis), smith_factors=list(map(int, im.smith_factors)))
            if representation == "minimal":
                same_lattice(matrix(QQ, old["source_basis_matrix"]), row.basis_matrix)
                prior = matrix(QQ, [old["image_basis"][old["schubert_words"].index(w)] for w in words])
                same_lattice(prior, im.basis)
                assert entry["smith_factors"] == old["smith_factors"]
            else:
                M = coefficient_matrix(source, row.monomials)
                B0 = row.basis_matrix/(24**d)
                entry["coefficient_rows"] = M.nrows()
                entry["completeness_certificate"] = certify_preimage(M, B0)
                a = sources[0][2].degree(d, max_matrix_entries=10000000).image()
                comparison = a.compare(im)
                entry["comparison"] = dict(
                    minimal_in_adjoint=comparison.left_in_right,
                    adjoint_in_minimal=comparison.right_in_left,
                    sum_basis=rows(comparison.sum.basis),
                    intersection_basis=rows(comparison.intersection.basis),
                    sum_mod_minimal=comparison.sum.quotient_by(a).as_dict(),
                    sum_mod_adjoint=comparison.sum.quotient_by(im).as_dict(),
                    minimal_mod_intersection=a.quotient_by(comparison.intersection).as_dict(),
                    adjoint_mod_intersection=im.quotient_by(comparison.intersection).as_dict())
                comparisons.append((d, entry["smith_factors"], comparison.left_in_right, comparison.right_in_left))
            case["degrees"].append(entry)
            print(kind, representation, d, "rank", len(powers), "elapsed",
                  round(perf_counter()-started, 2), flush=True)
        report["cases"].append(case)
        output.write_text(json.dumps(report, indent=2)+"\n")
    files += [ROOT/"experiments/e6_character_form.py",
              ROOT/"experiments/exceptional_published_rings.py",
              ROOT/"experiments/exceptional_minuscule_images.py",
              ROOT/"src/sage_borel/data/exceptional_generators.json"]
    files += sorted((ROOT/"src/sage_borel").glob("*.py"))
    report.update(result="passed", fresh_degree_rows=62,
                  independently_checked_borel_columns=2731,
                  seconds=round(perf_counter()-started, 3),
                  source_hash_semantics="SHA256 with CRLF normalized to LF.",
                  source_sha256={p.relative_to(ROOT).as_posix():sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest() for p in files})
    assert sum(len(c["degrees"]) for c in report["cases"]) == 62
    assert sum(r["source_rank"] for c in report["cases"] for r in c["degrees"]) == 2731
    output.write_text(json.dumps(report, indent=2)+"\n")
    print("PASSED", report["seconds"], comparisons, flush=True)


if __name__ == "__main__":
    main()
