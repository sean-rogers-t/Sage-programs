"""Exact G2 quotient coordinates and critical determinants; run with Sage."""
from pathlib import Path
import json
from sage.all import QQ, matrix, identity_matrix
from sage_borel import SpringerGeometry, springer_character_map, g2_adjoint_quotient


def main():
    reports = []
    for name, weight in (("minimal", (1, 0)), ("adjoint", (0, 1))):
        s = springer_character_map("G2", weight)
        g = SpringerGeometry(s)
        q = g2_adjoint_quotient(representation=name)
        print(name, "quotient", q.character_polynomials, flush=True)
        print(name, "degree-six pullback factored", q.character_polynomials[1].factor(), flush=True)
        J = g.torus_jacobian()
        assert matrix(QQ, [[f(1, 1) for f in row] for row in J]) == identity_matrix(2)
        psi = g.cayley_determinant()
        assert psi(1, 1) == 1
        torus = g.character_coordinates(J.det())
        full = g.character_coordinates(psi)
        print(name, "torus Jacobian", torus.factor(), flush=True)
        print(name, "full determinant", full.factor(), flush=True)
        reports.append(dict(representation=name, dimension=s.dimension,
                            target_invariants=list(map(str, q.invariants)),
                            quotient_coordinates=list(map(str, q.character_polynomials)),
                            quotient_coordinates_factored=[str(f.factor()) for f in q.character_polynomials],
                            quotient_jacobian=str(q.character_jacobian().det().factor()),
                            generic_degree=q.generic_degree(),
                            torus_jacobian=str(torus),
                            full_cayley_determinant=str(full),
                            factored_cayley_determinant=str(full.factor()),
                            laurent_terms=len(psi.dict())))
    destination = Path(__file__).with_name("springer-geometry-results.json")
    destination.write_text(json.dumps(dict(result="passed", cases=reports), indent=2)+"\n")


if __name__ == "__main__":
    main()
