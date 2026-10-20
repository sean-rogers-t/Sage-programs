"""First quotient coordinate for E6/E7 minimal and adjoint trace maps.

This is not a full adjoint-quotient or full-determinant computation.
"""
from pathlib import Path
import json
from time import perf_counter
from sage.all import QQ, matrix, vector, WeylCharacterRing
from sage_borel import springer_character_map, SpringerGeometry


def main():
    reports = []
    for kind, minimal, adjoint in (("E6", (1,0,0,0,0,0), (0,1,0,0,0,0)),
                                    ("E7", (0,0,0,0,0,0,1), (1,0,0,0,0,0,0))):
        for name, weight in (("minimal", minimal), ("adjoint", adjoint)):
            start = perf_counter()
            s = springer_character_map(kind, weight)
            g = SpringerGeometry(s)
            w = vector(s.ring.gens())
            q = (w*matrix(QQ, s.cartan_type.cartan_matrix())*w)/2
            f = s.pullback(q)
            answer = g.character_coordinates(f)
            # Independent representation-theoretic check. The radial Casimir
            # product rule gives |grad chi|^2 from the tensor-square
            # decomposition, without differentiating Laurent polynomials.
            WCR = WeylCharacterRing(kind, base_ring=QQ, style="coroots")
            character = WCR(weight)
            Cinv = matrix(QQ, s.cartan_type.cartan_matrix()).inverse()
            rho_twice = vector(QQ, [2]*s.rank)
            casimir = lambda a: vector(QQ, a)*Cinv*(vector(QQ, a)+rho_twice)
            c = casimir(weight)
            index = s.trace_gram[0, 0]/2
            alternative = WCR.zero()
            for nu, multiplicity in (character*character).monomial_coefficients().items():
                coroots = WCR.space().simple_coroots()
                coords = tuple(int(nu.inner_product(coroots[j+1])) for j in range(s.rank))
                alternative += multiplicity*(casimir(coords)-2*c)/(4*index**2)*WCR(coords)
            fundamentals = [WCR(tuple(int(i == j) for j in range(s.rank))) for i in range(s.rank)]
            assert answer(*fundamentals) == alternative
            print(kind, name, answer, "seconds", perf_counter()-start, flush=True)
            reports.append(dict(cartan_type=kind, representation=name, dimension=s.dimension,
                                target_invariant=str(q), character_pullback=str(answer),
                                laurent_terms=len(f.dict()),
                                character_nodes_used=[i+1 for i in sorted(g._fundamentals)],
                                tensor_casimir_check=True,
                                seconds=perf_counter()-start))
    Path(__file__).with_name("exceptional-quadratic-geometry.json").write_text(
        json.dumps(dict(result="passed", scope="Quadratic coordinate only; exact Weyl invariance, Laurent reconstruction, and an independent tensor-Casimir identity. Not a full matrix audit.", cases=reports), indent=2)+"\n")


if __name__ == "__main__":
    main()
