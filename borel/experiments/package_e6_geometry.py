"""Build runtime coefficient data from the exact interpolation certificates."""
from pathlib import Path
import json
from e6_quotient_interpolation import WEIGHTS, ROOTS


def main():
    root=Path(__file__).resolve().parents[1]
    data=dict(cartan_type='E6', degrees=[2,5,6,8,9,12],
              target='p2/12,p5,p6,p8,p9,p12; pd=trace_27(h^d)',
              weights_27=WEIGHTS,roots=ROOTS,representations={})
    for name in ('minimal','adjoint'):
        quotient=json.loads((root/('experiments/e6-%s-quotient-interpolation.json'%name)).read_text())
        determinant=json.loads((root/('experiments/e6-%s-determinant-interpolation.json'%name)).read_text())
        coordinates=sorted([c for b in quotient['cases'] for c in b['coordinates']],key=lambda c:c['degree'])
        data['representations'][name]=dict(coordinates=coordinates,
            torus_determinant=determinant['cases'][0]['coordinates'][0])
    destination=root/'src/sage_borel/data/e6_geometry.json'
    destination.write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(destination, destination.stat().st_size)
    from sage_borel import e6_adjoint_quotient, e6_cayley_determinant
    readable={}
    for name in ('minimal','adjoint'):
        q=e6_adjoint_quotient(representation=name)
        d=e6_cayley_determinant(representation=name)
        readable[name]=dict(coordinate_strings={str(k):str(f) for k,f in zip(q.degrees,q.character_polynomials)},
                            torus_determinant=str(d.torus_character_polynomial),
                            factored_torus_determinant=str(d.torus_character_polynomial.factor()))
    (root/'experiments/e6-geometry-formulas.json').write_text(json.dumps(readable,indent=2)+'\n')


if __name__=='__main__':
    main()
