"""Package certified E7 coefficients and readable character formulas."""
from pathlib import Path
import json
from e7_geometry_support import WEIGHTS, ROOTS, weights_912, DEGREES


def main():
    root=Path(__file__).resolve().parents[1]
    data=dict(cartan_type='E7',degrees=DEGREES,
              target='p2/24,p6,p8,p10,p12,p14,p18; pd=trace_56(h^d)',
              weights_56=WEIGHTS,roots=ROOTS,weights_912=weights_912(),representations={})
    for name in ('minimal','adjoint'):
        records=[]
        for block in ('low','high'):
            certificate=json.loads((root/('experiments/e7-%s-%s-interpolation.json'%(name,block))).read_text())
            assert certificate['result']=='passed'
            records.extend(certificate['cases'][0]['coordinates'])
        determinant=json.loads((root/('experiments/e7-%s-determinant-interpolation.json'%name)).read_text())
        assert determinant['result']=='passed'
        data['representations'][name]=dict(coordinates=sorted(records,key=lambda r:r['degree']),
            torus_determinant=determinant['cases'][0]['coordinates'][0])
    destination=root/'src/sage_borel/data/e7_geometry.json'
    destination.write_text(json.dumps(data,separators=(',',':'))+'\n')
    print(destination,destination.stat().st_size,flush=True)
    from sage_borel import e7_adjoint_quotient, e7_cayley_determinant
    readable={}
    for name in ('minimal','adjoint'):
        q=e7_adjoint_quotient(representation=name)
        d=e7_cayley_determinant(representation=name)
        readable[name]=dict(coordinate_strings={str(k):str(f) for k,f in zip(q.degrees,q.character_polynomials)},
                            torus_determinant=str(d.torus_character_polynomial),
                            factored_torus_determinant=str(d.torus_character_polynomial.factor()))
    (root/'experiments/e7-geometry-formulas.json').write_text(json.dumps(readable,indent=2)+'\n')


if __name__=='__main__':
    main()
