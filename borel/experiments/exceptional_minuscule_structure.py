"""Exact integral generator changes for the primitive E6/E7 image forms."""
from hashlib import sha256
import json
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing

ROOT = Path(__file__).resolve().parents[1]


def main():
    path = ROOT/'experiments/exceptional-minuscule-results.json'
    data = json.loads(path.read_text())
    assert data['result'] == 'passed'
    identities = 0
    rows = []
    for kind in ('E6','E7'):
        case = next(c for c in data['cases'] if c['kind']==kind and c['form']=='individually_primitive_images')
        R = PolynomialRing(QQ,['h','s'] if kind=='E6' else ['h','s','t'])
        A = PolynomialRing(ZZ,['H','S'] if kind=='E6' else ['H','S','T'])
        if kind == 'E6':
            h,s = R.gens(); H,S = A.gens()
            relations = [2*h**9+3*h*s**2-6*h**5*s,s**3-6*h**4*s**2+h**12]
            scaled = (h,16*s)
            formulas = {'y1':H,'e1':-H**2,'e2':7*H**4-S,
                        'e3':47*H**6-7*H**2*S,
                        'e4':-457*H**8+74*H**4*S-S**2,
                        'e5':H**2*(7*H**4-S)**2,'euler':-7*H**5+H*S}
            reverse = {'H':'y1','S':'7*y1^4-e2'}
            prime = 2
        else:
            h,s,t = R.gens(); H,S,T = A.gens()
            relations = [s*s-2*h*t,2*s*t-9*h**4*s*s+6*h**9*s-h**14,
                         t*t+10*h**3*s**3-9*h**8*s*s+2*h**13*s]
            scaled = (h,81*s,19683*t)
            formulas = {'omega7':H,'psi2':-H**2,'psi5':35*H**5-S,
                        'psi6':683*H**6-20*H*S,
                        'psi8':-67537*H**8+1960*H**3*S,
                        'psi9':-1666*H**9-3*H**4*S+T,
                        'psi12':2311007*H**12-222640*H**7*S+9020*H**3*T-9020*H**2*S**2}
            reverse = {'H':'omega7','S':'35*omega7^5-psi5',
                       'T':'psi9+1666*omega7^9+3*omega7^4*(35*omega7^5-psi5)'}
            prime = 3
        G = R.ideal(relations).groebner_basis()
        originals = {g['name']:R(g['published_ring_expression'])*QQ(case['multipliers'][g['name']])
                     for g in case['generators']}
        for name,f in formulas.items():
            assert not (originals[name]-R(f(*scaled))).reduce(G), (kind,name)
            identities += 1
        D = PolynomialRing(ZZ,len(originals),names=list(originals))
        for name,p in reverse.items():
            expected = scaled[list(A.variable_names()).index(name)]
            assert not (R(D(p)(*originals.values()))-expected).reduce(G)
            identities += 1
        for row in case['degrees']:
            for v in row['smith_factors']:
                while v % prime == 0:
                    v //= prime
                assert v == 1
        rows.append({'kind':kind,'scaled_schubert_generators':list(map(str,scaled)),
                     'primitive_catalogue_in_scaled_generators':{n:str(p) for n,p in formulas.items()},
                     'scaled_generators_in_primitive_catalogue':reverse,
                     'only_possible_quotient_prime':prime})
    result = {'result':'passed','integer_generator_identities':identities,'cases':rows,
              'scope':'Exact polynomial identities modulo the published integral presentations. Both directions use integer expressions; equality over QQ implies equality in torsion-free integral cohomology.',
              'source_hash_semantics':'SHA256 with CRLF normalized to LF.',
              'source_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
                               for p in (path,Path(__file__).resolve())}}
    (ROOT/'experiments/exceptional-minuscule-structure.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASSED',identities,'integer identities; E6 primitive image is Z[h,16s]; E7 is Z[h,81s,19683t].')


if __name__ == '__main__':
    main()
