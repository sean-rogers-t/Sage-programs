"""Bind the degree derivation, finite corroboration, tests and installed wheel."""
from hashlib import sha256
from pathlib import Path
import json
import zipfile
from sage.all import GF, PolynomialRing
import sage_borel
from monomial_staircase import staircase_length

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def read(name):
    data=json.loads((ROOT/name).read_text())
    assert data['result']=='passed',name
    return data


def main():
    assert staircase_length([(2,0),(0,2)])==4
    assert staircase_length([(3,0),(1,1),(0,4)])==6
    assert staircase_length([(0,0)])==0
    version=sage_borel.__version__
    regression=read('character-verification.json')
    docs=read('documentation-verification.json')
    wheel=read('geometry-wheel-verification.json')
    geometry=read('experiments/exceptional-degree-geometry.json')
    assert regression['version']==wheel['version']==version=='0.20.0'
    for name,expected in regression['source_sha256'].items():
        assert digest(ROOT/name)==expected,name
    wheel_path=ROOT/'dist'/wheel['wheel']
    assert sha256(wheel_path.read_bytes()).hexdigest()==wheel['wheel_sha256']
    runtime=list((ROOT/'src/sage_borel').glob('*.py'))
    runtime+=list((ROOT/'src/sage_borel/data').glob('*.json'))
    with zipfile.ZipFile(wheel_path) as archive:
        for path in runtime:
            name=path.relative_to(ROOT/'src').as_posix()
            assert archive.read(name).replace(b'\r\n',b'\n')==path.read_bytes().replace(b'\r\n',b'\n'),name
    degrees={}
    for entry in geometry['degree_deductions']:
        key=(entry['cartan_type'],entry['representation'])
        c=sage_borel.exceptional_springer_degree(key[0],representation=key[1])
        assert (c.normalized_volume,c.boundary_contribution,c.central_cover,c.degree)==(
            entry['normalized_volume'],entry['boundary_correction'],entry['central_cover'],entry['degree'])
        degrees[key]=c.degree
    fibers=[]
    fiber_paths=sorted((ROOT/'experiments').glob('e[67]-*-fiber-*.json'))
    for path in fiber_paths:
        f=json.loads(path.read_text())
        key=(f['cartan_type'],f['representation'])
        assert f['status']=='specialized-fiber-only'
        assert f['fiber_dimension']==0 and f['fiber_length']==degrees[key],path
        # Recompute the length from the recorded leading ideal. This checks
        # the report's staircase count, not the Groebner-basis derivation.
        n=(f['base_rank'] if 'module_leading_terms' in f else
           len(f['leading_exponents'][0]) if 'leading_exponents' in f else len(f['targets']))
        R=PolynomialRing(GF(f['characteristic']),n,names=['x%d'%(i+1) for i in range(n)])
        if 'module_leading_terms' in f:
            ideals=[[] for _ in range(f['module_rank'])]
            for row in f['module_leading_terms']:
                ideals[row[0]-1].append(R({tuple(row[1:]):1}))
            components=[int(R.ideal(gens).vector_space_dimension()) for gens in ideals]
            assert components==f['module_component_lengths']
            for gens,component_length in zip(ideals,components):
                assert staircase_length([tuple(g.exponents()[0]) for g in gens])==component_length
            assert f['preflight_lengths']==[1,4]
            if 'tower_reduction' in f:
                assert f['rank_four_preflight_length']==4
                assert f['tower_reduction']['first_target_matrix_identity']
                assert f['tower_reduction']['second_target_matrix_identity']
                assert f['tower_reduction']['commuting_regular_representation']
            length=sum(components)
        else:
            monomials=([R({tuple(e):1}) for e in f['leading_exponents']]
                       if 'leading_exponents' in f else [R(m) for m in f['leading_monomials']])
            length=int(R.ideal(monomials).vector_space_dimension())
            assert staircase_length([tuple(g.exponents()[0]) for g in monomials])==length
        cover=(f.get('center_chart') or {}).get('cover_degree',1)
        if f.get('center_chart'):
            assert f['center_chart']['excluded_original_divisor_empty']
            assert f['center_chart']['excluded_quotient_divisor_empty']
        assert cover*length==f['fiber_length'],path
        fibers.append(dict(file=path.relative_to(ROOT).as_posix(),cartan_type=key[0],
            representation=key[1],prime=f['characteristic'],seed=f['seed'],
            length=f['fiber_length'],engine=f['engine']))
    paths=fiber_paths+runtime+[ROOT/name for name in (
        'EXCEPTIONAL_DEGREES.md','experiments/exceptional_degree_geometry.py',
        'experiments/exceptional-degree-geometry.json','experiments/exceptional_degree_fibers.py',
        'experiments/exceptional-degree-bounded-attempts.json',
        'experiments/e7_adjoint_module_fiber.py',
        'experiments/e7_adjoint_shifted_fiber.py',
        'experiments/e7_adjoint_tower_fiber.py',
        'experiments/e6_matrix_oracle.py','experiments/e7_matrix_oracle.py',
        'README.md','../README.md','E6_GEOMETRY.md','E7_GEOMETRY.md',
        'RESEARCH_LEDGER.md','AGENTS.md',
        '../articles/research-program-proposal/main.tex',
        '../articles/research-program-proposal/README.md',
        'PROOF_STATUS.md','EXTENSIONS.md','SPRINGER_GEOMETRY.md',
        'tools/verify_documentation.py','tools/verify_geometry_wheel.py',
        'tools/monomial_staircase.py',
        'tests/test_exceptional_degrees.py','character-verification.json',
        'documentation-verification.json','geometry-wheel-verification.json')]+[Path(__file__).resolve()]
    report=dict(result='passed',version=version,
        characteristic_zero_proof='EXCEPTIONAL_DEGREES.md; toric localized top Chern and explicit boundary Segre classes',
        degree_deductions=geometry['degree_deductions'],specialized_fibers=fibers,
        full_fiber_cases_not_yet_corroborated=[list(key) for key in degrees
            if key not in {(f['cartan_type'],f['representation']) for f in fibers}],
        finite_check_scope='Specialized lengths corroborate the separate geometric proof; no inference of characteristic-zero generic degree from specialization.',
        standard_monomial_counts_rechecked='Singular vdim and independent integer exponent-slice recursion',
        regression_tests=regression['test_methods_run'],documentation_examples=docs['python_blocks_executed'],
        local_links_and_anchors=docs['local_links_and_anchors_checked'],isolated_wheel_tests=wheel['tests'],
        wheel_sha256=wheel['wheel_sha256'],wheel_sources_and_data_match=True,
        hash_semantics='SHA256 after CRLF to LF normalization',
        audit_sha256={p.relative_to(ROOT).as_posix():digest(p) for p in sorted(set(paths))})
    (ROOT/'release-0.20-verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='audit_sha256'},indent=2))


if __name__=='__main__':
    main()
