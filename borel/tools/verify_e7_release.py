"""Bind the E7 mathematical audit, regressions, docs and wheel to this tree."""
from hashlib import sha256
from pathlib import Path
import json
import zipfile
import sage_borel

ROOT=Path(__file__).resolve().parents[1]


def normalized_hash(path):
    return sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def read(name):
    data=json.loads((ROOT/name).read_text())
    assert data['result']=='passed',name
    return data


def main():
    version=sage_borel.__version__
    regression=read('character-verification.json')
    docs=read('documentation-verification.json')
    wheel=read('geometry-wheel-verification.json')
    audit=read('experiments/e7-geometry-audit.json')
    assert regression['version']==wheel['version']==version=='0.19.0'
    for name,expected in regression['source_sha256'].items():
        assert normalized_hash(ROOT/name)==expected,name
    wheel_path=ROOT/'dist'/wheel['wheel']
    assert sha256(wheel_path.read_bytes()).hexdigest()==wheel['wheel_sha256']
    runtime=list((ROOT/'src/sage_borel').glob('*.py'))
    runtime+=list((ROOT/'src/sage_borel/data').glob('*.json'))
    with zipfile.ZipFile(wheel_path) as archive:
        for path in runtime:
            name=path.relative_to(ROOT/'src').as_posix()
            assert archive.read(name).replace(b'\r\n',b'\n')==path.read_bytes().replace(b'\r\n',b'\n'),name
    for name in ('minimal','adjoint'):
        for block in ('low','high','determinant'):
            certificate=read('experiments/e7-%s-%s-interpolation.json'%(name,block))['cases'][0]
            n=certificate['candidate_count']
            assert n==len(certificate['points'])==certificate['modular_evaluation_rank']==certificate['exact_equations_checked']
    assert audit['simple_generator_relations']==231
    assert audit['independently_verified_commutators']==17689
    assert audit['independently_reconstructed_912_weights']==632
    coordinates=sum(c['exact_scalar_identities'] for case in audit['cases']
                    for c in case['unisolvent_certificates'] if c['block']!='determinant')
    determinants=sum(c['exact_scalar_identities'] for case in audit['cases']
                     for c in case['unisolvent_certificates'] if c['block']=='determinant')
    assert coordinates==32727 and determinants==196
    paths=list((ROOT/'experiments').glob('e7*'))
    paths+=[ROOT/'experiments/audit_e7_geometry.py',ROOT/'experiments/package_e7_geometry.py',
            ROOT/'experiments/audit_trace_source_survey.py',ROOT/'E7_GEOMETRY.md',
            ROOT/'E6_GEOMETRY.md',ROOT/'README.md',ROOT/'../README.md',
            ROOT/'PROOF_STATUS.md',ROOT/'EXTENSIONS.md',ROOT/'SPRINGER_GEOMETRY.md',
            ROOT/'tools/verify_documentation.py',ROOT/'tools/verify_geometry_wheel.py',
            Path(__file__).resolve()]
    paths=[p for p in paths if p.is_file()]
    report=dict(result='passed',version=version,regression_tests=regression['test_methods_run'],
                failures=regression['failures'],errors=regression['errors'],
                documentation_examples=docs['python_blocks_executed'],
                local_links_and_anchors=docs['local_links_and_anchors_checked'],
                isolated_wheel_tests=wheel['tests'],wheel_sha256=wheel['wheel_sha256'],
                wheel_sources_and_data_match=True,
                current_regression_dependency_hashes=len(regression['source_sha256']),
                representations=['E7 minimal 56','E7 adjoint 133'],
                quotient_degrees=[2,6,8,10,12,14,18],exact_coordinate_checks=coordinates,
                exact_cartan_determinant_checks=determinants,coefficientwise_root_blocks=252,
                independent_commutators=17689,
                scalar_normalization='d theta at identity is identity; both full determinants equal 1 there',
                determinant_format='Complete exact products: 7-dimensional Cartan determinant and all 63 opposite-root pairs; not expanded into one character polynomial.',
                hash_semantics='SHA256 after CRLF to LF normalization',
                audit_sha256={p.relative_to(ROOT).as_posix():normalized_hash(p) for p in sorted(set(paths))})
    destination=ROOT/'release-0.19-verification.json'
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='audit_sha256'},indent=2))


if __name__=='__main__':
    main()
