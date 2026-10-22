"""Reproducibly package audited E7 theorem data; this is not a new proof."""
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def packaged_data():
    source = ROOT/'experiments/e7-all-parabolics-summary.json'
    data = json.loads(source.read_text())
    assert data['result'] == 'passed' and data['complete']
    assert len(data['parabolics']) == 127
    rows = []
    for row in data['parabolics']:
        family = row['proof_family']
        note = ('E7_D4_REFINEMENT.md' if family == 'uniform-112' else
                'E7_REMAINING_PARABOLICS.md' if family == 'D6' else
                'E7_D5_A1_COMPARISON.md' if family == 'D5+A1' else
                'E7_D4_A1_COMPARISON.md' if family == 'D4+A1' else 'E7_A5_FAMILY.md')
        rows.append({**{key: row[key] for key in (
            'retained_nodes', 'levi_type', 'dimension', 'cohomology_ranks',
            'multiplicity_coefficients', 'proof_family')},
            'proof_references': ['E7_ALL_PARABOLICS.md', note]})
    return dict(schema_version=1, cartan_type='E7',
                certificate='experiments/e7-all-parabolics-summary.json',
                certificate_sha256=sha256(source.read_bytes().replace(b'\r\n', b'\n')).hexdigest(),
                parabolics=rows)


if __name__ == '__main__':
    target = ROOT/'src/sage_borel/data/e7_relative_images.json'
    target.write_text(json.dumps(packaged_data(), separators=(',', ':'))+'\n')
    print('Packaged 127 E7 parabolics:', target)
