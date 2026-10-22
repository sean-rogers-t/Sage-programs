"""Reproducibly package E6 relative theorem records, retaining all 3-powers."""
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def packaged_data():
    source = ROOT/'experiments/e6-all-parabolics-summary.json'
    data = json.loads(source.read_text())
    assert data['result'] == 'passed' and data['complete']
    assert len(data['parabolics']) == 63
    rows = []
    for row in data['parabolics']:
        references = (['E6_ROOT_WEIGHT_COMPARISON.md'] if row['proof_family'] == 'good-Levi-60'
                      else ['E6_REMAINING_CASES.md', 'E6_REMAINING_TABLES.md'])
        rows.append({**row, 'proof_references': references})
    return dict(schema_version=2, cartan_type='E6',
                certificate='experiments/e6-all-parabolics-summary.json',
                certificate_sha256=sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
                parabolics=rows)


if __name__ == '__main__':
    target = ROOT/'src/sage_borel/data/e6_relative_images.json'
    target.write_text(json.dumps(packaged_data(), separators=(',', ':'))+'\n')
    print('Packaged 63 E6 parabolics:', target)
