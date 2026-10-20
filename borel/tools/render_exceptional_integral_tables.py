"""Render the complete, explicitly normalized G2/F4 Smith tables."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT/'experiments/exceptional-integral-results.json').read_text())
assert data['result'] == 'passed'
lines = ['# Complete G2/F4 integral image tables', '',
         'Generated from [the exact dataset](experiments/exceptional-integral-results.json).',
         'See [the mathematical guide](EXCEPTIONAL_INTEGRAL.md) for normalization,',
         'proofs, API usage and independent verification. These are selected source',
         'Z-forms, not an identified representation-ring intersection.', '',
         'Degree d means H^(2d). Each tuple lists every nonzero Smith factor,',
         'including 1. All degree images have full rational rank. Degrees above',
         'the listed complex dimension have zero cohomology.', '',
         '| Space | Complex dimension | Top factor: denominators cleared | Top factor: primitive generator images |',
         '| --- | --- | --- | --- |']
pairs = []
for index in range(0, len(data['cases']), 2):
    a, b = data['cases'][index:index+2]
    assert (a['kind'],a['parabolic']) == (b['kind'],b['parabolic'])
    pairs.append((a,b))
    lines.append('| {}/P{} | {} | {} | {} |'.format(
        a['kind'],a['parabolic'],a['complex_dimension'],
        a['degrees'][-1]['smith_factors'][0],b['degrees'][-1]['smith_factors'][0]))
for a,b in pairs:
    lines += ['', '## {}/P{}'.format(a['kind'],a['parabolic']), '',
              '| d | Image denominators cleared | Individually primitive images |',
              '| --- | --- | --- |']
    for x,y in zip(a['degrees'],b['degrees']):
        fmt = lambda values: '(' + ', '.join(map(str,values)) + ')'
        lines.append('| {} | {} | {} |'.format(x['degree'],fmt(x['smith_factors']),fmt(y['smith_factors'])))
lines += ['', 'Regenerate with `python tools/render_exceptional_integral_tables.py`.', '']
(ROOT/'EXCEPTIONAL_INTEGRAL_TABLES.md').write_text('\n'.join(lines), encoding='utf-8')
print('Rendered', data['degree_rows'], 'degree/form rows.')
