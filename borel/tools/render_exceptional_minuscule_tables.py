"""Render all E6/P1 and E7/P7 degreewise integral image factors."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
data = json.loads((ROOT/'experiments/exceptional-minuscule-results.json').read_text())
assert data['result'] == 'passed'
lines = ['# Complete E6/P1 and E7/P7 integral image tables', '',
         'Generated from [the exact dataset](experiments/exceptional-minuscule-results.json).',
         'See [the mathematical guide](EXCEPTIONAL_MINUSCULE.md) for the source forms,',
         'published presentations, integer generator identities and independent checks.', '',
         'Degree d means H^(2d). Each tuple includes all Smith factors, including 1.',
         'Every image has full rational rank. Original catalogue images are already',
         'integral; the second column of factors divides each generator image by its',
         'own Schubert content before forming products. Higher target degrees vanish.', '']
for index in range(0,len(data['cases']),2):
    a,b = data['cases'][index:index+2]
    assert (a['kind'],a['parabolic']) == (b['kind'],b['parabolic'])
    assert a['form'] == 'image_denominators_cleared'
    assert b['form'] == 'individually_primitive_images'
    lines += ['## {}/P{}'.format(a['kind'],a['parabolic']), '',
              'Complex dimension: {}. Total Schubert classes: {}.'.format(a['complex_dimension'],a['fixed_point_count']), '',
              '| d | Original catalogue | Individually primitive images |',
              '| --- | --- | --- |']
    for x,y in zip(a['degrees'],b['degrees']):
        assert x['degree'] == y['degree']
        fmt = lambda v: '(' + ', '.join(map(str,v)) + ')'
        lines.append('| {} | {} | {} |'.format(x['degree'],fmt(x['smith_factors']),fmt(y['smith_factors'])))
    lines.append('')
lines += ['Regenerate with `python tools/render_exceptional_minuscule_tables.py`.', '']
(ROOT/'EXCEPTIONAL_MINUSCULE_TABLES.md').write_text('\n'.join(lines),encoding='utf-8')
print('Rendered',data['degree_rows'],'degree/form rows.')
