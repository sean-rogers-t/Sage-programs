"""Render all E7/P7 literal-character image and embedded comparisons."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'experiments/e7-character-form-results.json').read_text())
assert data['result']=='passed'
fmt=lambda xs:'('+', '.join(map(str,xs))+')'
lines=['# E7/P7: complete integral character-image comparison','',
       'See the [integral-form proof and methods](E7_CHARACTER_FORM.md) and',
       '[exact dataset](experiments/e7-character-form-results.json).','',
       'Polynomial degree d means cohomological degree 2d. All images have full',
       'rational rank. In the character column, multiply **each** displayed',
       'factor by 12^d. This removes the uniform trace-index scaling only for',
       'readability; the actual character image retains that factor. Degrees above 27 vanish.','',
       '| d | Character factors / 12^d | Original catalogue factors | Primitive-image factors |',
       '| --- | --- | --- | --- |']
for row in data['degrees']:
    scale=12**row['degree'];assert all(x%scale==0 for x in row['smith_factors'])
    a,b=row['comparisons']
    lines.append('| {} | {} | {} | {} |'.format(row['degree'],fmt([x//scale for x in row['smith_factors']]),
        fmt(a['old_smith_factors']),fmt(b['old_smith_factors'])))
lines+=['','## Embedded comparisons','',
        'Here C is the character image and J is the indicated earlier image.',
        'The tuple for (C+J)/J also describes C/(C intersection J).',
        'The tuple for (C+J)/C also describes J/(C intersection J).',
        'All tuples are computed using sums/intersections in shared Schubert coordinates.','',
        '| d | Earlier form J | C contained in J? | J contained in C? | (C+J)/J | (C+J)/C |',
        '| --- | --- | --- | --- | --- | --- |']
for row in data['degrees']:
    for c in row['comparisons']:
        name='Catalogue' if c['form']=='image_denominators_cleared' else 'Primitive'
        lines.append('| {} | {} | {} | {} | {} | {} |'.format(row['degree'],name,
            'Yes' if c['old_contains_character'] else 'No','Yes' if c['character_contains_old'] else 'No',
            fmt(c['sum_mod_old']),fmt(c['sum_mod_character'])))
lines+=['','## Unfactored character Smith lists','',
        '| d | Actual character-image factors |','| --- | --- |']
for row in data['degrees']:
    lines.append('| {} | {} |'.format(row['degree'],fmt(row['smith_factors'])))
lines+=['','Regenerate with `python tools/render_e7_character_tables.py`.','']
(ROOT/'E7_CHARACTER_TABLES.md').write_text('\n'.join(lines),encoding='utf-8')
print('Rendered',len(data['degrees']),'degree rows and',2*len(data['degrees']),'embedded comparisons.')
