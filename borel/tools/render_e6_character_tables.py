"""Render the complete literal-character/catalogue/primitive comparison."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'experiments/e6-character-form-results.json').read_text())
assert data['result']=='passed'
fmt=lambda xs:'('+', '.join(map(str,xs))+')'
lines=['# E6/P1: complete integral character-image comparison','',
       'See the [proof and methods](E6_CHARACTER_FORM.md) and',
       '[exact dataset](experiments/e6-character-form-results.json).', '',
       'Polynomial degree d means cohomological degree 2d. All image lattices',
       'have full rank in their target degree. Factors 1 are retained.',
       'The literal character image is contained in each of the two prior images;',
       'containments are checked in common Schubert coordinates. Degrees above 16 vanish.', '',
       '| d | Literal Levi characters | Original catalogue | Individually primitive images |',
       '| --- | --- | --- | --- |']
for row in data['degrees']:
    a,b=row['comparisons']
    assert a['form']=='image_denominators_cleared' and b['form']=='individually_primitive_images'
    lines.append('| {} | {} | {} | {} |'.format(row['degree'],fmt(row['smith_factors']),
                 fmt(a['old_smith_factors']),fmt(b['old_smith_factors'])))
lines+=['','## Quotients of the previous images by the character image','',
        'These are quotients between embedded lattices, not ratios of individual',
        'ambient Smith factors. They use the integer containment matrices.', '',
        '| d | Catalogue / character | Primitive / character |','| --- | --- | --- |']
for row in data['degrees']:
    a,b=row['comparisons']
    lines.append('| {} | {} | {} |'.format(row['degree'],fmt(a['quotient_factors']),fmt(b['quotient_factors'])))
lines+=['','## Removing the uniform trace-index scaling','',
        'The character form in degree d is 6^d times the invariant root-lattice',
        'form Sym_Z(Q)^(W_L). The following factors are therefore the character',
        'Smith factors divided by 6^d, not an alternative character normalization.', '',
        '| d | Root-lattice image factors |','| --- | --- |']
for row in data['degrees']:
    scale=6**row['degree']
    assert all(v%scale==0 for v in row['smith_factors'])
    lines.append('| {} | {} |'.format(row['degree'],fmt([v//scale for v in row['smith_factors']])))
lines+=['','Regenerate with `python tools/render_e6_character_tables.py`.','']
(ROOT/'E6_CHARACTER_TABLES.md').write_text('\n'.join(lines),encoding='utf-8')
print('Rendered',len(data['degrees']),'complete degree comparisons.')
