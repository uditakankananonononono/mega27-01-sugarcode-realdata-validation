"""structural_biophysics.vina_score vs AutoDock Vina 1.2.7 score_only on carbon-only ring stacks (all atoms C_H, so element typing is exact)."""
import math, sys, os, json
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from vina import Vina
from sugarcode.modules.structural_biophysics.core import Atom, vina_score
os.makedirs('/tmp/vt', exist_ok=True); os.chdir('/tmp/vt')
ring = [(1.39 * math.cos(k * math.pi / 3), 1.39 * math.sin(k * math.pi / 3), 0.0) for k in range(6)]
def line(i, x, y, z, res, ch, t): return f"ATOM  {i:5d}  C   {res} {ch}   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 {t:<2s}\n"
rec = ring + [(x + 4.2, y, z) for x, y, z in ring]
open('rec.pdbqt', 'w').write(''.join(line(i + 1, *c, 'BEN', 'A', 'A') for i, c in enumerate(rec)))
rows = []
for dz in (3.4, 3.6, 3.8, 4.0, 4.3, 4.6, 5.0, 6.0):
    lig = [(x + 2.1, y + 1.2, z + dz) for x, y, z in ring]
    open('lig.pdbqt', 'w').write('ROOT\n' + ''.join(line(i + 1, *c, 'LIG', 'L', 'A') for i, c in enumerate(lig)) + 'ENDROOT\nTORSDOF 0\n')
    v = Vina(sf_name='vina', verbosity=0); v.set_receptor('rec.pdbqt'); v.set_ligand_from_file('lig.pdbqt')
    v.compute_vina_maps(center=[2.1, 1.2, 2.0], box_size=[24, 24, 24]); e = v.score()
    o = vina_score([Atom(i, 'C', 'BEN', 'A', 1, '', 'C', c) for i, c in enumerate(rec)], [Atom(i, 'C', 'LIG', 'L', 1, '', 'C', c) for i, c in enumerate(lig)])
    rows.append({'separation_A': dz, 'vina_inter': round(float(e[1]), 4), 'ours': o['score_kcal_mol']})
res = {'reference': 'AutoDock Vina 1.2.7 (Python bindings) score_only', 'rows': rows, 'max_abs_diff': max(abs(r['vina_inter'] - r['ours']) for r in rows),
       'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else '',
       'pre_fix': {'rows_3.6_4.0_4.6': [-1.5646, -0.8557, -0.3662], 'vina': [-2.31, -1.65, -0.634], 'causes': ['Bondi radii instead of Vina XS radii in the surface distance', 'sulfur counted hydrophobic', 'rotor term added (w*Nrot) instead of dividing by 1+w*Nrot']},
       'sources': ['https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/src/lib/atom_constants.h', 'https://raw.githubusercontent.com/ccsb-scripps/AutoDock-Vina/develop/src/lib/conf_independent.cpp'],
       'scope': 'rotor division verified from source, not numerically (no torsion tree in the test ligands); element typing cannot separate C_H from C_P or assign donors without hydrogens'}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_vina_score.json'), 'w'), indent=1); print(res['rows'], res['max_abs_diff'])
