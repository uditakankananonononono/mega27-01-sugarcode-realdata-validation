"""structural_biophysics.shrake_rupley vs FreeSASA (Lee-Richards, same atoms and radii) on real PDB files."""
import sys, os, glob, json, time, numpy as np
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import freesasa
from scipy.stats import pearsonr
from sugarcode.modules.structural_biophysics.core import parse_pdb, shrake_rupley, VDW_RADII
files = sorted(glob.glob(os.path.expanduser('~/mega/m01/data/pdb_files/*.pdb')), key=os.path.getsize)[:int(sys.argv[1]) if len(sys.argv) > 1 else 6]
rows = []
for f in files:
    atoms = parse_pdb(open(f).read(), include_hetero=False)
    t = time.time(); o = shrake_rupley(atoms); dt = time.time() - t
    xyz = np.array([a.xyz for a in atoms], dtype=float).ravel().tolist()
    radii = [VDW_RADII.get(a.element, 1.70) for a in atoms]
    p = freesasa.Parameters({'algorithm': freesasa.LeeRichards, 'probe-radius': 1.4, 'n-slices': 100})
    r = freesasa.calcCoord(xyz, radii, p)
    fs = np.array([r.atomArea(i) for i in range(len(atoms))]); ours = np.array(o['atom_A2'])
    res_f = {}; [res_f.__setitem__(a.residue_id, res_f.get(a.residue_id, 0) + fs[i]) for i, a in enumerate(atoms)]
    rk = list(o['residue_A2']); rf = np.array([res_f[k] for k in rk]); ro = np.array([o['residue_A2'][k] for k in rk])
    rows.append({'pdb': os.path.basename(f)[:4], 'n_atoms': len(atoms), 'total_ours': o['total_A2'], 'total_freesasa': round(r.totalArea(), 3),
                 'total_rel_diff': round(abs(o['total_A2'] - r.totalArea()) / r.totalArea(), 5), 'atom_r': round(float(pearsonr(ours, fs)[0]), 5),
                 'residue_r': round(float(pearsonr(ro, rf)[0]), 5), 'residue_mean_abs_diff_A2': round(float(np.mean(np.abs(ro - rf))), 3), 'sec': round(dt, 1)})
    print(rows[-1], flush=True)
res = {'reference': 'FreeSASA ' + getattr(freesasa, '__version__', '') + ' Lee-Richards (100 slices), probe 1.4 A, module radii', 'rows': rows,
       'max_total_rel_diff': max(r['total_rel_diff'] for r in rows), 'min_residue_r': min(r['residue_r'] for r in rows)}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_sasa_freesasa_recheck14.json'), 'w'), indent=1)
print({k: v for k, v in res.items() if k != 'rows'})
