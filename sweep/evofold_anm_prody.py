"""evofold_4d.anm_modes vs ProDy ANM Hessian (eigen-solved with numpy) on 30 PDB structures; B-factor correlation."""
import sys, os, glob, json, numpy as np
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import prody; prody.confProDy(verbosity='none')
from scipy.stats import pearsonr
from sugarcode.modules.evofold_4d.core import anm_modes
rows = []
for f in sorted(glob.glob(os.path.expanduser('~/mega/m01/data/pdb_files/*.pdb'))):
    st = prody.parsePDB(f); s = st.select('protein and name CA')
    if s is None: continue
    ch = s.getChids()[0]; s = s.select(f'chain {ch}')
    if s.numAtoms() < 30 or s.numAtoms() > 400: rows.append({'pdb': os.path.basename(f)[:4], 'skipped': int(s.numAtoms())}); continue
    C = s.getCoords(); a = prody.ANM(); a.buildHessian(C, cutoff=10.0, gamma=1.0)
    H = a.getHessian(); H = H.toarray() if hasattr(H, 'toarray') else np.asarray(H)
    w, v = np.linalg.eigh(H); nzi = [i for i in range(len(w)) if w[i] > 1e-6][:6]; w6 = w[nzi]; V = v[:, nzi]; nzero = int(sum(w <= 1e-6))
    pf = sum(((V[:, k].reshape(-1, 3)) ** 2).sum(1) / w6[k] for k in range(6)) / 6
    o = anm_modes(C.tolist(), 6, 10.0); ow = np.array([m['frequency'] ** 2 for m in o['modes']])
    b = s.getBetas()
    rows.append({'pdb': os.path.basename(f)[:4], 'chain': ch, 'n_ca': len(C), 'n_zero_modes': nzero, 'max_rel_eig_diff': float(np.max(np.abs(ow - w6) / w6)), 'max_freq_diff': float(np.max(np.abs(np.sqrt(w6) - np.array([m['frequency'] for m in o['modes']])))),
                 'fluct_r_vs_prody': float(pearsonr(pf, o['fluctuation_profile'])[0]),
                 'bfactor_r_prody': float(pearsonr(pf, b)[0]) if b.std() > 0 else None, 'bfactor_r_ours': float(pearsonr(o['fluctuation_profile'], b)[0]) if b.std() > 0 else None})
sc = [r for r in rows if 'skipped' not in r]
res = {'reference': 'ProDy ' + prody.__version__ + ' ANM Hessian, cutoff 10 A, gamma 1; eigenvectors via numpy.linalg.eigh (ProDy calcModes fails with this SciPy)',
       'n_structures': len(sc), 'n_freq_match_at_4dp': sum(r['max_freq_diff'] <= 5.001e-5 for r in sc), 'n_with_extra_zero_modes': sum(r['n_zero_modes'] > 6 for r in sc), 'min_fluct_r': min(r['fluct_r_vs_prody'] for r in sc),
       'median_bfactor_r': float(np.median([r['bfactor_r_ours'] for r in sc if r['bfactor_r_ours'] is not None])),
       'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else '', 'rows': rows}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_evofold_anm.json'), 'w'), indent=1)
print({k: v for k, v in res.items() if k != 'rows'})
