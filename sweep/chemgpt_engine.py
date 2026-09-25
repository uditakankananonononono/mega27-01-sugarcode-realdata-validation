"""Sweep: chemgpt_engine vs PubChem molecular weights/identities + live ChEMBL similarity + recomputation."""
import json, os, sys, time, urllib.parse, urllib.request
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
from sugarcode.modules.chemgpt_engine import core as C
CD = os.path.join(os.path.dirname(__file__), '..', 'data', 'pubchem')
os.makedirs(CD, exist_ok=True)
res = {'module': 'chemgpt_engine', 'sources': ['PubChem PUG-REST', 'ChEMBL similarity API (live)', 'recomputation']}
def pug(smiles):
    fn = os.path.join(CD, 'cg_' + urllib.parse.quote(smiles, safe='') + '.json')
    if not os.path.exists(fn):
        u = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'
             + urllib.parse.quote(smiles, safe='') + '/property/MolecularWeight,CanonicalSMILES/JSON')
        for k in range(4):
            try:
                urllib.request.urlretrieve(u, fn); break
            except Exception:
                time.sleep(2 * (k + 1))
        time.sleep(0.3)
    return json.load(open(fn))['PropertyTable']['Properties'][0]
CASES = [(['benzene', 'hydroxyl'], 'phenol', 94.11), (['benzene', 'amine'], 'aniline', 93.13),
         (['benzene', 'carboxyl'], 'benzoic acid', 122.12), (['benzene', 'fluorine'], 'fluorobenzene', 96.10),
         (['pyridine', 'amine'], 'aminopyridine', 94.11), (['imidazole', 'methyl'], 'methylimidazole', 82.10),
         (['benzene', 'ethyl_link', 'amine'], 'phenethylamine', 121.18)]
rows = []
for frags, name, pubmw in CASES:
    m = C.score_molecule(frags)
    p = pug(m['smiles'])
    rows.append({'name': name, 'fragments': frags, 'smiles': m['smiles'], 'mw_module': m['mw'],
                 'mw_pubchem': p['MolecularWeight'], 'mw_diff': round(abs(m['mw'] - float(p['MolecularWeight'])), 3),
                 'pubchem_canonical': p.get('CanonicalSMILES') or p.get('ConnectivitySMILES'), 'logP_module': m['logP'],
                 'lipinski_violations': m['lipinski_violations']})
res['pubchem_mw'] = rows
# score recomputation
sc = C.score_molecule(['benzene', 'carboxyl'])
exp = {'atoms': 9, 'logp': round(1.7 - 0.3, 2), 'hbd': 1, 'hba': 2, 'rot': 1,
       'logS': round(0.5 - 1.4 - 0.01 * (122.1 - 200) / 50, 2)}
res['score_recompute'] = {'module': {k: sc[k] for k in ('logP', 'logS', 'hbd', 'hba', 'rotatable')},
                          'expected': exp, 'ok': sc['logP'] == exp['logp'] and sc['logS'] == exp['logS']
                          and sc['hbd'] == 1 and sc['hba'] == 2 and sc['rotatable'] == 1}
# pareto front vs brute force
gen = C.generate(n=24, target_logp=2.5, seed=42)
import itertools
def nd(cands):
    out = []
    for a in cands:
        dom = False
        for b in cands:
            if b is a: continue
            if all(b['objectives'][k] >= a['objectives'][k] for k in a['objectives']) and any(b['objectives'][k] > a['objectives'][k] for k in a['objectives']):
                dom = True; break
        if not dom: out.append(a)
    return out
# brute-force check: replicate the candidate list with the same RNG, filter non-dominated, compare
import random as _r
rng = _r.Random(42)
rings = [f for f in C.FRAGMENTS if C.FRAGMENTS[f][5] == 0 and C.FRAGMENTS[f][1] >= 5]
groups = [f for f in C.FRAGMENTS if f not in rings]
cands = []
for _ in range(24):
    while True:
        frag = [rng.choice(rings)] + [rng.choice(groups) for _ in range(rng.randint(1, 3))]
        try:
            m = C.score_molecule(frag); break
        except ValueError:
            continue
    m['objectives'] = C._objectives(m, 2.5); cands.append(m)
bf = nd(cands)
res['pareto'] = {'front_size': gen['pareto_size'], 'brute_force_size': len(bf),
                 'matches_brute_force': sorted(id(x) for x in []) == [] and {tuple(sorted(c['fragments'])) for c in gen['pareto_front']} == {tuple(sorted(c['fragments'])) for c in bf},
                 'deterministic': C.generate(n=24, seed=42)['pareto_front'] == gen['pareto_front'],
                 'front_sorted_by_composite': all(gen['pareto_front'][i]['composite'] >= gen['pareto_front'][i+1]['composite'] for i in range(len(gen['pareto_front'])-1))}
# live ChEMBL similarity for phenol and for the top designed molecule
sim1 = C.similarity_check('c1ccccc1O', cutoff=85, limit=5)
res['chembl_phenol'] = {'status': sim1['status'], 'closest': sim1.get('closest'), 'novelty': sim1['novelty']}
top = gen['pareto_front'][0]
sim2 = C.similarity_check(top['smiles'], cutoff=80, limit=5)
res['chembl_designed'] = {'smiles': top['smiles'], 'status': sim2['status'], 'closest': sim2.get('closest'), 'novelty': sim2['novelty']}
# molecular graph: rings closed, edges match bonds
g = C.molecular_graph(['imidazole', 'fluorine'])
res['graph'] = {'smiles': g['smiles'], 'n_nodes': len(g['nodes']), 'n_edges': len(g['edges']),
                'imidazole_ring_atoms': 5, 'note': 'imidazole core has 5 ring atoms; ring closure must add ring edges'}
# conformer ensemble boltzmann weights
ce = C.conformer_ensemble(['benzene', 'hydroxyl'], n=10, seed=0)
import numpy as _np
if isinstance(ce, dict):
    res['conformer'] = {k: ce[k] for k in ce if k != 'conformers'}
else:
    ens = ce[-1] if ce and isinstance(ce[-1], dict) else {}
    ws = [x.get('boltzmann_weight') for x in ce] if isinstance(ce, list) else None
    res['conformer'] = {'type': 'list', 'n': len(ce), 'keys': sorted(ce[0]) if ce else [],
                        'weights_sum': round(sum(ws), 8) if ws and all(w is not None for w in ws) else None}
fn = os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_chemgpt_engine.json')
json.dump(res, open(fn, 'w'), indent=1, default=str)
print(json.dumps(res, default=str)[:5000])
