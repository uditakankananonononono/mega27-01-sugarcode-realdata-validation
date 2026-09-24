"""molecule_eval: graph-key dedup vs RDKit canonical SMILES (stereo removed), randomized-SMILES invariance, validity vs RDKit, internal diversity vs RDKit Morgan, on ChEMBL approved drugs."""
import sys, os, csv, json, random, numpy as np
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from rdkit import Chem, DataStructs, RDLogger; RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import rdFingerprintGenerator
from sugarcode.modules.molecule_eval.core import _canonical_graph_key, evaluate_library
rows = list(csv.DictReader(open(os.path.expanduser('~/mega/m01/data/chembl_approved_1000.tsv')), delimiter='\t'))
random.seed(1)
ok = []; parse_fail = []
for r in rows:
    m = Chem.MolFromSmiles(r['smiles'])
    try: k = _canonical_graph_key(r['smiles'])
    except Exception as e: parse_fail.append({'id': r.get('chembl_id'), 'err': str(e)[:80], 'rdkit_ok': m is not None}); continue
    if m is None: continue
    Chem.RemoveStereochemistry(m); ok.append((r.get('chembl_id'), r['smiles'], k, Chem.MolToSmiles(m)))
# invariance: 3 random SMILES per molecule
var_fail = []
for cid, s, k, can in ok:
    m = Chem.MolFromSmiles(s)
    for _ in range(3):
        rs = Chem.MolToSmiles(m, doRandom=True, isomericSmiles=False)
        try: k2 = _canonical_graph_key(rs)
        except Exception as e: var_fail.append({'id': cid, 'random': rs, 'err': str(e)[:60]}); break
        if k2 != k: var_fail.append({'id': cid, 'random': rs, 'orig': s}); break
# collisions / splits vs RDKit canonical
from collections import defaultdict
byk = defaultdict(set); byc = defaultdict(set)
for cid, s, k, can in ok: byk[k].add(can); byc[can].add(k)
coll = [sorted(v) for v in byk.values() if len(v) > 1]; split = [c for c, v in byc.items() if len(v) > 1]
# diversity vs RDKit on first 300 unique
sub = [s for _, s, _, _ in ok[:300]]
ev = evaluate_library(sub, n_bits=2048, max_pairs=10**9)
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fr = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in sub]
d = [1 - DataStructs.TanimotoSimilarity(fr[i], fr[j]) for i in range(len(fr)) for j in range(i + 1, len(fr))]
res = {'reference': 'RDKit ' + Chem.rdBase.rdkitVersion + ' canonical SMILES (stereo removed) and Morgan r2/2048', 'n_input': len(rows),
  'sugarcode_parse_fail': len(parse_fail), 'parse_fail_rdkit_ok': sum(x['rdkit_ok'] for x in parse_fail), 'parse_fail_examples': parse_fail[:8],
  'n_compared': len(ok), 'random_smiles_tested': 3, 'invariance_failures': len(var_fail), 'invariance_examples': var_fail[:6],
  'key_collisions_distinct_rdkit': len(coll), 'collision_examples': coll[:6], 'rdkit_same_key_split': len(split), 'split_examples': split[:6],
  'unique_rdkit': len(byc), 'unique_keys': len(byk),
  'diversity_n': len(sub), 'internal_diversity_sugarcode': round(ev['internal_diversity']['estimate'], 4), 'internal_diversity_rdkit': round(float(np.mean(d)), 4),
  'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else ''}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_molecule_eval.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if 'examples' not in k}, indent=0)); print(res['parse_fail_examples'][:4]); print(res['invariance_examples'][:3]); print(res['collision_examples'][:3]); print(res['split_examples'][:3])
