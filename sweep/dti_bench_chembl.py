"""Validate dti_bench against live ChEMBL: ingestion correctness + model behavior."""
import json, os, sys, hashlib, urllib.request
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import numpy as np
from sugarcode.modules.dti_bench.chembl_pairs import ChEMBLPairClient
from sugarcode.modules.dti_bench.core import (protein_features, pair_features, fit_dti,
    predict_dti, validate_dti, _tan, _seq_identity)
from sugarcode.modules.qsar_bench.chemistry import morgan_fingerprint
T = ['CHEMBL203','CHEMBL279','CHEMBL1862','CHEMBL2971','CHEMBL4282','CHEMBL2842']
ACC = {'CHEMBL203':'P00533','CHEMBL279':'P35968','CHEMBL1862':'P00519','CHEMBL2971':'O60674','CHEMBL4282':'P31749','CHEMBL2842':'P42345'}
c = ChEMBLPairClient('data/dti_bench/cache', offline=True)
data = {t: c.pairs([t], per_target=300) for t in T}
res = {'module': 'dti_bench', 'targets': {t: len(d['records']) for t, d in data.items()}}
# 1) sequence spot-check vs UniProt (independent live source)
def uniprot_seq(acc):
    f = f'data/dti_bench/{acc}.fasta'
    if not os.path.exists(f):
        with urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.fasta', timeout=60) as r:
            open(f, 'wb').write(r.read())
    return ''.join(l.strip() for l in open(f) if not l.startswith('>'))
seqchk = {}
for t, acc in ACC.items():
    got = c.target_sequence(t); ref = uniprot_seq(acc)
    seqchk[t] = {'accession': acc, 'module_len': len(got), 'uniprot_len': len(ref), 'exact': got == ref}
res['target_sequences_vs_uniprot'] = seqchk
# 2) median aggregation independently recomputed from the raw cached activity page
def raw_page(tid):
    return c._get("activity.json", {"target_chembl_id": tid, "standard_type__in": ",".join(("IC50","Ki","Kd")),
                  "pchembl_value__isnull": "false", "limit": 300})
agree = mism = 0
for t in T[:2]:
    acts = raw_page(t)['activities']
    vals = {}
    for a in acts:
        if str(a.get('standard_relation') or '').strip() == '=' and a.get('pchembl_value'):
            vals.setdefault(a['molecule_chembl_id'], []).append(float(a['pchembl_value']))
    mine = {m: (sorted(v)[len(v)//2] if len(v) % 2 else sum(sorted(v)[len(v)//2-1:len(v)//2+1])/2)
            for m, v in vals.items()}
    # module groups by canonical SMILES; compare via molecule ids in each record
    for rec in data[t]['records']:
        for m in rec['molecule_chembl_ids']:
            if m in mine and abs(mine[m] - rec['pchembl_median']) < 1e-9: agree += 1
            elif m in mine: mism += 1
res['median_aggregation'] = {'checked_molecules': agree + mism, 'agree': agree, 'mismatch': mism}
# 3) assemble pooled dataset, warm split model vs naive-mean baseline
rows = [r for t in T for r in data[t]['records']]
S = [r['smiles'] for r in rows]; Q = [r['sequence'] for r in rows]
Y = [r['pchembl_median'] for r in rows]; TID = [r['target_chembl_id'] for r in rows]
YRS = [r['year'] for r in rows]
v1 = validate_dti(S, Q, Y, TID, strategy='warm_random', seed=23)
v2 = validate_dti(S, Q, Y, TID, strategy='warm_random', seed=23)
det = v1['test_indices'] == v2['test_indices']
te = np.array(v1['test_indices']); tr = np.array(v1['train_indices'])
ym = float(np.mean(np.array(Y)[tr])); base_rmse = float(np.sqrt(np.mean((np.array(Y)[te] - ym)**2)))
res['warm_random'] = {'n_train': len(tr), 'n_test': len(te), 'model_rmse': v1['metrics']['rmse'],
    'model_r2': v1['metrics']['r2'], 'naive_mean_rmse': base_rmse, 'deterministic_split': det,
    'inside_ad': v1['support']['inside_ad'], 'outside_ad': v1['support']['outside_ad']}
# 4) cold-target split (one held-out kinase)
vc = validate_dti(S, Q, Y, TID, strategy='cold_target', seed=23)
res['cold_target'] = {'n_train': len(vc['train_indices']), 'n_test': len(vc['test_indices']),
    'rmse': vc['metrics']['rmse'], 'r2': vc['metrics']['r2'],
    'inside_ad': vc['support']['inside_ad'], 'outside_ad': vc['support']['outside_ad'],
    'test_targets': sorted({TID[i] for i in vc['test_indices']})}
# 5) time split
vt = validate_dti(S, Q, Y, TID, strategy='time', years=YRS, seed=23)
res['time_split'] = {'n_train': len(vt['train_indices']), 'n_test': len(vt['test_indices']),
    'rmse': vt['metrics']['rmse'], 'r2': vt['metrics']['r2'],
    'train_max_year': max(YRS[i] for i in vt['train_indices']),
    'test_min_year': min(YRS[i] for i in vt['test_indices'])}
# 6) Tanimoto rank agreement vs RDKit Morgan generator (bit identities differ by design;
#    the validated criterion from sweep_qsar_descriptors is rank agreement)
from rdkit import Chem
from rdkit.Chem import DataStructs, rdFingerprintGenerator
from scipy.stats import spearmanr
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=128)
uniq = sorted(set(S))[:60]
fr = [gen.GetFingerprint(Chem.MolFromSmiles(x)) for x in uniq]
fo = [morgan_fingerprint(x, n_bits=128, radius=2) for x in uniq]
pairs = [(i, j) for i in range(len(uniq)) for j in range(i+1, len(uniq))]
a = [DataStructs.TanimotoSimilarity(fr[i], fr[j]) for i, j in pairs]
b = [_tan(fo[i], fo[j]) for i, j in pairs]
rho = float(spearmanr(a, b).correlation)
res['tanimoto_vs_rdkit'] = {'n_molecules': len(uniq), 'n_pairs': len(pairs), 'spearman': round(rho, 4),
    'criterion': 'rank agreement (sugarcode hash is not bit-identical to RDKit by design)'}
# 7) protein_features independent recomputation (composition part)
s = c.target_sequence('CHEMBL203'); pf = protein_features(s)
comp_ok = all(abs(pf[i] - s.upper().count(a)/len(s)) < 1e-12 for i, a in enumerate('ACDEFGHIKLMNPQRSTVWY'))
res['protein_features_composition'] = {'match_20aa': comp_ok}
json.dump(res, open('benchmarks/sweep_dti_bench.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
