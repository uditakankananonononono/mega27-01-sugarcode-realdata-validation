"""qsar_bench.morgan_fingerprint vs RDKit Morgan r=2/2048: pairwise Tanimoto rank agreement and nearest neighbours on 300 ChEMBL drugs."""
import sys, os, csv, random, json, numpy as np
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from rdkit import Chem, DataStructs, RDLogger; RDLogger.DisableLog('rdApp.*')
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
from sugarcode.modules.qsar_bench.chemistry import morgan_fingerprint
rows = list(csv.DictReader(open(os.path.expanduser('~/mega/m01/data/chembl_approved_1000.tsv')), delimiter='\t'))[:300]
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fr = [gen.GetFingerprint(Chem.MolFromSmiles(r['smiles'])) for r in rows]
fo = [morgan_fingerprint(r['smiles']) for r in rows]
def tan(x, y): inter = (x * y).sum(); return float(inter / (x.sum() + y.sum() - inter))
random.seed(0); pairs = random.sample([(i, j) for i in range(300) for j in range(i + 1, 300)], 5000)
a = [DataStructs.TanimotoSimilarity(fr[i], fr[j]) for i, j in pairs]; b = [tan(fo[i], fo[j]) for i, j in pairs]
nn = sum(max((DataStructs.TanimotoSimilarity(fr[i], fr[j]), j) for j in range(300) if j != i)[1] == max((tan(fo[i], fo[j]), j) for j in range(300) if j != i)[1] for i in range(300))
res = {'reference': 'RDKit Morgan radius 2, 2048 bits', 'n_molecules': 300, 'n_pairs': 5000, 'spearman_tanimoto': round(float(spearmanr(a, b).correlation), 4),
       'mean_abs_tanimoto_diff': round(float(np.mean(np.abs(np.array(a) - np.array(b)))), 4), 'nearest_neighbour_agree': nn,
       'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else ''}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_qsar_morgan.json'), 'w'), indent=1); print(res)
