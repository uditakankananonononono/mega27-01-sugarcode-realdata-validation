"""microbiome_exp alpha/beta diversity, PCoA and permutation differential abundance vs scipy/independent recomputation on synthetic cohorts."""
import sys, os, json, math
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import numpy as np
from scipy.stats import entropy
from scipy.spatial.distance import braycurtis
from sugarcode.modules.microbiome_exp.core import analyze_16s, cohort_analysis, differential_abundance
GENERA = ['Faecalibacterium','Bacteroides','Lactobacillus','Akkermansia','Escherichia','Bifidobacterium','Clostridium','Prevotella','Roseburia','Ruminococcus']
rng = np.random.default_rng(11)
table = {}
meta = {}
for i in range(12):
    a = rng.dirichlet(np.ones(len(GENERA)) * (0.4 if i % 2 else 1.2))
    n = int(rng.integers(5000, 20000))
    table[f'S{i}'] = {g: int(c) for g, c in zip(GENERA, rng.multinomial(n, a)) if c > 0}
    meta[f'S{i}'] = {'group': 'control' if i % 2 else 'patient'}
res = {}
# single-sample analyze_16s vs skbio
single_bad = 0
for s, counts in table.items():
    ids = sorted(counts); c = np.array([counts[g] for g in ids])
    got = analyze_16s(counts)['alpha_diversity']
    p1 = c / c.sum()
    if abs(got['shannon'] - round(float(entropy(p1)), 3)) > 1e-9: single_bad += 1
    if abs(got['simpson'] - round(1 - float(np.sum(p1 * p1)), 3)) > 1e-9: single_bad += 1
    if got['richness'] != int(np.sum(c > 0)): single_bad += 1
res['analyze16s_mismatches'] = single_bad
# cohort alpha + bray-curtis + PCoA vs skbio
co = cohort_analysis(table, meta)
taxa = co['taxa']; samples = co['samples']
C = np.array([[table[s].get(t, 0) for t in taxa] for s in samples])
alpha_bad = 0
for row, c in zip(co['alpha_diversity'], C):
    p1 = c / c.sum(); S = int(np.sum(c > 0))
    if abs(row['shannon'] - float(entropy(p1))) > 1e-9: alpha_bad += 1
    if abs(row['simpson'] - (1 - float(np.sum(p1 * p1)))) > 1e-9: alpha_bad += 1
    if row['richness'] != S: alpha_bad += 1
    if abs(row['pielou_evenness'] - (float(entropy(p1)) / math.log(S) if S > 1 else 0.0)) > 1e-9: alpha_bad += 1
res['cohort_alpha_mismatches'] = alpha_bad
P = C / C.sum(1, keepdims=True)
bc_ref = np.array([[0.0 if i == j else braycurtis(P[i], P[j]) for j in range(len(samples))] for i in range(len(samples))])
bc_got = np.array(co['bray_curtis'])
res['bray_curtis_max_abs_diff'] = float(np.max(np.abs(bc_ref - bc_got)))
# independent PCoA: Gower centering via explicit centering matrix, scipy eigh
import scipy.linalg
n_ = len(samples); J = np.eye(n_) - np.ones((n_, n_)) / n_
G = -0.5 * J @ (bc_ref ** 2) @ J
w, v = scipy.linalg.eigh(G)
res['pcoa_eig0_diff'] = float(abs(w[-1] - co['ordination']['eigenvalues'][0]))
res['pcoa_eig1_diff'] = float(abs(w[-2] - co['ordination']['eigenvalues'][1]))
# differential abundance: planted effect, convention and permutation-p sanity
table2 = dict(table); meta2 = dict(meta)
for s in samples[:6]:
    table2[s] = dict(table2[s]); table2[s]['Roseburia'] = table2[s].get('Roseburia', 0) + (3000 if meta2[s]['group'] == 'control' else 0)
d = differential_abundance(table2, meta2, 'group')
ros = next(x for x in d['results'] if x['taxon'] == 'Roseburia')
res['planted_effect'] = {'log2fc': round(ros['log2_fold_change'], 3), 'p': ros['permutation_p'],
 'numerator_group': ros['group_comparison'], 'sign_correct': bool(ros['log2_fold_change'] > 0 and ros['permutation_p'] < 0.05)}
# null case: no planted effect -> median p should not be tiny across taxa
d0 = differential_abundance(table, meta, 'group')
ps = [x['permutation_p'] for x in d0['results']]
res['null_median_p'] = float(np.median(ps)); res['null_min_p'] = min(ps)
res['sugarcode_commit'] = sys.argv[1] if len(sys.argv) > 1 else ''
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_microbiome.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
