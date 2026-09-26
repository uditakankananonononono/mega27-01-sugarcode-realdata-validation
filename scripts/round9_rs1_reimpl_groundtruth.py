#!/usr/bin/env python3
"""RS1 independent reimplementation v3 - against AUTHORITATIVE ground truth.
Activity + published scores: Doench 2014 Supplementary Table 7 (MOESM8, downloaded
from nature.com media.springernature.com, 2026-09-26). Activity = within-gene
percent-rank in marker-negative population. n=1,841 sgRNAs, 30nt Expanded Sequence.
Same paper-only pipeline as v2 (586 one-hot features, L1-SVM selection w/ gene-grouped
nested CV, logistic on top-quintile labels, 9-gene holdout). Also: shipped
crispr_opt.doench2014_ontarget verified against the published final-model scores.
"""
import sys, json
sys.path.insert(0, '/home/sandbox/sugarcode-ai/src')
from sugarcode.modules.crispr_opt import core as co
import openpyxl, numpy as np
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from scipy import stats as sps

wb = openpyxl.load_workbook('/home/sandbox/wave1/doench2014_suppl/MOESM8_Table7.xlsx', read_only=True)
rows = [r for r in list(wb.active.iter_rows(values_only=True))[2:] if r and r[0]]
mers  = np.array([str(r[1]).strip().upper() for r in rows])
genes = np.array([str(r[4]) for r in rows])
act   = np.array([float(r[7]) for r in rows])
pub   = np.array([float(r[8]) for r in rows])
shipped = np.array([co.doench2014_ontarget(m) for m in mers])
n = len(mers)

BASE = {'A':0,'C':1,'G':2,'T':3}
def featurize(m):
    x = np.zeros(586)
    b = [BASE[c] for c in m]
    for i in range(30): x[i*4 + b[i]] = 1.0
    for i in range(29): x[120 + i*16 + b[i]*4 + b[i+1]] = 1.0
    gc = sum(1 for c in m[4:24] if c in 'GC')
    x[584] = max(0, 10-gc); x[585] = max(0, gc-10)
    return x
X = np.stack([featurize(m) for m in mers])
y = np.zeros(n)
for g in set(genes):
    idx = np.where(genes==g)[0]
    y[idx] = (act[idx] >= np.quantile(act[idx], 0.8)).astype(float)

Cs = [0.001,0.0025,0.005,0.01,0.025,0.05,0.1,0.25]
gkf = GroupKFold(n_splits=9)
best_c, best_auc, cv_perf = None, -1, {}
for C in Cs:
    aucs = [roc_auc_score(y[te], LinearSVC(C=C,penalty='l1',dual=False,max_iter=8000).fit(X[tr],y[tr]).decision_function(X[te]))
            for tr,te in gkf.split(X,y,genes) if y[tr].sum()>0]
    cv_perf[str(C)] = round(float(np.mean(aucs)),4)
    if cv_perf[str(C)] > best_auc: best_auc, best_c = cv_perf[str(C)], C
svm = LinearSVC(C=best_c,penalty='l1',dual=False,max_iter=8000).fit(X,y)
sel = np.where(np.abs(svm.coef_[0])>1e-9)[0]
oof = np.zeros(n)
for tr,te in gkf.split(X[:,sel],y,genes):
    oof[te] = LogisticRegression(max_iter=2000).fit(X[np.ix_(tr,sel)],y[tr]).predict_proba(X[np.ix_(te,sel)])[:,1]

per_gene, per_gene_pub, per_gene_ship = {}, {}, {}
for g in sorted(set(genes)):
    idx = np.where(genes==g)[0]
    if len(idx)<30: continue
    per_gene[g] = round(float(sps.spearmanr(oof[idx],act[idx])[0]),3)
    per_gene_pub[g] = round(float(sps.spearmanr(pub[idx],act[idx])[0]),3)
    per_gene_ship[g] = round(float(sps.spearmanr(shipped[idx],act[idx])[0]),3)

result = {
 'experiment':'round9_rs1_reimpl_groundtruth',
 'source':'Doench 2014 Suppl Table 7 (MOESM8, nature.com media.springernature.com)',
 'n': int(n),
 'shipped_vs_published_model': {'overall_spearman': round(float(sps.spearmanr(shipped,pub)[0]),4),
                                'overall_pearson': round(float(sps.pearsonr(shipped,pub)[0]),4)},
 'published_vs_activity_overall': round(float(sps.spearmanr(pub,act)[0]),4),
 'shipped_vs_activity_overall': round(float(sps.spearmanr(shipped,act)[0]),4),
 'reimpl_heldout_vs_activity_overall': round(float(sps.spearmanr(oof,act)[0]),4),
 'reimpl_heldout_vs_published_overall': round(float(sps.spearmanr(oof,pub)[0]),4),
 'feature_selection_cv_auc': cv_perf, 'chosen_C': best_c, 'n_features_selected': int(len(sel)),
 'per_gene_spearman_vs_true_activity': {'reimpl_heldout': per_gene, 'published_model': per_gene_pub, 'shipped_rs1': per_gene_ship},
 'mean_per_gene': {'reimpl_heldout': round(float(np.mean(list(per_gene.values()))),3),
                   'published_model': round(float(np.mean(list(per_gene_pub.values()))),3),
                   'shipped_rs1': round(float(np.mean(list(per_gene_ship.values()))),3)},
}
print(json.dumps(result, indent=1))
open('/home/sandbox/wave1/logs/round9_rs1_reimpl_groundtruth.json','w').write(json.dumps(result, indent=1))
