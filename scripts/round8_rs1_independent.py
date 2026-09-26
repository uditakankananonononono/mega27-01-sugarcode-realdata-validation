#!/usr/bin/env python3
"""RS1 INDEPENDENT mathematical reimplementation - paper-only derivation (v2, corrected).

Source: Doench 2014 Nat Biotechnol PMC4262738 (methods text fetched from PMC; no source code).
- Activity direction confirmed from paper text: cells were FACS-sorted for the NEGATIVE
  population (lost the surface marker). Percent Peptide = retained signal, so
  HIGH activity = LOW Percent Peptide. (v1 had this backwards.)
- Assay restriction: 9 surface-marker genes, nodrug rows only (drug-screen rows use a
  different assay where Percent Peptide is not meaningful).
- Features (586): 120 one-hot single-nt (30x4) + 464 one-hot dinucleotide (29x16)
  + 2 GC-count deviation features (guide GC count vs 10, below/above).
- Feature selection: L1-regularized linear SVM, penalty chosen by gene-grouped nested CV.
- Classifier: logistic regression on top-quintile (lowest Percent Peptide) labels.
- Evaluation: gene-grouped holdout; per-gene top-quintile overlap and Spearman vs activity.
- Reference: sugarcode crispr_opt.doench2014_ontarget (CRISPOR-parameter RS1) as the
  reference implementation; 'predictions' column (RS2-generation, ships with the 2016
  data bundle) shown only as secondary context.
"""
import csv, json, sys
sys.path.insert(0, '/home/sandbox/sugarcode-ai/src')
from sugarcode.modules.crispr_opt import core as co
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from scipy import stats as sps

SURF = {'CD5','CD13','CD15','CD28','CD33','CD43','CD45','H2-K','THY1'}
rows = list(csv.DictReader(open('/home/sandbox/mega27-01-sugarcode-realdata-validation/data/doench2016/FC_plus_RES_withPredictions.csv')))
sub = [r for r in rows if r['Target gene'] in SURF and r['drug'] == 'nodrug']
mers = np.array([r['30mer'].strip().upper() for r in sub])
genes = np.array([r['Target gene'] for r in sub])
pep = np.array([float(r['Percent Peptide']) for r in sub])   # retained signal; LOW = active
act = -pep                                                   # activity, higher = better
pub = np.array([float(r['predictions']) for r in sub])
ref = np.array([co.doench2014_ontarget(m) for m in mers])    # reference RS1 implementation
n = len(mers)

BASE = {'A':0,'C':1,'G':2,'T':3}
def featurize(m):
    x = np.zeros(586)
    b = [BASE[c] for c in m]
    for i in range(30):
        x[i*4 + b[i]] = 1.0
    for i in range(29):
        x[120 + i*16 + b[i]*4 + b[i+1]] = 1.0
    gc = sum(1 for c in m[4:24] if c in 'GC')
    x[584] = max(0, 10 - gc); x[585] = max(0, gc - 10)
    return x
X = np.stack([featurize(m) for m in mers])

y = np.zeros(n)
for g in set(genes):
    idx = np.where(genes == g)[0]
    y[idx] = (act[idx] >= np.quantile(act[idx], 0.8)).astype(float)

Cs = [0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25]
gkf = GroupKFold(n_splits=9)
best_c, best_auc, cv_perf = None, -1, {}
for C in Cs:
    aucs = [roc_auc_score(y[te], LinearSVC(C=C, penalty='l1', dual=False, max_iter=8000).fit(X[tr], y[tr]).decision_function(X[te]))
            for tr, te in gkf.split(X, y, genes) if y[tr].sum() > 0]
    cv_perf[str(C)] = float(np.mean(aucs))
    if cv_perf[str(C)] > best_auc: best_auc, best_c = cv_perf[str(C)], C

svm = LinearSVC(C=best_c, penalty='l1', dual=False, max_iter=8000).fit(X, y)
sel = np.where(np.abs(svm.coef_[0]) > 1e-9)[0]
Xs = X[:, sel]

# held-out scores: logistic retrained per fold on selected features
oof = np.zeros(n)
for tr, te in gkf.split(Xs, y, genes):
    lr = LogisticRegression(max_iter=2000).fit(Xs[tr], y[tr])
    oof[te] = lr.predict_proba(Xs[te])[:, 1]

def quintile_overlap(score, activity, genes):
    ovs, sprs = {}, {}
    for g in sorted(set(genes)):
        idx = np.where(genes == g)[0]
        k = max(1, len(idx)//5)
        top_s = set(idx[np.argsort(-score[idx])[:k]])
        top_a = set(idx[np.argsort(-activity[idx])[:k]])
        ovs[g] = round(len(top_s & top_a)/k, 3)
        sprs[g] = round(float(sps.spearmanr(score[idx], activity[idx])[0]), 3)
    return ovs, sprs

ov_mine, sp_mine = quintile_overlap(oof, act, genes)
ov_ref,  sp_ref  = quintile_overlap(ref, act, genes)
ov_pub,  sp_pub  = quintile_overlap(pub, act, genes)

result = {
 'experiment': 'round8_rs1_independent_reimplementation_v2',
 'n': int(n), 'genes': sorted(set(genes)),
 'direction_fix': 'Percent Peptide = retained signal; LOW = active (paper: negative population isolated by FACS)',
 'feature_selection_cv_auc': cv_perf, 'chosen_C': best_c,
 'n_features_selected': int(len(sel)),
 'heldout_spearman_vs_activity_overall': round(float(sps.spearmanr(oof, act)[0]), 4),
 'ref_rs1_spearman_vs_activity_overall': round(float(sps.spearmanr(ref, act)[0]), 4),
 'heldout_spearman_vs_ref_rs1': round(float(sps.spearmanr(oof, ref)[0]), 4),
 'per_gene_topquintile_overlap': {'mine': ov_mine, 'ref_rs1': ov_ref, 'pub_col': ov_pub},
 'per_gene_spearman_vs_activity': {'mine': sp_mine, 'ref_rs1': sp_ref, 'pub_col': sp_pub},
 'mean_overlap': {'mine': round(float(np.mean(list(ov_mine.values()))),3),
                  'ref_rs1': round(float(np.mean(list(ov_ref.values()))),3),
                  'pub_col': round(float(np.mean(list(ov_pub.values()))),3),
                  'chance': 0.2, 'paper_claim': 0.8},
 'mean_pergene_spearman': {'mine': round(float(np.mean(list(sp_mine.values()))),3),
                           'ref_rs1': round(float(np.mean(list(sp_ref.values()))),3),
                           'pub_col': round(float(np.mean(list(sp_pub.values()))),3)},
}
print(json.dumps(result, indent=1))
open('/home/sandbox/wave1/logs/round8_rs1_independent.json','w').write(json.dumps(result, indent=1))
