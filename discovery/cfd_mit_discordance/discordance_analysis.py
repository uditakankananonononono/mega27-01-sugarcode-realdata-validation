#!/usr/bin/env python3
"""F1 discovery: does CFD-vs-MIT discordance carry biological information about
measured off-target cleavage? Data: 7-study experimental set (crisporPaper offtargets.tsv).
Method: per-row oracle CFD + MIT scores; labels from measured read fractions;
leave-one-guide-out comparison of single scores vs a 2-feature model."""
import csv, pickle, math, random
from collections import defaultdict

ORACLE = '/home/sandbox/wave1/oracle'
mm_scores = pickle.load(open(f'{ORACLE}/mismatch_score.pkl','rb'))
pam_scores = pickle.load(open(f'{ORACLE}/pam_scores.pkl','rb'))
_comp = {'A':'T','C':'G','G':'C','T':'A','U':'A'}
def calc_cfd(wt, sg, pam2):
    score = 1.0
    sg = sg.upper().replace('T','U'); wt = wt.upper().replace('T','U')
    for i in range(len(sg)):
        if wt[i] != sg[i]:
            k = 'r'+wt[i]+':d'+_comp[sg[i]]+','+str(i+1)
            if k not in mm_scores: return None
            score *= mm_scores[k]
    if pam2 not in pam_scores: return None
    return score * pam_scores[pam2]
hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]
def calc_mit(s1, s2):
    dists, mmCount, lastMmPos, score1 = [], 0, None, 1.0
    for pos in range(20):
        if s1[pos] != s2[pos]:
            mmCount += 1
            if lastMmPos is not None: dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    score2 = 1.0 if mmCount < 2 else 1.0/(((19-sum(dists)/len(dists))/19.0)*4+1)
    score3 = 1.0 if mmCount == 0 else 1.0/(mmCount**2)
    return score1*score2*score3*100

rows = list(csv.DictReader(open('/home/sandbox/mega27-01-sugarcode-realdata-validation/data/crispor_offtargets_7studies.tsv'), delimiter='\t'))
groups = defaultdict(list)
guides = {}
for r in rows:
    if r['type'] == 'on-target':
        guides[r['name']] = r['seq'].upper()
    groups[r['name']].append(r)

data = []
skipped = defaultdict(int)
for name, rs in groups.items():
    g = guides.get(name)
    if not g or len(g) != 23: skipped['no_guide23'] += 1; continue
    g20 = g[:20]
    for r in rs:
        if r['type'] != 'off-target': continue
        s = r['seq'].upper()
        if len(s) != 23 or any(c not in 'ACGT' for c in s): skipped['bad_seq'] += 1; continue
        nmm = sum(1 for a,b in zip(g20, s[:20]) if a != b)
        if nmm > 4: skipped['>4mm'] += 1; continue
        cfd = calc_cfd(g20, s[:20], s[-2:])
        if cfd is None: skipped['no_pam'] += 1; continue
        mit = calc_mit(g20, s[:20])
        data.append(dict(name=name, cfd=cfd, mit=mit, score=float(r['score']), nmm=nmm))
print(f'usable off-targets: {len(data)} across {len(set(d["name"] for d in data))} guides; skipped {dict(skipped)}')

POS = 0.001  # read fraction > 0.1% = detected cleavage (CRISPOR-paper style threshold)
for d in data:
    d['pos'] = d['score'] > POS
npos = sum(d['pos'] for d in data)
print(f'positives (readfrac>{POS}): {npos}/{len(data)}')

def auc(scores, labels):
    pos = [s for s,l in zip(scores,labels) if l]; neg = [s for s,l in zip(scores,labels) if not l]
    if not pos or not neg: return None
    w = sum(1.0 if p>n else 0.5 if p==n else 0.0 for p in pos for n in neg)
    return w/(len(pos)*len(neg))

# single-score AUROCs (global, then LOGO mean)
labels = [d['pos'] for d in data]
print('global AUROC  CFD:', round(auc([d['cfd'] for d in data], labels), 4),
      ' MIT:', round(auc([d['mit'] for d in data], labels), 4))

# 2-feature logistic, leave-one-guide-out
import numpy as np
def fit_logreg(X, y, epochs=400, lr=0.5):
    Xb = np.hstack([np.ones((len(X),1)), X])
    w = np.zeros(Xb.shape[1])
    for _ in range(epochs):
        p = 1/(1+np.exp(-Xb@w))
        w += lr * (Xb.T @ (y - p)) / len(y)
    return w
names = sorted(set(d['name'] for d in data))
X = np.array([[math.log10(max(d['cfd'],1e-6)), math.log10(max(d['mit'],1e-6))] for d in data])
y = np.array(labels, dtype=float)
preds = np.zeros(len(data))
for gname in names:
    tr = [i for i,d in enumerate(data) if d['name'] != gname]
    te = [i for i,d in enumerate(data) if d['name'] == gname]
    if not te: continue
    w = fit_logreg(X[tr], y[tr])
    Xb = np.hstack([np.ones((len(te),1)), X[te]])
    preds[te] = 1/(1+np.exp(-Xb@w))
print('LOGO 2-feat logistic AUROC:', round(auc(list(preds), labels), 4))
# discordance-only feature: residual of log cfd on log mit (fit globally, simple)
b = np.polyfit(X[:,1], X[:,0], 1)
resid = X[:,0] - (b[0]*X[:,1] + b[1])
print('discordance residual |r|: pos mean', round(float(np.mean([abs(resid[i]) for i in range(len(data)) if labels[i]])),4),
      'neg mean', round(float(np.mean([abs(resid[i]) for i in range(len(data)) if not labels[i]])),4))
