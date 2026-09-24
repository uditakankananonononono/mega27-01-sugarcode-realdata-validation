"""Gene-grouped 5-fold CV: splice-region SNV P vs B.
Models on identical folds: (a) PWM delta (deepsplice), (b) logistic on
PWM ref/alt/delta + position one-hot, (c) CNN on ref+alt one-hot 81-nt
context + the (b) features. Train genes of the PWM harvest are excluded."""
import csv, json, random, hashlib, sys
import numpy as np, torch, torch.nn as nn
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
torch.manual_seed(7); np.random.seed(7); random.seed(7)
torch.set_num_threads(2)
def auc(y, s):
    y = np.asarray(y); s = np.asarray(s); rk = rankdata(s); p = y.sum(); n = len(y) - p
    return float((rk[y == 1].sum() - p * (p + 1) / 2) / (p * n))
win = {r['vid']: r for r in csv.DictReader(open('windows.tsv'), delimiter='\t')}
rows = [r for r in csv.DictReader(open('pwm_scores.tsv'), delimiter='\t') if r['train_gene'] == '0' and r['cls'] in ('P', 'B')]
rng = random.Random(7)
data = []
for site in ('donor', 'acceptor'):
    P = [r for r in rows if r['site'] == site and r['cls'] == 'P']
    B = [r for r in rows if r['site'] == site and r['cls'] == 'B']
    B = rng.sample(B, min(len(B), 4 * len(P)))
    data += P + B
KS = [3, 4, 5, 6] + list(range(-14, -2))
def feats(r):
    f = [float(r['ref_s']), float(r['alt_s']), float(r['delta']), 1.0 if r['site'] == 'donor' else 0.0]
    f += [1.0 if int(r['k']) == k else 0.0 for k in KS]
    return f
M = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
def onehot(s):
    a = np.zeros((4, len(s)), np.float32)
    for i, c in enumerate(s):
        if c in M: a[M[c], i] = 1
    return a
X = np.stack([np.concatenate([onehot(win[r['vid']]['ref_ctx']), onehot(win[r['vid']]['alt_ctx'])]) for r in data])
Fz = np.array([feats(r) for r in data], np.float32)
y = np.array([1 if r['cls'] == 'P' else 0 for r in data])
genes = [r['gene'] for r in data]
fold = np.array([int(hashlib.md5(g.encode()).hexdigest(), 16) % 5 for g in genes])
class Net(nn.Module):
    def __init__(s, nf):
        super().__init__()
        s.c = nn.Sequential(nn.Conv1d(8, 64, 9, padding=4), nn.ReLU(), nn.Dropout(0.2),
                            nn.Conv1d(64, 64, 7, padding=3), nn.ReLU(), nn.AdaptiveMaxPool1d(1))
        s.h = nn.Sequential(nn.Linear(64 + nf, 32), nn.ReLU(), nn.Linear(32, 1))
    def forward(s, x, f): return s.h(torch.cat([s.c(x).squeeze(-1), f], 1)).squeeze(-1)
pred = {m: np.zeros(len(y)) for m in ('pwm', 'logit', 'cnn')}
for k in range(5):
    tr, te = fold != k, fold == k
    pred['pwm'][te] = -Fz[te, 2]
    lr = LogisticRegression(max_iter=2000).fit(Fz[tr], y[tr]); pred['logit'][te] = lr.predict_proba(Fz[te])[:, 1]
    net = Net(Fz.shape[1]); opt = torch.optim.Adam(net.parameters(), 1e-3, weight_decay=1e-4)
    Xt, Ft, yt = torch.tensor(X[tr]), torch.tensor(Fz[tr]), torch.tensor(y[tr], dtype=torch.float32)
    pw = torch.tensor((1 - y[tr]).sum() / max(1, y[tr].sum()), dtype=torch.float32)
    lossf = nn.BCEWithLogitsLoss(pos_weight=pw)
    for ep in range(12):
        net.train(); perm = torch.randperm(len(yt))
        for i in range(0, len(perm), 128):
            b = perm[i:i + 128]; opt.zero_grad(); l = lossf(net(Xt[b], Ft[b]), yt[b]); l.backward(); opt.step()
    net.eval()
    with torch.no_grad(): pred['cnn'][te] = torch.sigmoid(net(torch.tensor(X[te]), torch.tensor(Fz[te]))).numpy()
    print('fold', k, {m: round(auc(y[te], pred[m][te]), 4) for m in pred}, flush=True)
site = np.array([r['site'] for r in data]); kk = np.array([int(r['k']) for r in data])
res = {'n': int(len(y)), 'n_pos': int(y.sum()), 'overall': {m: round(auc(y, pred[m]), 4) for m in pred}}
for s in ('donor', 'acceptor'):
    res[s] = {m: round(auc(y[site == s], pred[m][site == s]), 4) for m in pred}
# within-position AUC (removes the position prior), weighted by positives
wp = {m: [] for m in pred}
for k in sorted(set(kk)):
    sel = kk == k
    if y[sel].sum() >= 10 and (1 - y[sel]).sum() >= 10:
        for m in pred: wp[m].append((auc(y[sel], pred[m][sel]), y[sel].sum()))
res['within_position_weighted'] = {m: round(sum(a * w for a, w in v) / sum(w for _, w in v), 4) for m, v in wp.items()}
print(json.dumps(res, indent=1)); json.dump(res, open('cv_results.json', 'w'), indent=1)
np.save('cv_pred.npy', np.stack([pred['pwm'], pred['logit'], pred['cnn'], y]))
