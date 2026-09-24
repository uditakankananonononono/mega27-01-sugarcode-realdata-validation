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

# ---- final model on all labeled data, then score every VUS (non-train genes)
net = Net(Fz.shape[1]); opt = torch.optim.Adam(net.parameters(), 1e-3, weight_decay=1e-4)
Xt, Ft, yt = torch.tensor(X), torch.tensor(Fz), torch.tensor(y, dtype=torch.float32)
lossf = nn.BCEWithLogitsLoss(pos_weight=torch.tensor((1 - y).sum() / y.sum(), dtype=torch.float32))
for ep in range(12):
    net.train(); perm = torch.randperm(len(yt))
    for i in range(0, len(perm), 128):
        b = perm[i:i + 128]; opt.zero_grad(); lossf(net(Xt[b], Ft[b]), yt[b]).backward(); opt.step()
net.eval(); torch.save(net.state_dict(), 'cnn_final.pt')
V = [r for r in csv.DictReader(open('pwm_scores.tsv'), delimiter='\t') if r['train_gene'] == '0' and r['cls'] == 'VUS']
out = open('vus_cnn.tsv', 'w'); out.write('vid\tgene\tname\tsite\tk\tdelta\tcnn\n')
for i in range(0, len(V), 512):
    ch = V[i:i + 512]
    xb = torch.tensor(np.stack([np.concatenate([onehot(win[r['vid']]['ref_ctx']), onehot(win[r['vid']]['alt_ctx'])]) for r in ch]))
    fb = torch.tensor(np.array([feats(r) for r in ch], np.float32))
    with torch.no_grad(): pr = torch.sigmoid(net(xb, fb)).numpy()
    for r, q in zip(ch, pr):
        out.write(f"{r['vid']}\t{r['gene']}\t{win[r['vid']]['name']}\t{r['site']}\t{r['k']}\t{r['delta']}\t{q:.5f}\n")
out.close(); print('scored', len(V))
