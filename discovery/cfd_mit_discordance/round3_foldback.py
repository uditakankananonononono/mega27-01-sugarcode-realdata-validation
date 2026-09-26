#!/usr/bin/env python3
"""F1 round-3 novelty fold-back (judge redirections):
Exp A: leave-guide-group-out generalization (5-fold grouped by guide) - does |CFD-MIT|
       help predict log(readfrac) on completely unseen guides?
Exp B: within-guide permutation test (1000x) - is the mixed-model signal real or
       inflated by within-guide structure? Statistic: t of disc in guide-demeaned OLS.
Exp C: assay interaction - does disc x study explain the Tsai leave-one-study-out harm?
"""
import csv, pickle, math, json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats as sps

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
guides, per = {}, {}
for r in rows:
    if r['type'] == 'on-target': guides[r['name']] = r['seq'].upper()
    per.setdefault(r['name'], []).append(r)
data = []
for name, rs in per.items():
    g = guides.get(name)
    if not g or len(g) != 23: continue
    g20 = g[:20]; study = name.split('_')[0].split('/')[0]
    for r in rs:
        if r['type'] != 'off-target': continue
        s = r['seq'].upper()
        if len(s) != 23 or any(c not in 'ACGT' for c in s): continue
        nmm = sum(1 for a,b in zip(g20, s[:20]) if a != b)
        if nmm > 4 or nmm == 0: continue
        cfd = calc_cfd(g20, s[:20], s[-2:])
        if cfd is None: continue
        mit_s = calc_mit(g20, s[:20])/100.0
        data.append(dict(guide=name, study=study, cfd=cfd, mit=mit_s,
                         disc=abs(cfd-mit_s), y=math.log10(float(r['score'])+1e-6)))
df = pd.DataFrame(data)
out = {'n_rows': len(df), 'n_guides': int(df.guide.nunique())}
rng = np.random.default_rng(42)

# Exp A: 5-fold grouped-by-guide CV
guides_u = np.array(sorted(df.guide.unique()))
rng.shuffle(guides_u)
folds = np.array_split(guides_u, 5)
pb_all = np.zeros(len(df)); pd_all = np.zeros(len(df))
fold_rows = []
for k, fg in enumerate(folds):
    te = df.guide.isin(fg).values; tr = ~te
    b = smf.ols('y ~ cfd + mit', data=df[tr]).fit()
    d_ = smf.ols('y ~ cfd + mit + disc', data=df[tr]).fit()
    pb = b.predict(df[te]); pdx = d_.predict(df[te])
    pb_all[te] = pb; pd_all[te] = pdx
    fold_rows.append({'fold': k, 'n_test': int(te.sum()),
        'spearman_base': float(sps.spearmanr(df.y[te], pb)[0]),
        'spearman_disc': float(sps.spearmanr(df.y[te], pdx)[0]),
        'rmse_base': float(np.sqrt(np.mean((df.y[te]-pb)**2))),
        'rmse_disc': float(np.sqrt(np.mean((df.y[te]-pdx)**2)))})
out['expA_unseen_guide_cv'] = {
  'folds': fold_rows,
  'pooled_spearman_base': float(sps.spearmanr(df.y, pb_all)[0]),
  'pooled_spearman_disc': float(sps.spearmanr(df.y, pd_all)[0]),
  'pooled_rmse_base': float(np.sqrt(np.mean((df.y-pb_all)**2))),
  'pooled_rmse_disc': float(np.sqrt(np.mean((df.y-pd_all)**2)))}

# Exp B: within-guide permutation of disc, guide-demeaned (within) OLS t-stat
def within_t(d, disc_vals):
    X = d[['cfd','mit']].copy(); X['disc'] = disc_vals
    X['y'] = d.y.values; X['g'] = d.guide.values
    Xd = X.groupby('g').transform(lambda s: s - s.mean())
    yd = Xd['y'].values
    A = np.column_stack([Xd['cfd'], Xd['mit'], Xd['disc']])
    beta, res, *_ = np.linalg.lstsq(A, yd, rcond=None)
    resid = yd - A@beta
    dof = len(yd) - A.shape[1] - d.guide.nunique()
    sigma2 = (resid@resid)/dof
    cov = sigma2*np.linalg.inv(A.T@A)
    return beta[2]/np.sqrt(cov[2,2])
obs_t = within_t(df, df['disc'].values)
null_ts = []
disc_arr = df['disc'].values.copy()
garr = df['guide'].values
for i in range(1000):
    perm = disc_arr.copy()
    for g in guides_u:
        m = garr == g
        perm[m] = rng.permutation(perm[m])
    null_ts.append(within_t(df, perm))
null_ts = np.array(null_ts)
out['expB_within_guide_permutation'] = {
  'observed_t': float(obs_t),
  'null_mean': float(null_ts.mean()), 'null_sd': float(null_ts.std()),
  'permutation_p_2sided': float((np.abs(null_ts) >= abs(obs_t)).mean()),
  'n_perm': 1000}

# Exp C: disc x study interaction (guide random intercept)
try:
    mi = smf.mixedlm('y ~ cfd + mit + disc*C(study)', data=df, groups=df['guide']).fit(reml=False, method='nm', maxiter=1500)
    mn = smf.mixedlm('y ~ cfd + mit + disc + C(study)', data=df, groups=df['guide']).fit(reml=False, method='nm', maxiter=1500)
    lr = 2*(mi.llf-mn.llf); k = len(mi.params)-len(mn.params)
    tsai_term = [p for p in mi.params.index if 'disc:C(study)' in p and 'Tsai' in p]
    out['expC_assay_interaction'] = {
      'lrt_chi2': float(max(lr,0)), 'lrt_df': int(k),
      'lrt_p': float(1-sps.chi2.cdf(max(lr,0),k)),
      'tsai_interaction_coef': float(mi.params[tsai_term[0]]) if tsai_term else None,
      'tsai_interaction_p': float(mi.pvalues[tsai_term[0]]) if tsai_term else None,
      'converged': bool(mi.converged)}
except Exception as e:
    out['expC_assay_interaction'] = {'error': str(e)}

print(json.dumps(out, indent=1))
json.dump(out, open('/home/sandbox/wave1/logs/round3_foldback_result.json','w'), indent=1)
