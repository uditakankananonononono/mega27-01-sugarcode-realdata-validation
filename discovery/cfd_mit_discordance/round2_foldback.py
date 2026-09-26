#!/usr/bin/env python3
"""F1 round-2 novelty fold-back: continuous cleavage modeling (judge redirection).
Exp1: OLS log(readfrac) ~ CFD + MIT + (CFD-MIT)
Exp2: mixed-effects, random intercepts for study + guide (crossed via vcomp)
Exp3: discordance predictive value stratified by mismatch class
Exp4: leave-one-study-out (unseen-study) validation, baseline vs discordance model
Data + scores: same oracle pipeline as discordance_analysis.py."""
import csv, pickle, math, json
import numpy as np
import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy import stats as sps
import pandas as pd

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
    if r['type'] == 'on-target':
        guides[r['name']] = r['seq'].upper()
    per.setdefault(r['name'], []).append(r)
data = []
for name, rs in per.items():
    g = guides.get(name)
    if not g or len(g) != 23: continue
    g20 = g[:20]
    study = name.split('_')[0].split('/')[0]
    for r in rs:
        if r['type'] != 'off-target': continue
        s = r['seq'].upper()
        if len(s) != 23 or any(c not in 'ACGT' for c in s): continue
        nmm = sum(1 for a,b in zip(g20, s[:20]) if a != b)
        if nmm > 4 or nmm == 0: continue
        cfd = calc_cfd(g20, s[:20], s[-2:])
        if cfd is None: continue
        mit = calc_mit(g20, s[:20])
        rf = float(r['score'])
        mit_s = mit/100.0
        data.append(dict(guide=name, study=study, cfd=cfd, mit=mit_s,
                         disc=abs(cfd - mit_s), disc_sign=cfd - mit_s,
                         inter=cfd*mit_s, y=math.log10(rf + 1e-6), nmm=nmm))
df = pd.DataFrame(data)
print(f'rows={len(df)} guides={df.guide.nunique()} studies={df.study.nunique()}')
out = {'n_rows': len(df), 'n_guides': int(df.guide.nunique()), 'n_studies': int(df.study.nunique())}

# Exp1: OLS baseline vs discordance
m0 = smf.ols('y ~ cfd + mit', data=df).fit()
m1 = smf.ols('y ~ cfd + mit + disc', data=df).fit()
out['exp1_ols'] = {
  'baseline_r2': float(m0.rsquared), 'disc_r2': float(m1.rsquared),
  'disc_coef': float(m1.params['disc']), 'disc_p': float(m1.pvalues['disc']),
  'lr_test_p': float(m1.compare_f_test(m0)[1])}

# Exp2: crossed mixed effects (study groups + guide variance component)
try:
    vc = {'guide': '0 + C(guide)'}
    try:
        mm2 = smf.mixedlm('y ~ cfd + mit + disc', data=df, groups=df['study'], vc_formula=vc).fit(reml=True, method='lbfgs', maxiter=400)
        mm1 = smf.mixedlm('y ~ cfd + mit', data=df, groups=df['study'], vc_formula=vc).fit(reml=False, method='lbfgs', maxiter=400)
        mm2n = smf.mixedlm('y ~ cfd + mit + disc', data=df, groups=df['study'], vc_formula=vc).fit(reml=False, method='lbfgs', maxiter=400)
    except Exception:
        # lbfgs hits singular matrices on this small crossed design; nm is robust
        mm2 = smf.mixedlm('y ~ cfd + mit + disc', data=df, groups=df['study'], vc_formula=vc).fit(reml=True, method='nm', maxiter=1000)
        mm1 = smf.mixedlm('y ~ cfd + mit', data=df, groups=df['study'], vc_formula=vc).fit(reml=False, method='nm', maxiter=1000)
        mm2n = smf.mixedlm('y ~ cfd + mit + disc', data=df, groups=df['study'], vc_formula=vc).fit(reml=False, method='nm', maxiter=1000)
    # lbfgs can silently drop the guide variance component (empty vcomp); refit with nm
    if len(getattr(mm2n, 'vcomp', [])) == 0:
        mm1 = smf.mixedlm('y ~ cfd + mit', data=df, groups=df['study'], vc_formula=vc).fit(reml=False, method='nm', maxiter=1000)
        mm2n = smf.mixedlm('y ~ cfd + mit + disc', data=df, groups=df['study'], vc_formula=vc).fit(reml=False, method='nm', maxiter=1000)
    lr = 2*(mm2n.llf - mm1.llf)
    p_lr = float(1 - sps.chi2.cdf(max(lr,0), 1))
    out['exp2_mixed'] = {
      'disc_coef': float(mm2n.params['disc']), 'disc_p': float(mm2n.pvalues['disc']),
      'lrt_disc_chi2': float(max(lr,0)), 'lrt_disc_p': p_lr,
      'study_var': float(mm2n.cov_re.iloc[0,0]) if mm2n.cov_re.size else None,
      'guide_var': float(mm2n.vcomp[0]) if len(getattr(mm2n,'vcomp',[])) else None,
      'resid_var': float(mm2n.scale),
      'converged': bool(mm2n.converged)}
except Exception as e:
    out['exp2_mixed'] = {'error': str(e)}

# Exp3: stratify by mismatch class - Spearman(residual of baseline fit, disc) per class
df['resid0'] = m0.resid
exp3 = {}
for c in [1,2,3,4]:
    d = df[df.nmm == c]
    if len(d) < 10: exp3[str(c)] = {'n': len(d), 'note': 'too few'}; continue
    rho, p = sps.spearmanr(d['resid0'], d['disc'])
    exp3[str(c)] = {'n': int(len(d)), 'spearman_disc_resid': float(rho), 'p': float(p)}
out['exp3_by_mismatch_class'] = exp3

# Exp4: leave-one-study-out
studies = sorted(df.study.unique())
exp4 = {}
for st in studies:
    tr, te = df[df.study != st], df[df.study == st]
    if len(te) < 10: continue
    b = smf.ols('y ~ cfd + mit', data=tr).fit()
    d_ = smf.ols('y ~ cfd + mit + disc', data=tr).fit()
    pb, pd_ = b.predict(te), d_.predict(te)
    exp4[st] = {'n_test': int(len(te)),
                'spearman_baseline': float(sps.spearmanr(te.y, pb)[0]),
                'spearman_disc': float(sps.spearmanr(te.y, pd_)[0])}
out['exp4_leave_one_study_out'] = exp4
wins = sum(1 for v in exp4.values() if v['spearman_disc'] > v['spearman_baseline'])
out['exp4_disc_wins_studies'] = f'{wins}/{len(exp4)}'

print(json.dumps(out, indent=1))
json.dump(out, open('/home/sandbox/wave1/logs/round2_foldback_result.json','w'), indent=1)
