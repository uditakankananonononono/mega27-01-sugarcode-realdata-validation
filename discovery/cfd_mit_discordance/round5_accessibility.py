#!/usr/bin/env python3
"""F1 precommitted final rescue: does matched-cell-type chromatin accessibility
(ENCODE DNase, hg19, UCSC API) collapse the guide-level residual variance in the
round-2 mixed model? Precommitted success metric: % of guide random-intercept
variance explained. Also: does |CFD-MIT| coefficient change?"""
import csv, pickle, math, json
import numpy as np, pandas as pd
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

dnase = json.load(open('logs/dnase_signals.json'))
rows = list(csv.DictReader(open('/home/sandbox/mega27-01-sugarcode-realdata-validation/data/crispor_offtargets_7studies.tsv'), delimiter='\t'))
guides, per = {}, {}
for r in rows:
    if r['type'] == 'on-target': guides[r['name']] = r['seq'].upper()
    per.setdefault(r['name'], []).append(r)
data, dropped = [], {'no_coords': 0, 'no_cfd': 0}
for name, rs in per.items():
    g = guides.get(name)
    if not g or len(g) != 23: continue
    g20 = g[:20]
    for r in rs:
        if r['type'] != 'off-target': continue
        s = r['seq'].upper()
        if len(s) != 23 or any(c not in 'ACGT' for c in s): continue
        nmm = sum(1 for a,b in zip(g20, s[:20]) if a != b)
        if nmm > 4 or nmm == 0: continue
        cfd = calc_cfd(g20, s[:20], s[-2:])
        if cfd is None: dropped['no_cfd'] += 1; continue
        d = dnase.get(s)
        if not d or 'mean' not in d: dropped['no_coords'] += 1; continue
        mit_s = calc_mit(g20, s[:20])/100.0
        data.append(dict(guide=name, cfd=cfd, mit=mit_s, disc=abs(cfd-mit_s),
                         dnase=math.log1p(d['mean']), y=math.log10(float(r['score'])+1e-6)))
df = pd.DataFrame(data)
out = {'n_rows': len(df), 'n_guides': int(df.guide.nunique()), 'dropped': dropped,
       'dnase_descriptives': {'mean': float(df.dnase.mean()), 'sd': float(df.dnase.std()),
                              'frac_zero': float((df.dnase==0).mean())}}
# baseline (round-2) model on THIS subset, then + accessibility
m0 = smf.mixedlm('y ~ cfd + mit + disc', data=df, groups=df['guide']).fit(reml=False)
m1 = smf.mixedlm('y ~ cfd + mit + disc + dnase', data=df, groups=df['guide']).fit(reml=False)
v0, v1 = float(m0.cov_re.iloc[0,0]), float(m1.cov_re.iloc[0,0])
lr = 2*(m1.llf - m0.llf)
out['result'] = {
  'guide_var_baseline': v0, 'guide_var_with_dnase': v1,
  'guide_var_explained_pct': float(100*(1 - v1/v0)),
  'dnase_coef': float(m1.params['dnase']), 'dnase_p': float(m1.pvalues['dnase']),
  'disc_coef_baseline': float(m0.params['disc']), 'disc_p_baseline': float(m0.pvalues['disc']),
  'disc_coef_with_dnase': float(m1.params['disc']), 'disc_p_with_dnase': float(m1.pvalues['disc']),
  'lrt_dnase_chi2': float(max(lr,0)), 'lrt_dnase_p': float(1-sps.chi2.cdf(max(lr,0),1))}
print(json.dumps(out, indent=1))
json.dump(out, open('logs/round5_accessibility_result.json','w'), indent=1)
