#!/usr/bin/env python3
"""F1 round-4 bench: guide-level feature decomposition.
Round-3 pivot established: |CFD-MIT| within-guide signal is real (permutation p<0.001)
but does NOT transfer across guides -> mismatch models omit GUIDE-LEVEL determinants.
This experiment: compute per-guide sequence features; ask (A) do guide features explain
the guide-level random-intercept variance (0.96 in round 2)? (B) does disc's
within-guide signal survive controlling guide features? (C) which features matter?"""
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

def guide_features(g20):
    gc = (g20.count('G')+g20.count('C'))/20.0
    seed = g20[12:]  # PAM-proximal 8nt
    gc_seed = (seed.count('G')+seed.count('C'))/8.0
    # dinucleotide Shannon entropy (sequence complexity)
    di = [g20[i:i+2] for i in range(19)]
    from collections import Counter
    c = Counter(di); n = len(di)
    ent = -sum((v/n)*math.log2(v/n) for v in c.values())
    # longest homopolymer run
    run = mx = 1
    for i in range(1,20):
        run = run+1 if g20[i]==g20[i-1] else 1
        mx = max(mx, run)
    # simple self-complementarity proxy: max contiguous Watson-Crick hairpin stem
    comp = {'A':'T','T':'A','G':'C','C':'G'}
    rc = ''.join(comp[b] for b in reversed(g20))
    best = 0
    for off in range(-19, 20):
        cur = 0
        for i in range(20):
            j = i+off
            if 0 <= j < 20 and g20[i] == rc[j]:
                cur += 1; best = max(best, cur)
            else: cur = 0
    return dict(gc=gc, gc_seed=gc_seed, dinuc_entropy=ent, homopolymer_max=mx, selfcomp_max=best)

rows = list(csv.DictReader(open('/home/sandbox/mega27-01-sugarcode-realdata-validation/data/crispor_offtargets_7studies.tsv'), delimiter='\t'))
guides, per = {}, {}
for r in rows:
    if r['type'] == 'on-target': guides[r['name']] = r['seq'].upper()
    per.setdefault(r['name'], []).append(r)
data = []
for name, rs in per.items():
    g = guides.get(name)
    if not g or len(g) != 23: continue
    g20 = g[:20]; gf = guide_features(g20)
    for r in rs:
        if r['type'] != 'off-target': continue
        s = r['seq'].upper()
        if len(s) != 23 or any(c not in 'ACGT' for c in s): continue
        nmm = sum(1 for a,b in zip(g20, s[:20]) if a != b)
        if nmm > 4 or nmm == 0: continue
        cfd = calc_cfd(g20, s[:20], s[-2:])
        if cfd is None: continue
        mit_s = calc_mit(g20, s[:20])/100.0
        data.append(dict(guide=name, cfd=cfd, mit=mit_s, disc=abs(cfd-mit_s),
                         y=math.log10(float(r['score'])+1e-6), nmm=nmm, **gf))
df = pd.DataFrame(data)
out = {'n_rows': len(df), 'n_guides': int(df.guide.nunique())}
gfeats = ['gc','gc_seed','dinuc_entropy','homopolymer_max','selfcomp_max']

# A: guide random-intercept variance before vs after adding guide features
m_pre = smf.mixedlm('y ~ cfd + mit + disc', data=df, groups=df['guide']).fit(reml=False)
f_full = 'y ~ cfd + mit + disc + ' + ' + '.join(gfeats)
m_post = smf.mixedlm(f_full, data=df, groups=df['guide']).fit(reml=False)
out['A_guide_variance'] = {
  'var_without_features': float(m_pre.cov_re.iloc[0,0]),
  'var_with_features': float(m_post.cov_re.iloc[0,0]),
  'var_explained_pct': float(100*(1 - m_post.cov_re.iloc[0,0]/m_pre.cov_re.iloc[0,0]))}

# B: disc within-guide signal after guide features
out['B_disc_after_features'] = {
  'disc_coef': float(m_post.params['disc']), 'disc_p': float(m_post.pvalues['disc'])}

# C: which guide features matter (population-level)
out['C_feature_effects'] = {f: {'coef': float(m_post.params[f]), 'p': float(m_post.pvalues[f])} for f in gfeats}
lr = 2*(m_post.llf - m_pre.llf); k = len(gfeats)
out['C_features_lrt'] = {'chi2': float(max(lr,0)), 'df': k, 'p': float(1-sps.chi2.cdf(max(lr,0),k))}

print(json.dumps(out, indent=1))
json.dump(out, open('/home/sandbox/wave1/logs/round4_guide_features_result.json','w'), indent=1)
