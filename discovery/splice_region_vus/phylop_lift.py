"""Does phyloP100way conservation at the variant position add to logit+MaxEnt for splice-region SNV P vs B?
Stratified subsample of 1,501 of the 9,235 modelled variants (fetch_phylop_sub.py). Gene-grouped 5-fold CV repeated over 20 seeds;
paired bootstrap on out-of-fold predictions for the AUC difference."""
import csv, json, numpy as np
from maxentpy import maxent
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
F=40; M5=maxent.load_matrix5(); M3=maxent.load_matrix3()
win={r['vid']:r for r in csv.DictReader(open('windows.tsv'),delimiter='\t')}
pwm={r['vid']:r for r in csv.DictReader(open('pwm_scores.tsv'),delimiter='\t')}
ph={}
for l in open('phylop_sub.tsv'):
    v,c,p,val=l.rstrip('\n').split('\t')
    if val not in ('NA','ERR'): ph[v]=float(val)
ids=[v for v in (l.strip() for l in open('phylop_sub_ids.txt')) if v in ph]
KS=[3,4,5,6]+list(range(-14,-2))
X0=[];XP=[];y=[];g=[];site=[]
for v in ids:
    r=win[v]; k=int(r['k']); R,A=r['ref_ctx'],r['alt_ctx']
    if k>0: b=F-k; rs,as_=maxent.score5(R[b-2:b+7],matrix=M5),maxent.score5(A[b-2:b+7],matrix=M5)
    else: e=F-k; rs,as_=maxent.score3(R[e-20:e+3],matrix=M3),maxent.score3(A[e-20:e+3],matrix=M3)
    p=pwm[v]; base=[float(p['ref_s']),float(p['alt_s']),float(p['delta']),rs,as_,rs-as_]+[int(k==kk) for kk in KS]
    X0.append(base); XP.append(base+[ph[v],ph[v]*(k>0)]); y.append(int(r['cls']=='P')); g.append(r['gene']); site.append('donor' if k>0 else 'acceptor')
X0=np.array(X0);XP=np.array(XP);y=np.array(y);g=np.array(g);site=np.array(site);phv=XP[:,-2]
def oof(X,seed):
    rng=np.random.default_rng(seed); ug=np.unique(g); perm={u:i for i,u in enumerate(rng.permutation(ug))}; gi=np.array([perm[x] for x in g])
    pr=np.zeros(len(y))
    for tr,te in GroupKFold(5).split(X,y,gi):
        m=make_pipeline(StandardScaler(),LogisticRegression(max_iter=2000,C=1.0)).fit(X[tr],y[tr]); pr[te]=m.predict_proba(X[te])[:,1]
    return pr
a0=[];aP=[];P0=[];PP=[]
for s in range(20):
    p0=oof(X0,s); pp=oof(XP,s); a0.append(roc_auc_score(y,p0)); aP.append(roc_auc_score(y,pp)); P0.append(p0); PP.append(pp)
p0=np.mean(P0,0); pp=np.mean(PP,0); rng=np.random.default_rng(7); d=[]
for _ in range(2000):
    i=rng.integers(0,len(y),len(y))
    if y[i].min()!=y[i].max(): d.append(roc_auc_score(y[i],pp[i])-roc_auc_score(y[i],p0[i]))
out={'n':int(len(y)),'n_pathogenic':int(y.sum()),'n_genes':int(len(np.unique(g))),'phylop_missing':1501-len(ids),'tool':'UCSC Genome Browser REST API getData/track phyloP100way (hg38)',
 'auc_phylop_alone':round(roc_auc_score(y,phv),4),'auc_logit_maxent_mean20':round(float(np.mean(a0)),4),'auc_logit_maxent_phylop_mean20':round(float(np.mean(aP)),4),
 'auc_diff_ci95_paired_bootstrap':[round(float(np.percentile(d,2.5)),4),round(float(np.percentile(d,97.5)),4)],
 'per_site':{s:{'n':int((site==s).sum()),'phylop_alone':round(roc_auc_score(y[site==s],phv[site==s]),4),'base':round(roc_auc_score(y[site==s],p0[site==s]),4),'with_phylop':round(roc_auc_score(y[site==s],pp[site==s]),4)} for s in ('donor','acceptor')},
 'phylop_median_P':round(float(np.median(phv[y==1])),3),'phylop_median_B':round(float(np.median(phv[y==0])),3)}
json.dump(out,open('phylop_lift.json','w'),indent=1); print(json.dumps(out,indent=1))
