import csv,json,numpy as np
from maxentpy import maxent
from sklearn.metrics import roc_auc_score
F=40
M5=maxent.load_matrix5(); M3=maxent.load_matrix3()
win={r['vid']:r for r in csv.DictReader(open('windows.tsv'),delimiter='\t')}
ids=[l.strip() for l in open('cv_order.txt') if l.strip()]
P=np.load('cv_pred.npy'); assert P.shape[1]==len(ids),(P.shape,len(ids))
y=P[3]; me=[]; site=[]
for v in ids:
    r=win[v]; k=int(r['k']); R,A=r['ref_ctx'],r['alt_ctx']
    if k>0: b=F-k; rs,as_=maxent.score5(R[b-2:b+7],matrix=M5),maxent.score5(A[b-2:b+7],matrix=M5)
    else: e=F-k; rs,as_=maxent.score3(R[e-20:e+3],matrix=M3),maxent.score3(A[e-20:e+3],matrix=M3)
    me.append(rs-as_); site.append('donor' if k>0 else 'acceptor')
me=np.array(me); site=np.array(site)
out={'n':len(ids),'maxentscan_delta':round(roc_auc_score(y,me),4),'pwm':round(roc_auc_score(y,P[0]),4),'logit':round(roc_auc_score(y,P[1]),4),'cnn':round(roc_auc_score(y,P[2]),4)}
for s in ('donor','acceptor'):
    m=site==s; out[s]={'n':int(m.sum()),'maxentscan':round(roc_auc_score(y[m],me[m]),4),'cnn':round(roc_auc_score(y[m],P[2][m]),4),'pwm':round(roc_auc_score(y[m],P[0][m]),4)}
b=np.random.default_rng(7); diffs=[]
for _ in range(1000):
    i=b.integers(0,len(y),len(y))
    if y[i].min()==y[i].max(): continue
    diffs.append(roc_auc_score(y[i],P[2][i])-roc_auc_score(y[i],me[i]))
out['cnn_minus_maxent_ci95']=[round(float(np.percentile(diffs,2.5)),4),round(float(np.percentile(diffs,97.5)),4)]
out['tool']='MaxEntScan (Yeo & Burge 2004) via maxentpy (github.com/kepbod/maxentpy)'
json.dump(out,open('maxent_cmp.json','w'),indent=1); print(json.dumps(out))
