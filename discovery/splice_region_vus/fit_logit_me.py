import json,numpy as np
exec(open('cnn_cv_me.py').read().split('class Net')[0])
from sklearn.linear_model import LogisticRegression
lr=LogisticRegression(max_iter=2000).fit(Fz,y)
P=np.load('cv_pred_me.npy'); s,yy=P[1],P[3]; cal={}
for t in (0.5,0.8,0.9,0.95,0.98):
    m=s>=t; tp=int((m&(yy==1)).sum()); fp=int((m&(yy==0)).sum()); cal[str(t)]={'tp':tp,'fp':fp,'sensitivity':round(tp/(yy==1).sum(),4),'fpr':round(fp/(yy==0).sum(),5)}
json.dump({'model':'logistic regression on deepsplice PWM ref/alt/delta, donor flag, MaxEntScan ref/alt/delta (/10), offset one-hot','feature_order':'[pwm_ref,pwm_alt,pwm_delta,is_donor,me_ref/10,me_alt/10,me_delta/10]+onehot(k in KS)','KS':KS,
 'coef':[round(float(c),6) for c in lr.coef_[0]],'intercept':round(float(lr.intercept_[0]),6),'n_train':int(len(y)),
 'cv_auroc':0.9674,'calibration':cal,'source':'mega27-01 discovery/splice_region_vus/cv_results_me.json'},open('/home/sandbox/mega/sc-ai/src/sugarcode/tools/data/splice_logit_me_v1.json','w'),indent=1)
print(cal['0.9'],cal['0.98'])
