import json,numpy as np,collections,sys
from sklearn.linear_model import LogisticRegression
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.modules.neohunter.core import HLA_SUPERTYPES
AA='ACDEFGHIKLMNPQRSTVWY'
data=collections.defaultdict(list)
for l in open('bdata.20130222.mhci.txt'):
    s=l.rstrip('\n').split('\t')
    if s[0]=='human' and s[2]=='9': data[s[1]].append((s[3],s[4],float(s[5])))
W={}
for a in HLA_SUPERTYPES:
    X=[];y=[]
    for p,ineq,m in data['HLA-'+a]:
        if not set(p)<=set(AA): continue
        if m<500 and ineq in('=','<'): y.append(1)
        elif m>=500 and ineq in('=','>'): y.append(0)
        else: continue
        v=np.zeros(180)
        for i,c in enumerate(p): v[i*20+AA.index(c)]=1
        X.append(v)
    m=LogisticRegression(C=0.5,max_iter=3000).fit(np.array(X),np.array(y))
    W[a]={'bias':round(float(m.intercept_[0]),5),'weights':[[round(float(m.coef_[0][i*20+j]),5) for j in range(20)] for i in range(9)],'n_train':len(y)}
json.dump({'source':'IEDB MHC-I binding data 2013 (bdata.20130222.mhci.txt), 9-mers, binder IC50<500nM','alphabet':AA,'model':'L2 logistic regression on one-hot 9x20 (C=0.5)','alleles':W},open('/home/sandbox/mega/sc-ai/src/sugarcode/modules/neohunter/iedb_pssm_9mer.json','w'))
print({a:W[a]['n_train'] for a in W})
