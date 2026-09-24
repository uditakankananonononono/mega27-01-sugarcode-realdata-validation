import sys,json,math,numpy as np,collections
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.modules.neohunter.core import hla_binding, HLA_SUPERTYPES
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
AA='ACDEFGHIKLMNPQRSTVWY'
def iedb(a): return 'HLA-'+a
data=collections.defaultdict(list)
for l in open('bdata.20130222.mhci.txt'):
    s=l.rstrip('\n').split('\t')
    if s[0]!='human' or s[2]!='9': continue
    data[s[1]].append((s[3],s[4],float(s[5])))
def onehot(p): 
    v=np.zeros(9*20)
    for i,c in enumerate(p):
        if c in AA: v[i*20+AA.index(c)]=1
    return v
R={'source':'IEDB MHC-I binding data 2013 (Kim et al. 2014, bdata.20130222.mhci.txt)','binder_def':'IC50<500nM with = or < inequality; non-binder IC50>=500 with = or >','length':9,'alleles':{}}
for a in HLA_SUPERTYPES:
    rows=[(p,ineq,m) for p,ineq,m in data[iedb(a)] if set(p)<=set(AA)]
    y=[];pep=[]
    for p,ineq,m in rows:
        if m<500 and ineq in ('=','<'): y.append(1)
        elif m>=500 and ineq in ('=','>'): y.append(0)
        else: continue
        pep.append(p)
    y=np.array(y)
    if len(pep)<50 or y.sum()<10: R['alleles'][a]={'n':len(pep),'note':'too few'}; continue
    heur=np.array([hla_binding(p,a)['score'] for p in pep])
    X=np.array([onehot(p) for p in pep]); rng=np.random.default_rng(7); f=rng.integers(0,5,len(pep)); pr=np.zeros(len(pep))
    for k in range(5):
        m=LogisticRegression(C=0.5,max_iter=2000).fit(X[f!=k],y[f!=k]); pr[f==k]=m.predict_proba(X[f==k])[:,1]
    R['alleles'][a]={'n':len(pep),'n_binders':int(y.sum()),'auc_module_heuristic':round(roc_auc_score(y,heur),4),'auc_pssm_logit_cv5':round(roc_auc_score(y,pr),4)}
    print(a,R['alleles'][a],flush=True)
json.dump(R,open('neo_bench.json','w'),indent=1)
