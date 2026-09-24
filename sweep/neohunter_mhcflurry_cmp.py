import json,collections,numpy as np
from sklearn.metrics import roc_auc_score
from mhcflurry import Class1AffinityPredictor
P=Class1AffinityPredictor.load()
AA='ACDEFGHIKLMNPQRSTVWY'
data=collections.defaultdict(list)
for l in open('bdata.20130222.mhci.txt'):
    s=l.rstrip('\n').split('\t')
    if s[0]=='human' and s[2]=='9': data[s[1]].append((s[3],s[4],float(s[5])))
R=json.load(open('neo_bench.json')); R['mhcflurry_note']='MHCflurry pan-allele models were trained on data that includes IEDB 2013; its AUROC here is in-sample (upper reference), not a fair held-out comparison.'
for a in R['alleles']:
    pep=[];y=[]
    for p,ineq,m in data['HLA-'+a]:
        if not set(p)<=set(AA): continue
        if m<500 and ineq in('=','<'): y.append(1)
        elif m>=500 and ineq in('=','>'): y.append(0)
        else: continue
        pep.append(p)
    pr=P.predict(peptides=pep,alleles=['HLA-'+a]*len(pep))
    R['alleles'][a]['auc_mhcflurry_insample']=round(roc_auc_score(y,-np.log(pr)),4)
    print(a,R['alleles'][a],flush=True)
R['mhcflurry_version']=__import__('mhcflurry').__version__
json.dump(R,open('neo_bench.json','w'),indent=1)
