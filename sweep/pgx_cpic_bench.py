import json,sys,re,collections
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.modules.pgx_guidelines.core import translate_phenotype
R={'source':'CPIC API https://api.cpicpgx.org/v1/diplotype (fetched 2026-09-24)'}
for g in ('CYP2C19','CYP2D6'):
    d=json.load(open(g+'.json')); n=agree=cov=0; mism=[]; unsup=collections.Counter()
    for r in d:
        n+=1; t=translate_phenotype(g,r['diplotype'])
        if t['phenotype']=='Missing':
            for a in re.split('/',r['diplotype']): unsup[a]+=1
            continue
        cov+=1; ok=t['phenotype'].lower()==r['generesult'].lower()
        agree+=ok
        if not ok: mism.append([r['diplotype'],r['generesult'],t['phenotype'],r['totalactivityscore'],t['activity_score']])
    R[g]={'cpic_diplotypes':n,'module_covered':cov,'agree':agree,'mismatch':len(mism),'mismatch_examples':mism[:15],'coverage':round(cov/n,4),'top_unsupported_alleles':unsup.most_common(10)}
json.dump(R,open('pgx_bench.json','w'),indent=1,ensure_ascii=False)
for g in ('CYP2C19','CYP2D6'): print(g,{k:v for k,v in R[g].items() if k!='top_unsupported_alleles'})
