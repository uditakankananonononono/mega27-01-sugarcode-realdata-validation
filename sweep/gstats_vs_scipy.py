"""bio.gstats on real 1000 Genomes phase 3 genotype counts (Ensembl GRCh37 REST, population_genotypes, 20 rsIDs fetched individually, 2026-09-25).
HWE exact p vs an independent full enumeration of the Levene/Haldane distribution; allelic/genotypic chi-square vs scipy.stats.chi2_contingency;
odds ratio + Woolf CI vs statsmodels Table2x2; BH/Bonferroni vs statsmodels multipletests."""
import sys, glob, json, math, os
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import gstats as G
import numpy as np, scipy, statsmodels
from scipy.stats import chi2_contingency
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import multipletests
def lgam(x): return math.lgamma(x+1)
def hwe_enum(aa,ab,bb):
    n=aa+ab+bb; nA=2*aa+ab; nB=2*n-nA
    def lp(h):
        ra=(nA-h)//2; rb=(nB-h)//2
        return lgam(n)-lgam(ra)-lgam(h)-lgam(rb)+h*math.log(2)+lgam(nA)+lgam(nB)-lgam(2*n)
    hs=[h for h in range(nA%2,min(nA,nB)+1,2)]
    ps=np.array([lp(h) for h in hs]); m=ps.max(); pr=np.exp(ps-m); pr/=pr.sum()
    po=pr[hs.index(ab)]; return float(min(1.0,pr[pr<=po*(1+1e-7)].sum()))
POPS=['AFR','AMR','EAS','EUR','SAS']
rows=[]
for f in sorted(glob.glob('data/ensembl_var/rs*.json')):
    d=json.load(open(f))
    if 'population_genotypes' not in d: continue
    rs=os.path.basename(f)[:-5]
    anc=d.get('ancestral_allele'); m=d['mappings'][0]['allele_string'].split('/') if d.get('mappings') else []
    for pop in [f'1000GENOMES:phase_3:{p}' for p in POPS]:
        g={x['genotype'].replace('/','|'):x['count'] for x in d['population_genotypes'] if x['population']==pop}
        if not g: continue
        al=sorted({a for k in g for a in k.split('|')})
        if len(al)>2: continue
        A=m[0] if m and m[0] in al else al[0]; B=[a for a in al if a!=A]; B=B[0] if B else ('N' if A!='N' else 'X')
        aa=g.get(f'{A}|{A}',0); bb=g.get(f'{B}|{B}',0); ab=g.get(f'{A}|{B}',0)+g.get(f'{B}|{A}',0)
        rows.append({'rs':rs,'pop':pop.split(':')[-1],'A':A,'B':B,'counts':[aa,ab,bb]})
h=[]; 
for r in rows:
    ours=G.hwe_exact(*r['counts']); ref=hwe_enum(*r['counts']); r['hwe_sugarcode']=ours; r['hwe_enum']=ref; h.append(abs(ours-ref)/max(ref,1e-300))
hwe={'n_tests':len(rows),'max_rel_err':max(h),'n_rel_err_lt_1e-6':sum(x<1e-6 for x in h),'n_p_lt_0.05':sum(r['hwe_sugarcode']<0.05 for r in rows)}
# association EUR vs EAS
assoc=[]
by={(r['rs'],r['pop']):r for r in rows}
for rs in sorted({r['rs'] for r in rows}):
    if (rs,'EUR') not in by or (rs,'EAS') not in by: continue
    ca=by[(rs,'EUR')]['counts']; co=by[(rs,'EAS')]['counts']
    if by[(rs,'EUR')]['A']!=by[(rs,'EAS')]['A']: continue
    tab=np.array([[2*ca[0]+ca[1],ca[1]+2*ca[2]],[2*co[0]+co[1],co[1]+2*co[2]]],float)
    if (tab.sum(0)==0).any(): continue
    rec={'rs':rs}
    for y in (False,True):
        o=G.allelic_test(tuple(ca),tuple(co),yates=y); st,p,_,_=chi2_contingency(tab,correction=y)
        rec['yates' if y else 'pearson']={'stat_sc':o['statistic'],'stat_scipy':float(st),'p_sc':o['p'],'p_scipy':float(p)}
    o=G.allelic_test(tuple(ca),tuple(co)); t=Table2x2(tab+ (0.5 if (tab==0).any() else 0)); lo,hi=t.oddsratio_confint()
    rec['or']={'sc':[o['odds_ratio'],o['ci_low'],o['ci_high']],'statsmodels':[float(t.oddsratio),float(lo),float(hi)]}
    gt=np.array([ca,co],float); keep=gt.sum(0)>0
    if keep.all():
        o=G.genotypic_test(tuple(ca),tuple(co)); st,p,_,_=chi2_contingency(gt,correction=False); rec['genotypic']={'stat_sc':o['statistic'],'stat_scipy':float(st),'p_sc':o['p'],'p_scipy':float(p)}
    assoc.append(rec)
def rel(a,b): return abs(a-b)/max(abs(b),1e-300)
mx=lambda key,sub: max(rel(r[key][sub+'_sc'],r[key][sub+'_scipy']) for r in assoc if key in r)
ps=[r['pearson']['p_sc'] for r in assoc]
bh=G.benjamini_hochberg(ps); bh_ref=multipletests(ps,method='fdr_bh')[1]; bo=G.bonferroni(ps); bo_ref=multipletests(ps,method='bonferroni')[1]
out={'reference':f'scipy {scipy.__version__}, statsmodels {statsmodels.__version__}, independent HWE enumeration','hwe':hwe,
 'association_EUR_vs_EAS':{'n_snps':len(assoc),'pearson_stat_max_rel':mx('pearson','stat'),'pearson_p_max_rel':mx('pearson','p'),'yates_stat_max_rel':mx('yates','stat'),'yates_p_max_rel':mx('yates','p'),
   'genotypic_n':sum('genotypic' in r for r in assoc),'genotypic_stat_max_rel':mx('genotypic','stat'),'genotypic_p_max_rel':mx('genotypic','p'),
   'or_max_rel':max(rel(a,b) for r in assoc for a,b in zip(r['or']['sc'],r['or']['statsmodels']))},
 'multiple_testing':{'bh_max_abs':float(max(abs(a-b) for a,b in zip(bh,bh_ref))),'bonferroni_max_abs':float(max(abs(a-b) for a,b in zip(bo,bo_ref)))},
 'rows':rows,'assoc':assoc}
json.dump(out,open('benchmarks/sweep_gstats_scipy.json','w'),indent=1,default=float)
print(json.dumps({k:out[k] for k in ('hwe','association_EUR_vs_EAS','multiple_testing')},indent=1))
worst=sorted(rows,key=lambda r:-abs(r['hwe_sugarcode']-r['hwe_enum'])/max(r['hwe_enum'],1e-300))[:3]; print([(r['rs'],r['pop'],r['counts'],r['hwe_sugarcode'],r['hwe_enum']) for r in worst])
