import gzip,re,json,math,sys,numpy as np
from scipy.stats import spearmanr
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.codon import cai, load_published_table, STANDARD_CODE
from sugarcode.modules.codon_opt import core
# CDS by locus tag
seqs={};gene={};cur=None
for l in gzip.open('cds.fna.gz','rt'):
    if l[0]=='>':
        m=re.search(r'\[locus_tag=(b\d+)\]',l); g=re.search(r'\[gene=([^\]]+)\]',l)
        cur=m.group(1) if m and 'pseudo=true' not in l else None
        if cur: seqs[cur]='';gene[cur]=g.group(1) if g else cur
    elif cur: seqs[cur]+=l.strip()
seqs={k:v for k,v in seqs.items() if len(v)%3==0 and len(v)>=300 and re.fullmatch('[ACGT]+',v)}
ab={}
for l in open('pax.txt'):
    if l[0]=='#':continue
    n,i,a=l.split('\t'); ab[i.split('.')[1]]=float(a)
keys=[k for k in seqs if k in ab and ab[k]>0]
y=np.log10([ab[k] for k in keys])
def counts(ks):
    c={}
    for k in ks:
        s=seqs[k]
        for i in range(0,len(s)-3,3): c[s[i:i+3]]=c.get(s[i:i+3],0)+1
    return c
def freq_table(c):
    t={}
    for cod,aa in STANDARD_CODE.items():
        tot=sum(c.get(x,0) for x,a in STANDARD_CODE.items() if a==aa)
        t[cod]=(c.get(cod,0)+0.5)/(tot+0.5) if tot else 0
    return t
R={'n_genes':len(keys),'abundance':'PaxDb 511145 WHOLE_ORGANISM integrated','cds':'NCBI GCF_000005845.2 cds_from_genomic'}
genome_tab=load_published_table('e_coli_316407') if 'e_coli_316407' else None
x=[cai(seqs[k],genome_tab) for k in keys]; R['module_cai_genome_table']=round(spearmanr(x,y)[0],4)
# Sharp&Li reference: ribosomal protein genes (rpl/rps/rpm), evaluate excluding them
ribo=[k for k in seqs if re.match(r'rp[lsm][A-Z]',gene[k])]
ev=[i for i,k in enumerate(keys) if k not in ribo]
t_ribo=freq_table(counts(ribo)); x2=[cai(seqs[keys[i]],t_ribo) for i in ev]
R['n_ribo_ref']=len(ribo); R['n_eval_excl_ribo']=len(ev)
R['module_cai_genome_table_excl_ribo']=round(spearmanr([x[i] for i in ev],y[ev])[0],4)
R['cai_ribo_ref_excl_ribo']=round(spearmanr(x2,y[ev])[0],4)
# 5-fold CV: reference = top 5% abundance in training folds
rng=np.random.default_rng(7); f=rng.integers(0,5,len(keys)); pr=np.zeros(len(keys))
for k5 in range(5):
    tr=[keys[i] for i in range(len(keys)) if f[i]!=k5]; ytr=y[f!=k5]
    top=[tr[i] for i in np.argsort(-ytr)[:int(.05*len(tr))]]
    t=freq_table(counts(top))
    for i in np.where(f==k5)[0]: pr[i]=cai(seqs[keys[i]],t)
R['cai_top5pct_ref_cv']=round(spearmanr(pr,y)[0],4)
# Biopython cross-check of module CAI implementation on same ribo table
try:
    from Bio.SeqUtils import CodonAdaptationIndex as C
    bc=C([seqs[k] for k in ribo]); xb=[bc.calculate(seqs[keys[i]]) for i in ev[:500]]
    xs=[cai(seqs[keys[i]],{c:v for c,v in freq_table(counts(ribo)).items()}) for i in ev[:500]]
    R['biopython_vs_module_pearson_500']=round(float(np.corrcoef(xb,xs)[0,1]),5)
    R['biopython_cai_ribo_ref_excl_ribo_500']=round(spearmanr(xb,y[ev[:500]])[0],4)
except Exception as e: R['biopython_err']=str(e)[:200]
json.dump(R,open('cai_bench.json','w'),indent=1); print(json.dumps(R,indent=1))
