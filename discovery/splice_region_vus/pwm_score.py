import sys,csv,json,random
import os; sys.path[:0]=[os.environ.get('SUGARCODE_SRC','../../../sugarcode-ai/src'), os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../src')]
from sugarcode.modules.deepsplice import core as ds
from sc_validate.auc import auc
meta=json.load(open(os.path.join(sys.path[0],'sugarcode/bio/data/splice_sites/harvest_meta.json')))
train_genes=set(meta['genes'].keys()) if isinstance(meta['genes'],dict) else set(g if isinstance(g,str) else g.get('gene') for g in meta['genes'])
F=40
out=open('pwm_scores.tsv','w'); out.write('vid\tgene\tcls\tsite\tk\tref_s\talt_s\tdelta\ttrain_gene\n')
for r in csv.DictReader(open('windows.tsv'),delimiter='\t'):
    k=int(r['k']); R,A=r['ref_ctx'],r['alt_ctx']
    if r['site']=='donor':
        b=F-k  # last exon base index; 9-mer = 3 exon + 6 intron
        rw,aw=R[b-2:b+7],A[b-2:b+7]
    else:
        e=F-k; rw,aw=R[e-14:e+1],A[e-14:e+1]
    try: v=ds.variant_effect(rw,aw,r['site'])
    except Exception as ex: continue
    out.write(f"{r['vid']}\t{r['gene']}\t{r['cls']}\t{r['site']}\t{k}\t{v.get('ref_score',v.get('ref'))}\t{v.get('alt_score',v.get('alt'))}\t{v['delta']}\t{int(r['gene'] in train_genes)}\n")
out.close()
print(len(train_genes), sorted(train_genes)[:5])
