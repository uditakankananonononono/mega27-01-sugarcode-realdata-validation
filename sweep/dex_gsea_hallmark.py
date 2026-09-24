"""Preranked GSEA (gseapy) of the PyDESeq2 dex results on MSigDB Hallmark (50 sets):
can a rank-based test recover the glucocorticoid biology that over-representation missed?"""
import gzip, csv, json, math, gseapy, pandas as pd, sys
D=sys.argv[1] if len(sys.argv)>1 else 'benchmarks/'
sym={}
for r in csv.reader(gzip.open('data/Homo_sapiens.gene_info.gz' if __import__('os').path.exists('data/Homo_sapiens.gene_info.gz') else 'gene_info.gz','rt'),delimiter='\t'):
    if not r[0].startswith('#'): sym[r[1]]=r[2]
rows=[]
for r in csv.DictReader(gzip.open(D+'deseq2_dex_paired.tsv.gz','rt'),delimiter='\t'):
    if r['pvalue'] in ('','nan','NA') or r['gene'] not in sym: continue
    rows.append((sym[r['gene']], math.copysign(-math.log10(max(float(r['pvalue']),1e-300)),float(r['log2fc']))))
rnk=pd.DataFrame(rows,columns=['gene','score']).drop_duplicates('gene').sort_values('score',ascending=False)
lib='MSigDB_Hallmark_2020'; gs=json.load(open(lib+'.json' if __import__('os').path.exists(lib+'.json') else 'benchmarks/'+lib+'.json'))
res=gseapy.prerank(rnk=rnk,gene_sets=gs,permutation_num=1000,seed=7,threads=1,min_size=10,max_size=500,outdir=None,verbose=False).res2d
res['FDR q-val']=res['FDR q-val'].astype(float); res=res.sort_values('FDR q-val')
top=[{'term':t,'nes':float(n),'fdr':float(q),'nom_p':float(p)} for t,n,q,p in zip(res.Term,res.NES,res['FDR q-val'],res['NOM p-val'])]
out={'tool':f'gseapy {gseapy.__version__} prerank (1000 permutations, seed 7)','gene_sets':lib+' via gseapy.get_library (Enrichr)',
 'rank_metric':'sign(log2FC)*-log10(p), PyDESeq2 paired ~cell+dex (benchmarks/deseq2_dex_paired.tsv.gz)','symbol_source':'NCBI Homo_sapiens.gene_info (FTP, 2026-09-25)',
 'n_ranked':len(rnk),'n_tested':len(res),'n_fdr_lt_0.05':int((res['FDR q-val']<0.05).sum()),'n_fdr_lt_0.25':int((res['FDR q-val']<0.25).sum()),'all':top}
json.dump(out,open(D+'sweep_dex_gsea_hallmark.json','w'),indent=1)
print(out['n_ranked'],out['n_tested'],out['n_fdr_lt_0.05'],out['n_fdr_lt_0.25'])
for t in top[:12]: print(round(t['nes'],2),round(t['fdr'],4),t['term'])
