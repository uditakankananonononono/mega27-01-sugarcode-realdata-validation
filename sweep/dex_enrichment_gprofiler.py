"""Enrichment of dex-induced genes (PyDESeq2 ~cell+dex, padj<0.05, log2FC>1) with g:Profiler g:GOSt
(custom background = all 16,802 tested genes) and the Reactome AnalysisService."""
import json, urllib.request, pandas as pd
d=pd.read_csv('benchmarks/deseq2_dex_paired.tsv.gz',sep='\t')
up=[str(g) for g in d[(d.padj<0.05)&(d.log2fc>1)].gene]; bg=[str(g) for g in d.gene]
q={'organism':'hsapiens','query':up,'numeric_ns':'ENTREZGENE_ACC','background':bg,'domain_scope':'custom','sources':['GO:BP','REAC','KEGG'],'user_threshold':0.05,'significance_threshold_method':'g_SCS','no_evidences':True}
r=json.load(urllib.request.urlopen(urllib.request.Request('https://biit.cs.ut.ee/gprofiler/api/gost/profile/',data=json.dumps(q).encode(),headers={'Content-Type':'application/json','User-Agent':'mega27-01'}),timeout=100))
res=r['result']; top=[{'source':x['source'],'native':x['native'],'name':x['name'],'p':x['p_value'],'intersection':x['intersection_size'],'term_size':x['term_size']} for x in sorted(res,key=lambda x:x['p_value'])[:25]]
out={'query':'dex-induced genes, PyDESeq2 ~cell+dex padj<0.05 and log2FC>1 (GSE52778)','n_query':len(up),'n_background':len(bg),'tool':'g:Profiler g:GOSt API (g:SCS correction, custom background)','gprofiler_meta_version':r.get('meta',{}).get('version'),'n_significant_terms':len(res),'top_terms':top}
json.dump(out,open('sweep_dex_enrichment.json','w'),indent=1)
for t in top[:15]: print(t['source'],t['native'],t['name'][:60],'%.2g'%t['p'],t['intersection'])
