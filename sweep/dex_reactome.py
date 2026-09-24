"""Reactome AnalysisService over-representation for dex-induced genes (symbols from NCBI esummary),
projected to human. Compared with g:Profiler's REAC results in the paper."""
import json, urllib.request, pandas as pd
d=pd.read_csv('benchmarks/deseq2_dex_paired.tsv.gz',sep='\t')
up=[str(g) for g in d[(d.padj<0.05)&(d.log2fc>1)].gene]
sym={}
for i in range(0,len(up),200):
    j=json.load(urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json&id='+','.join(up[i:i+200]),timeout=60))['result']
    for u in j.get('uids',[]): sym[u]=j[u]['name']
body='\n'.join(sym[g] for g in up if g in sym).encode()
r=json.load(urllib.request.urlopen(urllib.request.Request('https://reactome.org/AnalysisService/identifiers/projection?pageSize=25&page=1&sortBy=ENTITIES_FDR&order=ASC&resource=TOTAL',data=body,headers={'Content-Type':'text/plain','User-Agent':'Mozilla/5.0 (mega27-01 research)'}),timeout=100))
top=[{'stId':p['stId'],'name':p['name'],'fdr':p['entities']['fdr'],'p':p['entities']['pValue'],'found':p['entities']['found'],'total':p['entities']['total']} for p in r['pathways']]
out={'tool':'Reactome AnalysisService /identifiers/projection','summary':{k:r['summary'].get(k) for k in ('token','type','sampleName')},'identifiers_not_found':r.get('identifiersNotFound'),'n_query':len(sym),'pathways_found':r.get('pathwaysFound'),'top25':top,
 'note':'Reactome uses its whole human gene universe as background, not the expressed genes.'}
json.dump(out,open('sweep_dex_reactome.json','w'),indent=1)
for t in top[:10]: print(t['stId'],t['name'][:60],'fdr=%.2g'%t['fdr'],t['found'],'/',t['total'])
