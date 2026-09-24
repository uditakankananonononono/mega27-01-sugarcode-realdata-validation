"""STRING PPI-enrichment test for dex-induced genes (PyDESeq2 ~cell+dex, padj<0.05, log2FC>1, GSE52778).
Entrez IDs mapped to symbols with NCBI E-utilities esummary; STRING v12 API /ppi_enrichment, species 9606."""
import json, urllib.request, urllib.parse, pandas as pd
d=pd.read_csv('benchmarks/deseq2_dex_paired.tsv.gz',sep='\t')
up=[str(g) for g in d[(d.padj<0.05)&(d.log2fc>1)].sort_values('padj').gene]
sym={}
for i in range(0,len(up),200):
    j=json.load(urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json&id='+','.join(up[i:i+200]),timeout=60))['result']
    for u in j.get('uids',[]): sym[u]=j[u]['name']
syms=[sym[g] for g in up if g in sym]
def post(ep,params):
    return urllib.request.urlopen(urllib.request.Request('https://string-db.org/api/json/'+ep,data=urllib.parse.urlencode(params).encode()),timeout=100).read()
ver=json.loads(urllib.request.urlopen('https://string-db.org/api/json/version',timeout=30).read())
top=syms[:300]
e=json.loads(post('ppi_enrichment',{'identifiers':'\r'.join(top),'species':9606,'caller_identity':'mega27-01'}))
out={'tool':'STRING API /ppi_enrichment','string_version':ver,'input':'top 300 dex-induced genes by padj (of %d with symbols)'%len(syms),'n_mapped_symbols':len(syms),'result':e}
json.dump(out,open('sweep_dex_string_ppi.json','w'),indent=1); print(ver, e)
