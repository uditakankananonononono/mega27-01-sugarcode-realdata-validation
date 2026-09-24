"""Hub analysis of the STRING network among the top 300 dex-induced genes (networkx).
Asks whether known glucocorticoid-response genes sit at network hubs more than chance (degree permutation)."""
import json, gzip, csv, random, urllib.request, urllib.parse, networkx as nx, pandas as pd
d=pd.read_csv('benchmarks/deseq2_dex_paired.tsv.gz',sep='\t')
up=[str(g) for g in d[(d.padj<0.05)&(d.log2fc>1)].sort_values('padj').gene]
sym={}
for i in range(0,len(up),200):
    j=json.load(urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json&id='+','.join(up[i:i+200]),timeout=60))['result']
    for u in j.get('uids',[]): sym[u]=j[u]['name']
top=[sym[g] for g in up if g in sym][:300]
req=urllib.request.Request('https://string-db.org/api/json/network',data=urllib.parse.urlencode({'identifiers':'\r'.join(top),'species':9606,'caller_identity':'mega27-01'}).encode())
edges=json.loads(urllib.request.urlopen(req,timeout=100).read())
G=nx.Graph()
for e in edges: G.add_edge(e['preferredName_A'],e['preferredName_B'],w=e['score'])
deg=dict(G.degree()); btw=nx.betweenness_centrality(G)
comps=sorted(nx.connected_components(G),key=len,reverse=True)
# canonical GR targets named in the GSE52778 paper (Himes et al. 2014) and GR literature
gr=['FKBP5','TSC22D3','PER1','DUSP1','KLF15','ZBTB16','CRISPLD2','ERRFI1','SGK1','KLF9']
present=[g for g in gr if g in G]
obs=sum(deg[g] for g in present)/max(1,len(present))
nodes=list(G); random.seed(7)
null=[sum(deg[n] for n in random.sample(nodes,len(present)))/len(present) for _ in range(10000)] if present else []
p=(1+sum(x>=obs for x in null))/(1+len(null)) if present else None
out={'tool':'networkx %s on STRING API /network (v12)'%nx.__version__,'n_input':len(top),'n_nodes':G.number_of_nodes(),'n_edges':G.number_of_edges(),
 'largest_component':len(comps[0]) if comps else 0,'n_components':len(comps),
 'top_degree':sorted(deg.items(),key=lambda x:-x[1])[:15],'top_betweenness':[(k,round(v,4)) for k,v in sorted(btw.items(),key=lambda x:-x[1])[:15]],
 'gr_targets_checked':gr,'gr_targets_in_network':present,'gr_targets_not_in_network':[g for g in gr if g not in G],
 'gr_mean_degree':obs,'network_mean_degree':sum(deg.values())/len(deg),'perm_p_gr_degree':p,'n_perm':10000}
json.dump(out,open('benchmarks/sweep_dex_string_network.json','w'),indent=1)
print({k:out[k] for k in ['n_nodes','n_edges','largest_component','gr_targets_in_network','gr_mean_degree','network_mean_degree','perm_p_gr_degree']}); print(out['top_degree'][:8])
