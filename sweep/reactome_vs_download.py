"""bio.reactome.pathways_for_uniprot (live ContentService) vs Reactome's own bulk mapping file UniProt2Reactome.txt (lowest-level pathways, downloaded 2026-09-25; the 60-gene human subset is committed as data/reactome/UniProt2Reactome_60genes.tsv).
Accessions are HGNC's curated UniProt IDs for the 60 symbols in data/hgnc."""
import sys, glob, json, collections, time
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import reactome as R
ref=collections.defaultdict(set)
for l in open('data/reactome/UniProt2Reactome_60genes.tsv'): f=l.split('\t'); ref[f[0].split('-')[0]].add(f[1])  # canonical + isoform rows (P01116-1, -2)
accs=sorted({a for f in glob.glob('data/hgnc/*.json') for a in json.load(open(f))['response']['docs'][0].get('uniprot_ids',[])})
rows=[]
for a in accs:
    live={p['id'] for p in R.pathways_for_uniprot(a)}
    rows.append({'accession':a,'n_live':len(live),'n_file':len(ref[a]),'equal':live==ref[a],'only_live':sorted(live-ref[a])[:10],'only_file':sorted(ref[a]-live)[:10]})
tl=sum(r['n_live'] for r in rows); tf=sum(r['n_file'] for r in rows)
inter=sum(len(set()) for r in rows)
out={'reference':'Reactome UniProt2Reactome.txt (lowest-level, Homo sapiens; canonical plus isoform rows pooled per accession)','n_accessions':len(rows),'equal_sets':sum(r['equal'] for r in rows),'pairs_live':tl,'pairs_file':tf,'rows':rows}
json.dump(out,open('benchmarks/sweep_reactome_download.json','w'),indent=1)
print(out['n_accessions'],out['equal_sets'],tl,tf); print([r for r in rows if not r['equal']][:5])
