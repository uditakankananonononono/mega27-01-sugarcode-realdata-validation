"""bio.entrez.gene_id (live NCBI esearch) vs HGNC's curated NCBI Gene ID (entrez_id) for the same 60 human symbols (data/hgnc)."""
import sys, glob, json, os
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import entrez as E
rows=[]
for f in sorted(glob.glob('data/hgnc/*.json')):
    g=os.path.basename(f)[:-5]; ref=json.load(open(f))['response']['docs'][0].get('entrez_id')
    ours=E.gene_id(g); rows.append({'gene':g,'hgnc_entrez':ref,'sugarcode':ours,'match':ours==ref})
out={'reference':'HGNC REST entrez_id','n':len(rows),'match':sum(r['match'] for r in rows),'mismatches':[r for r in rows if not r['match']],'rows':rows}
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else 'benchmarks/sweep_entrez_hgnc.json','w'),indent=1)
print(out['n'],out['match'],out['mismatches'])
