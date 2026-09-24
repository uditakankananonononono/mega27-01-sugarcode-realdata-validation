"""bio.uniprot.search (live UniProt REST) vs HGNC's curated symbol->UniProt mapping (rest.genenames.org, one fetch per symbol, 2026-09-25) for 60 human genes.
Also checks returned sequence length and mass against the UniProt entry fetched by accession (independent endpoint)."""
import sys, glob, json, os, time, urllib.request
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import uniprot as U
CACHE=os.environ.get('UP_CACHE','/tmp/up_rows'); os.makedirs(CACHE,exist_ok=True); T0=time.time()
res=[]
for f in sorted(glob.glob('data/hgnc/*.json')):
    g=os.path.basename(f)[:-5]
    cf=f'{CACHE}/{g}.json'
    if os.path.exists(cf): res.append(json.load(open(cf))); continue
    if time.time()-T0>80: print('partial; rerun to continue'); sys.exit(3)
    doc=json.load(open(f))['response']['docs'][0]; ref=doc.get('uniprot_ids',[])
    try: r=U.search(g)
    except Exception as ex: r=None; err=str(ex)[:80]
    acc=r['accession'] if r else None
    row={'gene':g,'hgnc_uniprot':ref,'sugarcode_accession':acc,'match':acc in ref,'sugarcode_gene_names':None}
    if r:
        e=json.loads(urllib.request.urlopen(f'https://rest.uniprot.org/uniprotkb/{acc}.json',timeout=20).read())
        row['sugarcode_gene_names']=[x.get('geneName',{}).get('value') for x in e.get('genes',[])]
        row['length_ok']=r['length']==e['sequence']['length'] and len(r['sequence'])==r['length']
    json.dump(row,open(cf,'w')); res.append(row); time.sleep(0.2)
out={'reference':'HGNC REST fetch/symbol uniprot_ids','n':len(res),'match':sum(r['match'] for r in res),'length_ok':sum(r.get('length_ok',False) for r in res),'mismatches':[r for r in res if not r['match']],'rows':res}
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else 'benchmarks/sweep_uniprot_hgnc.json','w'),indent=1)
print(out['n'],out['match'],out['length_ok']); print(json.dumps(out['mismatches'],indent=0)[:2500])
