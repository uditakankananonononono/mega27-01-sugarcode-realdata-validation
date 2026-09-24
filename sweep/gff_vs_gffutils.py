"""bio.gff vs gffutils on the NCBI RefSeq E. coli K-12 MG1655 annotation (GCF_000005845.2_ASM584v2_genomic.gff.gz).
Checks every record's type/coords/strand/phase and attributes, parent->children links, and write->parse round trip."""
import gzip, json, sys, collections, gffutils
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import gff as G
P='data/refseq/GCF_000005845.2_genomic.gff.gz'
text=gzip.open(P,'rt').read()
a=G.parse_gff(text)
db=gffutils.create_db(text,':memory:',from_string=True,merge_strategy='create_unique',keep_order=True,sort_attribute_values=False)
fb=list(db.all_features(order_by=('file_order',)))  # creation order
n_ok=0; bad=[]
for r,f in zip(a['records'],fb):
    same=(r['seqid'],r['type'],r['start'],r['end'],r['strand'])==(f.seqid,f.featuretype,f.start,f.end,f.strand) and \
         (r['phase'] if r['phase'] is not None else '.')==(int(f.frame) if f.frame!='.' else '.') and \
         {k:(v if isinstance(v,list) else [v]) for k,v in r['attributes'].items()}=={k:list(v) for k,v in f.attributes.items()}
    n_ok+=same
    if not same and len(bad)<5: bad.append({'sc':[r['type'],r['start'],r['end'],dict(list(r['attributes'].items())[:3])],'gffutils':[f.featuretype,f.start,f.end,dict(list(f.attributes.items())[:3])]})
# parent links
ids=[r['attributes'].get('ID') for r in a['records']]; ids=[i[0] if isinstance(i,list) else i for i in ids if i]
genes=[i for i in ids if i.startswith('gene-')][:300]; link_ok=0
for g in genes:
    import re
    sc=sorted((c['attributes'].get('ID')[0] if isinstance(c['attributes'].get('ID'),list) else c['attributes'].get('ID')) for c in G.children_of(a,g))
    gu=sorted(re.sub(r'_\d+$','',c.id) if c.id not in sc else c.id for c in db.children(g,level=1))  # gffutils create_unique renames repeated IDs of multi-line (split) features X -> X_1
    link_ok+= sc==gu
rt=G.parse_gff(G.write_gff(a)); rt_ok=rt['records']==a['records']
out={'reference':'gffutils %s (sqlite, in-memory)'%gffutils.__version__,'file':P,'n_records_sugarcode':len(a['records']),'n_features_gffutils':len(fb),
 'records_identical':n_ok,'examples_differ':bad,'genes_checked_children':len(genes),'children_identical_multiset':link_ok,'roundtrip_lossless':rt_ok,
 'types':collections.Counter(r['type'] for r in a['records']).most_common(8)}
json.dump(out,open('benchmarks/sweep_gff_gffutils.json','w'),indent=1)
print({k:out[k] for k in ['n_records_sugarcode','n_features_gffutils','records_identical','genes_checked_children','children_identical_multiset','roundtrip_lossless']}); print(bad[:2])
