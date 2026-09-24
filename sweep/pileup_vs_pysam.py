"""bio.pileup vs pysam/htslib pileup on the same 10,943 NA12878 chr20 reads (data/sam/na12878_chr20_1M.sam.gz).
Same input, same filters (only unmapped excluded, no base/map-quality filter, overlaps not de-duplicated): per-position depth (bases + deletions) and A/C/G/T/N counts."""
import gzip, json, sys, collections, pysam
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import sam as S, pileup as PU
t=gzip.open('data/sam/na12878_chr20_1M.sam.gz','rt').read()
rows=PU.pileup(S.parse_sam(t))
open('/tmp/_p.sam','w').write(t)
pysam.sort('-o','/tmp/_p.bam','/tmp/_p.sam'); pysam.index('/tmp/_p.bam')
bam=pysam.AlignmentFile('/tmp/_p.bam')
ref={}
for col in bam.pileup('20',stepper='nofilter',flag_filter=4,ignore_overlaps=False,ignore_orphans=False,min_base_quality=0,min_mapping_quality=0,max_depth=100000,truncate=False):
    c=collections.Counter(); dels=0
    for pr in col.pileups:
        if pr.is_del or pr.is_refskip: dels+=pr.is_del; continue
        c[pr.alignment.query_sequence[pr.query_position].upper()]+=1
    ref[col.reference_pos+1]={'counts':c,'dels':dels}
ours={}
for r in rows:
    c=collections.Counter()
    for k,v in r['counts'].items():
        if k!='*': c[k.upper()]+=v  # '*' = deletion placeholder (mpileup convention), compared via del_depth
    ours[r['pos']]={'counts':c,'dels':r['del_depth'],'depth':r['depth']}
pos=sorted(set(ref)|set(ours)); cnt_eq=del_eq=0; ex=[]
for p in pos:
    a=ours.get(p,{'counts':collections.Counter(),'dels':0}); b=ref.get(p,{'counts':collections.Counter(),'dels':0})
    ce=+a['counts']==+b['counts']; de=a['dels']==b['dels']; cnt_eq+=ce; del_eq+=de
    if not (ce and de) and len(ex)<5: ex.append({'pos':p,'sc':dict(a['counts']),'sc_del':a['dels'],'pysam':dict(b['counts']),'pysam_del':b['dels']})
out={'reference':'pysam %s pileup (stepper nofilter, flag_filter=4, no quality filters, overlaps kept)'%pysam.__version__,'n_positions_union':len(pos),'n_positions_sugarcode':len(ours),'n_positions_pysam':len(ref),
 'base_counts_equal':cnt_eq,'deletion_counts_equal':del_eq,'examples':ex,'depth_semantics':'sugarcode also lists deletions as "*" in counts (samtools mpileup convention); compared via del_depth instead'}
json.dump(out,open('benchmarks/sweep_pileup_pysam.json','w'),indent=1); print({k:v for k,v in out.items() if k!='examples'}); print(ex[:3])
