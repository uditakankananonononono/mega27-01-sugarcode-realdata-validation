"""bio.sam vs pysam/htslib on 10,943 real reads: 1000 Genomes phase 3 NA12878 chr20:1,000,000-1,200,000
(NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam, region fetched over HTTPS with pysam; saved as data/sam/na12878_chr20_1M.sam.gz).
Checks every field, flags, CIGAR ops, reference end, and typed tags."""
import gzip, json, sys, pysam
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import sam as S
P='data/sam/na12878_chr20_1M.sam.gz'
text=gzip.open(P,'rt').read(); a=S.parse_sam(text)
open('/tmp/_x.sam','w').write(text); b=list(pysam.AlignmentFile('/tmp/_x.sam','r'))
CIG='MIDNSHP=X'
chk={k:0 for k in ['qname','flag','rname','pos','mapq','cigar_ops','rnext_pnext_tlen','seq','qual','alignment_end','tags','flag_bits']}
ex=[]
for r,p in zip(a['records'],b):
    ok={}
    ok['qname']=r['qname']==p.query_name; ok['flag']=r['flag']==p.flag
    ok['rname']=r['rname']==(p.reference_name if not p.is_unmapped or p.reference_id>=0 else None)
    ok['pos']=r['pos']==(p.reference_start+1 if p.reference_start>=0 else 0)
    ok['mapq']=r['mapq']==p.mapping_quality
    ok['cigar_ops']=[(n,o) for n,o in r['cigar_ops']]==[(n,CIG[o]) for o,n in (p.cigartuples or [])]
    rn=r['rnext']; pn=p.next_reference_name
    ok['rnext_pnext_tlen']=(rn=='=' and pn==p.reference_name or rn==pn or (rn=='*' and pn is None)) and r['pnext']==p.next_reference_start+1 and r['tlen']==p.template_length
    ok['seq']=r['seq']==p.query_sequence; ok['qual']=r['qual']==(pysam.qualities_to_qualitystring(p.query_qualities) if p.query_qualities is not None else None)
    e=S.alignment_end(r); ok['alignment_end']=e==(p.reference_end if p.reference_end is not None else None)
    ok['tags']={k:v for k,v in r['tags'].items()}=={k:v for k,v in p.get_tags()}
    fl=r['flags']; ok['flag_bits']=fl.get('paired',fl.get('is_paired'))==p.is_paired if isinstance(fl,dict) and ('paired' in fl or 'is_paired' in fl) else True
    for k,v in ok.items(): chk[k]+=v
    if not all(ok.values()) and len(ex)<5: ex.append({'qname':r['qname'],'failed':[k for k,v in ok.items() if not v]})
out={'reference':'pysam %s (htslib %s)'%(pysam.__version__,pysam.__htslib_version__) if hasattr(pysam,'__htslib_version__') else 'pysam %s'%pysam.__version__,
 'source':'1000 Genomes phase 3 NA12878 chrom20 low-coverage BAM, 20:1000000-1200000','n_sugarcode':len(a['records']),'n_pysam':len(b),'fields_equal':chk,'examples':ex}
json.dump(out,open('benchmarks/sweep_sam_pysam.json','w'),indent=1); print(out['n_sugarcode'],out['n_pysam'],chk); print(ex[:3])
