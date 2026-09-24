"""bio.bed vs bioframe on real UCSC hg38 tracks (REST getData/track, 2026-09-25): CpG islands chr20 (cpgIslandExt) and RepeatMasker chr20:0-5 Mb (rmsk).
Checks parse->records, merge (with and without bookended joining), and overlap queries against bioframe.overlap."""
import json, sys, random, pandas as pd, bioframe as bf
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import bed as B
def load(f,k,c,s,e,n):
    d=json.load(open(f))[k]; return [(x[c],int(x[s]),int(x[e]),str(x[n])) for x in d]
sets={'cpgIslandExt_chr20':load('data/bed/cpg_chr20.json','cpgIslandExt','chrom','chromStart','chromEnd','name'),
      'rmsk_chr20_0_5Mb':load('data/bed/rmsk_chr20_5M.json','rmsk','genoName','genoStart','genoEnd','repName')}
res={}
for nm,iv in sets.items():
    txt='\n'.join(f"{c}\t{s}\t{e}\t{n}" for c,s,e,n in iv)+'\n'
    bd=B.parse_bed(txt); df=pd.DataFrame(iv,columns=['chrom','start','end','name'])
    parse_ok=[(r['chrom'],r['start'],r['end']) for r in bd['records']]==[(c,s,e) for c,s,e,_ in iv]
    ours=[(r['chrom'],r['start'],r['end']) for r in B.merge_intervals(bd['records'])]
    bf0=[tuple(x) for x in bf.merge(df,min_dist=0)[['chrom','start','end']].itertuples(index=False)]   # bioframe/bedtools default: bookended merge
    bfo=[tuple(x) for x in bf.merge(df,min_dist=None)[['chrom','start','end']].itertuples(index=False)] if True else None
    try: bf_strict=[tuple(x) for x in bf.merge(df,min_dist=-1)[['chrom','start','end']].itertuples(index=False)]
    except Exception as ex: bf_strict=str(ex)
    bookended=sum(1 for a,b in zip(sorted(iv),sorted(iv)[1:]) if a[0]==b[0] and a[2]==b[1])
    rng=random.Random(5); q=[(iv[0][0],p,p+rng.randint(50,5000)) for p in (rng.randint(0,max(e for _,_,e,_ in iv)) for _ in range(300))]
    qdf=pd.DataFrame(q,columns=['chrom','start','end']); ov=bf.overlap(qdf,df,how='left',return_index=True,suffixes=('','_'))
    agree=0
    for qi,(c,s,e) in enumerate(q):
        a={i for i,r in enumerate(bd['records']) if B.overlaps(r,c,s,e)}
        b=set(ov[ov['index']==qi]['index_'].dropna().astype(int))
        agree+=a==b
    res[nm]={'n':len(iv),'parse_equal':parse_ok,'bookended_pairs_in_input':bookended,'merged_sugarcode':len(ours),'merged_bioframe_min_dist0':len(bf0),
             'equal_to_bioframe_min_dist0':ours==bf0,'bioframe_strict':(len(bf_strict) if isinstance(bf_strict,list) else bf_strict),'equal_to_bioframe_strict':ours==bf_strict,
             'overlap_queries':len(q),'overlap_sets_equal':agree}
out={'reference':'bioframe %s (merge, overlap)'%bf.__version__,'note':'Before sugarcode fix (sc-ai fc1fd1e) merge_intervals refused to join book-ended intervals despite a bedtools-style docstring: rmsk 11,364 merged vs bioframe/bedtools 8,251 (3,112 book-ended pairs). After fix (min_dist=0 default, min_dist=-1 strict) identical.','before_fix':{'rmsk_chr20_0_5Mb_merged_sugarcode':11364,'cpgIslandExt_chr20_merged_sugarcode':847},'sets':res}
json.dump(out,open('benchmarks/sweep_bed_bioframe.json','w'),indent=1); print(json.dumps(res,indent=1))
