"""bio.genbank vs Biopython SeqIO on 25 individually fetched NCBI nuccore GenBank records (efetch gbwithparts, 2026-09-25).
Compares sequence, feature count, per-feature spans/strand, multi-line qualifiers (/translation) and CDS extraction."""
import sys, glob, json, os
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import genbank as G
from Bio import SeqIO, __version__ as biov
def bio_spans(f):
    parts=f.location.parts
    return sorted((int(p.start)+1,int(p.end)) for p in parts), (f.location.strand or 1)
tot={'records':0,'seq_equal':0,'features_bio':0,'features_sc':0,'feature_count_equal':0,'spans_equal':0,'strand_equal':0,'translation_q_equal':0,'translation_q_total':0,'cds_extract_equal':0,'cds_total':0}
per=[]
for fn in sorted(glob.glob('data/genbank/*.gb')):
    txt=open(fn).read(); sc=G.parse_genbank(txt); rec=SeqIO.read(fn,'genbank')
    bf=[f for f in rec.features]; sf=sc['features']
    r={'acc':os.path.basename(fn)[:-3],'len_bio':len(rec.seq),'len_sc':len(sc['sequence']),'seq_equal':sc['sequence']==str(rec.seq).upper(),'n_bio':len(bf),'n_sc':len(sf)}
    tot['records']+=1; tot['seq_equal']+=r['seq_equal']; tot['features_bio']+=len(bf); tot['features_sc']+=len(sf); tot['feature_count_equal']+=len(bf)==len(sf)
    sp=st=tq=tqn=ce=cn=0; bad=[]
    for b,s in zip(bf,sf):
        bs,bstr=bio_spans(b); ok=sorted(s['spans'])==bs; sp+=ok; st+=(s['strand']==bstr)
        if not ok and len(bad)<3: bad.append({'key':b.type,'bio':str(b.location),'sc_loc':s['loc'][:80],'sc_spans':s['spans'][:4]})
        if b.type=='CDS' and 'translation' in b.qualifiers:
            tqn+=1; tq+=s['qualifiers'].get('translation','')==b.qualifiers['translation'][0]
            cn+=1
            ex=''.join(sc['sequence'][a-1:e] for a,e in s['spans'])  # file order (origin-spanning joins)
            
            if s['strand']==-1: ex=G.revcomp(ex)
            okx=ex==str(b.extract(rec.seq)).upper(); ce+=okx
            if not okx and len(bad)<3: bad.append({'key':'CDS-extract','bio':str(b.location),'sc_loc':s['loc'][:120]})
    r.update(spans_equal=sp,strand_equal=st,translation_q_equal=tq,translation_q_total=tqn,cds_extract_equal=ce,cds_total=cn,examples_mismatch=bad)
    for k in ('spans_equal','strand_equal','translation_q_equal','translation_q_total','cds_extract_equal','cds_total'): tot[k]+=r[k]
    per.append(r)
out={'reference':'Biopython %s SeqIO genbank'%biov,'totals':tot,'per_record':per}
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else 'benchmarks/sweep_genbank_biopython.json','w'),indent=1)
print(json.dumps(tot,indent=1))
for r in per:
    if not r['seq_equal'] or r['n_bio']!=r['n_sc'] or r['spans_equal']<min(r['n_bio'],r['n_sc']) or r['translation_q_equal']<r['translation_q_total']: print(r['acc'],r['len_bio'],r['len_sc'],r['n_bio'],r['n_sc'],r['spans_equal'],r['translation_q_equal'],r['translation_q_total'],r['cds_extract_equal'],r['examples_mismatch'][:2])
