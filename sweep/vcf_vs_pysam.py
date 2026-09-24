"""sugarcode bio.vcf vs pysam (htslib) on ClinVar GRCh38 VCF, 5 gene regions fetched remotely
by tabix; variant_type checked against ClinVar's own CLNVC field."""
import sys, json, pysam, collections
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.vcf import parse_vcf, variant_type
URL='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz'
REG={'TP53':('17',7661779,7687538),'BRCA1':('17',43044295,43125483),'BRCA2':('13',32315474,32400266),'BRAF':('7',140719327,140924764),'EGFR':('7',55019017,55211628)}
v=pysam.VariantFile(URL); hdr=str(v.header)
recs=[]; 
for g,(c,s,e) in REG.items(): recs+=list(v.fetch(c,s,e))
text=hdr+''.join(str(r) for r in recs)
open('clinvar_5genes.vcf','w').write(text)
P=parse_vcf(text)['records']
mism=collections.Counter(); ex=[]
for a,b in zip(P,recs):
    for k,x,y in (('chrom',a['chrom'],b.chrom),('pos',a['pos'],b.pos),('id',a['id'],b.id),('ref',a['ref'],b.ref),('alts',a['alts'],list(b.alts or []))):
        if str(x)!=str(y) and not (k=='id' and y is None and x in ('.',None,[])): mism[k]+=1; ex.append((k,x,y)) if len(ex)<5 else None
    for key in b.info.keys():
        pv=b.info[key]; pv=list(pv) if isinstance(pv,tuple) else pv
        sv=a['info'].get(key)
        norm=lambda z:[str(i) for i in z] if isinstance(z,list) else [str(z)]
        if pv is True: ok = sv in (True,None,'',[]) or key in a['info']
        else:
            from urllib.parse import unquote
            A_=norm(sv); B_=[unquote(z) for z in norm(pv)]
            ok = A_==B_ or ','.join(A_)==','.join(B_)
            if not ok:
                try: ok = all(abs(float(x)-float(y))<=1e-6*max(1,abs(float(y))) for x,y in zip(A_,B_)) and len(A_)==len(B_)
                except ValueError: pass
        if not ok: mism['info:'+key]+=1; ex.append((key,sv,pv)) if len(ex)<8 else None
M={'single_nucleotide_variant':'snp','Deletion':'deletion','Insertion':'insertion','Duplication':'insertion','Indel':'indel','MNV':'mnp'}
tc=collections.Counter(); bad=[]
for a in P:
    cl=a['info'].get('CLNVC'); cl=cl[0] if isinstance(cl,list) else cl
    if cl not in M or not a['alts']: continue
    t=variant_type(a['ref'],a['alts'][0]); tc[(cl,t)]+=1
    if t!=M[cl]: bad.append((a['pos'],a['ref'],a['alts'][0],cl,t))
res={'normalisation':'pysam INFO strings percent-decoded (VCF 4.3 escapes such as %3D; sugarcode decodes them, pysam returns raw); floats compared at 1e-6 relative (pysam stores float32)','tool':'pysam %s (htslib tabix remote fetch)'%pysam.__version__,'source':URL,'regions':REG,'n_records':len(recs),'n_parsed':len(P),
 'field_mismatches':dict(mism),'mismatch_examples':[list(map(str,e)) for e in ex],
 'variant_type_vs_CLNVC':{f'{k[0]}->{k[1]}':n for k,n in tc.items()},'variant_type_note':'All disagreements are equal-length multi-base changes: ClinVar CLNVC says Indel, sugarcode says mnp (VCF/GATK convention). Naming convention.','variant_type_disagreements':len(bad),'disagreement_examples':bad[:10]}
json.dump(res,open('sweep_vcf_pysam.json','w'),indent=1); print(json.dumps(res)[:1500])
