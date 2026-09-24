"""bio.sequence vs Biopython on real sequences: 25 NCBI GenBank records (data/genbank) and 500 UniProt reviewed human proteins (data/uniprot, bulk; not counted as accessions).
Checks translate (every table-1 CDS vs its /translation qualifier and vs Bio.Seq.translate), reverse_complement, gc_content, find_motif (IUPAC, vs Bio.SeqUtils.nt_search), molecular_weight (vs ProteinAnalysis), tm_wallace (short oligos vs MeltingTemp.Tm_Wallace) and orfs (vs an independent Biopython frame scan)."""
import sys, glob, json, random, re
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import sequence as S
from Bio import SeqIO, __version__ as biov
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, nt_search, MeltingTemp as MT
from Bio.SeqUtils.ProtParam import ProteinAnalysis
out={'reference':'Biopython %s'%biov}
recs=[SeqIO.read(f,'genbank') for f in sorted(glob.glob('data/genbank/*.gb'))]
# translate
t=dict(cds=0,skipped_other_table_or_frame=0,eq_qualifier=0,eq_bio=0,iupac_in_cds=0,eq_qualifier_iupac=0)
bad=[]
for r in recs:
    for f in r.features:
        if f.type!='CDS' or 'translation' not in f.qualifiers: continue
        if f.qualifiers.get('transl_table',['1'])[0] not in ('1','11') or f.qualifiers.get('codon_start',['1'])[0]!='1' or 'transl_except' in f.qualifiers or 'ribosomal_slippage' in f.qualifiers:
            t['skipped_other_table_or_frame']+=1; continue
        nt=str(f.extract(r.seq)).upper(); q=f.qualifiers['translation'][0]; t['cds']+=1
        iu=bool(re.search('[^ACGT]',nt)); t['iupac_in_cds']+=iu
        ours=S.translate(nt,to_stop=True)
        # first codon: alternative starts (GTG/TTG) are M in /translation
        ok=ours[1:]==q[1:]; t['eq_qualifier']+=ok; t['eq_qualifier_iupac']+=(ok and iu)
        t['eq_bio']+=S.translate(nt)==str(Seq(nt).translate())
        if not ok and len(bad)<4: bad.append({'acc':r.id,'gene':f.qualifiers.get('gene',['?'])[0],'len_nt':len(nt),'ours':ours[:30],'q':q[:30],'iupac':iu})
t['mismatch_examples']=bad; out['translate']=t
# revcomp + gc
out['revcomp_equal']=sum(S.reverse_complement(str(r.seq).upper())==str(r.seq.reverse_complement()).upper() for r in recs)
out['gc']=[{'acc':r.id,'ours':S.gc_content(str(r.seq)),'bio':gc_fraction(r.seq,ambiguous='ignore'),'nonACGT':len(re.findall('[^ACGT]',str(r.seq).upper()))} for r in recs]
out['gc_equal_1e12']=sum(abs(g['ours']-g['bio'])<1e-12 for g in out['gc'])
out['n_records']=len(recs)
# motifs
mot=['GAATTC','GGATCC','RGATCY','CANNTG','TATAWAW','GCCRCCATGG','AATAAA','YGCGY']; m=dict(queries=0,equal=0)
for r in recs[:25]:
    s=str(r.seq).upper()
    if re.search('[^ACGT]',s): continue
    for mo in mot:
        m['queries']+=1; b=nt_search(s,mo)[1:]; m['equal']+=S.find_motif(s,mo)==b
out['find_motif']=m
# proteins
prots=[str(p.seq) for p in SeqIO.parse('data/uniprot/human_reviewed_500.fasta','fasta')]
std=[p for p in prots if set(p)<=set('ACDEFGHIKLMNPQRSTVWY')]
d=[abs(S.molecular_weight(p)-ProteinAnalysis(p).molecular_weight()) for p in std]
rel=[x/ProteinAnalysis(p).molecular_weight() for x,p in zip(d,std)]
out['molecular_weight']={'n':len(std),'n_nonstandard_skipped':len(prots)-len(std),'max_abs_da':max(d),'max_rel':max(rel),'median_abs_da':sorted(d)[len(d)//2]}
# Tm short oligos from lambda
rng=random.Random(1); lam=str([r for r in recs if r.id.startswith('NC_001416')][0].seq)
ol=[lam[i:i+L] for i,L in ((rng.randint(0,40000),rng.randint(8,13)) for _ in range(500))]
out['tm_wallace_short']={'n':len(ol),'equal':sum(S.tm_wallace(o)==MT.Tm_Wallace(o) for o in ol)}
# ORFs
def ref_orfs(s,min_aa):
    res=[]
    for strand,x in (('+',s),('-',str(Seq(s).reverse_complement()))):
        for fr in range(3):
            i=fr
            while i<len(x)-2:
                if x[i:i+3]=='ATG':
                    p=str(Seq(x[i:i+ (len(x)-i)//3*3]).translate())
                    k=p.find('*')
                    if k<0: break
                    if k>=min_aa: res.append((strand,fr,i,i+3*k+3))
                    i+=3*k+3
                else: i+=3
    return res
o=dict(records=0,equal=0,n_orfs_bio=0,n_orfs_ours=0)
for r in recs:
    s=str(r.seq).upper()
    if re.search('[^ACGT]',s) or len(s)>60000: continue
    a=sorted((x['strand'],x['frame'],x['start'],x['end']) for x in S.orfs(s,min_aa=100)); b=sorted(ref_orfs(s,100))
    o['records']+=1; o['equal']+=a==b; o['n_orfs_bio']+=len(b); o['n_orfs_ours']+=len(a)
out['orfs_min100']=o
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else 'benchmarks/sweep_sequence_biopython.json','w'),indent=1)
print(json.dumps({k:v for k,v in out.items() if k!='gc'},indent=1)); print([g for g in out['gc'] if abs(g['ours']-g['bio'])>=1e-12])
