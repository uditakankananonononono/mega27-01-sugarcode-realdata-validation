"""sugarcode bio.orf.find_orfs vs orfipy and vs the NCBI-annotated CDS of 5 RefSeq mRNAs.
Checks: (1) longest ATG ORF on + strand == annotated CDS, (2) translation == annotated /translation,
(3) full six-frame ORF set (min 30 aa, longest-nested) agrees with orfipy."""
import sys, json, orfipy_core as oc
from Bio import SeqIO
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.orf import find_orfs, translate
res=[]
for rec in SeqIO.parse('/home/sandbox/mega/m01/data/refseq/five_refseq_mrna.gb','genbank'):
    s=str(rec.seq).upper(); cds=[f for f in rec.features if f.type=='CDS'][0]
    st,en=int(cds.location.start),int(cds.location.end); prot=cds.qualifiers['translation'][0]
    ours=find_orfs(s,min_aa=30)
    plus=[o for o in ours if o['strand'] in ('+',1)]; L=max(plus,key=lambda o:o['length_nt'])
    oo=oc.orfs(s,minlen=93,maxlen=10**7,starts=['ATG'],strand='b',include_stop=True,partial3=False,partial5=False,between_stops=False)
    ours_set={(o['start'],o['end'],'+' if o['strand'] in ('+',1) else '-') for o in ours}
    orfipy_set={(a,b,c) for a,b,c,_ in oo}
    res.append({'accession':rec.id,'annotated_cds':[st,en],'our_longest_plus':[L['start'],L['end']],'longest_equals_cds':(L['start'],L['end'])==(st,en),
      'translation_equals_annotation':translate(s[st:en],cds=True).rstrip('*')==prot,
      'n_ours':len(ours_set),'n_orfipy':len(orfipy_set),'n_shared':len(ours_set&orfipy_set),
      'only_ours_example':sorted(ours_set-orfipy_set)[:3],'only_orfipy_example':sorted(orfipy_set-ours_set)[:3]})
    print(res[-1],flush=True)
json.dump({'tool':'orfipy (orfipy_core.orfs, ATG starts, both strands, min 93 nt incl stop)','data':'NCBI RefSeq GenBank records (E-utilities)','results':res},open('sweep_orf_orfipy.json','w'),indent=1)
