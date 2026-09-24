"""sugarcode bio.primer.tm_nn vs primer3-py calc_tm (SantaLucia 1998, Na 50 mM, no Mg/dNTP)
on every 20-mer (step 7) of 5 RefSeq transcripts. Matching concentrations: ours uses
C = primer - template/2 = 12.5 nM; primer3 uses dna_conc/4, so dna_conc = 50 nM."""
import json, primer3, sys, numpy as np
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.primer import tm_nn, hairpin_max_stem
import importlib.util as iu; sp=iu.spec_from_file_location('oldp','/tmp/oldsc/old_primer.py'); oldp=iu.module_from_spec(sp); sp.loader.exec_module(oldp)
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
recs={}; cur=None
for l in open('tx.fa'):
    l=l.strip()
    if l.startswith('>'): cur=l[1:].split()[0]; recs[cur]=''
    else: recs[cur]+=l.upper()
rows=[]
for acc,s in recs.items():
    for i in range(0,len(s)-20,7):
        o=s[i:i+20]
        if set(o)-set('ACGT'): continue
        a=tm_nn(o); b=primer3.calc_tm(o,mv_conc=50,dv_conc=0,dntp_conc=0,dna_conc=50,tm_method='santalucia',salt_corrections_method='santalucia')
        c=mt.Tm_NN(Seq(o),dnac1=25,dnac2=25,Na=50,saltcorr=5,nn_table=mt.DNA_NN4); c3=mt.Tm_NN(Seq(o),dnac1=25,dnac2=25,Na=50,saltcorr=5); old=oldp.tm_nn(o)
        hp=primer3.calc_hairpin(o,mv_conc=50,dv_conc=0,dntp_conc=0,dna_conc=50)
        rows.append((acc,i,a,b,c,hairpin_max_stem(o),hp.dg/1000 if hp.structure_found else 0.0,c3,old))
A=np.array([r[2] for r in rows]); P=np.array([r[3] for r in rows]); Bp=np.array([r[4] for r in rows])
from scipy.stats import spearmanr
res={'n_oligos':len(rows),'transcripts':list(recs),'tools':'primer3-py %s calc_tm/calc_hairpin; Biopython Tm_NN'%primer3.__version__,
 'tm_vs_primer3':{'mean_diff':round(float((A-P).mean()),4),'max_abs_diff':round(float(abs(A-P).max()),4),'pearson':round(float(np.corrcoef(A,P)[0,1]),6)},
 'primer3_vs_biopython_NN3_max_abs':round(float(max(abs(r[3]-r[7]) for r in rows)),4),'old_code_vs_biopython_NN4':{'n_affected':int(sum(abs(r[8]-r[4])>1e-6 for r in rows)),'max_abs_diff':round(float(max(abs(r[8]-r[4]) for r in rows)),4)},'note':'sugarcode vendors DNA_NN4 (SantaLucia & Hicks 2004 unified); primer3 and Biopython default use SantaLucia 1998 (DNA_NN3); the ~0.17 C mean offset vs primer3 is a table choice, not a bug.','tm_vs_biopython_NN4':{'mean_diff':round(float((A-Bp).mean()),4),'max_abs_diff':round(float(abs(A-Bp).max()),4)},
 'hairpin_stem_vs_primer3_dG_spearman':round(float(spearmanr([r[5] for r in rows],[-r[6] for r in rows])[0]),4)}
json.dump(res,open('sweep_primer_primer3.json','w'),indent=1); print(res)
