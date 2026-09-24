"""modules.cfd_offtarget vs the official Doench 2016 CFD implementation logic (cfd-score-calculator.py) applied to the published
mismatch_score.pkl / pam_scores.pkl (data/CFD_Scoring), on the 600 BRCA1 guide/off-target pairs (data/crispr_fresh_pairs.tsv)
crossed with all 16 PAM dinucleotides (NNN with positions 2-3 varied)."""
import sys, csv, pickle, json, itertools
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.modules.cfd_offtarget import core as C
mm=pickle.load(open('data/CFD_Scoring/mismatch_score.pkl','rb')); pam=pickle.load(open('data/CFD_Scoring/pam_scores.pkl','rb'))
def revcom(s): return s.translate(str.maketrans('ACGTU','TGCAA'))[::-1]
def official(wt,off,p):  # Doench et al. 2016 cfd-score-calculator.py
    wt=wt.replace('T','U'); off=off.replace('T','U'); score=1
    for i,sl in enumerate(off):
        if wt[i]==sl: continue
        key='r'+wt[i]+':d'+revcom(sl)+','+str(i+1); score*=mm[key]
    return score*pam[p]
rows=list(csv.DictReader(open('data/crispr_fresh_pairs.tsv'),delimiter='\t'))
n=0; mx=0.0; exact=0; bad=[]
for r in rows:
    g=r['guide'][:20]; o=r['offtarget'][:20]
    for a,b in itertools.product('ACGT',repeat=2):
        ref=official(g,o,a+b); ours=C.score(g,o,pam='N'+a+b).score; n+=1
        d=abs(ours-ref); mx=max(mx,d); exact+=d<1e-12
        if d>=1e-12 and len(bad)<3: bad.append((g,o,a+b,ours,ref))
out={'reference':'Doench 2016 CFD pickles + official scoring logic','n_scores':n,'n_pairs':len(rows),'pams':16,'max_abs_diff':mx,'exact_1e-12':exact,'examples_bad':bad}
json.dump(out,open('benchmarks/sweep_cfd_doench.json','w'),indent=1); print(out)
