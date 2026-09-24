import sys,time,math,json; sys.path.insert(0,'.')
import phmm as P
from pyhmmer.easel import MSAFile
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1]/'src'))
from sc_validate.auc import auc
D=str(__import__('pathlib').Path(__file__).resolve().parents[1]/'data')+'/'
msa=MSAFile(D+'globins25_hmmalign.sto').read()
rows=list(msa.alignment); names=[s.name.decode() if isinstance(s.name,bytes) else s.name for s in msa.sequences]
def fa(p):
    out={};k=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; out[k]=''
        else: out[k]+=l
    return out
G=fa(D+'globins25.fa'); N=fa(D+'nonglobins25.fa')
Nseq=[s for s in N.values() if set(s)<=set("ACDEFGHIKLMNPQRSTVWY")] 
t=time.time(); pos=[];neg=[[] for _ in Nseq]; posn=[];negn=[[] for _ in Nseq]
for f in range(5):
    tr=[r for i,r in enumerate(rows) if i%5!=f]; te=[names[i] for i in range(25) if i%5==f]
    m=P.build_profile_hmm(tr, background="alignment")
    for n in te:
        s=''.join(c for c in G[n] if c.isalpha()); r=P.forward(m,s); pos.append(r['log_odds_bits']); posn.append(r['log_odds_bits']/len(s))
    for j,s in enumerate(Nseq):
        r=P.forward(m,s); neg[j].append(r['log_odds_bits']); negn[j].append(r['log_odds_bits']/len(s))
    print("fold",f,round(time.time()-t,1),flush=True)
neg=[sum(v)/len(v) for v in neg]; negn=[sum(v)/len(v) for v in negn]
res={"cv":"5-fold on hmmalign MSA","n_pos":len(pos),"n_neg":len(neg),"auc_log_odds_bits":auc(pos,neg),"auc_per_residue":auc(posn,negn),"min_pos":min(pos),"max_neg":max(neg)}
print(json.dumps(res)); json.dump(res,open('/tmp/phm/res.json','w'))
