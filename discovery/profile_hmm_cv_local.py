import sys,json; sys.path.insert(0,__import__('os').environ.get('SUGARCODE_SRC','../sc-ai/src'))
from sugarcode.modules.profile_hmm import build_profile_hmm, forward_local
from pyhmmer.easel import MSAFile
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1]/'src')); from sc_validate.auc import auc
D=str(__import__('pathlib').Path(__file__).resolve().parents[1]/'data')+'/'
msa=MSAFile(D+'globins25_hmmalign.sto').read(); rows=list(msa.alignment); names=[s.name.decode() if isinstance(s.name,bytes) else s.name for s in msa.sequences]
def fa(p):
    out={};k=None
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; out[k]=''
        else: out[k]+=l
    return out
G=fa(D+'globins25.fa'); N=[s for s in fa(D+'nonglobins25.fa').values() if set(s)<=set("ACDEFGHIKLMNPQRSTVWY")]
pos=[];neg=[[] for _ in N]
for f in range(5):
    m=build_profile_hmm([r for i,r in enumerate(rows) if i%5!=f], background="alignment")
    for i in range(25):
        if i%5==f: pos.append(forward_local(m,''.join(c for c in G[names[i]] if c.isalpha()))['log_odds_bits'])
    for j,s in enumerate(N): neg[j].append(forward_local(m,s)['log_odds_bits'])
    print('fold',f,flush=True)
neg=[sum(v)/len(v) for v in neg]
r={'auc_forward_local':auc(pos,neg),'min_pos':min(pos),'max_neg':max(neg),'n_pos':len(pos),'n_neg':len(neg)}
print(r); json.dump(r,open('res_local.json','w'))
