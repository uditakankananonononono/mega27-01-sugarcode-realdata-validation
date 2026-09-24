"""Replay UPGMA merges on RF00005 and record whether any merge step had a tied minimum."""
import sys,random; sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio import phylo
from Bio import AlignIO
aln=AlignIO.read('/home/sandbox/mega/m01/data/rfam/RF00005.sto','stockholm'); recs=[r for r in aln]; random.Random(1).shuffle(recs); recs=recs[:30]
seqs={f"s{i}":str(r.seq).upper().replace('U','T').replace('.','-') for i,r in enumerate(recs)}; names=list(seqs)
dm=phylo.distance_matrix(seqs,model='pdistance'); M=dm['matrix'] if isinstance(dm,dict) and 'matrix' in dm else dm
if isinstance(M,dict): M=[[M[a][b] for b in names] for a in names]
cl={i:[i] for i in range(len(names))}; D={(i,j):M[i][j] for i in range(len(names)) for j in range(i+1,len(names))}; ties=0; nxt=len(names)
while len(cl)>1:
    m=min(D.values()); tied=[k for k,v in D.items() if abs(v-m)<1e-12]
    if len(tied)>1: ties+=1
    i,j=tied[0]; ni,nj=len(cl[i]),len(cl[j]); new=cl.pop(i)+cl.pop(j)
    for k in list(cl):
        a=D.pop((min(i,k),max(i,k))); b=D.pop((min(j,k),max(j,k))); D[(k,nxt)]=(a*ni+b*nj)/(ni+nj)
    D.pop((i,j)); cl[nxt]=new; nxt+=1
print('merge steps with tied minimum:',ties)
