"""sugarcode bio.align nw/sw (BLOSUM62, gap open -11 / extend -1) vs parasail nw_scan/sw_scan
on all pairs of the 25 Pfam PF00042 globins (m01 data/globins25.fa). Tests both gap conventions."""
import sys, json, itertools, parasail
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.align import nw_align, sw_align
d={};k=None
for l in open('/home/sandbox/mega/m01/data/globins25.fa'):
    l=l.strip()
    if l.startswith('>'): k=l[1:].split()[0]; d[k]=''
    else: d[k]+=l.upper()
S={k:v for k,v in d.items() if set(v)<=set('ARNDCQEGHILKMFPSTWYVBZX')}
rows=[]
for a,b in itertools.combinations(list(S),2):
    for mode,f,pf in (('global',nw_align,parasail.nw_scan),('local',sw_align,parasail.sw_scan)):
        o=f(S[a],S[b],matrix='BLOSUM62',gap_open=-11,gap_extend=-1)['score']
        rows.append((mode,o,pf(S[a],S[b],11,1,parasail.blosum62).score,pf(S[a],S[b],12,1,parasail.blosum62).score))
res={'n_seqs':len(S),'n_pairs':len(rows)//2,'tool':'parasail %s'%'.'.join(map(str,parasail.version())),'params':'sugarcode gap_open=-11 gap_extend=-1; parasail open=11 or 12, extend=1'}
for mode in ('global','local'):
    R=[r for r in rows if r[0]==mode]
    res[mode]={'n':len(R),'exact_vs_parasail_open11':sum(abs(r[1]-r[2])<1e-9 for r in R),'exact_vs_parasail_open12':sum(abs(r[1]-r[3])<1e-9 for r in R),
      'max_abs_diff_open11':float(max(abs(r[1]-r[2]) for r in R)),'max_abs_diff_open12':float(max(abs(r[1]-r[3]) for r in R))}
json.dump(res,open('sweep_align_parasail.json','w'),indent=1); print(res)
