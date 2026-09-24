"""bio.stockholm vs Biopython AlignIO on 15 Rfam seed alignments (12 fetched 2026-09-25 from rfam.org/family/<acc>/alignment/stockholm + 3 earlier).
Checks: sequence ids/order, aligned sequences, SS_cons, round-trip write->parse, pairwise identity vs an independent count."""
import json, sys, glob, os, Bio
from Bio import AlignIO
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import stockholm as st
res=[]
for f in sorted(glob.glob('data/rfam_sto/*.sto')+glob.glob('data/rfam/*.sto')):
    t=open(f).read(); a=st.parse_stockholm(t); b=AlignIO.read(f,'stockholm')
    ids_ok=[n for n,_ in a['seqs']]==[r.id for r in b]
    seq_ok=all(s==str(r.seq) for (_,s),r in zip(a['seqs'],b))
    ss_ok=a['gc'].get('SS_cons')==b.column_annotations.get('secondary_structure')
    rt=st.parse_stockholm(st.write_stockholm(a)); rt_ok=rt['seqs']==a['seqs'] and rt['gc']==a['gc'] and rt['gf']==a['gf']
    x,y=a['seqs'][0][1],a['seqs'][1][1]
    cols=[(p,q) for p,q in zip(x,y) if p not in '-.' and q not in '-.']
    ind=sum(p.upper()==q.upper() for p,q in cols)/len(cols) if cols else None
    pid=st.pairwise_identity(x,y)
    res.append({'file':f,'ac':a['gf'].get('AC'),'n_seqs':len(a['seqs']),'width':len(a['seqs'][0][1]),'ids_match':ids_ok,'seqs_match':seq_ok,'ss_cons_match':ss_ok,'roundtrip':rt_ok,
                'pid_sugarcode':pid,'pid_independent':ind})
out={'reference':'Biopython %s AlignIO stockholm'%Bio.__version__,'n_alignments':len(res),'n_all_match':sum(r['ids_match'] and r['seqs_match'] and r['ss_cons_match'] and r['roundtrip'] for r in res),
     'total_sequences':sum(r['n_seqs'] for r in res),'per_alignment':res}
json.dump(out,open('benchmarks/sweep_stockholm_biopython.json','w'),indent=1)
print(out['n_alignments'],out['n_all_match'],out['total_sequences'])
for r in res:
    if not (r['ids_match'] and r['seqs_match'] and r['ss_cons_match'] and r['roundtrip']) or (r['pid_independent'] is not None and abs((r['pid_sugarcode'] or 0)-r['pid_independent'])>1e-9): print(r)
