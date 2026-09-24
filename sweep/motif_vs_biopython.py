"""bio.motif (JASPAR I/O + log-odds scanning) vs Biopython Bio.motifs on 10 JASPAR 2024 profiles (jaspar.elixir.no REST, 2026-09-25)
scanned over BRCA1 10 kb (data/brca1_10kb.fa) and phage lambda (NC_001416.1): every window score on both strands, and the hit set at 80% of the score range."""
import json, sys, glob, os, numpy as np, Bio
from Bio import motifs, SeqIO
from Bio.Seq import Seq
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import motif as M
seqs={os.path.basename(p):str(next(SeqIO.parse(p,'fasta')).seq).upper() for p in ('data/brca1_10kb.fa','data/refseq/NC_001416.1.fa')}
res=[]
for f in sorted(glob.glob('data/jaspar/*.jaspar')):
    ours=M.read_jaspar(open(f).read())[0]
    with open(f) as h: bm=motifs.read(h,'jaspar')
    pssm=bm.counts.normalize(pseudocounts=0.5).log_odds({'A':.25,'C':.25,'G':.25,'T':.25})
    lod_diff=max(abs(ours['lod'][i][b]-pssm[b][i]) for i in range(ours['length']) for b in 'ACGT')
    row={'motif':ours['name'],'length':ours['length'],'max_lod_diff':float(lod_diff),'min_eq':abs(ours['min_score']-pssm.min)<1e-9,'max_eq':abs(ours['max_score']-pssm.max)<1e-9}
    for nm,s in seqs.items():
        L=ours['length']; fw=pssm.calculate(Seq(s)); rv=pssm.reverse_complement().calculate(Seq(s))
        sc_f=np.array([M.score_site(s[i:i+L],ours) for i in range(len(s)-L+1)])
        thr=M.threshold_score(ours,0.8)
        bio_hits={(i,'+') for i,v in enumerate(fw) if v>=thr-1e-9}|{(i,'-') for i,v in enumerate(rv) if v>=thr-1e-9}
        our_hits={(h['start'],h['strand']) for h in M.scan(s,ours,threshold_fraction=0.8)}
        row[nm]={'windows':len(sc_f),'max_abs_score_diff_fwd':float(np.nanmax(np.abs(sc_f-fw))),'hits_sugarcode':len(our_hits),'hits_biopython':len(bio_hits),'hit_sets_equal':our_hits==bio_hits,
                 'n_diff':len(our_hits^bio_hits)}
    res.append(row)
out={'reference':'Biopython %s Bio.motifs (jaspar reader, pseudocount 0.5, uniform background, log2)'%Bio.__version__,'n_motifs':len(res),
 'all_hit_sets_equal':all(r[k]['hit_sets_equal'] for r in res for k in seqs),'per_motif':res}
json.dump(out,open('benchmarks/sweep_motif_biopython.json','w'),indent=1)
for r in res: print(r['motif'][:22],r['max_lod_diff'],r['min_eq'],r['max_eq'],[(r[k]['hits_sugarcode'],r[k]['hits_biopython'],r[k]['n_diff'],'%.1e'%r[k]['max_abs_score_diff_fwd']) for k in seqs])
