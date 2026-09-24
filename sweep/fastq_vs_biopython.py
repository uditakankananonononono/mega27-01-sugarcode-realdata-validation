"""bio.fastq vs Biopython SeqIO on the first 20,000 reads of ENA SRR622461_1.fastq.gz (NA12878, 1000 Genomes low coverage).
Checks ids, sequences, per-base Phred scores, mean Phred, GC fraction, and write->parse round trip."""
import gzip, json, sys, numpy as np, Bio
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import fastq as F
P='data/fastq/SRR622461_1_first20k.fastq.gz'
t=gzip.open(P,'rt').read(); a=F.parse_fastq(t); b=list(SeqIO.parse(gzip.open(P,'rt'),'fastq'))
idk=[k for k in ('id','name','header') if k in a[0]][0]
ids=sum((r[idk]==x.description or r[idk]==x.id) for r,x in zip(a,b)); seq=sum(r['sequence']==str(x.seq) for r,x in zip(a,b))
ph=sum(F.phred_scores(r['quality'])==x.letter_annotations['phred_quality'] for r,x in zip(a,b))
mp=max(abs(F.mean_phred(r)-np.mean(x.letter_annotations['phred_quality'])) for r,x in zip(a,b))
s=F.stats(a); allseq=''.join(str(x.seq) for x in b)
gc_ref=gc_fraction(allseq,ambiguous='remove'); gc_ign=gc_fraction(allseq,ambiguous='ignore')
rt=F.parse_fastq(F.write_fastq(a))==a
qd=np.bincount(np.concatenate([x.letter_annotations['phred_quality'] for x in b])).tolist()
out={'reference':'Biopython %s SeqIO fastq-sanger + gc_fraction'%Bio.__version__,'source':'ENA SRR622461_1.fastq.gz first 20,000 reads','n_sugarcode':len(a),'n_biopython':len(b),
 'ids_equal':ids,'seqs_equal':seq,'phred_equal':ph,'max_mean_phred_diff':float(mp),'gc_sugarcode':s['gc'],'gc_biopython_remove_N':round(gc_ref,4),'gc_biopython_ignore_N':round(gc_ign,4),'gc_note':'sugarcode divides by A+C+G+T (Biopython ambiguous=remove); ambiguous=ignore divides by all bases incl. N','detected_offset':s['offset'],'roundtrip':rt,'phred_histogram':qd}
json.dump(out,open('benchmarks/sweep_fastq_biopython.json','w'),indent=1); print({k:v for k,v in out.items() if k!='phred_histogram'})
