import sys, numpy as np, pandas as pd
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.de import de_analysis_deseq2
C=pd.read_csv('https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE52778&format=file&file=GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz',sep='\t',index_col=0)
U=['GSM1275862','GSM1275866','GSM1275870','GSM1275874']; D=['GSM1275863','GSM1275867','GSM1275871','GSM1275875']
X=C[U+D]; X=X[(X>=10).sum(1)>=4]
t={'genes':[str(g) for g in X.index],'samples':U+D,'counts':X.values.tolist()}
r=de_analysis_deseq2(t,['untrt']*4+['trt']*4,blocks=['c1','c2','c3','c4']*2)
pd.DataFrame(r['results']).to_csv('deseq2_dex_paired.tsv',sep='\t',index=False); print(r['design'],len(r['results']))
