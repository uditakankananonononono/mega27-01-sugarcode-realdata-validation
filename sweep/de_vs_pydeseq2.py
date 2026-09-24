"""sugarcode bio.rnaseq/bio.de vs PyDESeq2 on GEO GSE52778 (airway smooth muscle, Himes et al. 2014):
untreated vs dexamethasone, 4 cell lines. NCBI-generated raw counts (GRCh38.p13)."""
import json, sys, numpy as np, pandas as pd
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.rnaseq import size_factors
from sugarcode.bio.de import de_analysis
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
C=pd.read_csv('https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE52778&format=file&file=GSE52778_raw_counts_GRCh38.p13_NCBI.tsv.gz',sep='\t',index_col=0)
U=['GSM1275862','GSM1275866','GSM1275870','GSM1275874']; D=['GSM1275863','GSM1275867','GSM1275871','GSM1275875']
X=C[U+D]; X=X[(X>=10).sum(1)>=4]
meta=pd.DataFrame({'dex':['untrt']*4+['trt']*4,'cell':['c1','c2','c3','c4']*2},index=U+D)
sf_ours=np.array(size_factors(X.values.tolist()))
def run(design):
    dds=DeseqDataSet(counts=X.T,metadata=meta,design=design,quiet=True,n_cpus=1,low_memory=True); dds.deseq2()
    st=DeseqStats(dds,contrast=['dex','trt','untrt'],quiet=True,n_cpus=1); st.summary()
    sf=np.asarray(dds.obs['size_factors'] if 'size_factors' in dds.obs else dds.obsm['size_factors']).ravel()
    return st.results_df, sf
R,sf_py=run('~cell + dex'); R0,_=run('~dex')
table={'genes':[str(g) for g in X.index],'samples':U+D,'counts':X.values.tolist()}
genes=np.array([str(g) for g in X.index])
def sig(p): return set(genes[np.nan_to_num(np.asarray(p,float),nan=1.0)<0.05])
py=sig(R['padj']); py0=sig(R0['padj'])
out={'dataset':'GEO GSE52778 NCBI raw counts GRCh38.p13; untreated GSM1275862/66/70/74 vs Dex GSM1275863/67/71/75; genes with >=10 counts in >=4 samples',
 'tool':'PyDESeq2 (NB GLM, Wald test, n_cpus=1)','n_genes_tested':int(len(X)),
 'size_factor_max_rel_diff':float(np.max(np.abs(sf_ours-sf_py)/sf_py)),
 'pydeseq2_paired_n_sig':len(py),'pydeseq2_unpaired_n_sig':len(py0)}
W={}
for m in ('welch','wilcoxon'):
    o=de_analysis(table,['untrt']*4+['trt']*4,method=m); rows=o['results']
    pad=np.array([np.nan if r.get('padj') is None else r['padj'] for r in rows],float); lfc=np.array([r['log2fc'] for r in rows],float)
    s=sig(pad); W[m]=pad
    out[m]={'n_sig_padj05':len(s),'overlap_with_pydeseq2_paired':len(s&py),'overlap_with_pydeseq2_unpaired':len(s&py0),
      'lfc_pearson_vs_pydeseq2':round(float(np.corrcoef(lfc,R['log2FoldChange'].values)[0,1]),4)}
canon={'2289':'FKBP5','1831':'TSC22D3','5187':'PER1','1843':'DUSP1','28999':'KLF15','7704':'ZBTB16','83716':'CRISPLD2'}
idx={g:i for i,g in enumerate(genes)}
out['canonical_dex_genes_note']='Gene symbols verified via NCBI E-utilities esummary (db=gene)'
out['canonical_dex_genes']={canon[g]:({'pydeseq2_paired_padj':float(R['padj'].values[idx[g]]),'pydeseq2_lfc':round(float(R['log2FoldChange'].values[idx[g]]),3),'welch_padj':float(W['welch'][idx[g]]),'wilcoxon_padj':float(W['wilcoxon'][idx[g]])} if g in idx else 'filtered') for g in canon}
json.dump(out,open('sweep_de_pydeseq2.json','w'),indent=1); print(json.dumps(out,indent=1))
