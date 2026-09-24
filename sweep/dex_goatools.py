"""Local GO:BP over-representation of dex-induced genes with goatools (Fisher + BH), as an independent cross-check of the g:Profiler result."""
import gzip, csv, json, collections, sys, pandas as pd, goatools
import statsmodels.sandbox.stats.multicomp as _mc, statsmodels.stats.multitest as _mt; _mc.multipletests=_mt.multipletests  # goatools vs statsmodels 0.15 import path
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
B=sys.argv[1] if len(sys.argv)>1 else 'benchmarks/'
sym={}
for r in csv.reader(gzip.open(sys.argv[2] if len(sys.argv)>2 else 'gene_info.gz','rt'),delimiter='\t'):
    if not r[0].startswith('#'): sym[r[1]]=r[2]
d=pd.read_csv(B+'deseq2_dex_paired.tsv.gz',sep='\t'); d['gene']=d.gene.astype(str)
pop={sym[g] for g in d.gene if g in sym}
study={sym[g] for g in d[(d.padj<0.05)&(d.log2fc>1)].gene if g in sym}
dag=GODag('go-basic.obo',optional_attrs=[],load_obsolete=False,prt=None)
assoc=collections.defaultdict(set)
for l in gzip.open('goa_human.gaf.gz','rt'):
    if l[0]=='!': continue
    f=l.split('\t')
    if f[8]!='P' or 'NOT' in f[3] or f[4] not in dag: continue
    assoc[f[2]].add(f[4])
assoc={k:v for k,v in assoc.items() if k in pop}
g=GOEnrichmentStudy(pop,assoc,dag,propagate_counts=True,alpha=0.05,methods=['fdr_bh'],log=None)
res=[r for r in g.run_study(study,prt=None) if r.enrichment=='e']
res.sort(key=lambda r:r.p_fdr_bh)
gp=json.load(open(B+'sweep_dex_enrichment.json'))
gp_bp={t['native'] for t in gp['top_terms'] if t['source']=='GO:BP'}
sig=[r for r in res if r.p_fdr_bh<0.05]
gc=[r for r in res if 'glucocorticoid' in r.name]
out={'tool':'goatools %s GOEnrichmentStudy (Fisher exact, BH FDR, propagated counts)'%goatools.__version__,'go_release':dag.version if hasattr(dag,'version') else 'go-basic.obo 2026-07-26',
 'annotations':'GOA human GAF (current.geneontology.org, 2026-09-25), BP, NOT excluded','n_population':len(pop),'n_study':len(study),'n_annotated_pop':len(assoc),
 'n_sig_bp_fdr05':len(sig),'top25':[{'go':r.GO,'name':r.name,'study':r.ratio_in_study[0],'pop':r.ratio_in_pop[0],'p':r.p_uncorrected,'fdr':r.p_fdr_bh} for r in res[:25]],
 'glucocorticoid_terms':[{'go':r.GO,'name':r.name,'study':r.ratio_in_study[0],'pop':r.ratio_in_pop[0],'p':r.p_uncorrected,'fdr':r.p_fdr_bh} for r in gc],
 'gprofiler_top_bp_in_goatools_sig':sorted(gp_bp & {r.GO for r in sig}),'gprofiler_top_bp_total':len(gp_bp)}
json.dump(out,open(B+'sweep_dex_goatools.json','w'),indent=1)
print(len(pop),len(study),len(sig),len(out['gprofiler_top_bp_in_goatools_sig']),'/',len(gp_bp))
for t in out['top25'][:8]: print(t['go'],t['name'][:50],t['study'],'%.2g'%t['fdr'])
for t in out['glucocorticoid_terms']: print('GC',t['name'],t['study'],t['pop'],'%.2g'%t['fdr'])
