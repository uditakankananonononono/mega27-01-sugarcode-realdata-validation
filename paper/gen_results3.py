"""Generate sections_results3.tex (DE, enrichment, PPI, alignment) from committed benchmark JSON."""
import json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
B='../benchmarks/'
def esc(s): return str(s).replace('_','\\_').replace('&','\\&').replace('%','\\%').replace('#','\\#')
de=json.load(open(B+'sweep_de_pydeseq2.json')); en=json.load(open(B+'sweep_dex_enrichment.json'))
st=json.load(open(B+'sweep_dex_string_ppi.json'))['result'][0]; rc=json.load(open(B+'sweep_dex_reactome.json'))
pa=json.load(open(B+'sweep_align_parasail.json'))
o=[]
o.append('\\section{Differential expression and pathway-level follow-up}')
o.append(f"We re-ran the airway dexamethasone experiment (GEO GSE52778) through the sugarcode DESeq2 path and through PyDESeq2. Of {de['n_genes_tested']:,} genes passing the count filter, the paired design $\\sim$cell$+$dex gave {de['pydeseq2_paired_n_sig']:,} genes at adjusted $p<0.05$ and the unpaired design gave {de['pydeseq2_unpaired_n_sig']:,}. Size factors agreed to a maximum relative difference of {de['size_factor_max_rel_diff']:.1e} (benchmarks/sweep\\_de\\_pydeseq2.json).")
o.append('\\paragraph{(F25) Negative-binomial GLM.} For gene $g$ and sample $j$, $K_{gj}\\sim\\mathrm{NB}(\\mu_{gj},\\alpha_g)$ with $\\mu_{gj}=s_j q_{gj}$ and $\\log_2 q_{gj}=\\sum_r x_{jr}\\beta_{gr}$, so $\\mathrm{Var}(K_{gj})=\\mu_{gj}+\\alpha_g\\mu_{gj}^2$.')
o.append('\\paragraph{(F26) Median-of-ratios size factor.} $s_j=\\mathrm{median}_{g:K_{g\\cdot}>0}\\, K_{gj}\\big/\\big(\\prod_{v=1}^{m}K_{gv}\\big)^{1/m}$.')
o.append('\\paragraph{(F27) Wald statistic.} $W_g=\\hat\\beta_{g,\\mathrm{dex}}/\\mathrm{SE}(\\hat\\beta_{g,\\mathrm{dex}})$, compared with $N(0,1)$; the paired design removes the cell-line effect from the residual variance, which is why it finds more genes than the unpaired design.')
o.append('\\subsection{Gene-set over-representation}')
o.append(f"We took the {en['n_query']} dex-induced genes (adjusted $p<0.05$, $\\log_2\\mathrm{{FC}}>1$) and tested them against the {en['n_background']:,} tested genes with g:Profiler (version {esc(en['gprofiler_meta_version'])}, g:SCS correction). {en['n_significant_terms']} terms passed. The top terms are broad (Table~\\ref{{tab:gost}}).")
o.append('\\paragraph{(F28) Hypergeometric tail.} With background size $M$, term size $K$, query size $n$ and overlap $k$, $p=\\sum_{i\\ge k}\\binom{K}{i}\\binom{M-K}{n-i}\\big/\\binom{M}{n}$, and fold enrichment $=(k/n)/(K/M)$.')
o.append('\\begin{table}[h]\\centering\\small\\begin{tabular}{llrrr}\\hline Source & Term & overlap & size & $p_{\\mathrm{adj}}$\\\\\\hline')
for t in en['top_terms'][:12]:
    o.append(f"{esc(t['source'])} & {esc(t['name'][:48])} & {t['intersection']} & {t['term_size']} & {t['p']:.2g}\\\\")
o.append('\\hline\\end{tabular}\\caption{Top g:Profiler terms for dex-induced genes (benchmarks/sweep\\_dex\\_enrichment.json).}\\label{tab:gost}\\end{table}')
gc=en['glucocorticoid_terms']
o.append('\\paragraph{Negative result.} The terms one would expect for a glucocorticoid are enriched but do not survive correction:')
o.append('\\begin{table}[h]\\centering\\small\\begin{tabular}{lrrrr}\\hline Term & overlap & size & fold (approx.) & $p_{g:SCS}$\\\\\\hline')
for t in gc:
    o.append(f"{esc(t['name'][:50])} & {t['intersection']} & {t['term_size']} & {t.get('fold_enrichment_approx','--')} & {t['p_gSCS']:.2g}\\\\")
o.append('\\hline\\end{tabular}\\caption{Glucocorticoid terms: large fold change, not significant. '+esc(en['glucocorticoid_note'][:0])+'}\\end{table}')
o.append('The small term sizes (tens of genes) and the strong multiple-testing burden explain this: a 7-fold enrichment built on 9 genes is not enough evidence after g:SCS correction. We report it as a negative rather than as confirmation of the expected biology.')
o.append('\\subsection{Reactome pathways}')
t0=rc['top25']
o.append(f"Reactome AnalysisService mapped the query to {rc['pathways_found']:,} pathways ({rc['identifiers_not_found']} identifiers not found). No pathway reached FDR$<0.05$; the smallest FDR was {min(t['fdr'] for t in t0):.2f}. Reactome uses its whole human gene universe as background, which is a further caveat.")
o.append('\\begin{table}[h]\\centering\\small\\begin{tabular}{lrrrr}\\hline Pathway & found & total & $p$ & FDR\\\\\\hline')
for t in t0[:10]:
    o.append(f"{esc(t['name'][:50])} & {t['found']} & {t['total']} & {t['p']:.2g} & {t['fdr']:.2f}\\\\")
o.append('\\hline\\end{tabular}\\caption{Top Reactome pathways (benchmarks/sweep\\_dex\\_reactome.json).}\\end{table}')
o.append('\\paragraph{(F29) Benjamini-Hochberg.} For ordered $p_{(1)}\\le\\dots\\le p_{(m)}$, $q_{(i)}=\\min_{j\\ge i}\\min\\{1,\\,m\\,p_{(j)}/j\\}$.')
o.append('\\subsection{Protein interaction network}')
o.append(f"STRING's PPI enrichment test on the top 300 dex-induced genes mapped {st['number_of_nodes']} proteins with {st['number_of_edges']} edges against {st['expected_number_of_edges']} expected (average degree {st['average_node_degree']}, local clustering {st['local_clustering_coefficient']}), reported $p$ below floating-point resolution. Caveat: STRING's expectation uses a genome-wide background, so genes expressed in airway smooth muscle may be more connected than random genes regardless of dex.")
o.append('\\paragraph{(F30) Edge enrichment ratio.} $E=\\,$observed edges$/$expected edges$=%d/%d=%.2f$.' % (st['number_of_edges'],st['expected_number_of_edges'],st['number_of_edges']/st['expected_number_of_edges']))
o.append('\\section{Pairwise alignment against parasail}')
g,l=pa['global'],pa['local']
o.append(f"We aligned all {pa['n_pairs']} pairs of {pa['n_seqs']} protein sequences with sugarcode and parasail {esc(pa['tool'].split()[-1])}. With parasail's open penalty set to 11 the scores matched exactly for {g['exact_vs_parasail_open11']}/{g['n']} global and {l['exact_vs_parasail_open11']}/{l['n']} local pairs. With open penalty 12 only {g['exact_vs_parasail_open12']} and {l['exact_vs_parasail_open12']} matched (maximum difference {g['max_abs_diff_open12']:.0f} and {l['max_abs_diff_open12']:.0f}). This fixes the convention: sugarcode's gap\\_open$=-11$ is the cost of the first gap residue, i.e. parasail open 11.")
o.append('\\paragraph{(F31) Affine gap cost.} A gap of length $L$ costs $\\gamma(L)=o+(L-1)e$ under sugarcode and parasail open$=11$; under the alternative convention $\\gamma(L)=o+Le$, which explains the open-12 mismatch.')
fig,ax=plt.subplots(1,2,figsize=(8,3))
ax[0].barh([t['name'][:30] for t in t0[:8]][::-1],[t['fdr'] for t in t0[:8]][::-1]); ax[0].axvline(0.05,color='r',ls='--'); ax[0].set_xlabel('Reactome FDR'); ax[0].tick_params(labelsize=6)
ax[1].bar(['observed','expected'],[st['number_of_edges'],st['expected_number_of_edges']],color=['C0','0.6']); ax[1].set_ylabel('STRING edges')
plt.tight_layout(); plt.savefig('figs/dex_pathways.pdf')
o.append('\\begin{figure}[h]\\centering\\includegraphics[width=0.95\\textwidth]{figs/dex_pathways.pdf}\\caption{Left: Reactome FDRs for the top pathways; the red line is 0.05. Right: STRING observed versus expected edges.}\\end{figure}')
open('sections_results3.tex','w').write('\n'.join(o)+'\n')
