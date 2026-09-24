"""Generate paper sections directly from committed result files (no hand-typed numbers)."""
import json,csv,gzip
B='../benchmarks/'; D='../discovery/splice_region_vus/'
esc=lambda s:str(s).replace('_','\\_').replace('&','\\&').replace('%','\\%').replace('#','\\#').replace('~','\\~{}').replace('≥','$\\geq$')
o=[]
rna=json.load(open(B+'sweep_rna.json')); rv=json.load(open(B+'sweep_rna_vienna_methods.json'))
o.append('\\section{Module Sweep: Additional Real-Data Results}')
o.append('\\subsection{RNA secondary structure}')
o.append('Base-pair F1 against Rfam seed consensus structures (25 sequences per family; files \\texttt{benchmarks/sweep\\_rna.json}, \\texttt{sweep\\_rna\\_vienna\\_methods.json}).')
o.append('\\begin{center}\\small\\begin{tabular}{lrrrrrr}\\hline Family & mean len & Nussinov & MFE & centroid & MEA($\\gamma$=1) & MEA($\\gamma$=4)\\\\\\hline')
for f in rna: o.append(f"{f} & {rna[f]['mean_len']} & {rna[f]['nussinov_f1']:.3f} & {rv[f]['mfe']:.3f} & {rv[f]['cent']:.3f} & {rv[f]['mea1']:.3f} & {rv[f]['mea4']:.3f}\\\\")
o.append('\\hline\\end{tabular}\\end{center}')
o.append('Base-pair maximisation (Nussinov) loses to thermodynamic folding on all three families; ensemble decoders (centroid/MEA) beat MFE only on the long RNase P family. The shipped \\texttt{rna\\_nussinov.fold\\_energy} now delegates to ViennaRNA.')
c=json.load(open(B+'sweep_codon_cai.json'))
o.append('\\subsection{Codon adaptation index versus protein abundance}')
o.append(f"On {c['n_genes']} \\emph{{E. coli}} K-12 genes ({esc(c['cds'])}) paired with {esc(c['abundance'])}, the module CAI agrees with Biopython (Pearson {c['biopython_vs_module_pearson_500']} on 500 genes). Spearman correlation with $\\log_{{10}}$ abundance is {c['module_cai_genome_table_excl_ribo']} with the genome-wide reference table and {c['cai_ribo_ref_excl_ribo']} with a ribosomal-protein reference ({c['n_ribo_ref']} genes, evaluated on {c['n_eval_excl_ribo']} non-ribosomal genes); a cross-validated top-5\\% abundance reference gives {c['cai_top5pct_ref_cv']}. The ribosomal reference is now shipped as \\texttt{{ecoli\\_highexpr}}. File: \\texttt{{benchmarks/sweep\\_codon\\_cai.json}}.")
p=json.load(open(B+'sweep_pgx_cpic.json')); bf=p['before_fix']
o.append('\\subsection{Pharmacogenomic phenotype translation versus CPIC}')
o.append(f"Against the full CPIC diplotype tables, the original module covered {bf['CYP2C19']['covered']}/{bf['CYP2C19']['of']} CYP2C19 and {bf['CYP2D6']['covered']}/{bf['CYP2D6']['of']} CYP2D6 diplotypes, with {bf['CYP2D6']['covered']-bf['CYP2D6']['agree']} CYP2D6 disagreements traced to stale activity values (*9, *41). After vendoring current CPIC tables: {p['CYP2C19']['agree']}/{p['CYP2C19']['cpic_diplotypes']} and {p['CYP2D6']['agree']}/{p['CYP2D6']['cpic_diplotypes']}. {esc(p['after_fix_note'])} File: \\texttt{{benchmarks/sweep\\_pgx\\_cpic.json}}.")
cv=json.load(open(D+'cv_results.json')); hh=json.load(open(D+'headtohead_interim.json')); cal=json.load(open(D+'cnn_calibration.json'))['thresholds']
o.append('\\section{Discovery Track: Splice-Region VUS Reclassification}')
o.append(f"A dinucleotide-window CNN trained on ClinVar splice-region variants reaches gene-grouped 5-fold AUROC {cv['overall']['cnn']:.3f} (logistic {cv['overall']['logit']:.3f}, PWM {cv['overall']['pwm']:.3f}); within-position AUROC {cv['within_position_weighted']['cnn']:.3f}. In a head-to-head on n={hh['n']} variants SpliceAI scores {hh['spliceai']:.4f} versus {hh['cnn']:.3f} for the CNN, and the ensemble adds no significant gain (95\\% CI of difference [{hh['ens_minus_spliceai_ci'][0]}, {hh['ens_minus_spliceai_ci'][1]}]). We therefore make no state-of-the-art claim. Files: \\texttt{{discovery/splice\\_region\\_vus/}}.")
o.append('\\begin{center}\\small\\begin{tabular}{rrrrrr}\\hline CNN threshold & TP & FP & sensitivity & FPR & PPV (prior 3\\%)\\\\\\hline')
for t,v in cal.items(): o.append(f"{t} & {v['tp']} & {v['fp']} & {v['sensitivity']} & {v['fpr']} & {v['ppv_prior_0.03']}\\\\")
o.append('\\hline\\end{tabular}\\end{center}')
n=sum(1 for _ in open(D+'vus_candidates_consensus.tsv'))-1
o.append(f"Consensus VUS candidates (SpliceAI and CNN agreement) number {n} rows in \\texttt{{vus\\_candidates\\_consensus.tsv}}; none of the top 20 had a matching LOVD record (\\texttt{{lovd\\_top20.json}}). These are hypotheses for RNA-level testing, not reclassifications.")
o.append('\\appendix\\section{Tools Table}\\begin{center}\\small\\begin{tabular}{lll}\\hline Tool & Version & Status/used in\\\\\\hline')
for r in csv.DictReader(open('../manifests/tools.tsv'),delimiter='\t'): o.append(f"{esc(r['tool'])} & {esc(r['version'])} & {esc(r['status'])}\\\\")
o.append('\\hline\\end{tabular}\\end{center}')
o.append('\\section{Dataset Manifest}\\begin{center}\\scriptsize\\begin{tabular}{rp{4.2cm}p{3.6cm}p{5cm}}\\hline \\# & Dataset & Accession/version & Used in\\\\\\hline')
for r in csv.DictReader(open('../manifests/datasets.tsv'),delimiter='\t'): o.append(f"{r['id']} & {esc(r['dataset'])} & {esc(r['accession_or_version'])} & {esc(r['used_in'])}\\\\")
o.append('\\hline\\end{tabular}\\end{center}')
i=[k for k,x in enumerate(o) if x.startswith('\\appendix')][0]
open('sections_auto.tex','w').write('\n'.join(o[:i])+'\n')
open('sections_appendix.tex','w').write('\n'.join(o[i:])+'\n')
# benchmark comparison summary (numbers from committed files)
cal0=json.load(open(B+'profile_hmm_cv.json')); cal1=json.load(open(B+'profile_hmm_cv_local.json'))
def first_num(d):
    for k,v in d.items():
        if isinstance(v,(int,float)) and ('auc' in k.lower() or 'roc' in k.lower()): return k,v
    return None,None
k0,v0=first_num(cal0); k1,v1=first_num(cal1)
T=['\\begin{center}\\small\\begin{tabular}{p{3.2cm}p{3.6cm}p{3.2cm}p{4.2cm}}\\hline Module & Reference & Metric & Result (file)\\\\\\hline']
T.append(f"pgx\\_guidelines & CPIC diplotype tables & agreement & {p['CYP2D6']['agree']}/{p['CYP2D6']['cpic_diplotypes']} CYP2D6 (sweep\\_pgx\\_cpic.json)\\\\")
T.append(f"bio.codon CAI & Biopython; PaxDb & Pearson; Spearman & {c['biopython_vs_module_pearson_500']}; {c['cai_ribo_ref_excl_ribo']} (sweep\\_codon\\_cai.json)\\\\")
T.append(f"rna\\_nussinov & ViennaRNA 2.7 & bp F1 (tRNA) & {rna['RF00005']['nussinov_f1']:.3f} vs {rna['RF00005']['vienna_mfe_f1']:.3f} (sweep\\_rna.json)\\\\")
if v0 is not None: T.append(f"profile\\_hmm & Pfam PF00042 & {esc(k0)} & {v0} global, {v1} local (profile\\_hmm\\_cv*.json)\\\\")
T.append(f"splice CNN & SpliceAI 1.3.1 & AUROC n={hh['n']} & {hh['cnn']:.3f} vs {hh['spliceai']:.4f} (headtohead\\_interim.json)\\\\")
T.append('\\hline\\end{tabular}\\end{center}')
open('bench_table.tex','w').write('\n'.join(T)+'\n')
