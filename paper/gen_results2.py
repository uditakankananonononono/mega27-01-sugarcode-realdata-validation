"""Sections for the MaxEnt/logit, tAI/ENC, cross-species CAI and docking results.
Every number is read from a committed JSON/JSONL; nothing is typed by hand."""
import json, csv, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np
B='../benchmarks/'; D='../discovery/splice_region_vus/'
def esc(s): return str(s).replace('_','\\_').replace('&','\\&').replace('%','\\%').replace('#','\\#')
ms=json.load(open(B+'sweep_codon_cai_multispecies.json')); te=json.load(open(B+'sweep_codon_tai_enc.json'))
dk=json.load(open(B+'sweep_docking.json')); dr=[json.loads(l) for l in open(B+'sweep_docking_raw.jsonl')]
cv=json.load(open(D+'cv_results_me.json')); mx=json.load(open(D+'maxent_cmp.json')); lm=json.load(open(D+'splice_logit_me_v1.json'))
o=[]
# ---- splice
o.append('\\section{Splice Track II: A 2004 Baseline Beats Our CNN, and What Survives}')
o.append(f"We compared the ClinVar-trained CNN against MaxEntScan \\cite{{yeo2004}} on the same gene-grouped folds (n={mx['n']} variants). "
 f"MaxEntScan's ref$-$alt score difference reached AUROC {mx['maxentscan_delta']}, above the CNN ({mx['cnn']}), the logistic PWM model ({mx['logit']}) and the raw PWM delta ({mx['pwm']}). "
 f"The paired bootstrap 95\\% CI for CNN minus MaxEntScan is [{mx['cnn_minus_maxent_ci95'][0]}, {mx['cnn_minus_maxent_ci95'][1]}], so this is a clear negative result for the CNN. "
 f"The gap is mostly at acceptors (CNN {mx['acceptor']['cnn']} vs MaxEntScan {mx['acceptor']['maxentscan']}, n={mx['acceptor']['n']}); at donors the two are close ({mx['donor']['cnn']} vs {mx['donor']['maxentscan']}, n={mx['donor']['n']}).")
o.append(f"As a pivot we added the three MaxEntScan features to both models. The logistic model then reached {cv['overall']['logit_me']} (CI of gain over MaxEntScan [{cv['logit_me_minus_maxent_ci95'][0]}, {cv['logit_me_minus_maxent_ci95'][1]}]), while the CNN with the same features reached {cv['overall']['cnn_me']} (CI [{cv['cnn_me_minus_maxent_ci95'][0]}, {cv['cnn_me_minus_maxent_ci95'][1]}], not significant). "
 f"Within-position weighted AUROC, which removes the easy signal from offset alone, is nearly flat: MaxEntScan {cv['within_position_weighted']['maxent']}, logistic {cv['within_position_weighted']['logit_me']}, CNN {cv['within_position_weighted']['cnn_me']}. "
 "The gain is real but small. The deployed tool (splice-vus-triage v2) therefore uses the logistic model and reports its calibrated operating points (Table~\\ref{tab:lmcal}).")
o.append('\\begin{table}[h]\\centering\\small\\caption{Logistic PWM+MaxEntScan model: operating points from out-of-fold scores.}\\label{tab:lmcal}\\begin{tabular}{rrrrr}\\hline Threshold & TP & FP & Sensitivity & FPR\\\\\\hline')
for t,c in lm['calibration'].items(): o.append(f"{t} & {c['tp']} & {c['fp']} & {c['sensitivity']} & {c['fpr']}\\\\")
o.append('\\hline\\end{tabular}\\end{table}')
o.append('\\begin{figure}[h]\\centering\\includegraphics[width=0.85\\textwidth]{figs/splice_models.png}\\caption{AUROC by model and site type, gene-grouped 5-fold CV.}\\end{figure}')
fig,ax=plt.subplots(figsize=(6,3)); ms_=['maxent','logit_me','cnn_me']; lab=['MaxEntScan','logistic+ME','CNN+ME']; x=np.arange(3)
for i,s in enumerate(['overall','donor','acceptor']): ax.bar(x+i*0.25-0.25,[cv[s][m] for m in ms_],0.25,label=s)
ax.set_xticks(x); ax.set_xticklabels(lab); ax.set_ylim(0.9,0.98); ax.set_ylabel('AUROC'); ax.legend(fontsize=7); fig.tight_layout(); fig.savefig('figs/splice_models.png',dpi=200)
# ---- cross-species CAI
bh=ms['bh_fdr_rho_ribo_gt0']; bq={r['organism']:r['q'] for r in json.load(open(B+'sweep_codon_bh_fdr.json'))['per_taxon']}
o.append('\\section{Codon Track II: Cross-Species Test of the Ribosomal Reference Set}')
o.append(f"Question: {esc(ms['question'])} Data: {esc(ms['data'])}. Of {ms['n_taxa_scanned']} PaxDb taxa scanned, {ms['n_usable']} had a matched RefSeq assembly and enough genes, and {ms['n_with_ribo_ref']} had a usable ribosomal-protein reference set. "
 f"The ribosomal reference gave the higher Spearman correlation with measured abundance in {ms['ribo_better']}/{ms['n_with_ribo_ref']} taxa (median gain {ms['median_rho_gain']}; Wilcoxon signed-rank $p={ms['wilcoxon_p']:.2g}$, sign test $p={ms['sign_test_p']:.2g}$). "
 f"The per-taxon correlation was significant after Benjamini-Hochberg correction in {bh['n_significant_q05']}/{bh['n']} taxa. Exceptions: {esc(', '.join(ms['exceptions']))}. "
 "The idea goes back to Sharp and Li \\cite{sharpli1987}; what is new here is the systematic multi-taxon test with a pre-stated falsification criterion (the claim fails if the ribosomal reference is not better in a majority of taxa at $p<0.01$).")
o.append('\\begin{figure}[h]\\centering\\includegraphics[width=0.9\\textwidth]{figs/cai_species.png}\\caption{Spearman $\\rho$ (CAI vs PaxDb abundance) with the genome-wide versus ribosomal reference, per taxon.}\\end{figure}')
pt=[t for t in ms['per_taxon'] if t.get('rho_ribo_ref') is not None]
fig,ax=plt.subplots(figsize=(6,3.2)); g=[t['rho_genome_table'] for t in pt]; r=[t['rho_ribo_ref'] for t in pt]
ax.scatter(g,r,s=14); lo=min(g+r)-0.05; hi=max(g+r)+0.05; ax.plot([lo,hi],[lo,hi],'k--',lw=0.8)
ax.set_xlabel('$\\rho$, genome-wide reference'); ax.set_ylabel('$\\rho$, ribosomal reference'); fig.tight_layout(); fig.savefig('figs/cai_species.png',dpi=200)
o.append('{\\scriptsize\\begin{longtable}{lp{4.3cm}lrrrrr}\\caption{Per-taxon results (all taxa with a matched assembly).}\\\\\\hline Taxid & Organism & Assembly & $n$ & $n_{ribo}$ & $\\rho_{gen}$ & $\\rho_{ribo}$ & $q$\\\\\\hline\\endhead')
for t in ms['per_taxon']:
    q=bq.get(t['organism'],'')
    q=f"{q:.2g}" if isinstance(q,(int,float)) else '--'
    o.append(f"{t['taxid']} & {esc(t['organism'])} & {esc(t['assembly'])} & {t['n_matched']} & {t['n_ribo'] if t['n_ribo'] is not None else '--'} & {t['rho_genome_table']} & {t['rho_ribo_ref'] if t['rho_ribo_ref'] is not None else '--'} & {q}\\\\")
o.append('\\hline\\end{longtable}}')
o.append('$q$: Fisher-$z$ two-sided test of $\\rho_{ribo}$, Benjamini-Hochberg adjusted at full precision (benchmarks/sweep\\_codon\\_bh\\_fdr.json).')
o.append('\\subsection{tRNA adaptation index and effective number of codons}')
o.append(f"Tools: {esc(te['tools'])}. Table~\\ref{{tab:tai}} gives Spearman $\\rho$ with abundance. No single index wins everywhere: tAI leads in yeast, the ribosomal-reference CAI leads in \\textit{{E. coli}}, and $-$ENC is weakest in all three organisms.")
o.append('\\begin{table}[h]\\centering\\small\\caption{Codon indices vs abundance (Spearman $\\rho$).}\\label{tab:tai}\\begin{tabular}{p{5cm}rrrrr}\\hline Organism & $n$ & CAI$_{gen}$ & CAI$_{ribo}$ & tAI & $-$ENC\\\\\\hline')
for t in te['results']: o.append(f"{esc(t['organism'])} & {t['n_genes']} & {t['rho_cai_genome']} & {t['rho_cai_ribo']} & {t['rho_tai']} & {t['rho_neg_enc']}\\\\")
o.append('\\hline\\end{tabular}\\end{table}')
# ---- docking
sp=dk['scoring_power']; rd=dk['redocking']; ok=[d for d in dr if 'rmsd_top1' in d and d.get('rmsd_top1') is not None]
o.append('\\section{Docking Track: Real AutoDock Vina Against LP-PDBBind}')
o.append(f"Source: {esc(dk['source'])}. Engine: {esc(dk['engine'])}. {dk['n_ok']}/{dk['n_attempted']} complexes ran; the {len(dk['errors'])} failures are listed in Table~\\ref{{tab:dockerr}}. "
 f"Redocking: top-1 pose RMSD $<2$\\,\\AA{{}} in {rd['top1_rmsd_lt2A']*100:.1f}\\% and best-of-9 in {rd['best_of_9_rmsd_lt2A']*100:.1f}\\% ({esc(rd['note'])}). "
 f"Scoring power is weak for every score, including heavy-atom count alone (Pearson {sp['n_heavy']['pearson']}): Vina top pose {sp['vina_top_E']['pearson']} (95\\% CI {sp['vina_top_E']['pearson_ci95']}), sugarcode's own Vina-like score on the crystal pose {sp['sugarcode_crystal']['pearson']} (CI {sp['sugarcode_crystal']['pearson_ci95']}). "
 "The honest reading: sugarcode's scoring term adds nothing beyond ligand size on this split, so the module now calls real Vina (sugarcode-ai commit e0564e4) instead of claiming affinity prediction.")
o.append('\\begin{figure}[h]\\centering\\includegraphics[width=0.9\\textwidth]{figs/docking.png}\\caption{Left: redocking RMSD distribution. Right: Vina top-pose energy vs experimental pK.}\\end{figure}')
fig,(a1,a2)=plt.subplots(1,2,figsize=(7,3)); a1.hist([d['rmsd_top1'] for d in ok],bins=20); a1.axvline(2,color='k',ls='--'); a1.set_xlabel('top-1 RMSD (\u00c5)'); a1.set_ylabel('complexes')
a2.scatter([d['vina_top_E'] for d in ok],[d['pK'] for d in ok],s=10); a2.set_xlabel('Vina top-pose energy (kcal/mol)'); a2.set_ylabel('experimental pK'); fig.tight_layout(); fig.savefig('figs/docking.png',dpi=200)
o.append('\\begin{table}[h]\\centering\\small\\caption{Docking failures.}\\label{tab:dockerr}\\begin{tabular}{ll}\\hline PDB & Error\\\\\\hline')
for e in dk['errors']: o.append(f"{e['pdb']} & {esc(e['error'])}\\\\")
o.append('\\hline\\end{tabular}\\end{table}')
o.append('{\\scriptsize\\begin{longtable}{lrrrrrrr}\\caption{Per-complex docking results.}\\\\\\hline PDB & pK & heavy & $E_{cryst}$ & $E_{min}$ & $E_{top}$ & RMSD$_1$ & RMSD$_{9}$\\\\\\hline\\endhead')
f=lambda v:('--' if v is None else f"{v:.2f}")
for d in ok: o.append(f"{d['pdb']} & {d['pK']} & {d['n_heavy']} & {f(d['vina_crystal'])} & {f(d['vina_min'])} & {f(d['vina_top_E'])} & {f(d['rmsd_top1'])} & {f(d.get('rmsd_best_of_9'))}\\\\")
o.append('\\hline\\end{longtable}}')
open('sections_results2.tex','w').write('\n'.join(o)+'\n')
# ---- accession manifest appendix
A=['\\section{Accession-Level Dataset Manifest}','Each row is one identifier-backed record fetched individually and used in a committed analysis (counting rule in manifests/README.md).','{\\scriptsize\\begin{longtable}{rllp{3.2cm}p{4.2cm}}\\hline \\# & Kind & Study & Accession & Used in\\\\\\hline\\endhead']
for i,r in enumerate(csv.DictReader(open('../manifests/accessions.tsv'),delimiter='\t'),1): A.append(f"{i} & {esc(r['kind'])} & {esc(r['study_id'])} & {esc(r['accession'])} & {esc(r['used_in'])}\\\\")
A.append('\\hline\\end{longtable}}'); open('sections_accessions.tex','w').write('\n'.join(A)+'\n')
print('ok', len(o), len(A))
# ---- SASA and primer (appended)
sa=json.load(open(B+'sweep_sasa_freesasa.json')); ar=json.load(open(B+'sasa_altloc_recheck.json'))['after_altloc_fix']; pr=json.load(open(B+'sweep_primer_primer3.json'))
okS=[x for x in sa['per_structure'] if 'error' not in x]; bad=[x for x in sa['per_structure'] if 'error' in x]
X=['\\section{Structure and Primer Utilities Against Reference Tools}',
 f"\\paragraph{{SASA.}} Against {esc(sa['tool'])} on {len(sa['per_structure'])} RCSB receptors, the {sa['n']} structures with identical atom sets agreed to a median relative difference of {sa['median_rel_diff']*100:.1f}\\% (max {sa['max_abs_rel_diff']*100:.1f}\\%) with median per-atom Pearson {sa['median_per_atom_pearson']:.3f}. "
 f"The other {len(bad)} exposed a parser bug: every alternate location was kept, so overlapping conformers were all counted. After the fix, atom counts match and per-atom Pearson is "+', '.join(f"{x['per_atom_pearson']} ({x['pdb']})" for x in ar)+".",
 f"\\paragraph{{Primers.}} On {pr['n_oligos']} 20-mers from {', '.join(pr['transcripts'])}, the melting temperature equals Biopython's DNA\\_NN4 result (max difference {pr['tm_vs_biopython_NN4']['max_abs_diff']}). Against primer3 (SantaLucia 1998 table) the mean difference is {pr['tm_vs_primer3']['mean_diff']}\\,$^\\circ$C (max {pr['tm_vs_primer3']['max_abs_diff']}), a table choice. "
 f"The stem-length hairpin heuristic has Spearman {pr['hairpin_stem_vs_primer3_dG_spearman']} with primer3's hairpin $\\Delta G$, a negative result; the module now calls primer3 when installed."]
open('sections_results2.tex','a').write('\n'.join(X)+'\n')
