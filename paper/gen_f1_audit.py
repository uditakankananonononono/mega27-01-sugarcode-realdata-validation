"""sections_f1_audit.tex: F1 reference-tool reproducibility audit chapter.
Sources: discovery/crisprscan_audit/*.json + verification/sweep_reverify_20260926.json.
Failure-framing rule: lead with verified positives; compact limitations block."""
import json
A = '../discovery/crisprscan_audit/'
enum = json.load(open(A+'round6_enumeration_audit.json'))
agree = json.load(open(A+'round6_crisprscan_agreement.json'))
rs2 = json.load(open(A+'round7_rs2_fixture_audit.json'))
rs1 = json.load(open(A+'round7_doench2014_param_diff.json'))
fb = json.load(open(A+'round6_foldback_bioimpact.json'))
d1m = json.load(open(A+'round7_d1_mapping.json'))
d1 = json.load(open(A+'round7_bioimpact_d1.json'))
doe = json.load(open(A+'round7_bioimpact_doench.json'))
doeg = json.load(open(A+'round7_doench_pergene.json'))
extm = json.load(open(A+'round8_extscreen_mapping.json'))
ext = json.load(open(A+'round8_bioimpact_extscreen.json'))
cfdv = json.load(open('../verification/sweep_reverify_20260926.json'))
locus = enum['locus'].replace('_', '\\_')
o = []
o.append('\\section{Reference-tool reproducibility audit: the CRISPR scoring stack}')
o.append('\\paragraph{Motivation.} CRISPR guide selection runs on a stack of reference implementations: an enumerator that lists candidate off-targets, scoring models that rank them (CFD, MIT), and efficiency models that rank the guides themselves (CRISPRscan, Rule Set 1, Rule Set 2/Azimuth). A silent error in any layer propagates into every downstream experiment. This chapter audits each layer of that stack against independent re-implementations and official test fixtures, quantifies the one defect found, and measures its biological consequence across three independent activity datasets.')
o.append('\\subsection{Enumeration layer}')
o.append(f"On the TP53 locus ({locus}, {enum['length']:,} bp) our independent enumerator reproduces the reference ({enum['crispor_guides']:,} candidate sites) exactly: {enum['agree_keys']:,}/{enum['crispor_guides']:,} sites agree on key and sequence, with zero sites unique to either side and {enum['seq_mismatch_same_key']} sequence mismatches at shared keys. Edge-of-locus guides were checked explicitly. The enumeration layer of the reference stack is reproducible.")
o.append('\\subsection{CFD off-target scoring layer}')
c = cfdv['execution_level_sample']['cfd_offtarget']
o.append(f"Against the official Doench 2016 mismatch and PAM weight tables (the authors' own pickle files, driven by a Python-3 port of the published scoring logic), our CFD implementation matches on {c['exact_1e-12']:,}/{c['n_scores_compared']:,} scores (600 guide--off-target pairs $\\times$ 16 PAM contexts) to within $1\\times10^{{-12}}$ (max abs diff {c['max_abs_diff']:.2e}). The CFD layer is exact.")
o.append('\\subsection{CRISPRscan efficiency layer: the AA19 defect}')
o.append(f"Scoring all {agree['n_guides_scored']:,} TP53 NGG guides with both the production reference implementation and our corrected re-implementation shows two distinct discrepancies. First, {agree['n_different_gt_0.005']:,} guides ({100*agree['frac_different']:.1f}\\%) differ at the third significant figure, but almost all of these trace to integer truncation of the production score (documented). Second, and consequential, {agree['n_defect_only_apples_to_apples']} guides ({100*agree['frac_defect_only']:.1f}\\%) carry a genuine scoring defect: the reference reads the AA dinucleotide feature one nucleotide off its published position (model term AA19; upstream issue crisporWebsite\\#76 with minimal reproducer). Because the defect is a single-feature misplacement, its effect on any guide is analytic: $s_{{buggy}} = s_{{correct}} + \\beta_{{AA19}}\\,(\\mathbb 1_{{bug}} - \\mathbb 1_{{correct}})$, an identity we verified exactly on all {fb['identity_check']['n']:,} TP53 guides against both real implementations ({fb['identity_check']['within_1_int']:,}/{fb['identity_check']['n']:,} within one integer point).")
o.append('\\paragraph{Rank consequence.} ')
a = fb['expA_rank_perturbation']
tk = a['topk_overlap']
o.append(f"On TP53, corrected versus defective rankings correlate at Spearman {a['spearman_buggy_vs_correct']}, Kendall $\\tau$ {a['kendall_tau']}, but the shortlists researchers actually synthesize change materially: top-10 overlap {tk['10']}/10, top-25 {tk['25']}/25, top-50 {tk['50']}/50, top-100 {tk['100']}/100, top-500 {tk['500']}/500. Roughly 1 in 10--25 shortlisted candidates is displaced by the defect.")
o.append('\\subsection{Rule Set layers}')
o.append(f"Rule Set 2 (Azimuth): our Python port matches Microsoft's official 948-guide test fixture to a maximum absolute error of {rs2['full_model']['max_abs_err']:.2e} on both the full and no-position models ({rs2['full_model']['n']}/948 guides; one guide rejected on input validation). Rule Set 1: our re-implementation agrees with the reference on all {rs1['sugarcode_n']}/70 published model parameters with zero weight mismatches and zero positional shifts; we describe this as reproduction with shared code lineage, not independent validation, and we do not claim the latter.")
o.append('\\subsection{Biological consequence of the defect}')
o.append('We tested whether correcting the defect improves prediction of experimentally measured guide activity on three independent datasets spanning three species, three laboratories, and three assay types. In every case both the production (defective) and corrected scores were computed from true genomic 35-nt contexts; for the Doench 2016 set, 30-mers were padded (2 nt on the 5-prime side, 3 nt on the 3-prime side) with the padding touching only edge features shared by both implementations, and two different paddings agreed.')
rows = [
 ('Zebrafish in-vitro library (Moreno-Mateos 2015, d1)', d1['n_uniquely_mapped'], d1['defect_only_no_truncation']['n_defect_affected'], d1['defect_only_no_truncation']['delta'], d1['defect_only_no_truncation']['bootstrap_ci95']),
 ('Human FACS peptide (Doench 2016)', doe['n'], doe['defect_only']['n_defect_affected'], doe['defect_only']['delta'], doe['defect_only']['bootstrap_ci95']),
 ('Human essentiality HL60 (Labuhn 2017 ext.)', ext['n'], ext['n_defect_affected'], ext['hl60']['delta_defect_only'], ext['hl60']['delta_defect_only_ci95']),
 ('Human essentiality KBM7 (Labuhn 2017 ext.)', ext['n'], ext['n_defect_affected'], ext['kbm7']['delta_defect_only'], ext['kbm7']['delta_defect_only_ci95']),
]
o.append('\\begin{table}[h]\\centering\\small\\begin{tabular}{lrrrl}\\hline dataset & $n$ & defect-affected & $\\Delta$ Spearman & 95\\% CI\\\\\\hline')
for name, n, na, dlt, ci in rows:
    o.append(f"{name} & {n:,} & {na} ({100*na/n:.1f}\\%) & {dlt:+.4f} & [{ci[0]:+.4f}, {ci[2]:+.4f}]\\\\")
o.append('\\hline\\end{tabular}\\caption{Defect-only consequence on bulk activity prediction: corrected minus defective Spearman correlation with measured activity. Positive = correction helps. All estimates use the identity-verified analytic defect delta; intervals are bootstrap 95\\% over guides.}\\end{table}')
o.append(f"Every point estimate is positive and every effect is small ($\\leq$0.007); the $n$-weighted mean delta is about +0.004. Only the Doench pooled interval excludes zero, and only barely; a within-gene refinement on the same data ({doeg['sign_test']['n_positive']}/{doeg['sign_test']['n_genes']} genes positive, sign test $p={doeg['sign_test']['two_sided_p']}$; weighted mean {doeg['weighted_mean_delta']:+}) has an interval that grazes zero. Our conclusion: the defect imposes a small, directionally consistent penalty on bulk prediction, but its material harm is at the candidate-selection layer, where 6--8\\% of guides are mis-scored and 1 in 10--25 shortlist picks changes.")
o.append('\\paragraph{Library-composition observation.} ')
o.append(f"While mapping the Labuhn external screen library to hg19, we found {extm['multi']} of {extm['n_guides']:,} guides ({100*extm['multi']/extm['n_guides']:.0f}\\%) target multi-copy or repeat loci with valid NGG sites on multiple canonical chromosomes (verified tandem gene-family repeats). These guides were excluded from the unique-locus analysis, but their prevalence means candidate-selection scoring errors compound for a substantial fraction of a real library.")
o.append('\\subsection{Limitations}')
o.append('(i) The zebrafish library activity comes from the same study that built the model, so in-domain correlations may be inflated. (ii) The Doench baseline correlation is near floor (Spearman $\\approx$0.03), which inflates the relative reading of its small delta; Percent Peptide is a proxy activity measure and duplicate 30-mers across genes induce row non-independence (addressed by the within-gene refinement). (iii) 1,339 of 3,139 Labuhn guides did not map to hg19 (older-assembly design or PAM loss) and 712 multi-copy guides were excluded rather than analysed. (iv) RS1 agreement is shared-lineage reproduction. (v) The in-vivo zebrafish validation set (35 guides, 1 affected) is structurally underpowered and reported as such.')
o.append('\\paragraph{(F49) Defect identity.} For guide $g$, $s_{buggy}(g)=s_{correct}(g)+\\beta_{AA19}(\\mathbb 1[g_{12:13}=AA]-\\mathbb 1[g_{13:14}=AA])$ with $\\beta_{AA19}=-0.097377$ from the published model table; verified on 4,824/4,824 TP53 guides within one integer point of the production implementation.')
o.append('\\paragraph{(F50) Consequence estimator.} Per dataset, $\\Delta=\\rho(s_{correct},y)-\\rho(s_{buggy},y)$ with $\\rho$ the Spearman correlation against measured activity $y$; 95\\% intervals from 5,000--10,000 guide-level bootstrap resamples; negative controls apply the same-magnitude perturbation at uniformly random dinucleotide offsets (500 trials).')
open('sections_f1_audit.tex','w').write('\n'.join(o)+'\n')
print('written', sum(len(x) for x in o), 'chars')
