# Verified success inventory - pass 1 (2026-10-01)

Scope: 56 science repos under uditakankananonononono (mega27-*, selective-transfer-certificate). Public repos cloned at HEAD; private repos (01, 18, 27, selective-transfer) read through the GitHub API, README only.
Rule: a success needs a measured positive in a committed result file. Negatives, FAILs, scoped behaviors, README-only claims and reproductions that match a published number do not count.

## Count
- Tier A, file-verified successes: **5**
- Tier B, real but partial or caveated (not counted in the 5): 2
- Tier C, README or memory claims not yet checked against a result file (not counted): 8 candidates
- Pass 1 covered about 10 of 56 repos at file level. This is a floor, not the program total.

## Tier A (counted)
| # | Repo | Unit | Result | Numbers | Evidence |
|---|---|---|---|---|---|
| 1 | mega27-05-yeast-metabolic-twin | GCN essentiality ranking vs FBA-rule benchmark | Beats the published-benchmark FBA rule on AUC | paired delta AUC +0.1347, CI95 [0.114, 0.155]; README: AUC 0.831 vs 0.697 | results/cv_auc_delta.json @5612baa |
| 2 | mega27-05-yeast-metabolic-twin | Isozyme over-rescue discovery | FBA-missed essentials carry isozyme backups far more often | 20.2% vs 1.5%, Fisher OR 16.2, p=3.75e-4, 94 missed vs 65 caught, 19 genes named | results/overrescue_audit.json @5612baa |
| 3 | mega27-12-3d-drug-discovery | Vina redocking gate, SARS-CoV-2 Mpro (7KX5) | Protocol validated against crystal pose | RMSD 0.970 A (gate <2.0 A), affinity -9.319 kcal/mol | results/redock_7KX5.json @fb717cd |
| 4 | mega27-02-virtual-cell | E. coli stacked ensemble vs single layers | Ensemble beats every single layer on AUROC | OOF AUROC 0.7225 CI95 [0.681, 0.761] vs best single iJO1366 FBA-min 0.666 | results/ensemble_results.json @cb63e26. Caveat: README notes McNemar loses hard calls 27:50 (p=0.012) |
| 5 | mega27-09c-peptide-solubility-anticancer | BLAST novelty screen of designed peptides | 32 of 32 candidates novel | 32 novel, 0 near, 0 not novel (rule: no hit >=90% id over >=80% length; blastp nr + mmseqs2 cross-check) | results/c6_blast_novelty.json @ffbf74e |

## Tier B (real, not counted)
- mega27-12: GNN rescoring of docked poses, LOO AUROC 0.833 (5 seeds, sd 0.010) vs raw Vina 0.534. n=59 ligands, LOO, in-house set. results/bigscreen_analysis.json.
- mega27-09c: pepx-GNN aggregation AUC 0.8391 vs FoldAmyloid 0.7480 on pep424 (+0.091). Strongest comparators (PASTA 2.0 reported 0.8573 on its own benchmark, AmyloGram) not yet run on the same partition, so the beat gate is open. docs/COMPARATOR_AGGREGATION.md.

## Tier C (unverified, need a result-file check)
- mega27-05: dsdh2 dach1 double knockout, obligatory succinate 0.0892 at 99.2% WT growth, +49% over Raab-2010 quadruple (README). A5 strain-design verdict file shows candidates "screened_out" with filter B pending, so check before counting.
- mega27-03 virtual organoid: size-aware attenuation 0.249 [0.165, 0.326], 14 donors (README). Same README says the stricter size test and specificity check FAILED, so likely exploratory only.
- mega27-07: codon expression CNN Spearman 0.61 held out (README); reproducible CLI checkpoint gives 0.564. Need codon_benchmark.json vs CAI-greedy.
- mega27-14 digital embryo: embryosim reproduces 0.781% fragility basin and 2.39% vs Liu 4.94% information threshold (README). Reproduction, plus the corrected clamp effects (memory). Check results.
- mega27-18 mirna (private): CNN top non-conserved candidates less dosage-sensitive than matched genes, survives controls (README). Memory says lane-18 replicated positive on GSE97056. Need result files.
- mega27-04 NASA Blautia 16S vs shotgun discrepancy (memory, candidate-new). Repo has more recent commits than I read.
- mega27-19 xenobot stage-1, 7 audit-clean designs (memory). The repos I read show near-zero forecasting gain and no discovery, so locate the stage-1 result.
- mega27-09a / 09b / 11a: need result-file pass.

## Checked and excluded
- mega27-06b transfer: README says CNN human to mouse Spearman 0.78 vs 0.63. results/transfer.json shows CNN 0.665 human held-out and about 0.358 mouse transfer. README does not match the file. Not counted; README needs fixing.
- mega27-02 Bernstein reproduction 0.839 vs published 0.843: reproduction only, verdict "negative so far".
- mega27-09c HEMOPI1 MCC 0.927 vs comparator 0.93: no beat. DeepSol arm: all arms beat=false.
- mega27-17s G-A1/G-A5 sanity passes, G-A3 benchmark FAIL (0.014 vs ST-Net 0.19). Sanity gates are not wins.
- mega27-20 R3: test CI 0.723 vs >0.878 gate, terminal negative.
- mega27-22 sleep EEG, selective-transfer-certificate: no measured result by their own READMEs.
- mega27-13b: simulations only, README says gates unmet.

## Method note
No tool here can push; this file is delivered to the parent for Key-ops to commit additively to mega27-01-sugarcode-realdata-validation (e.g. audits/verified_success_inventory_2026-10-01.md).
