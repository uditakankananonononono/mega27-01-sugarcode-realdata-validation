# Verified success inventory - audited v3 (2026-10-02)

Strict count: 10 strict of 100, 3 provisional. Flagged strict: 09a (marginal), 11b (baseline-only), 02 pathway concordance (confounded). Everything else: no.
Rule: committed result file shows a measured positive vs a baseline and/or a measured new finding. Memory, README and prereg claims are not evidence. The user's bar (baseline beaten AND something new) is met by few entries; legs are named per line.
Audit: every strict and provisional line below was checked at file level (file exists at the cited path/commit, number matches, flag present) on 2026-10-02. Items marked [via routing] cite commits I could not inspect in my clones (they were pushed after my clone); they come from the parent's routing and should be verified at push time.
Coverage: 56 repos. Deep file-level: 02, 03, 04, 05, 07, 08, 09a, 09b, 09c, 10, 11a, 11b, 12, 14, 18, 21, 25 family, 06b. 01 (private) README only. Grep screen only: 10b, 13, 13b, 15, 16, 17, 17s, 19 x2, 20, 22, 23a, 23b, 24 family, 26b, selective-transfer, 06/06a, 27.

## STRICT (10)
1. 05 yeast - GCN beats FBA-rule essentiality benchmark. delta AUC +0.1347 CI [0.1144, 0.1550]. results/cv_auc_delta.json @5612baa (identical at HEAD).
2. 05 yeast - isozyme over-rescue finding. 20.2% vs 1.5%, OR 16.2, p=3.75e-4. results/overrescue_audit.json. Replicated in independent audits: E. coli Keio OR 11.5, E. coli Gerdes 9.0, B. subtilis 22.7, B. subtilis Kobayashi 9.2, M. tuberculosis 7.5, Mtb Griffin 9.07, Mtb Zhang 8.15, S. aureus (p=0.021), S. typhimurium (p=0.0002). NOT replicated: H. pylori n.s. (10 missed / 24 caught, p=0.294). Finding leg; no named published baseline.
3. 12 Mpro - Vina redock RMSD 0.970 A (gate <2.0), -9.319 kcal/mol. results/redock_7KX5.json @fb717cd (identical at HEAD). Protocol validation, not a discovery.
4. 02 virtual cell - ensemble OOF AUROC 0.7225 [0.681, 0.761] vs best single (fba_min) 0.6657. results/ensemble_results.json @cb63e26 (identical at HEAD). McNemar loses hard calls 27:50.
5. 09c - 32/32 designed peptides novel (blastp nr + mmseqs2 Swiss-Prot cross-check). Novel = no hit >=90% identity over >=80% length. results/c6_blast_novelty.json @ffbf74e (identical at HEAD). Novelty leg only, no baseline leg.
6. 14 embryo - CNN CF decoder RMSE 2.38% EL (R2 0.854), disjoint-line 5-fold CV, 250 embryos / 29 lines. results/cf_decoder.json. Comparator caveats: (i) the "Liu 4.94% (R2 0.374)" figure is a REFIT of Liu's threshold model on this same data, not the published number; paired bootstrap gain 2.57 points, CI [1.99, 3.11]. (ii) The published Liu 2013 bar in the file is 10.5%. (iii) Stronger in-file baselines: Liu-form with fitted lnA 2.86%, ridge on lnA 2.55%; the margin over ridge is narrow (0.17 points).
7. 07 pancreatic - TP53 external validation, CPTAC-PDAC n=140: AUC 0.741 [0.643, 0.830], AP 0.896 vs prevalence 0.75. results/cptac_external_validation.json. TP53 only; KRAS 0.384, CDKN2A 0.497, SMAD4 0.347 do not count. The file's own limits: features contain mutation burden, so the association can reflect mutation load; distinct study identifiers are not a patient-identity overlap audit; the model is a simple balanced logistic regression, not validation of the fitted GNN/CNN; genomic inference, not diagnosis.
8. 09a AMP [FLAG: marginal; replicated on a second fold seed, margins within 1 SD] - beats published Feb2020 CV on 4/5 metrics: AUC 0.9636 vs 0.962, MCC 0.8076 vs 0.799, acc 0.903 vs 0.899, spec 0.907 vs 0.891; sens 0.899 vs 0.906 (below). Every delta is within 1 SD (AUC delta 0.0016, SD 0.0066). results/study14_feb2020_full.json. Repeat with fold seed 101 beats on 5/5 (AUC 0.9636, MCC 0.8129, acc 0.906, spec 0.902, sens 0.910): results/study14_repeat_foldseed101.json @ec1b2c2c [via routing]. Two-run average sens 0.9045, 0.0015 below published.
9. 11b protein [FLAG: baseline-only; baseline beat replicated on seed 1 (r_inv 0.592 and 0.660 vs 0.48)] - Ssym inverse r 0.592 CI [0.508, 0.674] vs PoPMuSiCsym 0.48; sigma 1.43 vs 1.62. results/ssym_benchmark_v2_e300.json (identity-only ridge baseline r_inv 0.25). Seed 1: r_inv 0.660, sigma 1.26, direct r 0.689 (committed @b2594b8 / 6da785c3 [via routing]). P53 transfer collapses (r 0.054 on seed 1), no discovery leg. Second benchmark (Broom, n=605): r 0.597 vs DDGun3D 0.62, no beat.
10. 02 pathway concordance [FLAG: confounded with supplement identity, not replicated in 8 other bacteria] - preregistered primary SUPPORTED: median -0.97 on-pathway vs -4.24 off-pathway (650 vs 272 pairs, 26 vs 52 genes), MWU p=5.7e-75, cluster-bootstrap CI of the difference [2.16, 3.82]. results/pathway_concordance.json (last change cb63e26). Finding leg; secondary strata also positive but share the same confound.

## PROVISIONAL (finding leg only or internal baseline) (3)
- 04 keystone: sulfate-reducer enrichment BH q=0.019 (NCBI), GTDB replicate q=0.020, rests on 3 genera (Desulfovibrio, Bilophila, Desulfobulbus), fragile and exploratory; anaerobe share rho 0.30 q=2.8e-5, within-quintile permutation p=2e-4. results/finding_keystone.md. No named baseline.
- 03 organoid: size-dependent attenuation replicates within cohort, later portion 12/12 donors positive, median 0.199, sign p=2.4e-4 (earliest 6/8, p=0.145); H1 PASS. results/withincohort_replication.json. README: stricter donor x plate x dose test FAILED (p=0.073) and specificity check failed. Our segmentation AP50 0.7387 is below the published 0.76. Exploratory assay analysis, no baseline.
- 04 NASA Blautia assay discrepancy: V4 16S higher than shotgun in 20/20 mouse pairs (assigned-only also 20/20; results/nasa_rr5_assay_20260930/blautia_sensitivity.json) and 43/48 fecal pairs in RR-6 (results/nasa_rr6_library_20260930/README.md). Caveat: the 20 pairs come from two cages, so they are not independent (pseudoreplication) [via routing]; no independent replication. Candidate-new finding.

## NOT COUNTED (verdict - gap)
- 12 GNN rescoring LOO AUROC 0.833 vs Vina 0.534 - n=59 in-house set, no published comparator.
- 09c aggregation GNN 0.839 vs FoldAmyloid 0.748 - PASTA 2.0 (0.857) and AmyloGram not run on pep424; gate open.
- 09c HEMOPI1 MCC 0.927 vs 0.93 - no beat. DeepSol arm - beat=false on all arms. CamSol arm - AUC 0.621, gate not shown met.
- 18 miRNA (post-hoc, preregistered hypotheses falsified 13/13): preregistered LOEUF hypotheses H1-H3 all lost 13/13. The reversed direction (CNN top-K less constrained) is post-hoc exploratory only and counts nowhere. Open falsifier run 2026-10-02 (artifacts rw_falsifier_18.json): untrained random-weight DuplexCNNs reproduce the depletion in 1 of 5 seeds (seed 0: +0.148 vs trained +0.186, 13/13 miRNAs), nothing in 3, and the opposite significant direction in 1 (seed 4). The direction is seed-dependent; trained weights are not in git. Also internal ridge baseline only (0.814 vs 0.779); miRDB beats it on miRTarBase per README.
- 11b second benchmark (Broom, n=605): r 0.597 [0.542, 0.645], RMSE 1.834 vs DDGun3D r 0.62 / RMSE 1.66-1.68 and DDGun r 0.52 / 1.77-1.78 (Montanucci 2019 Table 2). No beat. Recorded as a negative.
- 09b B3 EL benchmark 41 alleles AUROC 0.953 - self-set bar (>=0.90); NetMHCpan-4.1 ~0.98 is above it. Not a beat.
- 09b B1 unseen allele - >90% peptide overlap flagged, invalid. CPP classifier AUC 0.927 - no comparator.
- 07 codon optimizer - CNN-opt beats WT (p=8e-5) but CAI-greedy is higher, evaluator same family as optimizer.
- 05 double knockout 0.0892 - not in results files; A5 candidates screened out, filter B pending.
- 06b transfer 0.78 vs 0.63 README - results/transfer.json shows 0.665 human, ~0.358 mouse. README is wrong.
- 02 Bernstein 0.839 vs 0.843 published - reproduction, verdict negative so far.
- 10 malaria - R-M1 mean acc 0.9477 vs gate 0.957 FAIL; R-M2 0.9471 FAIL.
- 11a - audit/negatives only (myoglobin r 0.455 to 0.24-0.29 after decontamination).
- 21 - Cox clin 0.668/full 0.680 beat all deep models; recalibration pivot band checks false in rows read ("POSITIVE" label rejected); GSE2034 AUC 0.557.
- 17s - sanity gates pass, benchmark G-A3 FAIL (0.014 vs 0.19). 20 R3 - 0.723 vs >0.878 FAIL.
- 25 family - PPD methylation 0 FDR hits; PPD GNN 0.695 < RWR 0.803; TB challenger 0.734 < 0.784; no performance claim holds. Memory "45/100 passed" not in any repo.
- 19 xenobot - near-zero forecast gain, no discovery per STATUS_GATES. Stage-1 "7 designs" not found in repos read.
- 22, selective-transfer, 13b, 24 family - no measured result per own docs.
- 08 phage - hybrid beats own k-mer LR (+4.5 to +6.3 AUROC, 3 seeds), but genus-holdout (Salmonella) AUROC 0.67-0.70 vs baseline 0.49 mixed with choice accuracy 0.32-0.43; internal baseline only. Not counted; not fully read.


## Unchecked or partial
10b, 13, 15, 16, 17, 23a, 23b, 26b, 01 (private), 19 stage-1, 23a stage_b assignment_audit PASS (audit, not a result). 03 and 04 independent replications need the parked rented-box run (RR-8), not run.
