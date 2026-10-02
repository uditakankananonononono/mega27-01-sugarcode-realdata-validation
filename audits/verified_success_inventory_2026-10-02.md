# Verified success inventory - sweep 1 final (2026-10-02)

Strict count: 9 of 100 (baseline-beat leg). Two carry flags (09a marginal, 11b baseline-only). Provisional finding-leg-only: 3. Everything else: no.
Rule: committed result file shows a measured positive. Memory and README claims are not evidence.
Coverage: 56 repos. Deep file-level: 01 (README only, private), 02, 03, 04, 05, 07, 08, 09a, 09b, 09c, 10, 11a, 11b, 12, 14, 18 (partial), 21, 25 family, 06b. Grep screen only (no positive found, not exhaustively read): 10b, 13, 13b, 15, 16, 17, 17s, 19 x2, 20, 22, 23a, 23b, 24 family, 26b, selective-transfer, 06/06a, 27, 01.

## STRICT (9)
1. 05 yeast - GCN beats FBA-rule essentiality benchmark. delta AUC +0.1347 CI [0.114,0.155]. results/cv_auc_delta.json @5612baa
2. 05 yeast - isozyme over-rescue finding. 20.2% vs 1.5%, OR 16.2, p=3.75e-4. results/overrescue_audit.json
3. 12 Mpro - Vina redock RMSD 0.970 A (gate <2.0), -9.319 kcal/mol. results/redock_7KX5.json @fb717cd. Protocol validation, not a discovery.
4. 02 virtual cell - ensemble OOF AUROC 0.7225 [0.681,0.761] vs best single 0.666. results/ensemble_results.json @cb63e26. McNemar loses hard calls 27:50.
5. 09c - 32/32 designed peptides novel (blastp nr + mmseqs2). results/c6_blast_novelty.json @ffbf74e. Novelty leg, no baseline leg.
6. 14 embryo - CNN CF decoder RMSE 2.38% EL (R2 0.854) vs Liu 4.94% (R2 0.374), disjoint-line 5-fold CV. results/cf_decoder.json
7. 07 pancreatic - TP53 external CPTAC-PDAC n=140 AUC 0.741, AP 0.896 vs prevalence 0.75. results/cptac_external_validation.json. TP53 only; KRAS 0.384, CDKN2A 0.497, SMAD4 0.347 do not count.
8. 09a AMP [FLAG: marginal] - beats published Feb2020 CV on 4/5 metrics: AUC 0.9636 vs 0.962, MCC 0.8076 vs 0.799, acc 0.903 vs 0.899, spec 0.907 vs 0.891; sens 0.899 vs 0.906 (below). Within 1 SD. results/study14_feb2020_full.json. Cluster-held-out MCC 0.816 (study23); replicated 2026-10-02 on a second fold seed (fold 101, CNN seeds 101/102), margins within 1 SD of published (mega27-09a results/study14_repeat_foldseed101.json, commit ec1b2c2c)
9. 11b protein [FLAG: baseline-only] - Ssym inverse r 0.592 [0.508,0.674] vs PoPMuSiCsym 0.48; sigma 1.43 vs 1.62. results/ssym_benchmark_v2_e300.json. One seed; does not generalize to SOD1.; baseline beat replicated 2026-10-02 on independent seed 1: r_inv 0.660 (seed 0: 0.592) vs PoPMuSiCsym 0.48, sigma_inv 1.258 vs 1.62 (mega27-11b results/ssym_benchmark_v2_e300_seed1.json); stays baseline-only - discovery leg missing, P53 transfer weak; second benchmark (Broom, n=605): r 0.597 vs DDGun3D 0.62, no beat

## PROVISIONAL (finding leg only or internal baseline)
- 04 keystone: sulfate-reducer enrichment q=0.019, GTDB replicate q=0.020, rests on 3 genera; anaerobe share rho 0.30 q=2.8e-5, perm p=2e-4. results/finding_keystone.md. No named baseline. Exploratory.
- 03 organoid: size-dependent attenuation replicates within cohort, later portion 12/12 donors positive, median 0.199, sign p=2.4e-4 (earliest 6/8, p=0.14); H1 PASS. results/withincohort_replication.json. README: stricter donor x plate x dose test FAILED (p=0.073) and specificity check failed. Exploratory assay analysis, no baseline.
- 18 miRNA: DuplexCNN held-out-gene Pearson 0.814 vs internal ridge+context 0.779. results/context_pp_benchmark.json. miRDB beats it on miRTarBase per README. LOEUF finding file not located.; random-weight falsifier 2026-10-02: untrained CNN reproduces ~80% of the trained LOEUF effect in 1 of 5 seeds, nothing in 2, opposite significant sign in 1; direction seed-dependent; trained weights not in repo - stays post-hoc exploratory, not counted

## NOT COUNTED (verdict - gap)
- 12 GNN rescoring LOO AUROC 0.833 vs Vina 0.534 - n=59 in-house set, no published comparator.
- 09c aggregation GNN 0.839 vs FoldAmyloid 0.748 - PASTA 2.0 (0.857) and AmyloGram not run on pep424; gate open.
- 09c HEMOPI1 MCC 0.927 vs 0.93 - no beat. DeepSol arm - beat=false on all arms. CamSol arm - AUC 0.621, gate not shown met.
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
10b, 13, 15, 16, 17, 23a, 23b, 26b, 01 (private), 18 LOEUF, 04 NASA Blautia, 19 stage-1, 23a stage_b assignment_audit PASS (audit, not a result).
