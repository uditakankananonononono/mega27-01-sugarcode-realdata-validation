# RS1 verification verdict (2026-09-26, supersedes round8_rs1_independent.json)

## Bottom line
sugarcode's `crispr_opt.doench2014_ontarget` (CRISPOR-lineage 70-param RS1) IS the
published Doench 2014 model: Spearman 0.9936 / Pearson 0.9912 against the authors'
own final-model scores (Supplementary Table 7, 1,841 sgRNAs, downloaded from
nature.com MOESM8 on 2026-09-26). The "second parameter-port defect" candidate
raised in commit 75f83a0 is REFUTED. No RS1 defect exists.

## Why the candidate was raised (void analysis, preserved honestly)
round8_rs1_independent (v2) used the Azimuth-bundle file
`FC_plus_RES_withPredictions.csv` column `Percent Peptide` as the activity
measure. Per-gene Spearman between that column and the true activity (Table 7
within-gene percent-rank in the marker-negative population) ranges from -0.395
(H2-K) to +0.366 (Thy1): it is noise, in neither direction. The file's
`predictions` column is a later-generation (RS2-era) score, per-gene only ~0.68
vs the RS1 published score. Benchmarking any RS1 implementation against those
columns yields garbage; that is what produced the apparent ~null reference
signal (+0.007) and the apparent reimplementation "win" (+0.152). Both numbers
are void.

## Ground-truth results (round9_rs1_reimpl_groundtruth.json)
- shipped RS1 vs published model scores: 0.9936 Spearman (per-gene 0.991-0.995)
- shipped RS1 vs true activity: 0.5417 overall, mean per-gene 0.545 (9/9 positive)
- published model vs true activity: 0.5047 overall, mean per-gene 0.508
- paper-only independent reimplementation (586 one-hot features, L1-SVM selection
  via gene-grouped nested CV, logistic, honest 9-gene holdout): 0.4772 overall,
  mean per-gene 0.472 (9/9 positive); held-out agreement with published model 0.8823

The judge's round-6 option (independent mathematical reimplementation of RS1) is
CLOSED: the reimplementation independently reproduces the published model's
behavior (0.88 held-out agreement) and confirms the shipped port.

## Data-hygiene finding (real, minor)
The widely-used `FC_plus_RES_withPredictions.csv` (Azimuth data bundle) is unfit
for RS1 benchmarking: its activity column does not track the published activity
and its score column is not RS1. One limitations paragraph in the F1 chapter.
