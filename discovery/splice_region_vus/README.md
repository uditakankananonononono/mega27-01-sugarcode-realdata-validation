# Splice-region VUS reclassification (item-1 discovery angle)

Pipeline (all free/public data, reproducible):
1. `pull.sh` - stream NCBI ClinVar variant_summary.txt.gz (GRCh38), keep intronic SNVs at c.N+/-3..20 -> 428,164 rows (2026-09-24 release).
2. `extract.py` - hg38.2bit (UCSC) 81-nt windows, strand from canonical GT/AG at the implied natural site; donor +3..+6 and acceptor -3..-14 -> 139,794 variants (2,883 P/LP, 97,034 B/LB, 39,877 VUS).
3. `pwm_score.py` - sugarcode deepsplice variant_effect delta; PWM training genes flagged and excluded.
4. `cnn_cv.py` - gene-grouped 5-fold CV (md5(gene) % 5), B subsampled to 4x P per site (seed 7): PWM vs logistic vs CNN.
5. `spliceai_run.py` - official SpliceAI weights (pip spliceai, tensorflow-cpu 2.15.1), D=50, on a seeded 1,500-variant subset of the same CV set.

CV result (cv_results.json): AUC PWM 0.889, logistic 0.912, CNN 0.936; within-position weighted PWM 0.809, CNN 0.900.
Caveat: ClinVar splice-region labels partly derive from in-silico evidence (often SpliceAI), which biases any ClinVar benchmark toward that tool.
Large inputs (ClinVar extract, hg38.2bit) are not committed; the scripts re-fetch them.

## Head-to-head vs SpliceAI (interim, n=757 of the 1,500 seeded subset, 185 P)
AUC: SpliceAI 0.978 (95% bootstrap CI 0.964-0.988) vs our CNN 0.929 (0.907-0.948). The rank-average ensemble reaches 0.978; the ensemble-minus-SpliceAI CI is -0.009 to +0.008, i.e. no gain. Verdict: SpliceAI beats our model on ClinVar splice-region variants. This is a preserved negative, with the circularity caveat above (ClinVar submitters use SpliceAI as evidence).

## Named candidate list (vus_candidates_consensus.tsv)
- Our final CNN scored 37,359 ClinVar VUS (non-training genes); the top 300 went to SpliceAI.
- 267/300 have a SpliceAI delta >= 0.5, and 206 (180 genes) have >= 0.8, SpliceAI's high-precision band.
- Each row gives VariationID, gene, HGVS, GRCh38 position, PWM delta, CNN probability and the SpliceAI DS_AG/AL/DG/DL.
- These are predictions, not discoveries. Each is falsifiable by minigene/RT-PCR or by future ClinVar reclassification.
- CV calibration at CNN >= 0.98: PPV 0.975 at the CV class mix, 0.78 at a 3% pathogenic prior.
- The top 300 are 95% donor +5 variants (a position prior), so SpliceAI agreement is what makes each candidate credible.
