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
