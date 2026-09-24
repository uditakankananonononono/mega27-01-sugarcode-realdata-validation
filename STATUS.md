# Independent external validation of SugarCode AI - honest status

## Verified (hermetic tests, green; 10/10 pytest, no-stub audit clean)
- CFD oracle reproduces CRISPOR reference doctests exactly (4/4 values).
- Frozen fixtures: 600 fresh CRISPR pairs + oracle CFD values; 150 Doench
  RS1 30-mer contexts + CRISPOR percentile scores; 500-variant ClinVar
  sample; 150-mutation SKEMPI subset; BiGG e_coli_core; Kazusa table;
  25 globin UniProt sequences; Pfam PF00042.29 globin HMM.
- Independent FBA reference on e_coli_core: 0.873922 h^-1 (published 0.8739).
- Pfam PF00042.29 (HMMER3) separates globin fixtures perfectly (AUC 1.0) -
  fixture labels are sound and the discrimination task is solvable.
- SKEMPI naive baselines: mean-predictor RMSE 1.629, |volume-change|
  Pearson 0.053, hydropathy-loss Pearson -0.129 (benchmarks/baselines.json).

## Module comparisons run (benchmarks/results.json, 2026-09-24, run 2)
All pass or honestly scoped:
- crispr_opt CFD vs oracle, 600/600 pairs exact to 5e-7 (max |diff| 5.0e-07 =
  the module's 6-dp rounding).
- crispr_opt RS1 raw score vs frozen CRISPOR percentiles: Spearman 1.000 (n=150).
- openclinvar vs frozen-500 ClinVar labels (offline, no label leakage): n=324
  called, accuracy 0.997 (150 TP, 173 TN, 1 FN, 0 FP, CI95 0.991-1.000),
  176 abstain by design (missense and +3..+20 splice-region calls stay VUS).
  clinvar_live_resolution: 40/40 sampled names resolve to their own
  VariationID via the transcript-qualified live query. (Run 1 abstained
  500/500 because ClinVar-style names parsed to 'unknown'; sugarcode-ai
  commit 7a80dd5 fixed the parser and the live query.)
- deepsplice PWM on held-out RefSeqGene junctions (NG_008245.1, disjoint
  from the module PROVENANCE): donor AUC 1.000 / acceptor 0.990 vs shuffled
  negatives; donor 0.964 / acceptor 0.911 vs 1,385 real intronic GT and
  1,786 intronic AG decoys >=30 nt from intron ends; consensus-seed matrix
  baseline 0.761. (Run 1's acceptor 0.547 was a harness off-by-one on the
  acceptor frame, not a module defect.)
- virtual_cell FBA parity: module 0.873922 == independent HiGHS 0.873922 ==
  published 0.873922 on BiGG e_coli_core.
- codon CAI vs independent Kazusa geometric-mean recompute: max |diff| 0.0
  over 40 random sequences.
- acmg_bayesian vs Tavtigian point math: posterior matches OP=350^(pts/8)
  exactly over -16..16; bands match the documented prior-0.10 bands; the
  paper's 6-point -> 0.900 identity reproduces.

Not applicable / reference only:
- mutdock ddG vs SKEMPI: no honest adapter exists (module ddG needs ligand
  SMILES + pocket string; SKEMPI rows have neither). Documented gap,
  contextualized by the naive baselines.
- profile_hmm module AUC: blocked - it trains from an MSA and no globin
  alignment ships. Pfam PF00042.29 (HMMER3) reference separates the fixtures
  perfectly (AUC 1.0): task solvable, labels sound, module verdict pending.
