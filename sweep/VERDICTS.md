# Module sweep verdicts (item 1)

| module | external reference | held-out data | result | verdict |
|---|---|---|---|---|
| chem_descriptors | RDKit 2026.03.6 | 993 ChEMBL phase-4 small molecules (ChEMBL API, 2026-09-24) | mol_wt, exact_mol_wt, hbd, rotatable_bonds, tpsa, fraction_csp3, heavy atoms, formula: 993/993 exact. hba: 993/993 vs the RDKit Lipinski HAcceptor SMARTS the module implements; 867/993 vs `Lipinski.NumHAcceptors` in RDKit 2026.03, which now calls the C++ CalcNumHBA and no longer equals that SMARTS | pass (version divergence documented) |
| chem_similarity | RDKit Morgan r=2, 2048 bits | 500 random pairs of the same 993 drugs | Tanimoto 500/500 exact (max abs diff 0.0) | pass |
| mit_offtarget | CRISPOR `calcHitScore` (verbatim, maximilianh/crisporWebsite) | 600 held-out BRCA1 guide/off-target pairs | 600/600 exact (max diff 1.1e-16) | pass |
| crisprscan_score | CRISPOR `calcCrisprScanScores` (verbatim) | 300 held-out BRCA1 35-nt contexts | 275/300 integer-equal, max diff 10 points | reference discrepancy found (see below) |
| rna_nussinov | ViennaRNA 2.7.2 (Turner 2004) | Rfam 15 seeds RF00005 tRNA, RF00001 5S, RF00010 RNase P; 25 seqs each (seed 7), consensus SS projected | base-pair F1 Nussinov 0.328 / 0.186 / 0.189 vs ViennaRNA MFE 0.714 / 0.531 / 0.545 | exact algorithm but biologically weak - closed by integrating `fold_energy` (MFE/centroid/MEA) into sugarcode-ai e00e90a |
| profile_hmm | HMMER3 (pyhmmer) + Pfam PF00042.29 | 25 globins (hmmalign MSA, 5-fold CV) vs 25 non-globins | AUC global forward 0.952 (benchmarks/profile_hmm_cv.json); after adding forward_local 0.979 (benchmarks/profile_hmm_cv_local.json); HMMER/Pfam 1.000 (benchmarks/results.json) | gap narrowed by sugarcode-ai 1d9091a; residual gap: no length-corrected null or E-values |

## Finding: two public CRISPRscan implementations disagree on one coefficient
Of the 91 CRISPRscan features, 90 are identical between crisprVerse/crisprScore (`inst/crisprscan/crisprscan_coefficients.csv`, devel and master, checked 2026-09-24) and CRISPOR (`crisporEffScores.py`, `paramsCRISPRscan`). The 91st (weight -0.0973770966031, the largest negative) sits on dinucleotide AA at position 19 in crisprScore but at position 18 in CRISPOR. sugarcode follows crisprScore. On 300 held-out BRCA1 contexts, this changes 25 scores (8.3%) by up to 10 points. The CRISPOR doctest (score 77) does not exercise that feature, so both implementations pass their own tests. Which position is correct has to be settled against the original Moreno-Mateos 2015 supplementary model; that is open.

## codon_opt / bio.codon CAI (2026-09-24)
Data: NCBI GCF_000005845.2 CDS (E. coli MG1655) x PaxDb 511145 integrated protein abundance, n=3489 genes.
- CAI implementation matches Biopython CodonAdaptationIndex (Pearson 0.9997 on 500 genes, same reference). VERIFIED.
- Reference table matters: CAI with vendored genome-wide table, Spearman vs log abundance 0.496 (non-ribosomal genes); with ribosomal-protein reference 0.580; top-5%-abundance reference under 5-fold CV 0.574.
- Fix shipped: sugarcode-ai ECOLI_HIGHEXPR / HOST_TABLES["ecoli_highexpr"]. Default left unchanged (optimizer choice of top codon is identical for most AAs; not re-benchmarked).
- Remaining gap: TASEP, metabolic_load, attention, evolutionary_robustness have no real-data validation - verdict THIN (heuristic, no scientific claim).
Numbers: benchmarks/sweep_codon_cai.json.

## pgx_guidelines vs CPIC (2026-09-24)
Reference: CPIC API diplotype tables (CYP2C19 666 diplotypes, CYP2D6 16,836).
- Before: covered 45/666 (C19) and 253/16,836 (D6); D6 agreement 217/253. 36 mismatches came from stale activity values (*9 and *41 at 0.5; CPIC now 0.25). BUG CONFIRMED.
- After (sugarcode-ai fix): 666/666 and 16,836/16,836 agree. C19 is a verbatim lookup (agreement by construction). D6 is computed from allele activity sums and thresholds, so it is a real check.
- Also fixed: "Likely poor/intermediate metabolizer" now triggers the clopidogrel alternative (was falling through to "no recommendation").
Numbers: benchmarks/sweep_pgx_cpic.json. Verdict: VERIFIED for phenotype translation; recommendation text covers only 2 gene-drug pairs (THIN scope).

## neohunter hla_binding vs IEDB (2026-09-24)
Data: IEDB MHC-I binding 2013 (Kim et al. 2014), human 9-mers, binder = IC50<500 nM.
- Anchor heuristic AUROC: A*02:01 0.830, A*03:01 0.801, A*24:02 0.815, B*07:02 0.848, B*44:03 0.827.
- One-hot logistic PSSM, 5-fold CV: 0.954, 0.942, 0.897, 0.962, 0.916.
- Fix shipped: sugarcode-ai neohunter.hla_binding_iedb (trained on all IEDB 2013 rows for the 5 alleles). Not yet compared with NetMHCpan/MHCflurry (next).
Numbers: benchmarks/sweep_neohunter_iedb.json. Verdict: heuristic THIN; PSSM VERIFIED (CV).
- MHCflurry 2.x pan-allele (in-sample reference, trained on data including IEDB 2013): 0.960, 0.948, 0.905, 0.953, 0.938. Our CV PSSM is within 0.02 on all five alleles and above MHCflurry's in-sample AUROC on B*07:02 (0.962 vs 0.953). This is not a SOTA claim: the training regimes differ.

## docking_studio vs AutoDock Vina on LP-PDBBind (2026-09-24)
73/80 complexes ran (7 prep errors listed in the JSON). Redocking with real Vina 1.2.7: top-1 RMSD<2A 52.1%, best-of-9 83.6%.
Scoring power (Pearson with pK, bootstrap 95% CI): Vina crystal 0.274 [0.04,0.47]; Vina top pose 0.285; sugarcode Vina-form scorer on crystal pose 0.298 [0.03,0.48]; heavy-atom count 0.245 [-0.04,0.49].
Verdict: sugarcode's pair scorer is on par with Vina for ranking crystal poses (the CIs overlap heavily), but its dock_vina_grid (straight-line ligand) is NOT a docking method. Fix shipped: real Vina in docking_studio/vina_real.py (sugarcode-ai e0564e4). No scoring-power SOTA claim: all scorers are weak on this split.
Numbers: benchmarks/sweep_docking.json.

## Cross-species CAI reference test (2026-09-24), a named, falsifiable claim
Claim: CAI computed against a ribosomal-protein reference set tracks measured protein abundance better than CAI against the organism's genome-wide codon table.
Test: all 405 PaxDb taxa were scanned; 40 were usable (at least 300 genes matched to NCBI RefSeq CDS); 29 have at least 15 ribosomal-protein genes.
Result: the ribosomal reference is better in 27/29 taxa (median Spearman gain +0.082, Wilcoxon p=3.3e-7, sign test p=1.6e-6). The gains are largest in fast-growing microbes (Klebsiella +0.32, S. cerevisiae +0.18, P. aeruginosa +0.17) and near zero in mammals. The exceptions are Plasmodium falciparum and rat.
Falsifiable prediction: for any new PaxDb organism with a strongly skewed codon usage, the ribosomal reference beats the genome table. For mammals and AT-rich parasites the difference is within 0.02.
Novelty caveat: Sharp & Li (1987) proposed highly expressed reference sets. This is a systematic 29-taxon abundance-based test of that choice, not a new theory.
Numbers: benchmarks/sweep_codon_cai_multispecies.json (per-taxon rows with assembly and PaxDb file IDs).
