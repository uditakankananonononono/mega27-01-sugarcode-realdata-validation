# Module sweep verdicts (item 1)

| module | external reference | held-out data | result | verdict |
|---|---|---|---|---|
| chem_descriptors | RDKit 2026.03.6 | 993 ChEMBL phase-4 small molecules (ChEMBL API, 2026-09-24) | mol_wt, exact_mol_wt, hbd, rotatable_bonds, tpsa, fraction_csp3, heavy atoms, formula: 993/993 exact. hba: 993/993 vs the RDKit Lipinski HAcceptor SMARTS the module implements; 867/993 vs `Lipinski.NumHAcceptors` in RDKit 2026.03, which now calls the C++ CalcNumHBA and no longer equals that SMARTS | pass (version divergence documented) |
| chem_similarity | RDKit Morgan r=2, 2048 bits | 500 random pairs of the same 993 drugs | Tanimoto 500/500 exact (max abs diff 0.0) | pass |
| mit_offtarget | CRISPOR `calcHitScore` (verbatim, maximilianh/crisporWebsite) | 600 held-out BRCA1 guide/off-target pairs | 600/600 exact (max diff 1.1e-16) | pass |
| crisprscan_score | CRISPOR `calcCrisprScanScores` (verbatim) | 300 held-out BRCA1 35-nt contexts | 275/300 integer-equal, max diff 10 points | reference discrepancy found (see below) |
| rna_nussinov | ViennaRNA 2.7.2 (Turner 2004) | Rfam 15 seeds RF00005 tRNA, RF00001 5S, RF00010 RNase P; 25 seqs each (seed 7), consensus SS projected | base-pair F1 Nussinov 0.328 / 0.186 / 0.189 vs ViennaRNA MFE 0.714 / 0.531 / 0.545 | exact algorithm but biologically weak - closed by integrating `fold_energy` (MFE/centroid/MEA) into sugarcode-ai e00e90a |

## Finding: two public CRISPRscan implementations disagree on one coefficient
Of the 91 CRISPRscan features, 90 are identical between crisprVerse/crisprScore (`inst/crisprscan/crisprscan_coefficients.csv`, devel and master, checked 2026-09-24) and CRISPOR (`crisporEffScores.py`, `paramsCRISPRscan`). The 91st (weight -0.0973770966031, the largest negative) sits on dinucleotide AA at position 19 in crisprScore but at position 18 in CRISPOR. sugarcode follows crisprScore. On 300 held-out BRCA1 contexts, this changes 25 scores (8.3%) by up to 10 points. The CRISPOR doctest (score 77) does not exercise that feature, so both implementations pass their own tests. Which position is correct has to be settled against the original Moreno-Mateos 2015 supplementary model; that is open.
