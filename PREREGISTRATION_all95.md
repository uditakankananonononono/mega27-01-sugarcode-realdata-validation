# Preregistration: real-data validation of all 95 SugarCode AI modules

Date: 2026-09-26. Owner repo: mega27-01-sugarcode-realdata-validation.
Companion platform repo: sugarcode-ai (95 registered modules: 77 spec + 18 beyond-spec;
source of truth for the inventory: src/omega/registry.py).

## Method gates (apply to every family)
1. Real public data only. Frozen fixtures with provenance (source + retrieval date) in data/README.md. No stubs, no simulated labels.
2. Every family has a named public benchmark/baseline the modules must BEAT (not tie). If a module loses, we improve it (sugarcode-ai side) and re-run until it beats, or honestly record the negative and pivot via ChatGPT redirection. Negatives are never terminal and are preserved in the paper.
3. Every family must produce at least one NEW DISCOVERY: a novel, verifiable finding preregistered below as a hypothesis target, so discoveries cannot be retrofitted.
4. Minimum 10 ChatGPT judge rounds per family on weaknesses/improvements, logged verbatim in JUDGE_ROUNDS/<family>.md (cloud browser, user's ChatGPT account). Negatives/failures trigger ChatGPT redirection rounds per the 2026-09-26 rules.
5. ISEF-winner archetype lens on writeup: verified facts only, honest limitations, clear novelty statement per family.
6. Deliverable: single 50+ page text-body paper (paper/paper.tex), per-family chapters, honest negatives preserved, free-first external tools (NCBI eutils, RCSB, UniProt, EBI, ChEMBL, BiGG, Kazusa, Pfam, CPIC, IEDB, DepMap, JUMP-CP, etc.).

## Status at preregistration (from mega27-01 STATUS.md / benchmarks/results.json run 2)
Already beaten/validated (7 modules): crispr_opt (CFD 600/600 exact; RS1 Spearman 1.000),
openclinvar (acc 0.997, CI95 0.991-1.000), deepsplice (donor AUC 0.964 / acceptor 0.911 vs hard
decoys; beats 0.761 consensus baseline), virtual_cell (FBA parity 0.873922), codon_opt (CAI exact),
acmg_bayesian (Tavtigian math exact). Open gaps: mutdock (no honest SKEMPI adapter), profile_hmm
(no globin MSA ships). 88 modules unvalidated externally.

## Family plan (16 families, 95 modules, each module in exactly one family)

F1 CRISPR on/off-target scoring (6): cfd_offtarget, mit_offtarget, crisprscan_score, crisprater, crispr_opt, crispr_muse.
  Data: CRISPOR oracle pairs (extend 600 -> 2000), GUIDE-seq / CHANGE-seq public off-target sets, Doench 2016, Moreno-Mateos 2015.
  Beat: CRISPOR percentile rank correlation; published CFD/MIT score concordance. Discovery: systematic discordance map between CFD and MIT scoring on endogenous sites (unreported at scale).

F2 Prime/epi editing & delivery design (4): prime_design, epi_edit, crispr_cargo, vector_opt.
  Data: published pegRNA efficiency datasets (prime editing screens), NCBI epigenome tracks, AAV serotype tropism literature tables.
  Beat: baseline pegRNA heuristics / pegIT-style rules. Discovery: PBS/RT-template feature interaction not in published rules.

F3 Variant interpretation & clinical genomics (7): openclinvar, acmg_bayesian, pgx_guidelines, rarenet_ai, str_scope, gene_analysis, infinite_diagnosis.
  Data: frozen ClinVar 500 (extend to 5000, held-out), CPIC/PharmGKB PGx guidelines, ClinGen, STR expansion loci literature.
  Beat: ClinVar concordance vs baseline pathogenicity predictors on held-out set; CPIC guideline rule coverage. Discovery: reproducible misclassification cluster in a specific variant class.

F4 Splicing & RNA biology (5): deepsplice, rna_decoder, rna_nussinov, riboswitch, dark_genome.
  Data: held-out RefSeqGene junctions (extend beyond NG_008245.1), SpliceAI benchmark junctions, RNAstructure/ViennaRNA as folding oracle, published m6A datasets.
  Beat: PWM vs SpliceAI-derived scores on disjoint junctions; Nussinov vs brute-force optima. Discovery: splice-region variant class where PWM and SpliceAI systematically disagree.

F5 Protein structure & mutation effects (7): alpha_fold_ui, structural_biophysics, mutdock, evofold_4d, protein_painter, stability_ai, profile_hmm.
  Data: AlphaFold DB (pLDDT parity), SKEMPI 2.0 (build honest adapter: use complex structures from PDB + ligand-free interface rows, or pivot to ProTherm/published ddG sets), Pfam families with shipped MSAs.
  Beat: naive ddG baselines (already computed: RMSE 1.629 mean predictor); Pfam HMMER reference. Discovery: mutation classes where simple biophysics beats ML-ddG tools (publishable negative for the field).

F6 Docking, DTI & cheminformatics (8): docking_studio, neodti_engine, dti_bench, qsar_bench, molecule_eval, chem_descriptors, chem_similarity, chemgpt_engine.
  Data: RDKit oracle (chem_descriptors/chem_similarity already exact), PDBbind refined set, DAVIS/KIBA DTI benchmarks, ChEMBL bioactivity, DUD-E.
  Beat: Vina-style scoring baseline on PDBbind; published DTI baseline AUCs on DAVIS/KIBA. Discovery: scaffold family where descriptor-only QSAR matches deep models.

F7 Codon/gene design & central dogma (4): codon_opt, gene_tx_opt, gene_explorer, promoter_lib.
  Data: Kazusa usage tables (done), GenBank annotations, published promoter strength libraries.
  Beat: CAI-only optimization baseline on expression proxy datasets. Discovery: codon-pair effects missed by single-codon optimization.

F8 Metabolic modeling & cell-free (3): virtual_cell, metabodesigner, cell_free_opt.
  Data: BiGG models (e_coli_core done; extend to iJO1366), published CFPS yield measurements.
  Beat: published growth-rate predictions; CFPS yield vs literature baseline. Discovery: gap-filled reaction set improving out-of-sample growth prediction.

F9 Synbio circuits & biosensors (5): synbio_studio, syn_bio_studio, living_computer, bio_switch, syn_stab_ai.
  Data: published repressilator/toggle-switch parameters, iGEM part registry measurements.
  Beat: ODE reference solutions vs published circuit dynamics. Discovery: parameter region predicting instability not flagged in original publications.

F10 Minimal genome & evolution (3): synthetic_life, synlife_evo, bioplayground.
  Data: JCVI-syn3.0 essential gene set, DEG database of essential genes.
  Beat: essentiality prediction vs DEG labels baseline. Discovery: cross-organism essentiality patterns (quasi-essential genes) quantified.

F11 Cell fate, organoids & tissue (8): cell_twin, fate_predictor, cellfatenet, organoid_ai, syndroid, tissue_eng, bio_material, bioprint_pro.
  Data: Mogrify/CellNet published TF sets, public scRNA-seq atlases (Human Cell Atlas), biomaterial property literature tables.
  Beat: Mogrify TF-recall on held-out cell-type pairs. Discovery: TF combinations predicted for conversions absent from current databases.

F12 Oncology & immunotherapy (6): neohunter, oncocircuit, car_t_designer, pdx_insight, liquid_biopsy, organoid_screen.
  Data: TESLA neoantigen consortium dataset, IEDB + MHCflurry/NetMHCpan baselines, CCLE/DepMap.
  Beat: MHCflurry binding-rank baseline on TESLA held-out. Discovery: immunogenicity features outside binding affinity that improve ranking.

F13 Microbiome, phage & living therapeutics (8): microaiverse, microbiome_exp, micro_tx, microbiome_rx, phage_tx, phage_designer, phageforge, living_tx.
  Data: HMP/MGnify public cohorts, phage-host databases, PHASTER prophage calls.
  Beat: taxonomy-classifier baseline on held-out HMP samples; phage host-prediction baseline. Discovery: phage-host pairs predicted but absent from databases (testable).

F14 Bioimaging & neuroimaging (4): bioimage_ai, cellpainter, cellpainter_4d, neuroplan_ai.
  Data: BBBC benchmark sets, JUMP-CP public plates, public MRI segmentation sets (e.g., BraTS).
  Beat: CellProfiler-style baseline features on BBBC; baseline segmentation Dice on a public MRI subset. Discovery: morphology-feature correlates of perturbation not annotated in the dataset.

F15 Neuro/literature intelligence (7): neuro_hub, neuro_pipeline, biogpt_lit, genomegpt, evidence_mining, bio_copilot, synbio_wizard.
  Data: BC5CDR and similar gold biomedical NLP corpora, PubMed baselines, regulatory-motif datasets for genomegpt.
  Beat: published entity-extraction F1 baselines on BC5CDR; motif-discovery baseline (MEME-style) on public ChIP-seq peaks. Discovery: literature-contradiction pairs surfaced by the temporal KG, verified against retracted/corrected papers.

F16 Platform & infrastructure (10): neuro_hub_dashboard, dna_to_code, biosimvr, biofactory_1_a, robotic_flow, omega_stats, enterprise_bio, ecosystem, nexus_support, report_studio.
  Honest scope: no scientific benchmark exists for UI/orchestration modules. Validation = integration/functional test matrix + omega_stats numerical checks vs scipy/R references. Beat: scipy stats parity within tolerance; report_studio round-trip integrity. Discovery: N/A by design - recorded as honest scoping, not failure.

## Sequencing
Wave 1 (harness exists, extend): F1, F3, F4, F5 (mutdock adapter), F6, F7, F8.
Wave 2 (new fixtures): F2, F9, F10, F12, F13, F14.
Wave 3: F11, F15, F16 + paper assembly + judge-round closeout.
Each family: preregister -> freeze fixtures -> run module vs benchmark -> 10 judge rounds -> improve/pivot -> discovery verification -> chapter.
