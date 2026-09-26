# Data provenance

Every dataset used by this project is listed here with source URL, accession, retrieval date, and license. No fabricated data.
- clinvar_sample_500.tsv: NCBI ClinVar variant_summary.txt.gz (ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/, retrieved 2026-09-24), GRCh38 SNVs, criteria-provided review status, reservoir sample 250 pathogenic-class + 250 benign-class from 1,517,551 eligible.
- skempi_v2.csv: SKEMPI 2.0 (Jankauskaite et al. 2019, PMID 30598530), life.bsc.es/pid/skempi2, retrieved 2026-09-24.
- e_coli_core.json: BiGG Models e_coli_core (bigg.ucsd.edu, Orth et al. 2011 PMID 21988831), retrieved 2026-09-24.
- kazusa_83333.html: Kazusa CUTG E. coli K12 codon usage (gbbct, 14 CDS / 5122 codons sample), retrieved 2026-09-24. Caveat: small-sample table; used only for triplet frequency sanity checks.
- crispr_fresh_pairs.tsv / crispr_cfd_oracle.tsv: 200 guides mined from BRCA1 RefSeqGene NG_005905.2:60000-70000 (NCBI efetch, 2026-09-24; disjoint from the sugarcode deepsplice training accessions), 600 guide/offtarget pairs (1/2/4 mismatches), oracle CFD values computed from the CRISPOR reference distribution (github.com/maximilianh/crisporWebsite, Concordet & Haeussler 2018 PMID 30341956).
- globins25.fa: 25 reviewed globin-family UniProtKB sequences (rest.uniprot.org, 2026-09-24) for profile-HMM validation; alignment via EBI Clustal Omega REST (job clustalo-R20260924-152926-0644-7040537-p1m).
- skempi_subset_150.tsv: frozen 150-mutation random subset (seed 7) of the 4,956 usable SKEMPI singles; mutation format = <chain><wt><pos><mut>; ddg_kcal = 0.596*ln(Kd_mut/Kd_wt) at 300K.
- doench_rs1_oracle.tsv: 150 30-mer guide contexts (4bp+20nt+NGG+3bp) from the same held-out BRCA1 segment, with Rule Set 1 percentile scores computed by the CRISPOR reference implementation (crisporEffScores.py, retrieved 2026-09-24).

## pfam_PF00042_globin.hmm
Pfam PF00042.29 (Globin) HMMER3 model, fetched 2026-09-24 from
https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00042?annotation=hmm . Real-world reference
for profile_hmm validation (scored with pyhmmer 0.12.3). Replaces the queued EBI Clustal route.
- crispr_tp53_pairs.tsv / tp53_NG_017013.2.fa: 350 NGG guides mined from TP53 RefSeqGene NG_017013.2 (NCBI efetch, 2026-09-26; disjoint from the BRCA1 fixtures and deepsplice training accessions), 1,400 guide/off-target pairs (2x1mm, 1x2mm, 1x4mm per guide, seed 20260926) with oracle CFD (CRISPOR CFD_Scoring, 6/6 doctests exact) and MIT Hsu scores. Generator: scripts/gen_f1_tp53_cfd_fixtures.py.
- clinvar_extension_4998.tsv: 4,998-variant extension sample (2,498 pathogenic-class + 2,500 benign-class after removing 2 overlapping the frozen 500), GRCh38 SNVs, criteria-provided review status, reservoir seed 20260926, from variant_summary.txt.gz retrieved 2026-09-26 (188,675 path / 1,328,876 ben eligible). Generator: scripts/gen_f3_clinvar_extension.py.
- crispor_offtargets_7studies.tsv: 764 experimentally measured on/off-target rows (read-fraction scores) from 7 studies (Cho2014, Frock2015, Hsu2013, Kim2015, Ran2015, Tsai2015, Wang2015), collected by the CRISPOR paper (github.com/maximilianh/crisporPaper, retrieved 2026-09-26). Experimental-validation arm for the F1 CFD-vs-MIT discordance discovery study (judge round 1 fold-back).
