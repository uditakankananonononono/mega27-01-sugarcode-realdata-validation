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

## Splice: MaxEntScan baseline and fix (2026-09-24)
Same gene-grouped folds, n=9235. MaxEntScan (Yeo & Burge 2004) delta has AUROC 0.9635, beating our CNN at 0.9355 (the CNN is worse by 0.023 to 0.033, 95% CI). NEGATIVE: the CNN lost to a 2004 method, mainly at acceptors (0.879 vs 0.965).
Pivot: add MaxEntScan ref/alt/delta as features. The logistic model with MaxEnt features reaches 0.9674 and beats MaxEnt alone by 0.0027 to 0.0053 (95% CI). CNN with MaxEnt features reaches 0.9622, which is no better than MaxEnt (-0.0042 to +0.0017).
Numbers: discovery/splice_region_vus/cv_results_me.json, maxent_cmp.json.

## Codon: tAI and ENC comparison (2026-09-24)
Spearman with log abundance (non-ribosomal genes):
- E. coli: CAI-genome 0.496, CAI-ribo 0.581, tAI 0.528, -ENC 0.384.
- B. subtilis: 0.330, 0.416, 0.415, 0.174.
- S. cerevisiae: 0.410, 0.592, 0.660, 0.447.
tAI beats CAI-ribo in yeast but not in E. coli. The tools used were codon-bias and GtRNAdb gene copy numbers. Numbers: benchmarks/sweep_codon_tai_enc.json.
Cross-species BH-FDR (statsmodels): ribosomal-reference CAI is significantly correlated with abundance (q<0.05) in 27/29 taxa.

## bio.structures SASA vs FreeSASA (2026-09-24)
- 28 RCSB receptors from the docking set (<=6000 atoms, H removed). 18 had identical atom sets: total SASA median relative difference +1.2% (max 2.6%), per-atom Pearson median 0.980 (benchmarks/sweep_sasa_freesasa.json). VERIFIED (radii differ from FreeSASA ProtOr, hence the small bias).
- 10/28 disagreed on atom count: sugarcode kept every alternate location, so overlapping conformers were all counted. BUG CONFIRMED, fixed in sugarcode-ai (keep first altloc). Recheck on 3 of them: atom counts now match, per-atom Pearson 0.975-0.983 (benchmarks/sasa_altloc_recheck.json).

## bio.primer vs primer3-py (2026-09-24)
- 5,404 20-mers tiled from 5 RefSeq mRNAs. tm_nn equals Biopython Tm_NN with DNA_NN4 exactly (max diff 0.0). Versus primer3 (SantaLucia 1998 table) mean -0.17 C, max 0.98 C, Pearson 0.99995: table choice, not a bug (benchmarks/sweep_primer_primer3.json).
- hairpin_max_stem (stem-length heuristic) vs primer3 hairpin dG: Spearman 0.04. NEGATIVE. Fix: sugarcode-ai bio.primer.hairpin_dg calls primer3 when installed.

## bio.phylo NJ/UPGMA vs DendroPy (2026-09-24)
- Same distance matrices (K80; p-distance for RF00005 where K80 saturates to infinity), 30 random sequences from each of 3 Rfam seed alignments.
- NJ: identical topology (RF 0) and patristic distances (max diff 1e-6) in 3/3 families. UPGMA: identical in RF00001 and RF00010; RF00005 has RF 2 (patristic Pearson 0.998), traced to 2 merge steps with tied minimum distance (sweep/phylo_upgma_ties.py), a tie-breaking difference, not a bug. VERIFIED (benchmarks/sweep_phylo_dendropy.json).
- p-distance differs from Biopython identity by design (pairwise gap deletion vs gaps counted): Pearson 0.79-0.98.

## bio.orf vs orfipy and NCBI CDS annotation (2026-09-24)
- 5 RefSeq mRNAs (TP53, BRCA1, BRCA2, BRAF, EGFR; GenBank records in data/refseq/). The longest + strand ATG ORF equals the annotated CDS in 5/5, and its translation equals the NCBI /translation in 5/5. VERIFIED.
- Six-frame ORF sets (min 30 aa) are identical to orfipy in 5/5 when orfipy minlen is 89-90; at 91-93 orfipy drops the 30-aa ORFs. Length-threshold convention, not a bug (benchmarks/sweep_orf_orfipy.json).

## bio.vcf vs pysam/htslib on ClinVar GRCh38 (2026-09-24)
- 46,603 ClinVar records across 5 gene regions, fetched by remote tabix. CHROM/POS/ID/REF/ALT: 0 mismatches. INFO: 0 mismatches after two normalisations (pysam returns percent-escapes raw where sugarcode decodes them; pysam stores floats as float32), except 2 CLNVI values where ClinVar uses %2B, which is not a VCF 4.3 escape and sugarcode correctly leaves as-is (benchmarks/sweep_vcf_pysam.json). VERIFIED.
- Small spec gap fixed: %3A (':') was not decoded; sugarcode-ai now decodes it.
- variant_type vs ClinVar CLNVC: agrees on 45,685 records. The other 302 are equal-length multi-base changes that ClinVar labels Indel and sugarcode labels mnp (VCF convention). Naming difference.

## bio.de / bio.rnaseq vs PyDESeq2 on GEO GSE52778 (2026-09-24)
- Airway smooth muscle, dexamethasone vs untreated, 4 cell lines; 16,802 genes after filtering.
- Median-of-ratios size factors equal PyDESeq2's (max relative diff 1.8e-9). VERIFIED.
- DE calls: PyDESeq2 ~cell + dex finds 4,451 genes at padj<0.05 (3,203 unpaired). sugarcode Welch finds 351 (331 shared with PyDESeq2) and Wilcoxon finds 0 (with 4 vs 4 the smallest possible p cannot survive BH). All 7 canonical dex-response genes (FKBP5, TSC22D3, PER1, DUSP1, KLF15, ZBTB16, CRISPLD2; PyDESeq2 padj 1e-19 to 1e-140) are missed by both simple tests. Log2FC agrees (Pearson 0.993). NEGATIVE (benchmarks/sweep_de_pydeseq2.json); the module's own scope note already warned about this.
- Fix: sugarcode-ai bio.de.de_analysis_deseq2 (PyDESeq2 NB GLM with an optional blocking factor for paired designs).

## bio.align NW/SW vs parasail (2026-09-24)
- All 300 pairs of 25 Pfam PF00042 globins, BLOSUM62, gap open 11 / extend 1: global and local scores equal parasail nw_scan/sw_scan in 300/300 each (max diff 0.0). Same gap convention as parasail (first gap residue costs the open penalty). VERIFIED (benchmarks/sweep_align_parasail.json).

## Dex-response enrichment with the new DESeq2 path (2026-09-24)
- sugarcode de_analysis_deseq2 (~block + group) reproduces the direct PyDESeq2 run exactly: 4,451 genes at padj<0.05 (536 up with log2FC>1, 529 down with log2FC<-1). Table in benchmarks/deseq2_dex_paired.tsv.gz.
- g:Profiler (custom background of 16,802 tested genes, g:SCS): top terms are generic (multicellular organismal process, signalling, cell adhesion; KEGG cytoskeleton in muscle cells). "Cellular response to glucocorticoid stimulus" is over-represented about 7-fold (9 of 43 genes; approximate raw hypergeometric p 3.4e-6) but does NOT pass g:SCS correction. Reported as is (benchmarks/sweep_dex_enrichment.json).
- STRING v12 PPI enrichment, top 300 dex-induced genes (284 mapped): 402 edges vs 184 expected (STRING reports p = 0.0, below its numeric floor). Caveat: STRING's expectation uses a genome-wide background, not the 16,802 expressed genes, which inflates the enrichment (benchmarks/sweep_dex_string_ppi.json).
- Reactome AnalysisService on the 536 dex-induced genes: no pathway reaches FDR<0.05; the best are anchoring fibril formation, FOXO-mediated transcription, laminin interactions and metallothioneins (all FDR 0.11). Negative at the pathway level, consistent with g:Profiler (benchmarks/sweep_dex_reactome.json).
- STRING network + networkx (dex top 300): 201 connected nodes, 402 edges, largest component 183. Hubs: FOXO1 (21), CEBPD, FOXO3. Eight of ten canonical GR targets (fixed list, set before the run) are in the network with mean degree 6.25 vs 4.0 overall; degree-permutation p=0.049 (10k perms). Marginal: 8 genes, one pre-set list, no multiple-test correction applies but treat as weak (benchmarks/sweep_dex_string_network.json).
- Preranked GSEA (gseapy, MSigDB Hallmark, 50 sets, 1000 perms) on the full PyDESeq2 dex ranking (16,525 genes): 1 set at FDR<0.05, Adipogenesis (NES 1.70, FDR 0.038); 5 sets at FDR<0.25 (UV response up, fatty acid metabolism, TNF-alpha/NF-kB, ROS). The rank-based test picks up a glucocorticoid-linked programme (adipogenesis) that the threshold-based ORA (g:Profiler, Reactome) missed; one set, modest FDR (benchmarks/sweep_dex_gsea_hallmark.json).
- goatools local GO:BP ORA (Fisher + BH, propagated, GOA human GAF, expressed-gene background of 16,525; 527 dex-induced genes mapped): 38 BP terms at FDR<0.05, led by cellular response to hormone stimulus and response to hormone (both FDR 3.4e-5). Cellular response to glucocorticoid stimulus: 7/34 genes, FDR 0.026 under BH, versus not significant under g:Profiler's stricter g:SCS. The glucocorticoid call depends on the correction method; report both. 9 of g:Profiler's 22 top BP terms are also significant here (benchmarks/sweep_dex_goatools.json).

## bio.proteinprops and bio.restriction vs Biopython (2026-09-25)
- proteinprops on 500 UniProt reviewed human proteins: MW, GRAVY, aromaticity, instability index and both extinction coefficients match ProtParam exactly (max rel diff 2e-16); pI within 6e-5. VERIFIED.
- restriction on pBR322 (J01749.1, circular) and lambda (NC_001416.1, linear), 610 enzymes: BUG FOUND. The reverse-strand pattern for degenerate sites used a plain DNA revcomp, so IUPAC codes were not complemented (AccI GTMKAC searched as GTKMAC on the other strand), giving extra false sites (AccI on pBR322: 5 vs the known 2). Also linear DNA reported Type IIS cuts outside the molecule. Before: 533/610 and 529/610 enzymes agree; after sugarcode-ai fe2185d: 608/610 and 608/610. The 2 left (LpnPI, MspJI) are modification-dependent enzymes that do not cut unmethylated DNA; open.
Numbers: benchmarks/sweep_proteinprops_restriction.json.

## bio.stockholm vs Biopython AlignIO (2026-09-25)
15 Rfam seed alignments (12 fetched individually today + RF00001/5/10), 4,471 sequences: ids and order, aligned sequences and SS_cons identical to Biopython in 15/15; write->parse round trip lossless in 15/15; pairwise identity equals an independent count. VERIFIED. Numbers: benchmarks/sweep_stockholm_biopython.json.

## bio.gff vs gffutils (2026-09-25)
NCBI RefSeq E. coli K-12 MG1655 annotation (GCF_000005845.2_ASM584v2_genomic.gff.gz), 9,523 records: seqid/type/coords/strand/phase/attributes identical to gffutils in 9,523/9,523; children of 300 genes identical in 300/300; write->parse lossless. The first comparison showed 2/300 child mismatches; both were split CDS features (b0240, b4587) with a repeated ID that gffutils renames X_1, a comparison artifact, not a sugarcode error. VERIFIED. Numbers: benchmarks/sweep_gff_gffutils.json.

## bio.newick vs Bio.Phylo and DendroPy (2026-09-25)
6 Open Tree of Life synthetic subtrees (Insecta, Laurasiatheria, Carnivora, Metazoa, Aves, Primates; 3,655 leaves incl. quoted labels and empty leaves from height truncation): leaf sets equal Biopython and DendroPy, internal-node counts and internal labels equal Biopython, write->parse lossless, 6/6. VERIFIED (topology/labels only; these trees carry no branch lengths). Numbers: benchmarks/sweep_newick_biophylo.json.

## bio.sam vs pysam/htslib (2026-09-25)
10,943 real reads (1000 Genomes NA12878 chr20:1.0-1.2 Mb, low-coverage bwa BAM): qname, flag, rname, pos, mapq, CIGAR ops, mate fields, seq, qual and typed tags equal pysam for 10,943/10,943. BUG FOUND in alignment_end: 38 unmapped mates (flag 0x4) that bwa places at the mate's position with a CIGAR got an end coordinate; htslib gives none. Fixed in sugarcode-ai 87cf78e: 10,905 -> 10,943/10,943. Numbers: benchmarks/sweep_sam_pysam.json.

## bio.fastq vs Biopython SeqIO (2026-09-25)
First 20,000 reads of ENA SRR622461_1 (NA12878): ids, sequences and per-base Phred equal Biopython for 20,000/20,000; mean Phred identical; offset detected as 33; round trip lossless. GC 0.4252 equals Biopython gc_fraction(ambiguous='remove'), i.e. N excluded from the denominator (0.4220 with 'ignore'); convention documented, not a bug. VERIFIED. Numbers: benchmarks/sweep_fastq_biopython.json.

## bio.pdb (PDB + mmCIF) vs gemmi (2026-09-25)
30 RCSB entries (first 30 docking IDs), both formats, 94,121 atom sites incl. all altlocs/models. PDB format: 30/30 identical to gemmi on every field. BUG FOUND in mmCIF: 0/30 real RCSB files parsed (POSIX shlex mis-split primes such as O5'; loops ended early on unquoted item-name values like _database_2.pdbx_DOI). Fixed in sugarcode-ai fc1fd1e with a CIF 1.1 tokenizer and row-boundary loop ends: 30/30 identical to gemmi. Numbers: benchmarks/sweep_pdb_gemmi.json.

## Splice: phyloP conservation adds to logit+MaxEnt (2026-09-25)
Stratified subsample of 1,501 of the 9,235 modelled splice-region SNVs (379 pathogenic, 1,054 genes), phyloP100way at each position from the UCSC REST API (0 missing). Gene-grouped 5-fold CV x 20 seeds: logit+MaxEnt 0.9630 -> 0.9678 with phyloP; paired bootstrap CI for the gain [+0.0008, +0.0091], excludes 0. The gain is at donors (0.9557 -> 0.9634); acceptors barely move (0.9702 -> 0.9723). phyloP alone: 0.838 (donor 0.911, acceptor 0.697). Caveat: ClinVar submitters may use conservation as supporting evidence (ACMG PP3/BP4), so part of the gain can be label circularity; also a subsample, not all 9,235. POSITIVE, modest. Numbers: discovery/splice_region_vus/phylop_lift.json.

## bio.pileup vs pysam pileup (2026-09-25)
Same 10,943 NA12878 chr20 reads, matched filters (unmapped excluded only): 197,437 positions; A/C/G/T/N counts and deletion counts equal pysam at 197,437/197,437. sugarcode also lists deletions as '*' in the counts (samtools mpileup convention), which the comparison maps to del_depth. VERIFIED. Numbers: benchmarks/sweep_pileup_pysam.json.

## bio.motif + bio.pwm vs Biopython Bio.motifs (2026-09-25)
10 JASPAR 2024 profiles (SP1, CTCF, TP53, TBP, NR3C1, ESR1, SPI1, FOS::JUN, NFKB1, TFAP2A) scanned over BRCA1 10 kb and lambda (both strands): log-odds matrices equal to 1.8e-15, min/max scores equal, per-window scores equal to 3.8e-6 (Biopython stores float32), and the 80%-of-range hit sets identical in 20/20 motif x sequence pairs. VERIFIED (bio.motif uses bio.pwm's log_odds/score/min/max, so both are covered). Numbers: benchmarks/sweep_motif_biopython.json.

## bio.bed vs bioframe (2026-09-25)
UCSC hg38 cpgIslandExt chr20 + rmsk chr20:0-5Mb. Parse and 600 overlap queries identical. BUG: merge_intervals did not join book-ended intervals (rmsk 11,364 vs 8,251). Fixed sc-ai 1fc4a9b. benchmarks/sweep_bed_bioframe.json

## bio.genbank vs Biopython (2026-09-25)
25 NCBI nuccore records, 1,867 features. BUG: multi-line qualifiers leaked into locations (/translation 5/302; spans 1,829/1,867); partial/single-base locations unparsed. Fixed sc-ai 2e21f86: 1,867/1,867 spans, 302/302 translations. benchmarks/sweep_genbank_biopython.json

## bio.sequence vs Biopython (2026-09-25)
All real-data checks agree (translate 285/285, ORFs 22/22, motifs 192/192). Latent IUPAC frame-shift bug (code review) fixed sc-ai 143d1c4. benchmarks/sweep_sequence_biopython.json

## bio.kmer vs sourmash (2026-09-25)
6 coronavirus genomes: exact k-mer sets/Jaccard identical on 15 pairs. NEGATIVE: minimizer-sketch Jaccard biased (MAE 0.0121 vs FracMinHash 0.0044). benchmarks/sweep_kmer_sourmash.json

## bio.gstats vs scipy/statsmodels (2026-09-25)
1000 Genomes phase 3, 22 rsIDs, 110 tables: HWE exact 110/110 vs enumeration; association/OR/BH match. No bug. benchmarks/sweep_gstats_scipy.json

## bio.uniprot vs HGNC (2026-09-25)
BUG: CDKN2A -> CDKN2A-AS1 (Q9UH64); 59->60/60 after gene_exact fix. BUG: salted-hash cache keys in 4 connectors. Fixed sc-ai 6c11d68. benchmarks/sweep_uniprot_hgnc.json

## bio.entrez vs HGNC (2026-09-25)
BUG: HTT -> SLC6A4, TTN -> TTR (alias hits). 58->60/60. Fixed sc-ai 93d5b12. benchmarks/sweep_entrez_hgnc.json

## bio.reactome vs UniProt2Reactome.txt (2026-09-25)
61/61 identical (isoform rows pooled). No bug; isoform pooling caveat. benchmarks/sweep_reactome_download.json

## bio.chembl activities ordering (2026-09-25)
BUG: "best first" was best of an arbitrary 50-row page (ABL1 1,200 nM vs 0.015 nM). Fixed sc-ai c6db085 (server-side pChEMBL order). benchmarks/sweep_chembl_order.json

## bio.gnomad constraint (2026-09-25)
BUG: v2 LOEUF cutoff 0.35 on r4 data; gnomAD v4 guidance 0.45. pLI agreement 51->57/60. Fixed sc-ai 0a265e7; rate-limit fail-fast ce87ccb. benchmarks/sweep_gnomad_constraint.json

## bio.splice junction_map vs ClinVar (2026-09-25)
159/159 MANE-named canonical-site SNVs (11 genes). No bug; CDH1 not in vendored list. benchmarks/sweep_splice_clinvar.json

## modules.cfd_offtarget vs Doench 2016 pickles (2026-09-25)
600 BRCA1 pairs x 16 PAMs = 9,600 scores, max diff 2.2e-16. No bug. benchmarks/sweep_cfd_doench.json
- str_scope (2026-09-25): BUG FIXED (sugarcode 727344a). find_strs vs pytrf 1.5.0 on 25 GenBank records: 208/211 before, 211/211 after, 0 extra calls; 209 identical coordinates, 2 same span in disease-motif phase (by design). Cause: greedy scan skipped past compound repeats that started inside a flank homopolymer. Also expansion_call: legacy delta rule mislabels HTT 27, 28, 36 and 154 FMR1 counts (all premutations 55-200 called pathogenic) vs GeneReviews NBK1305/NBK1384; optional locus bands added, 0 HTT disagreements. Evidence: benchmarks/sweep_str_pytrf.json.

## modules.acmg_bayesian vs ClinGen VCEP classifications (2026-09-25)
BUG FIXED (sugarcode cdf160b): ClinGen suffix notation (PM2_Supporting, PVS1_Very Strong, BS1_Stand Alone) was rejected; 5,699/7,929 rows lost a code, 4,583/7,929 panel calls reproduced. After: 0 rejected, 7,515/7,929 (94.8%) reproduced; ACMG-2015 rules baseline 7,267/7,929. Known difference kept: lone benign supporting (-1 pt) -> LB vs panels' VUS. benchmarks/sweep_acmg_erepo.json

## modules.openclinvar vs ClinVar records + review-status table (2026-09-25)
BUG FIXED (sugarcode 7067655). HGVS consequence on 682 ClinVar records (12 genes): before 646 right / 32 abstain / 4 wrong (exon-spanning ranges judged by start only), after 651 / 30 / 1 (older-transcript annotation, call kept). clinvar_stars: "criteria provided, multiple submitters" gave 1 star, official 2; now 9/9. benchmarks/sweep_openclinvar.json

## modules.crispr_opt on-target scorers vs Doench 2016 FC+RES (2026-09-25)
No code bug. 5,310 guides, 17 genes. Spearman vs measured: RS2+position 0.7106, RS2 sequence-only 0.6693, RS1 0.3948, edge heuristic 0.1125 (in-sample for RS1/RS2). Open: design_guides ranks with RS1; switching to RS2 needs a held-out screen. benchmarks/sweep_crispr_ontarget_fcres.json

## modules.gene_analysis.power_estimate vs statsmodels NormalIndPower (2026-09-25)
BUG FIXED (sugarcode 4defcdb). Two-entry z table matched 7/41 grid cases (power 0.5 gave 16 vs 8; power 0.99 gave 21 vs 37); exact quantiles now 41/41. Formula check only; rest of gene_analysis not validated (not counted in inventory). benchmarks/sweep_power_statsmodels.json

## modules.qsar_bench descriptors vs RDKit on 993 ChEMBL phase-4 drugs (2026-09-25)
BUG FIXED (sugarcode bf8665b). Before: 110/993 rejected (stereo bonds, salt ions); on 883 parsed MW 843, rotatable 96 (vs NonStrict), fraction_csp3 296. After: 993/993 parse; MW, heavy/hetero atoms, charge, ring rank, rotatable (NonStrict), fraction_csp3 all 993/993; HBD 992/993; HBA heuristic 264/993 (labelled). benchmarks/sweep_qsar_descriptors.json
- qsar_bench morgan_fingerprint (sugarcode 9a33eb4): invariants now total H + ring membership. Tanimoto Spearman vs RDKit Morgan r2/2048 on 5,000 pairs 0.7815 -> 0.9005; nearest neighbour 202 -> 236 /300. benchmarks/sweep_qsar_morgan.json
