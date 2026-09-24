# Manifests

tools.tsv: 42 tools actually imported or called by committed scripts (versions read from the environment 2026-09-24).
datasets.tsv: 24 study-level datasets (transparency count).
accessions.tsv: 327 accession-level records individually fetched and used (gate count under the uniform rule: each identifier-backed record fetched and used = 1). Study-level rows count once each; each RCSB PDB entry downloaded individually for docking counts once; each species PaxDb dataset and NCBI assembly fetched for the cross-species CAI benchmark counts once.
Excluded (bulk-fetched, not individually): 1000 ChEMBL compound IDs, 50 UniProt accessions, 500 UniProt reviewed human proteins (one search, sweep_proteinprops_restriction), 20 empty LOVD lookups.
Gate per spec: 40 tools / 120 datasets per project. Current: tools 53 listed incl. 3 infra = 50/40, datasets 327/120.
Tool-gate convention (program-wide ruling 2026-09-25): infrastructure (git/GitHub, Drive, pytest, generic file-format/transport plumbing, plotting) does not count; tools.tsv column "gate" marks each entry.
