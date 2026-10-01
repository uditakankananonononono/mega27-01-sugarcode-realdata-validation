# SugarCode fresh self-audit, October 1, 2026

Product source: https://github.com/uditakankananonononono/sugarcode-ai
Product commit: 1ccbd0305cd9b639b06324fa136654d686483b7a
Live remote HEAD rechecked after both final runs. Tracked working tree clean.
Validation source: https://github.com/uditakankananonononono/mega27-01-sugarcode-realdata-validation
No product changes, commits, or pushes performed.

## Reproduced results
- Two final dependency-complete full runs: 2,319 passed, zero skipped, zero failed each.
- 2,319 JUnit case identities and outcomes match between final runs.
- Run times: 90.92s and 90.45s. Both report two warnings: packaging license deprecation and DESeq2 low residual degrees of freedom on its small fixture. No assertion failure.
- 95 registry entries match 95 package directories, no missing or orphan packages.
- 95 imports succeed; 95 packages or core files have method docstrings.
- 77 spec documents, 18 beyond-spec registry entries.
- All 95 have passed tests attributed via test-file imports or in-package PGx test location. This attribution is non-exclusive and is not individual requirement verification.
- Seven independent behavioral checks reproduce byte-identically: FBA objective 15 versus separately assembled scipy LP; CAI 0.8532194221288482 versus geometric mean; seeded TASEP reproducibility and density bounds; ANM six zero modes and first frequency 0.1329 versus separately assembled Hessian; seeded chemical-Langevin reproducibility and nonnegative states; Hill at Kd 0.51 and broad SSA steady-state check; poly-Ala Chou-Fasman helix fraction 0.75.

## Optional dependency and real-fixture coverage
All initially skipped cases now run. Installed CPU torch, maxentpy, DESeq2, Biopython, primer3, RDKit, Vina/OpenBabel; set isolated-environment PATH for the OpenBabel executable. Downloaded the real BBBC001 fixture from Broad.

## Real-image benchmark negative
Source: https://data.broadinstitute.org/bbbc/BBBC001
Images: https://data.broadinstitute.org/bbbc/BBBC001/BBBC001_v1_images_tif.zip
Counts: https://data.broadinstitute.org/bbbc/BBBC001/BBBC001_v1_counts.txt
The six-image mean relative absolute cell-count deviation is 0.08182388535627162 (8.1824%). It passes the existing suite gate (<11%, described by Broad as inter-human variation), but does NOT beat Broad's listed published method (6.2%). This is benchmark-negative despite suite-positive. It is a small known fixture, not a fresh held-out validation or new discovery. Dataset license: CC BY-NC-SA 3.0; raw images are not included in this evidence patch. Download hashes are in provenance.json.

## Rechecked prior findings
Current STATUS identifies Hill kinetics with synbio_wizard/virtual_cell, not gene_analysis. It identifies living_computer as chemical-Langevin/Euler-Maruyama. It clearly labels the historic test count as historic and says a newer full run was missing. PGx's six tests are now in configured testpaths. Current module/core docstrings exist for all 95. None of these older findings should be reported as a fresh defect.

## Limits
Suite success is not 95 external benchmark wins, clinical validity, a new discovery, or every spec sentence independently verified. Functional fixtures can include synthetic inputs; this report does not reclassify them as real biological data. Stochastic checks use fixed seeds and broad expected-value bounds, not a new scientific validation. Some imported-file tests exercise other modules, so per-module attribution cannot be added up or treated as isolated module verification. The validation repo's earlier benchmark counts and paper are not freshly reproduced here. Specs label themselves aspirational; suite success validates shipped behavior rather than proving every original aspiration. Per-sentence spec verification remains open. No deployment occurred.

Reproduce from the pinned product source using the listed environment:
PATH=/tmp/sugarcode-venv/bin:$PATH /tmp/sugarcode-venv/bin/python -m pytest -q -rs --junitxml=results.xml
Run independent probes with SUGARCODE_SRC pointing to the pinned source src/ directory and PROBE_OUTPUT pointing to the output JSON file. Only seven selected algorithm checks are covered.
