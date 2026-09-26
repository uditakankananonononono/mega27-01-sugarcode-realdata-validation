# Sugarcode self-audit (2026-09-27, post-judge-route regime)

User direction (WhatsApp 2026-09-27 00:06:30 IST): "Atlas and sugarcode you conduct
your own audit" - sugarcode validation no longer routes through ChatGPT/DeepSeek
judges. This is the internal self-audit: full test suite + module/spec inventory.

## Test suite (HEAD 5981cfa, copied to /tmp, pytest 9.1.1, PYTHONPATH=src)
- 2,286 tests collected from tests/ (pyproject testpaths).
- **2,279 passed, 1 failed, 8 skipped in 86s.**
- The 1 failure: `tests/test_drop54.py::test_zero_empty_function_bodies` - the
  repo's own meta-gate finds 3 ellipsis protocol stubs in `self_improve`
  (planner.py:refine, gate.py:request, gate.py:decision). These are intentional
  protocol stubs with concrete implementations directly below (NullRefiner,
  ManualApprovalGate); functional behavior unaffected. Honest FAIL stands:
  the repo fails its own style gate.
- The 8 skips are all optional-dependency/environment gates: RNA (ViennaRNA),
  torch, pydeseq2, primer3, vina x2, BBBC001 cache, setuptools<61 wheel.

## Module coverage (95 modules in src/sugarcode/modules)
- 94/95 have import-linked tests inside the configured suite (per-module
  pass/fail attribution in inventory.json, parsed from test-file imports).
- pgx_guidelines: no suite-linked tests, but ships an in-module
  test_core.py that pytest's configured testpaths does NOT collect; run
  directly it passes 6/6. Gap: not wired into the suite.
- Spec coverage: 77/95 modules have spec/*.md files. The 18 without:
  acmg_bayesian, cfd_offtarget, chem_descriptors, chem_similarity, crisprater,
  crisprscan_score, dti_bench, evidence_mining, mit_offtarget, molecule_eval,
  neuro_hub, pgx_guidelines, profile_hmm, qsar_bench, report_studio,
  rna_nussinov, structural_biophysics, synbio_studio.

## Open gaps (per standing rules)
- Spec-vs-implementation verification for the 77 specced modules is the next
  work item (this commit is suite + inventory only).
- Per-module attribution is import-linked; meta-gate/CLI test files (91) are
  counted at suite level only.
- 18 unspecced modules need spec definitions or docstring-based verification.

Artifacts: pytest_full.log (verbatim), inventory.json, scripts/selfaudit_inventory.py.

## Leg 2: registry, Tier-1 headlines, Tier-2 algorithm probes (2026-09-27 00:15)
- Registry reconciliation (STATUS claim: 95 registered slugs == package dirs):
  PASS live - omega.registry 95 slugs, 95 package dirs, exact match both ways.
- Tier-1 headline regeneration (scripts/evidence_report.py, hermetic, re-run by
  auditor): PASS - all headlines regenerate exactly as STATUS.md claims
  (2,720 unique pathogenic + 86 benign; canonical 2,414/2,414 U2 + 16/16 AT-AC;
  benign specificity 85/86; GC-donor 30 cases zero benign; VUS 2,196/179;
  conflicting 977/106; exonic 1,117+216; cryptic-recall null; ESRseq null).
- Tier-2 behavioral probes (independent of the suite, scripts/selfaudit_tier2_probes.py):
  7/7 algorithm claims PASS:
  * virtual_cell FBA: module objective 15.0 == independent scipy.linprog/HiGHS 15.0.
  * codon_opt CAI: exact match to independent Sharp & Li 1987 recomputation
    (0.8532194221288482) under the published Met/Trp/stop exclusion convention
    (auditor's first probe used the wrong convention - probe corrected, module correct).
  * codon_opt TASEP: deterministic under fixed input, density within [0,1].
  * evofold_4d ANM: independent Hessian rebuild gives exactly 6 zero modes
    (24-atom helix) and first nonzero frequency 0.1329 == module's 0.1329.
  * living_computer: deterministic under seed, nonnegative counts. NOTE: the
    stochastic engine is chemical-Langevin Euler-Maruyama, a real published
    method - STATUS.md's "Gillespie" wording is imprecise for this module.
  * synbio_wizard: true Gillespie SSA (exponential waiting times); Hill
    promoter at regulator==Kd gives exactly the analytic 0.51; SSA steady-state
    tail means consistent with analytic values (mRNA 2.44 vs 2.0, protein 43.1
    vs 40.0 at 500h, small-system noise).
  * alpha_fold_ui Chou-Fasman: poly-Ala scores 75% helix, poly-Pro/Gly zero.
- DOCUMENTATION FINDINGS (module code correct; ledger stale/imprecise):
  1. STATUS.md test-count line is stale: "1831 passed, 0 failed, 8 skipped
     (2026-09-24)" vs live 2,279 passed, 1 failed (drop54), 8 skipped.
  2. STATUS.md Tier-2 credits gene_analysis with "Hill-kinetics ODE models" -
     no Hill code exists in gene_analysis; Hill kinetics live in
     synbio_wizard.promoter_kinetics and virtual_cell._transcription_rates.
  3. living_computer stochastic method wording ("Gillespie") should read
     chemical-Langevin Euler-Maruyama.

## Open gaps after leg 2
- Tier-3 (~55 heuristic modules): smoke-level probe (import, docstring, probe
  call) not yet run.
- Tier-1 rows beyond the regenerated headlines are suite-locked hermetic
  claims (2,279 passing) - not individually re-verified live (network).
- 18 unspecced modules: docstring-level verification pending.

## Leg 3: Tier-3 smoke (78 modules, 2026-09-27 00:33)
- Import: 78/78 PASS.
- Module-level docstring (>=20 chars, states method): 61/78 PASS.
- The 17 without module docstrings document method via top-of-file comments
  (verified on samples: sigma70 consensus in promoter_lib, DRACH in
  rna_decoder) plus partial function-level docstrings: coverage ranges from
  100% (liquid_biopsy 9/9) down to 8% (syn_bio_studio 1/12, bioplayground
  1/6, dna_to_code 1/6, car_t_designer 1/6). STATUS.md's "each module's
  docstring states its method" is inaccurate for these 17 - DOCUMENTATION
  FINDING #4 (severity low: method comments exist, Python docstrings thin).
- Entry-point probes: only 5 of 78 expose zero-arg demo/diagnostics entries;
  all 5 PASS. The other 73 require arguments - smoke level does not probe
  them (suite-level tests cover behavior; counted at suite level).

## Documentation findings so far (for Atlas lane; code correct)
1. STATUS.md stale test count (1,831/0-fail vs live 2,279/1-fail).
2. STATUS.md Tier-2 gene_analysis Hill misattribution.
3. STATUS.md living_computer "Gillespie" wording (actually CLE/Euler-Maruyama).
4. STATUS.md "each module's docstring states its method" - 17/78 Tier-3
   modules have no module docstring; function-level docs 8-100%.

## Open gaps after leg 3
- 73 Tier-3 modules untested beyond suite+import (no zero-arg probe) - deeper
  behavioral verification would need per-module fixtures; out of smoke scope.
- Leg 4: docstring-level verification of the 18 unspecced modules.
- Leg 5: final per-module pass/fail table (95 rows).

## Leg 4: unspecced-module docstring verification (18 modules, 00:33)
18/18 PASS: every unspecced module carries a module docstring naming its
method (Tavtigian 2018, Doench 2016 CFD, Labuhn 2018 CRISPRater, Moreno-Mateos
CRISPRscan, Hsu/MIT, Durbin ch.5 profile HMM, Nussinov-Jacobson 1980, CPIC);
function-level doc coverage 67-100%. Detail: leg4_unspecced_docstrings.json.

## Leg 5: FINAL per-module table (95 rows, FINAL_TABLE.json)
- Module level: 95/95 PASS at their verification level
  (Tier-1 9: headline regen + suite | Tier-2 8: independent behavioral probes |
  Tier-3 78: import + docstring + suite attribution).
- Claim level: 1 FAIL - gene_analysis's STATUS.md Tier-2 claim ("Hill-kinetics
  ODE") is unattributable; the module itself passes at smoke+docstring level.
- Module docstrings: 76/95 have module-level docstrings; the 19 without
  document method via header comments (finding #4, low severity).
- Suite: 2,279 passed / 1 failed (drop54 style gate) / 8 env skips.
- pgx_guidelines: 6 passing tests live outside configured testpaths (wiring gap).

## Audit verdict (honest, per standing rules)
The repo's own suite is green except one style-gate failure (intentional
protocol stubs); every tier claim we could verify independently checks out
behaviorally; the defects found are documentation-level: a stale STATUS.md
test count, one Tier-2 misattribution, one method-wording imprecision, one
docstring-coverage gap, one test-wiring gap. No functional module defect was
found in this audit beyond what the repo itself documents.
Open gaps: Tier-1 rows beyond regenerated headlines are suite-locked, not
individually live-reverified; 73 Tier-3 modules have no zero-arg probe
(behavioral coverage via the suite only); Ensembl-dependent routes remain
Missing per the repo's own ledger (not re-tested).

## Leg 6: live network spot-check of Tier-1 claims (2026-09-27 01:14)
- openclinvar BRCA1 c.5266dup: live ClinVar match "Pathogenic, reviewed by
  expert panel, 3-star" - STATUS claim VERIFIED live.
- gnomAD r4 SCN1A: LOEUF 0.1067 (claim 0.107), strongly LOF-constrained -
  VERIFIED live.
- gnomAD r4 TP53: LOEUF 0.4184 (claim 0.418) - value verified, but the
  STATUS row "TP53 0.418 honestly no [constrained]" is STALE: the module
  moved to the documented gnomAD v4 guidance threshold LOEUF < 0.45
  (code comment cites the gnomAD help page, read 2026-09-25), under which
  TP53 is now flagged constrained. DOCUMENTATION FINDING #6 (stale ledger
  row after a deliberate threshold change; the v4-guidance citation is as
  documented in code, not independently re-verified by the auditor).
