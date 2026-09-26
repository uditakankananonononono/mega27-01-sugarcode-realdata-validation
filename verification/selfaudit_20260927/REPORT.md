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
