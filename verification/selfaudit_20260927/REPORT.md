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
