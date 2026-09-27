# Atlas AI self-audit — 2026-09-27

Repo: uditakankananonononono/atlas-ai (public), audited HEAD `ea2755cdaf1263384e8f303aa594a1f94a9e540d` ("Wire free local Hermes and explicit owner OpenClaw bridge", 93-commit history). Audit is read-only; no writes to atlas-ai.

## Test suite (run, not claimed)
- Env: Python 3.12 venv via uv, `uv pip install -e '.[dev]'` (251 packages, exit 0) + `en_core_web_sm`.
- `python -m pytest tests -q`: **9,071 passed / 1 failed / 12 skipped in 373.7s** (log: pytest_full.log).
- The 1 failure is a stale assertion, not a behavior defect: `test_thin_module_wiring_wave.py::test_m22_ships_free_official_discovery_sources` expects collectors `['github','pypi','npm']`; `default_collectors()` actually ships 11 free collectors (github, pypi, npm, devto, wordpress, hackernews, medium, itunes, itunes-episodes, gitlab, codeberg). Code grew; test not updated. Finding A.

## Module inventory
- 28 module dirs under backend/app/modules (m00–m25 plus extra m17_advice_essay and m25_knowledge_copilot).
- registry.py IMPLEMENTED_SPECS: 26 specs; all mounted at /api/v1 behind require_tenant via the main.py loop (verified by import; app boots, 44 top-level route entries incl. mounts).
- m25_knowledge_copilot has its own router mounted directly in main.py. m17_advice_essay is not in the registry and has no direct mount in main.py (has routes.py + 19 test files; mounting path not confirmed). Finding D.
- 451 test files total under tests/ (368 in tests/modules).

## Ledger claims verified against HEAD
- audits/ledger-140.json: 140/140 rows status verified-pushed; all 140 referenced implementation paths and all referenced test files exist at HEAD.
- docs/IMPLEMENTATION_AUDIT.md: 2,010 rows parsed; statuses 1,719 verified-practitioner-depth + 291 upgraded-practitioner-depth; 15/15 randomly sampled named row-tests exist in the named test files.

## Behavioral spot-checks (live, this audit)
1. Free-first model routing VERIFIED: `paid_allowed()` returns False by default; with all keys/ATLAS_* scrubbed, `default_chain()` = [ollama (local), openai_compat (local)] only, zero paid entries. Paid routes require ATLAS_ALLOW_PAID=true (model_catalog.py:96,116).
2. Stripe test-mode enforced in code: stripe_client.py raises ValueError unless secret key starts `sk_test_`.
3. Approval center m00 spec/router import and mount verified; approval tests green in suite.
4. m22 collectors live-verified by direct call (11 free collectors, superset of test's stale 3).

## Findings (honest negatives)
- A. 1 stale failing test (m22 collector count, above).
- B. Provenance gap: the ledgers' generation commits are not in current public history — IMPLEMENTATION_AUDIT.md header commit a6e59c9d…, ledger audited_head 1a3b800, and all 7 distinct implementation_commit hashes are not git objects in the 93-commit clone (history appears rewritten/squashed). All ledger claims were therefore verified by file existence and behavior at HEAD, not by the cited commits.
- C. Stale audit-doc paths after module reorg: 115 rows' impl paths and 165 rows' test paths do not resolve verbatim (m18_project_builder → m14_project_builder/engineering_methods_510_574.py; m12 humanities_support_1860_1909 → m09_knowledge_workspace; env tests → test_m16_climate_environment_1610_1659.py). Every renamed equivalent exists and its tests pass; the doc table was not regenerated after the move.
- D. 2 module dirs outside the registry: m17_advice_essay (no confirmed mount), m25_knowledge_copilot (direct mount, no registry spec).

## Open gaps (not covered by this audit)
- 12 skipped tests are environment-dependent; suite ran on a sandbox, not project CI.
- Live external acceptance not exercised: Gmail OAuth/watch, Stripe live, Devpost endpoint, FreeTSA timestamps, RSS fetches — verified only via offline tests and code inspection.
- Frontend (Next.js dashboard) not audited.
- Docker/bubblewrap sandbox backends not run on a live host.
