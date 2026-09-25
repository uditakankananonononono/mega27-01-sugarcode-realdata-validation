# Independent external validation of SugarCode AI

Research question: do the SugarCode AI platform's internally verified claims
survive held-out public benchmarks its own tests never touched?

- Status ledger: STATUS.md (verified / pending / missing, updated each push)
- Validation package: `src/sc_validate` (pytest; hermetic fixtures in `data/`)
- Live re-runs: `scripts/` (NCBI eutils, RCSB, UniProt, EBI services)
- Results: `benchmarks/results.json` + `paper/` (Times, LaTeX)
- Numbered-formula audit: `audits/numbered_formula_gate.md` records 10 rendered, benchmark-used numbered equations; this does not add research services or biological accessions
- Provenance: `data/README.md` - every fixture lists source + retrieval date

Companion platform repo: github.com/uditakankananonononono/sugarcode-ai
