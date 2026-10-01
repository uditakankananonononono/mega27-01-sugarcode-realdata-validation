# Validation checkpoint, October 1, 2026

Product: https://github.com/uditakankananonononono/sugarcode-ai
Pinned product HEAD: 1ccbd0305cd9b639b06324fa136654d686483b7a.
Validation destination: https://github.com/uditakankananonononono/mega27-01-sugarcode-realdata-validation
Evidence-patch base: 754cb3455762f0f526aa6bcc685e11a9ed9d0c58.
No product edits or deployments. Evidence remains additive pending durable commit.

## Landed audit evidence
- Two final complete product-suite runs: 2,319 passed, zero failures/skips; identical test identity/outcome sets.
- 95 registered/package modules, 95 imports/method-docstrings/test-file attribution; file attribution is non-exclusive, not semantic spec compliance.
- Seven selected independent behavioral probes repeat identically.
- Six BBBC001 images: mean count deviation 8.1824%, passes 11% suite gate but worse than published 6.2% benchmark.
- All optional suite skips closed in local isolated environment.

## Confirmed findings
1. Promoter Lib homopolymer_max is numerically wrong on four independent controls; one positive control passes.
2. Generated Ecosystem POST SDK route absent from shipped app: in-process 404 twice, metadata GET positive control 200.
3. CellPainter nearest_mechanism cosine changes 1.0 -> 0.0 when only dictionary insertion order changes.
4. Promoter length and MetaboDesigner max_steps inputs ignored; repeated output probes.
5. RareNet inheritance input validated/echoed but not used in score/ranking.
6. NeuroPlan geometry-only planning routes return unsupported image-derived confidence/tractography labels.
7. Spec attribution contamination: 35 duplicate long paragraphs across six files; named subjects conflict with containing module in some files.
8. AlphaFold UI MSA MI changes do not affect sequence-generated geometry or confidence analogs; source explicitly disclaims learned folding model.

## Per-spec verification remaining
Preliminary punctuation inventory: 1,436 units from 77 documents. It is not a fixed requirement denominator. Five tranches record 49 semantic reviews and 15 attribution-blocked units. 1,372 units remain unreviewed. No full compliance percentage supported. Prior reviewed tranches: NeuroPlan/RareNet/Promoter, DNA-to-Code/Omega Stats/AlphaFold intro, remaining AlphaFold, Ecosystem/Nexus Support, CellPainter/4D.

Next bounded tranche: remaining short module specs, normalize ambiguous subject attribution without silently assigning foreign paragraphs, build independent checks for action-driving numerical claims. Recheck source HEAD before extending findings. Preserve aspirations versus shipped-scope distinctions. Do not count rendered visual claims without inspecting actual pixels. Product fixes require owner direction and are not included in this read-only audit.

The superseded functional-regression patch containing pyc must be discarded; use sugarcode-functional-regression-20261001-clean.patch. patch_manifest.json lists clean patches and exact destination files.
