# F1 reproducibility audit - experiment 1 (2026-09-26)
Preregistered arm (PREREGISTRATION_all95.md amendment #2): "Do CRISPR design tools agree
on the same biological reality?" First data: TP53 RefSeqGene locus (NG_017013.2, 32,772 bp).

## Result 1a - enumeration layer: perfectly reproducible
Two independent NGG enumerators (CRISPOR findPams/findPat semantics vs independent regex
reimplementation) agree on ALL 4,825 targetable guides: same positions, strands,
sequences, including 3 boundary guides within 23bp of a locus end. Existence/PAM layer
is stable. (round6_enumeration_audit.json)

## Result 1b - scoring layer: 8.5% of guides mis-scored by one displaced coefficient
CRISPOR's production CRISPRscan vs the supplement-verified implementation (sugarcode
crisprscan_score, cross-checked against Moreno-Mateos 2015 MOESM640):
- Parameter diff is EXACTLY one coefficient: dinucleotide AA, weight -0.097377097, at
  position 18 in CRISPOR vs position 19 in the original model (round6_param_diff.json;
  upstream issue maximilianh/crisporWebsite#76).
- Apples-to-apples (both sides int-truncated the way CRISPOR truncates): 411/4,824 TP53
  guides (8.5%) receive a different score solely from this defect.
- Rank impact: top-100 guides by each implementation overlap in only 92 - eight guides
  enter/leave the top-100 selection set because of the defect.
- Contexts carrying AA at pos 18: 338; at pos 19: 413 (affected universe).
- Honest noise separation: raw comparison shows 54.5% of guides differ, but that is
  CRISPOR's int(100*x) truncation (quantization up to 0.01); the defect-only rate is 8.5%.
  Both numbers preserved; the naive 54.5% must not be quoted as the defect rate.
- Raw-score Spearman 0.987; max single-guide delta 0.107 (on a 0-1 score).

Conclusion for the reproducibility map: on this locus, guide existence/PAM assignment is
reproducible across independent implementations, while the scoring layer carries a real,
rank-changing defect traceable to a single mislocated coefficient.
