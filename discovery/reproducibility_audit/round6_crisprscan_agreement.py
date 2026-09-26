#!/usr/bin/env python3
"""F1 reproducibility audit, experiment 1b: score-layer agreement.
CRISPOR's production CRISPRscan (crisporEffScores.calcCrisprScanScores - has the AA19
misplacement, crisporWebsite#76) vs sugarcode's corrected implementation on all TP53
NGG guides with available 35nt context. Quantifies real-world impact of the defect:
score deltas, rank instability, and the affected-context subset (AA at pos 18 vs 19)."""
import re, json, sys, math
sys.path.insert(0, '/home/sandbox/sugarcode-ai/src')
from sugarcode.modules.crisprscan_score.core import score as sg_score
import importlib.util
spec = importlib.util.spec_from_file_location('ce', '/home/sandbox/wave1/oracle/crisporEffScores.py')
ce = importlib.util.module_from_spec(spec); spec.loader.exec_module(ce)

seq = ''.join(l.strip() for l in open('/home/sandbox/wave1/fixtures/tp53_NG_017013.2.fa') if not l.startswith('>')).upper()
L = len(seq)
comp = {'A':'T','T':'A','G':'C','C':'G'}
def rc(s): return ''.join(comp[b] for b in reversed(s))

guides = []
for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', seq):
    guides.append((m.start(), '+'))
rcseq = rc(seq)
for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', rcseq):
    guides.append((L - (m.start() + 23), '-'))

contexts, keys = [], []
for pos, strand in guides:
    # 35-mer in guide-strand orientation: 6bp flank + 20 guide + 3 PAM + 6 flank
    if strand == '+':
        if pos < 6 or pos + 29 > L: continue
        ctx = seq[pos-6:pos+29]
    else:
        # minus-strand guide+PAM spans forward pos..pos+22; guide-oriented 35mer =
        # rc(forward[pos-6:pos+29]) = 6 flank + 20 guide + 3 PAM + 6 flank
        if pos < 6 or pos + 29 > L: continue
        ctx = rc(seq[pos-6:pos+29])
    if len(ctx) != 35: continue
    contexts.append(ctx); keys.append((pos, strand))

crispor_scores = ce.calcCrisprScanScores(contexts)  # int 0-100
rows = []
skipped_ref = 0
for k, ctx, cs in zip(keys, contexts, crispor_scores):
    try:
        sg = sg_score(ctx).score
    except Exception:
        skipped_ref += 1
        continue
    rows.append((k, ctx, cs/100.0, sg, cs, int(100*sg)))
n = len(rows)
diffs = [abs(r[2] - r[3]) for r in rows]
nz = [(r[0], r[1], r[2], r[3]) for r in rows if abs(r[2] - r[3]) > 0.005]
defect_only = [(r[0], r[1], r[4], r[5]) for r in rows if r[4] != r[5]]
aa18 = sum(1 for r in rows if r[1][17:19] == 'AA')
aa19 = sum(1 for r in rows if r[1][18:20] == 'AA')
# rank instability: top 100 by each
by_c = sorted(rows, key=lambda r: -r[2])[:100]
by_s = sorted(rows, key=lambda r: -r[3])[:100]
overlap = len(set(r[0] for r in by_c) & set(r[0] for r in by_s))
import statistics
out = {
 'locus': 'NG_017013.2 TP53', 'n_guides_scored': n,
 'crispor_impl': 'crisporWebsite crisporEffScores.calcCrisprScanScores (AA19 defect, issue #76)',
 'reference_impl': 'sugarcode crisprscan_score (verified vs Moreno-Mateos 2015 supplement)',
 'n_different_gt_0.005': len(nz),
 'frac_different': len(nz)/n,
 'n_skipped_ref': skipped_ref,
 'n_defect_only_apples_to_apples': len(defect_only),
 'frac_defect_only': len(defect_only)/n,
 'max_abs_diff': max(diffs), 'mean_abs_diff': statistics.mean(diffs),
 'contexts_with_AA_at_pos18': aa18, 'contexts_with_AA_at_pos19': aa19,
 'top100_overlap': overlap,
 'spearman': None}
from scipy import stats as sps
out['spearman'] = float(sps.spearmanr([r[2] for r in rows], [r[3] for r in rows])[0])
ex = sorted(nz, key=lambda r: -abs(r[2]-r[3]))[:10]
out['largest_delta_examples'] = [{'pos': k[0], 'strand': k[1], 'context': ctx, 'crispor': c, 'reference': s} for k, ctx, c, s in ex]
print(json.dumps(out, indent=1))
json.dump(out, open('/home/sandbox/wave1/logs/round6_crisprscan_agreement.json','w'), indent=1)
