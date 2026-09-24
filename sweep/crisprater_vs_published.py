"""crisprater.score vs author-computed CRISPRater scores for 3,141 sgRNAs from Labuhn 2018 Supplementary Table 3
(external screen / knockout / protein expression sheets; PMC5814880 supplementary, nar-02788-h-2017-File006.xlsx).
The pre-fix model used a GC4-14 window; published is GC4-13."""
import sys, os, json, openpyxl
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from sugarcode.modules.crisprater.core import score
wb = openpyxl.load_workbook(os.path.expanduser('~/mega/m01/data/crisprater/nar-02788-h-2017-File006.xlsx'), read_only=True)
rows = []
for sn in wb.sheetnames:
    hdr = None
    for r in wb[sn].iter_rows(values_only=True):
        if r and r[0] == 'target_number': hdr = r; continue
        if hdr and r and r[1]:
            d = dict(zip(hdr, r))
            if d.get('CRISPRater_score') is not None and len(str(d['sgRNA_sequence']).strip()) == 20:
                rows.append((sn, str(d['sgRNA_sequence']).strip().upper(), float(d['CRISPRater_score'])))
diffs = []
for sn, seq, pub in rows:
    got = score(seq).score
    diffs.append((abs(got - pub), sn, seq, pub, got))
diffs.sort(reverse=True)
mism = [d for d in diffs if d[0] > 5e-4]
# pre-fix comparison: GC over positions 4-14 (11 nt) instead of 4-13 (10 nt)
def pre_fix_score(seq):
    gc = (seq[3:14].count('G') + seq[3:14].count('C')) / 11
    vals = [gc, seq[19] == 'G', seq[2] in 'AT', seq[11] in 'AG', seq[5] == 'G',
            seq[3] in 'AT', seq[17] in 'AG', seq[4] in 'AC', seq[13] == 'G', seq[14] == 'A']
    w = [.14177385, .06966514, .04216254, .03303432, .02355430, -.04746424, -.04878001, -.06981921, -.07087756, -.081607]
    return .6505037 + sum(float(v) * x for v, x in zip(vals, w))
pre_mism = sum(abs(pre_fix_score(seq) - pub) > 5e-4 for _, seq, pub in rows)
res = {'reference': 'Labuhn 2018 NAR 46(3):1375 (PMC5814880) Supplementary Table 3 author-computed CRISPRater scores',
 'n_sgrnas': len(rows), 'sheets': sorted(set(sn for sn, _, _ in rows)),
 'post_fix': {'mismatches_gt_5e-4': len(mism), 'max_abs_diff': round(diffs[-1][0] if mism and False else min(d[0] for d in diffs), 12),
              'max_abs_diff_overall': round(max(d[0] for d in diffs), 12),
              'worst': [{'sheet': d[1], 'seq': d[2], 'published': d[3], 'got': d[4]} for d in mism[:5]]},
 'pre_fix_gc_window_4_14': {'mismatches_gt_5e-4': pre_mism},
 'class_thresholds_published': {'low_max_exclusive': 0.56, 'high_min_exclusive': 0.74},
 'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else ''}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_crisprater.json'), 'w'), indent=1)
print(json.dumps(res, indent=1)[:900])
