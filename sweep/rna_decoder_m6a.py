"""rna_decoder.predict_m6a vs SRAMP mature-mode benchmark labels (Zhou 2016, Genome Biol 17:138;
test CSV mirrored in deepSRAMP repo) on 40 Ensembl cDNA transcripts fetched live from
rest.ensembl.org, plus miCLIP motif-composition check vs GSE63753 (Linder 2015, Cell 163:1667)
single-nucleotide m6A sites (9,536 BED rows, hg38)."""
import sys, os, json, csv, re
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from collections import defaultdict
from sugarcode.modules.rna_decoder.core import predict_m6a
from sugarcode.bio.sequence import find_motif, gc_content

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data', 'm6a')
res = {}

# ---- load labels for the 40 selected transcripts
sel = [l.strip() for l in open(os.path.join(DATA, 'selected_transcripts.txt')) if l.strip()]
labs = defaultdict(list)
r = csv.reader(open(os.path.join(DATA, 'sramp_mature_test.csv'), encoding='utf-8-sig')); next(r)
selset = set(sel)
for tid, p, c in r:
    if tid in selset:
        labs[tid].append((int(p), 1 if c == '1' else 0))
res['n_transcripts'] = len(sel)
res['n_pos'] = sum(1 for t in sel for _, c in labs[t] if c == 1)
res['n_neg'] = sum(1 for t in sel for _, c in labs[t] if c == 0)

# ---- coordinate convention: positives should land on A
seqs = {t: open(os.path.join(DATA, 'ensembl', f'{t}.fa')).read().strip().upper() for t in sel}
a1 = a0 = tot = 0
for t in sel:
    s = seqs[t]
    for p, c in labs[t]:
        if c != 1: continue
        tot += 1
        if p - 1 < len(s) and s[p - 1] == 'A': a1 += 1
        if p < len(s) and s[p] == 'A': a0 += 1
res['pos_base_is_A_1based'] = a1
res['pos_base_is_A_0based'] = a0
res['n_pos_checked'] = tot
OFF = -1 if a1 >= a0 else 0  # choose convention

# ---- run predict_m6a, score each labeled site by the DRACH centered on it
scores, y, drift = [], [], 0
for t in sel:
    s = seqs[t]
    pred = {x['position']: x['m6a_probability'] for x in predict_m6a(s, threshold=0.0)}
    drach_centers = {m + 2 for m in find_motif(s, 'DRACH')}
    for p, c in labs[t]:
        idx = p + OFF
        if idx >= len(s) or s[idx] != 'A' or idx not in drach_centers:
            drift += 1  # annotation drift vs 2016 Ensembl: not a DRACH A in current sequence
            continue
        scores.append(pred.get(idx, 0.0))
        y.append(c)
res['excluded_annotation_drift'] = drift
res['evaluated'] = len(y)
res['evaluated_pos'] = sum(y)
res['evaluated_neg'] = len(y) - sum(y)

def auc(sc, yy):
    order = sorted(range(len(sc)), key=lambda i: sc[i])
    ranks = [0.0] * len(sc); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and sc[order[j + 1]] == sc[order[i]]: j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1): ranks[order[k]] = r
        i = j + 1
    P = sum(yy); N = len(yy) - P
    if P == 0 or N == 0: return None
    return (sum(r for r, c in zip(ranks, yy) if c == 1) - P * (P + 1) / 2) / (P * N)

res['auc_module_score'] = round(auc(scores, y), 4)  # module defaults (auto-ORF CDS after fix)
# reproducible pre-fix baseline: the old whole-sequence default, via explicit whole-sequence args
scores0 = []
for t in sel:
    s = seqs[t]
    drach_centers = {m + 2 for m in find_motif(s, 'DRACH')}
    pred = {x['position']: x['m6a_probability'] for x in predict_m6a(s, cds_start=0, cds_end=len(s), threshold=0.0)}
    for p, c in labs[t]:
        idx = p + OFF
        if idx >= len(s) or s[idx] != 'A' or idx not in drach_centers:
            continue
        scores0.append(pred.get(idx, 0.0))
res['auc_module_score_old_wholetx_default'] = round(auc(scores0, y), 4)
# exposure-only baseline (the module's only varying feature under default CDS-spanning args)
expo = []
for t in sel:
    s = seqs[t]
    drach_centers = {m + 2 for m in find_motif(s, 'DRACH')}
    for p, c in labs[t]:
        idx = p + OFF
        if idx >= len(s) or s[idx] != 'A' or idx not in drach_centers:
            continue
        w = s[max(0, idx - 15):idx + 16]
        expo.append((1 - gc_content(w), c))
keep = expo
res['auc_exposure_baseline'] = round(auc([e for e, _ in keep], [c for _, c in keep]), 4)
# DRACH recall of the benchmark labels
res['drach_recall_all_labels'] = round(1 - drift / max(1, res['n_pos'] + res['n_neg']), 4)
# score separation
mp = [s for s, c in zip(scores, y) if c == 1]; mn = [s for s, c in zip(scores, y) if c == 0]
res['mean_score_pos'] = round(sum(mp) / len(mp), 4)
res['mean_score_neg'] = round(sum(mn) / len(mn), 4)

# ---- CDS-aware run: real coding-region coordinates from Ensembl type=cds aligned into the cDNA
scores2, y2, reg_pos, reg_neg = [], [], defaultdict(int), defaultdict(int)
for t in sel:
    cdsf = os.path.join(DATA, 'ensembl', f'{t}.cds')
    if not os.path.exists(cdsf): continue
    cs, ce = map(int, open(cdsf).read().split())
    s = seqs[t]
    drach_centers = {m + 2 for m in find_motif(s, 'DRACH')}
    pred = {x['position']: x for x in predict_m6a(s, cds_start=cs, cds_end=ce, threshold=0.0)}
    for p, c in labs[t]:
        idx = p + OFF
        if idx >= len(s) or s[idx] != 'A' or idx not in drach_centers:
            continue
        x = pred.get(idx)
        if x is None: continue
        scores2.append(x['m6a_probability']); y2.append(c)
        (reg_pos if c == 1 else reg_neg)[x['region']] += 1
res['cds_aligned_transcripts'] = sum(1 for t in sel if os.path.exists(os.path.join(DATA, 'ensembl', f'{t}.cds')))
res['evaluated_cds'] = len(y2)
res['evaluated_cds_pos'] = sum(y2)
res['auc_module_score_with_cds'] = round(auc(scores2, y2), 4)
res['region_counts_pos'] = dict(reg_pos)
res['region_counts_neg'] = dict(reg_neg)

# ---- GSE63753 miCLIP: motif composition of 9,536 single-nucleotide m6A sites
IUPAC = {'D': 'AGT', 'R': 'AG', 'A': 'A', 'C': 'C', 'H': 'ACT'}
n = 0; drach = 0; motifs = defaultdict(int)
for line in open(os.path.join(DATA, 'GSE63753_hek293.abcam.CIMS.m6A.9536.hg38.bed')):
    f = line.rstrip('\n').split('\t')
    if len(f) < 4: continue
    mm = re.search(r'_([ACGTN]{5})$', f[3])
    if not mm: continue
    m5 = mm.group(1); n += 1; motifs[m5] += 1
    if all(b in IUPAC[c] for b, c in zip(m5, 'DRACH')): drach += 1
res['miclip_sites_parsed'] = n
res['miclip_drach_fraction'] = round(drach / n, 4)
res['miclip_top5'] = sorted(motifs.items(), key=lambda kv: -kv[1])[:5]
res['datasets'] = ['SRAMP_mature_test(Zhou2016)', 'GSE63753'] + [f'ensembl:{t}' for t in sel]
res['n_datasets'] = len(res['datasets'])

out = os.path.join(HERE, '..', 'benchmarks', 'sweep_rna_decoder.json')
json.dump(res, open(out, 'w'), indent=2)
print(json.dumps(res, indent=2))
