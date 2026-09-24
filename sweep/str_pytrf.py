"""Validate sugarcode str_scope.find_strs against pytrf 1.5.0 on real GenBank genomes."""
import json, sys, time, glob, os
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import pytrf
from Bio import SeqIO
from sugarcode.modules.str_scope.core import find_strs
TH = {1: 12, 2: 7, 3: 5, 4: 4, 5: 4, 6: 4}  # pytrf defaults
def ours(seq):
    out = []
    for h in find_strs(seq, 1, 6, 4):
        if h['repeats'] >= TH[h['unit_len']]:
            out.append((h['start'], h['start'] + h['length'], h['unit_len'], h['canonical_unit']))
    return out
def theirs(name, seq):
    out = []
    for s in pytrf.STRFinder(name, seq, *[TH[k] for k in range(1, 7)]):
        m = s.motif; can = min(m[i:] + m[:i] for i in range(len(m)))
        out.append((s.start - 1, s.end, len(m), can))
    return out
def match(a, b):
    # one-to-one by unit length + canonical unit + overlap
    used = set(); tp = []; exact = 0
    idx = {}
    for j, y in enumerate(b): idx.setdefault(y[3], []).append(j)
    for x in a:
        for j in idx.get(x[3], []):
            y = b[j]
            if j not in used and x[0] < y[1] and y[0] < x[1]:
                used.add(j); tp.append((x, y)); exact += (x[0] == y[0] and x[1] == y[1]); break
    return tp, exact, used
res = {'tool': 'pytrf 1.5.0', 'thresholds': TH, 'genomes': []}
for gb in sorted(glob.glob(os.path.expanduser('~/mega/m01/data/genbank/*.gb'))):
    rec = next(SeqIO.parse(gb, 'genbank'))
    seq = str(rec.seq).upper()
    if set(seq) - set('ACGT'):
        seq = ''.join(c if c in 'ACGT' else 'N' for c in seq)
    t = time.time(); a = ours(seq); ta = time.time() - t
    t = time.time(); b = theirs(rec.id, seq); tb = time.time() - t
    tp, exact, used = match(a, b)
    fn = [b[j] for j in range(len(b)) if j not in used]
    fp = [x for x in a if x not in [p[0] for p in tp]]
    res['genomes'].append({'accession': rec.id, 'length': len(seq), 'ours': len(a), 'pytrf': len(b),
        'matched': len(tp), 'exact_span': exact, 'only_ours': len(fp), 'only_pytrf': len(fn),
        'examples_only_ours': fp[:5], 'examples_only_pytrf': fn[:5], 'sec_ours': round(ta, 2), 'sec_pytrf': round(tb, 3)})
    print(res['genomes'][-1], flush=True)
tot = {k: sum(g[k] for g in res['genomes']) for k in ('ours', 'pytrf', 'matched', 'exact_span', 'only_ours', 'only_pytrf')}
res['totals'] = tot
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_str_pytrf.json'), 'w'), indent=1)
print(tot)
# Allele-band audit: legacy delta rule vs GeneReviews bands (NBK1305, NBK1384)
from sugarcode.modules.str_scope.core import expansion_call, LOCUS_BANDS
def gr_htt(n): return 'normal' if n<=26 else 'intermediate' if n<=35 else 'pathogenic'
def gr_fmr1(n): return 'not_pre' if n<=54 else 'premutation' if n<=200 else 'full'
def legacy(n, ref):
    c = expansion_call(n, ref)['classification']
    return 'normal' if c.startswith('normal') else 'intermediate' if c.startswith('intermediate') else 'pathogenic'
bands = {}
htt_bad = [n for n in range(5, 121) if legacy(n, 26) != gr_htt(n)]
fmr_map = {'normal': 'not_pre', 'intermediate': 'premutation', 'pathogenic': 'full'}
fmr_bad = [n for n in range(5, 301) if fmr_map[legacy(n, 44)] != gr_fmr1(n)]
bands['HTT_legacy_misclassified'] = htt_bad
bands['FMR1_legacy_misclassified_count_5_300'] = len(fmr_bad)
bands['FMR1_legacy_misclassified_ranges'] = [min(fmr_bad), max(fmr_bad)] if fmr_bad else []
bands['banded_errors'] = sum(1 for n in range(5, 121) if not expansion_call(n, 26, 'CAG', 'HTT')['classification'].startswith({'normal':'normal','intermediate':'intermediate','pathogenic':('reduced','full')}[gr_htt(n)] if gr_htt(n)!='pathogenic' else ('reduced','full')))
bands['sources'] = ['https://www.ncbi.nlm.nih.gov/sites/books/NBK1305/table/huntington.T.methods_to_characterize_htt/?report=objectonly', 'https://www.ncbi.nlm.nih.gov/books/NBK1384/table/fragilex.T.types_of_fmr1_repeat_expansio/']
res['allele_bands'] = bands
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_str_pytrf.json'), 'w'), indent=1)
print(bands)
