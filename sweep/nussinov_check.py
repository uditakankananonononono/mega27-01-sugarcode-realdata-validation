"""riboswitch.nussinov_fold vs (a) an independent Nussinov DP with a different decomposition, (b) exhaustive enumeration over all non-crossing matchings for n<=16, (c) ViennaRNA RNA.fold MFE structures, on real Rfam sequences."""
import sys, os, json, random, itertools
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import RNA
from functools import lru_cache
from sugarcode.modules.riboswitch.core import nussinov_fold
AL = {"AU","UA","GC","CG","GU","UG"}

def indep_nussinov(seq, min_loop=3):
    # independent formulation: for (i,j), either i unpaired, or i pairs with some k
    n = len(seq)
    @lru_cache(maxsize=None)
    def f(i, j):
        if j - i <= min_loop: return 0
        best = f(i + 1, j)
        for k in range(i + min_loop + 1, j + 1):
            if seq[i] + seq[k] in AL:
                best = max(best, 1 + f(i + 1, k - 1) + f(k + 1, j))
        return best
    return f(0, n - 1)

def exhaustive_max(seq, min_loop=3):
    # brute force over all non-crossing matchings; positions used at most once
    n = len(seq); best = [0]
    def rec(i, pairs, occ, cnt):
        if i >= n:
            best[0] = max(best[0], cnt); return
        if i in occ: rec(i + 1, pairs, occ, cnt); return
        if cnt + (n - i) // 2 <= best[0]: return  # cannot beat best
        rec(i + 1, pairs, occ, cnt)  # i unpaired
        for j in range(i + min_loop + 1, n):
            if j in occ or seq[i] + seq[j] not in AL: continue
            if any((a < i < b < j) or (i < a < j < b) for a, b in pairs): continue
            rec(i + 1, pairs + ((i, j),), occ | {i, j}, cnt + 1)
    rec(0, (), frozenset(), 0); return best[0]

def read_fam(p, n_take, rng):
    seqs = []; cur = None
    for line in open(p):
        if line.startswith('>'):
            if cur: seqs.append(cur)
            cur = ''
        elif cur is not None: cur += line.strip().replace('.', '').replace('-', '')
    if cur: seqs.append(cur)
    seqs = [s for s in seqs if 30 <= len(s) <= 150 and set(s) <= set('ACGU')]
    rng.shuffle(seqs); return seqs[:n_take]

rng = random.Random(7)
fams = {'RF00174': 'cobalamin riboswitch', 'RF00059': 'TPP riboswitch', 'RF00050': 'FMN riboswitch', 'RF00504': 'glycine riboswitch'}
real = []
for f, name in fams.items():
    p = os.path.expanduser(f'~/mega/m01/data/rfam/{f}.fa')
    got = read_fam(p, 50, rng); real += [(f, s) for s in got]
    print(f, len(got))
# exhaustive oracle: random + short real prefixes
oracle = []
short = [s[:16] for _, s in real[:40]] + [''.join(rng.choice('ACGU') for _ in range(rng.randint(8, 16))) for _ in range(60)]
mism = 0
for s in short:
    m = nussinov_fold(s)['pair_count']; e = exhaustive_max(s)
    if m != e: mism += 1; oracle.append({'seq': s, 'module': m, 'exhaustive': e})
# independent DP + ViennaRNA on real sequences
rows = []; dp_mism = 0
for fam, s in real:
    m = nussinov_fold(s)
    ind = indep_nussinov(s)
    if m['pair_count'] != ind: dp_mism += 1
    v_struct, v_mfe = RNA.fold(s)
    vp = set(); st = []
    for i, c in enumerate(v_struct):
        if c == '(': st.append(i)
        elif c == ')': vp.add((st.pop(), i))
    mp = set(map(tuple, m['pairs']))
    tp = len(mp & vp); f1 = 2 * tp / (len(mp) + len(vp)) if mp or vp else 1.0
    rows.append({'fam': fam, 'len': len(s), 'module_pairs': len(mp), 'vienna_pairs': len(vp), 'f1_vs_vienna': round(f1, 3)})
res = {'reference': 'independent Nussinov DP + exhaustive non-crossing-matching oracle (n<=16) + ViennaRNA ' + RNA.__version__ + ' RNA.fold',
 'n_oracle': len(short), 'oracle_mismatches': mism, 'oracle_examples': oracle[:5],
 'n_real': len(real), 'independent_dp_mismatches': dp_mism, 'families': fams,
 'f1_vs_vienna_median': round(sorted(r['f1_vs_vienna'] for r in rows)[len(rows) // 2], 3),
 'f1_vs_vienna_by_family': {f: round(sum(r['f1_vs_vienna'] for r in rows if r['fam'] == f) / max(1, sum(1 for r in rows if r['fam'] == f)), 3) for f in fams},
 'pair_ratio_module_over_vienna_median': round(sorted(r['module_pairs'] / max(1, r['vienna_pairs']) for r in rows)[len(rows) // 2], 3),
 'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else ''}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_nussinov.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in res.items() if k not in ('oracle_examples',)}, default=str, indent=0))
