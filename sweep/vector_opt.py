"""Sweep: vector_opt AAV capsid regions/conservation vs real AAV2(P03135)/AAV9(Q6JC40) VP1 sequences + exact scoring recomputation."""
import json, os, sys, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
from sugarcode.modules.vector_opt import core as V

CD = os.path.join(os.path.dirname(__file__), '..', 'data', 'aav')
res = {'module': 'vector_opt', 'sources': ['UniProtKB P03135 (AAV2 VP1)', 'UniProtKB Q6JC40 (AAV9 VP1)', 'Biopython PairwiseAligner', 'NCBI PubMed']}

# 1. real capsid sequences, global alignment, per-region identity vs module CONSERVATION
s2 = json.load(open(os.path.join(CD, 'P03135.json')))['sequence']['value']
s9 = json.load(open(os.path.join(CD, 'Q6JC40.json')))['sequence']['value']
from Bio import Align
al = Align.PairwiseAligner()
al.mode = 'global'; al.match_score = 1; al.mismatch_score = 0
al.open_gap_score = -0.5; al.extend_gap_score = -0.1
a = al.align(s2, s9)[0]
blocks2, blocks9 = a.aligned
pairs = [(s2[i:j], s9[k:l]) for (i, j), (k, l) in zip(blocks2, blocks9)]
aln2 = ''.join(p[0] for p in pairs); aln9 = ''.join(p[1] for p in pairs)
m = sum(1 for x, y in zip(aln2, aln9) if x == y and x != '-')
gaps = sum(1 for x, y in zip(aln2, aln9) if (x == '-') != (y == '-'))
res['alignment'] = {'aav2_len': len(s2), 'aav9_len': len(s9), 'aligned_len': len(aln2),
                    'identical': m, 'global_identity': round(m / max(len(aln2) - gaps, 1), 4)}
# AAV2 VP1 literature coordinates (1-based inclusive; Govindasamy 2006, Xie 2002, Lochrie 2006)
REG = {'VR-I': (262, 269), 'VR-II': (326, 332), 'VR-III': (381, 389), 'VR-IV': (450, 470),
       'VR-V': (488, 502), 'VR-VI': (531, 538), 'VR-VII': (546, 558), 'VR-VIII': (585, 593),
       'VR-IX': (712, 724), 'HI_loop': (660, 672), 'GH_loop': (590, 620)}
# map AAV2 1-based coords to aligned columns
pos2col = {}; col = 0
for i, ch in enumerate(aln2):
    if ch != '-':
        pos2col[i + 1] = col
    col += 1
def region_identity(lo, hi):
    cols = range(pos2col[lo], pos2col[hi] + 1)
    sub = [(aln2[c], aln9[c]) for c in cols]
    mm = sum(1 for x, y in sub if x == y and x != '-')
    return round(mm / max(len([s for s in sub if s[0] != '-' and s[1] != '-']), 1), 3)
rid = {name: region_identity(lo, hi) for name, (lo, hi) in REG.items()}
res['region_identity_aav2_vs_aav9'] = rid
res['module_conservation'] = V.CONSERVATION
res['conservation_notes'] = ('module CONSERVATION vs measured AAV2/AAV9 window identity: VR-IV module 0.3, VR-VIII 0.4, '
                             'GH_loop 0.35, HI_loop 0.7, threefold_spike 0.5 (spike proxied by VR-IV window)')
# 2. heparin receptor motif: AAV2 R585/R588 in VR-VIII (Summerford & Samulski 1998)
res['heparin_motif'] = {'AAV2_585': s2[584], 'AAV2_588': s2[587], 'AAV9_at_585_588': s9[583:589],
                        'aav2_has_R585_R588': s2[584] == 'R' and s2[587] == 'R',
                        'module_vrviii_receptor_contact': V.CAPSID_REGIONS['VR-VIII']['receptor_contact'],
                        'note': 'R585/R588 are the AAV2 heparan-sulfate receptor contact (Summerford 1998)'}
json.dump(res, open('/tmp/vector_opt_part1.json', 'w'), indent=1, default=str)
print(json.dumps({k: res[k] for k in ('alignment', 'region_identity_aav2_vs_aav9', 'module_conservation', 'heparin_motif')}, default=str))

# 3. exact scoring recomputation with the same RNG call sequence
def recompute(capsid, target, nab_escape, n_variants, seed):
    rng = random.Random(seed); out = []
    regs = list(V.CAPSID_REGIONS)
    for i in range(n_variants):
        region = rng.choice(regs); props = V.CAPSID_REGIONS[region]; n_mut = rng.randint(1, 4)
        gain = (0.15 * n_mut * (1 if props['receptor_contact'] else 0.2) * rng.uniform(0.7, 1.3))
        loss = (0.2 * n_mut if props['nab_epitope'] else 0.02 * n_mut) * rng.uniform(0.7, 1.3)
        integ = 1.0 - 0.08 * n_mut * (1 - V.CONSERVATION[region])
        trans = 0.5 + gain - (1 - integ)
        out.append({'variant_id': f'{capsid}-V{i+1}', 'region': region, 'n_mutations': n_mut,
                    'receptor_affinity_gain': round(gain, 3), 'nab_escape_score': round(min(loss, 1.0), 3),
                    'capsid_integrity': round(integ, 3), 'predicted_transduction': round(max(0.0, min(1.0, trans)), 3)})
    return out
score_checks = []
for seed in (1, 2, 3):
    got = V.engineer_capsid(capsid='AAV9', target_receptor='generic', nab_escape=True, n_variants=6, seed=seed)
    exp = recompute('AAV9', 'generic', True, 6, seed)
    score_checks.append({'seed': seed, 'variants_match': got['variants'] == sorted(exp, key=lambda v: -(v['nab_escape_score']*0.5 + v['predicted_transduction']*0.5)),
                         'lead_ok': got['lead'] == got['variants'][0],
                         'docking_gain_ok': got['docking_summary']['lead_receptor_gain'] == got['lead']['receptor_affinity_gain'],
                         'sera_ok': got['nab_panel_prediction']['sera_escape_fraction'] == round(got['lead']['nab_escape_score'] * 0.8, 3)})
res['scoring_recomputation'] = score_checks
# determinism: same seed twice
res['deterministic'] = V.engineer_capsid(seed=7)['variants'] == V.engineer_capsid(seed=7)['variants']
# directed evolution: round r must equal engineer_capsid(seed=1+r); lead chain
de = V.directed_evolution(capsid='AAV9', n_rounds=3, n_variants=6, seed=1)
de_ok = all(de['rounds'][r]['variants'] == V.engineer_capsid(capsid=de['lineage'][r], seed=1 + r)['variants'] for r in range(3))
res['directed_evolution'] = {'lineage': de['lineage'], 'round_seeds_ok': de_ok,
                             'dead_rng_line': 'directed_evolution creates rng=random.Random(seed) but never uses it (documented, cosmetic)'}
# 4. PubMed citations for receptor/region claims
from sugarcode.bio import entrez
import json as _j
def esearch(term):
    d = entrez._get('esearch.fcgi', {'db': 'pubmed', 'term': term, 'retmode': 'json', 'retmax': 3})
    return _j.loads(d)['esearchresult']['idlist']
summerford = esearch('Summerford Samulski heparan sulfate adeno-associated virus receptor')
govinda = esearch('Govindasamy adeno-associated virus capsid variable regions structure')
res['citations'] = {'summerford1998_candidates': summerford, 'govindasamy_candidates': govinda}
fn = os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_vector_opt.json')
json.dump(res, open(fn, 'w'), indent=1, default=str)
print(json.dumps({'scoring_recomputation': score_checks, 'deterministic': res['deterministic'], 'directed_evolution': res['directed_evolution'], 'citations': res['citations']}, default=str))
