"""Sweep: dark_genome vs JASPAR matrices + UCSC hg38 CpG islands + internal math."""
import json, math, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
import numpy as np
from sugarcode.modules.dark_genome import core as D

JD = os.path.join(os.path.dirname(__file__), '..', 'data', 'jaspar')
UD = os.path.join(os.path.dirname(__file__), '..', 'data', 'ucsc')
IUPAC = {'A':'A','C':'C','G':'G','T':'T','R':'AG','Y':'CT','W':'AT','S':'GC','K':'GT','M':'AC','N':'ACGT',
         'V':'ACG','H':'ACT','D':'AGT','B':'CGT'}
res = {'module': 'dark_genome', 'sources': ['https://jaspar.elixir.no/api/v1/matrix/<id>/', 'https://api.genome.ucsc.edu/getData/track?track=cpgIslandExt', 'https://api.genome.ucsc.edu/getData/sequence?genome=hg38']}

# 1. JASPAR consensus compatibility
MAT = {'AP-1': 'MA0099.3', 'NF-kB': 'MA0107.1', 'p53': 'MA0106.3', 'TATA-box': 'MA0108.2',
       'E-box': 'MA0093.3', 'GC-box': 'MA0079.4', 'CCAAT': 'MA0060.2', 'HNF4': 'MA0114.4'}
COMP = {}
for a, b in [('A','T'),('T','A'),('C','G'),('G','C'),('R','Y'),('Y','R'),('W','W'),('S','S'),('K','M'),('M','K'),('N','N'),('V','B'),('B','V'),('H','D'),('D','H')]:
    COMP.setdefault(a, b)
def topseq(pfm):
    return ''.join(max('ACGT', key=lambda b: pfm[b][i]) for i in range(len(pfm['A'])))
def consensus(pfm):
    out = []
    for i in range(len(pfm['A'])):
        tot = sum(pfm[b][i] for b in 'ACGT') or 1
        fr = {b: pfm[b][i] / tot for b in 'ACGT'}
        top = max(fr, key=fr.get)
        if fr[top] >= 0.5:
            out.append(top)
        else:
            allowed = ''.join(sorted(b for b in 'ACGT' if fr[b] >= 0.25))
            out.append(next(k for k, v in IUPAC.items() if set(v) == set(allowed) and len(k) == 1))
    return ''.join(out)
def rc(x):
    return ''.join(COMP[c] for c in reversed(x))
jas = []
for tf, mid in MAT.items():
    d = json.load(open(os.path.join(JD, mid + '.json')))
    cons = consensus(d['pfm']); top = topseq(d['pfm'])
    motif = D.TF_MOTIFS[tf]
    pat = re.compile(''.join('[' + IUPAC[c] + ']' for c in motif))
    hit = pat.search(top) or pat.search(rc(top))
    jas.append({'tf': tf, 'jaspar_id': mid, 'jaspar_name': d['name'], 'jaspar_consensus': cons,
                'jaspar_top_sequence': top, 'module_motif': motif,
                'module_motif_matches_jaspar_top_sequence': bool(hit),
                'match_offset': hit.start() if hit else None})
res['jaspar_check'] = jas

# 2. UCSC CpG islands vs module detector
def cpg_for(seqfile, trackfile, chrom, win_start):
    seq = json.load(open(os.path.join(UD, seqfile)))['dna'].upper()
    track = json.load(open(os.path.join(UD, trackfile)))['cpgIslandExt']
    inside = [{'start': x['chromStart'] - win_start, 'end': x['chromEnd'] - win_start, 'name': x['name']}
              for x in track if x['chromStart'] >= win_start and x['chromEnd'] <= win_start + len(seq)]
    det = D._cpg_islands(seq)
    rows = []
    for t in inside:
        ov = [d for d in det if d['start'] < t['end'] and d['end'] > t['start']]
        rows.append({'track_island': t, 'detected_overlaps': len(ov),
                     'detected_span': [min(o['start'] for o in ov), max(o['end'] for o in ov)] if ov else None})
    fp = [d for d in det if not any(d['start'] < t['end'] and d['end'] > t['start'] for t in inside)]
    return rows, fp, det
ml, ml_fp, _ = cpg_for('seq_mlh1.json', 'cpg_mlh1_wide.json', 'chr3', 36990000)
ck, ck_fp, _ = cpg_for('seq_cdkn2a.json', 'cpg_cdkn2a.json', 'chr9', 21990000)
ds, ds_fp, _ = cpg_for('seq_desert.json', 'cpg_cdkn2a.json', 'chr9', 21976000)
res['cpg_islands'] = {'mlh1_region': ml, 'cdkn2a_region': ck, 'desert_track_islands': ds,
                      'false_positives_mlh1': len(ml_fp), 'false_positives_cdkn2a': len(ck_fp),
                      'false_positives_desert': len(ds_fp),
                      'note': 'module criteria are Gardiner-Garden (GC>0.5, CpG O/E>0.6) on 200 bp windows without UCSC length/merge filtering; boundaries differ by algorithm'}

# 3. motif_binding_energy IUPAC bug demo
def correct_match(window, motif):
    return sum(a in IUPAC[b] for a, b in zip(window, motif)) / len(motif)
demo_seq = 'C' * 20
out_c = D.motif_binding_energy(demo_seq, 'RRR')
res['binding_energy_bug_demo'] = {'sequence': 'C x20', 'motif': 'RRR', 'module_match': out_c['match'],
                                  'correct_iupac_match': correct_match('CCC', 'RRR'),
                                  'note': 'module credits ANY base at degenerate motif positions (b in NRWYKMSVHD), so R matches C/T; match is inflated'}
real = json.load(open(os.path.join(UD, 'seq_mlh1.json')))['dna'].upper()[:2000]
diffs = []
for tf, m in D.TF_MOTIFS.items():
    mod = D.motif_binding_energy(real, m)
    scores = []
    for i in range(len(real) - len(m) + 1):
        w = real[i:i + len(m)]
        mt = correct_match(w, m)
        sh = .5 * (w.count('A') + w.count('T')) / len(m)
        scores.append((-8 * mt + .8 * sh, i, mt))
    e, i_c, mt = min(scores)
    diffs.append({'tf': tf, 'module_energy': mod['energy_kcal_mol'], 'correct_energy': round(e, 3),
                  'module_match': mod['match'], 'correct_match': round(mt, 3),
                  'inflated': mod['energy_kcal_mol'] < e - 1e-9})
res['binding_energy_real'] = diffs

# 4. decode on real sequence: verify hit positions independently
seq = json.load(open(os.path.join(UD, 'seq_mlh1.json')))['dna'].upper()
dec = D.decode(seq)
bad = 0
for h in dec['tf_motif_hits']:
    pat = ''.join('[' + IUPAC[c] + ']' for c in D.TF_MOTIFS[h['tf']])
    if not re.fullmatch(pat, seq[h['position']:h['position'] + len(D.TF_MOTIFS[h['tf']])]):
        bad += 1
res['decode_real'] = {'n_hits': len(dec['tf_motif_hits']), 'hits_failing_independent_regex': bad,
                      'n_clusters': len(dec['enhancer_clusters']), 'n_cpg': len(dec['cpg_islands']),
                      'n_lnc': len(dec['lncrna_candidates'])}

# 5. internal math recomputation
s200 = (seq[:200])
cf = D.coordinate_field(s200, {'atac': .8, 'h3k27ac': .6, 'methylation': .2})
w = s200[:50]
from sugarcode.bio.sequence import gc_content as _gc, clean_dna as _cd
from sugarcode.bio.sequence import find_motif as _fm
gc = _gc(_cd(w)); md = sum(len(_fm(_cd(w), m)) for m in D.TF_MOTIFS.values()) / len(w)
want = .3 * gc + .3 * .8 + .25 * .6 + .1 * (1 - .2) + .05 * min(1, md * 20)
res['coordinate_field'] = {'module': cf['field'][0]['regulatory_potential'], 'recompute': want,
                           'match': abs(cf['field'][0]['regulatory_potential'] - want) < 1e-9}
g = D.chromatin_graph(s200, bin_size=50)
prop = D.graph_propagate(g, {0: 1.0}, steps=3, damping=.6)
n = len(g['nodes']); A = np.zeros((n, n))
for e in g['edges']: A[e['source'], e['target']] = A[e['target'], e['source']] = e['contact']
rows = A.sum(1, keepdims=True); N = np.divide(A, rows, out=np.zeros_like(A), where=rows > 0)
st = np.array([x['regulatory_potential'] for x in g['nodes']]); st[0] += 1.0
init = st.copy()
for _ in range(3): st = (1 - .6) * init + .6 * N @ st
res['graph_propagate'] = {'match': bool(np.allclose(prop['propagated'], st, atol=1e-12)),
                          'edges': len(g['edges']), 'contact_decay': '(1+d)^-0.8, threshold 0.15'}
nuc = D.nucleosome_positioning(seq[:500])
ok = True
for e in nuc:
    w = seq[e['start']:e['end']]
    per = sum(w[j] in 'AT' and w[j + 10] in 'AT' for j in range(len(w) - 10)) / (len(w) - 10)
    ok &= abs(e['occupancy'] - min(1, .25 + .45 * _gc(w) + .3 * per)) < 1e-9
res['nucleosome_formula_match'] = bool(ok)
ca = D.enhancer_gene_causality([1, 2, 3, 4, 5], [2, 4, 5, 8, 11], [0, 0, 1, 1, 1])
res['causality'] = {'module_r': ca['observational_correlation'], 'numpy_r': float(np.corrcoef([1,2,3,4,5],[2,4,5,8,11])[0,1]),
                    'module_effect': ca['interventional_effect'], 'expected_effect': float(np.mean([5,8,11]) - np.mean([2,4]))}
res['logic'] = {'and': D.regulatory_logic([.5, .5])['activity'], 'or': D.regulatory_logic([.5, .5], 'OR')['activity'],
                'not': D.regulatory_logic([.3], 'NOT')['activity']}
tt = D.temporal_trajectory([.5, .5], [[.9, .1], [.2, .8]], steps=5)
stt = np.array([.5, .5])
for _ in range(5): stt = stt @ np.array([[.9, .1], [.2, .8]])
res['temporal_trajectory'] = {'match': bool(np.allclose(tt['final_state'], stt, atol=1e-12))}
pf = D.perturbation_feedback({'mean': 1.0, 'std': .5}, 1.4, .2)
p, q = 1 / .25, 1 / .04
res['perturbation_feedback'] = {'module': pf['mean'], 'expected': (p * 1.0 + q * 1.4) / (p + q),
                                'std_match': abs(pf['std'] - math.sqrt(1 / (p + q))) < 1e-12}
vc = D.variant_counterfactual(seq[:400], 100, 'T', samples=200, seed=7)
vc2 = D.variant_counterfactual(seq[:400], 100, 'T', samples=200, seed=7)
res['variant_counterfactual'] = {'seed_reproducible': vc == vc2, 'ref': vc['ref'], 'alt': vc['alt'],
                                 'delta_rp': vc['cell_states']['generic']['delta_regulatory_potential']}
out = sys.argv[1] if len(sys.argv) > 1 else 'benchmarks/sweep_dark_genome.json'
json.dump(res, open(out, 'w'), indent=1, default=str)
print(json.dumps({'jaspar': [(j['tf'], j['module_motif_matches_jaspar_top_sequence'], j['match_offset']) for j in jas],
                  'cpg_mlh1': [(r['track_island']['name'], r['detected_overlaps']) for r in ml],
                  'cpg_cdkn2a': [(r['track_island']['name'], r['detected_overlaps']) for r in ck],
                  'desert_fp': len(ds_fp), 'mlh1_fp': len(ml_fp), 'cdkn2a_fp': len(ck_fp),
                  'bug_demo': res['binding_energy_bug_demo'],
                  'inflated_tfs': [d['tf'] for d in diffs if d['inflated']],
                  'decode_bad': bad}, indent=1))
