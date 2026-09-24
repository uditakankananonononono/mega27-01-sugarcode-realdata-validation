"""Sweep: liquid_biopsy vs TCGA PanCancer mutation frequencies (cBioPortal), cfDNA fragment-length
literature (Underhill 2016), and independent recomputation (scipy)."""
import json, math, os, sys, urllib.request
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
import numpy as np
from sugarcode.modules.liquid_biopsy import core as L
from sugarcode.bio import entrez
API = 'https://www.cbioportal.org/api'
CD = os.path.join(os.path.dirname(__file__), '..', 'data', 'cbioportal')
res = {'module': 'liquid_biopsy', 'sources': ['cBioPortal REST API (TCGA PanCancer Atlas 2018 studies)', 'NCBI PubMed', 'scipy']}

def post(url, body):
    r = urllib.request.Request(url, data=json.dumps(body).encode(), headers={'Content-Type': 'application/json', 'Accept': 'application/json'})
    return json.load(urllib.request.urlopen(r, timeout=60))
def get(url): return json.load(urllib.request.urlopen(url, timeout=60))

# 1. marker panels vs TCGA mutation frequencies
STUDY = {'colorectal': 'coadread', 'lung': 'luad', 'breast': 'brca', 'pancreatic': 'paad', 'prostate': 'prad', 'melanoma': 'skcm'}
PRE = {'colorectal': ['KRAS', 'APC', 'TP53'], 'lung': ['EGFR', 'KRAS', 'ALK'], 'breast': ['PIK3CA', 'ESR1', 'TP53'],
       'pancreatic': ['KRAS', 'CDKN2A'], 'prostate': ['AR', 'TMPRSS2-ERG'], 'melanoma': ['BRAF', 'NRAS', 'TERT']}
EXTRA = {'lung': ['TP53', 'STK11', 'KEAP1', 'BRAF'], 'colorectal': ['PIK3CA', 'BRAF', 'SMAD4'], 'breast': ['GATA3', 'CDH1', 'MAP3K1'],
         'pancreatic': ['TP53', 'SMAD4'], 'prostate': ['ERG', 'SPOP', 'TP53', 'FOXA1'], 'melanoma': ['NF1', 'TP53']}
fn = os.path.join(CD, 'tcga_panel_mutations_sequenced.json')
if not os.path.exists(fn):
    syms = sorted({g for t in STUDY for g in PRE[t] + EXTRA[t] if '-' not in g})
    eid = {g['hugoGeneSymbol']: g['entrezGeneId'] for g in post(API + '/genes/fetch?geneIdType=HUGO_GENE_SYMBOL', syms)}
    raw = {'entrez': eid, 'studies': {}}
    for t, c in STUDY.items():
        st = f'{c}_tcga_pan_can_atlas_2018'
        genes = [g for g in PRE[t] + EXTRA[t] if g in eid]
        n = get(f'{API}/sample-lists/{st}_sequenced')['sampleCount']
        muts = post(f'{API}/molecular-profiles/{st}_mutations/mutations/fetch?projection=SUMMARY',
                    {'sampleListId': f'{st}_sequenced', 'entrezGeneIds': [eid[g] for g in genes]})
        raw['studies'][t] = {'study': st, 'n_sequenced': n, 'mutations': [[m['entrezGeneId'], m['sampleId'], m.get('mutationType')] for m in muts]}
    json.dump(raw, open(fn, 'w'))
raw = json.load(open(fn)); inv = {v: k for k, v in raw['entrez'].items()}
panels = {}
for t, d in raw['studies'].items():
    n = d['n_sequenced']
    samp = {}
    for g, s, ty in d['mutations']:
        if ty != 'Silent': samp.setdefault(inv[g], set()).add(s)
    freq = {g: round(100 * len(samp.get(g, ())) / n, 1) for g in PRE[t] + EXTRA[t] if '-' not in g}
    def cov(panel): return round(100 * len(set().union(*[samp.get(g, set()) for g in panel])) / n, 1)
    post_panel = L.CTDNA_MARKERS[t]
    top = max(freq, key=freq.get)
    panels[t] = {'study': d['study'], 'n_sequenced': n, 'mutation_pct': freq, 'most_frequent_checked': top,
                 'pre_panel': PRE[t], 'post_panel': post_panel, 'coverage_pre_pct': cov(PRE[t]), 'coverage_post_pct': cov(post_panel),
                 'top_gene_in_post_panel': top in post_panel}
res['marker_panels_vs_tcga'] = panels
res['marker_notes'] = ('TERT (melanoma) is driven by promoter mutations outside exome capture; AR and ESR1 are acquired under therapy '
                       '(mCRPC, endocrine-treated metastatic breast) so low TCGA primary frequencies are expected; TMPRSS2-ERG is a fusion '
                       'and not in mutation calls. TP53 was the single most frequent checked gene in LUAD and PAAD yet absent from the '
                       'lung and pancreatic panels -> curation fix.')
res['literature'] = {k: {'title': v.get('title'), 'source': v.get('source'), 'pubdate': v.get('pubdate')}
                     for k, v in entrez.esummary('pubmed', ['25079552', '28810144', '27428049', '26771485']).items()}

# 2. fragment-length anchors vs Underhill 2016 abstract (PMID 27428049): ctDNA 134-144 bp vs background 167 bp
ab = entrez.pubmed_abstracts(['27428049'])[0]['abstract']
fl = {tf: L.fragment_length_model(tf) for tf in (0.0, 0.1, 0.3, 0.6)}
res['fragment_model'] = {
    'underhill_abstract_has_134_144_and_167': ('134-144 bp' in ab and '167 bp' in ab),
    'module_anchor_healthy_bp': 166, 'module_anchor_ctdna': '134-144',
    'per_tf': {str(tf): {k: f[k] for k in ('modal_length_bp', 'shannon_entropy_bits', 'short_fraction_100_150bp', 'kl_divergence_to_healthy_ref')} for tf, f in fl.items()},
    'entropy_monotone': all(fl[a]['shannon_entropy_bits'] < fl[b]['shannon_entropy_bits'] for a, b in [(0.0, 0.1), (0.1, 0.3), (0.3, 0.6)]),
    'kl_monotone': all(fl[a]['kl_divergence_to_healthy_ref'] < fl[b]['kl_divergence_to_healthy_ref'] for a, b in [(0.0, 0.1), (0.1, 0.3), (0.3, 0.6)]),
    'short_frac_monotone': all(fl[a]['short_fraction_100_150bp'] < fl[b]['short_fraction_100_150bp'] for a, b in [(0.0, 0.1), (0.1, 0.3), (0.3, 0.6)]),
}
from scipy.stats import entropy as H, gamma as Gm
d = np.array(fl[0.3]['histogram']['density'])
res['fragment_model']['entropy_recomputed_scipy'] = round(float(H(d[d > 0], base=2)), 4)
k = 1 + (166 / 18) ** 2; th = 166 / (k - 1)
res['fragment_model']['gamma_mode_check'] = {'mode': round((k - 1) * th, 6), 'sd': round(math.sqrt(k) * th, 3), 'scipy_argmax_bp': int(np.argmax(Gm(k, scale=th).pdf(np.arange(50, 500))) + 50)}

# 3. haplotype EM: BUG 45 regression
fr = [[1, 1, -1]] * 30 + [[-1, 1, 1]] * 30 + [[0, 0, -1]] * 30 + [[-1, 0, 0]] * 30
h = L.bayesian_haplotype_inference(fr)['haplotypes']
res['haplotype_em'] = {'truth': ['111', '000'], 'post_fix_top2': [(x['haplotype'], x['posterior_mean']) for x in h[:2]],
                       'pre_fix': "top3 000 0.493, 011 0.253, 110 0.253 (111 never a candidate; sugarcode c685085, before the fix)"}

# 4. longitudinal slope vs scipy linregress; detect_ctdna threshold identity
from scipy.stats import linregress
samples = [{'time': t, 'burden': b} for t, b in [(0, .02), (28, .011), (56, .006), (84, .004), (112, .009)]]
lt = L.longitudinal_trajectory(samples)
res['longitudinal'] = {'slope_module': lt['slope'], 'slope_scipy': round(linregress([s['time'] for s in samples], [s['burden'] for s in samples]).slope, 8), 'trend': lt['trend']}
rng = np.random.default_rng(3); sig = np.abs(rng.normal(0, .0005, 60)); sig[[10, 30]] = [.02, .008]
dc = L.detect_ctdna(sig, 'lung', depth=30000)
thr = max(3 * 0.001, np.percentile(sig, 90) * 0.5)
res['detect_ctdna'] = {'called_loci': [c['locus_index'] for c in dc['candidates']], 'expected_loci': [int(i) for i in np.where(sig >= thr)[0]],
                       'markers_lung': dc['suggested_biomarkers'],
                       'note': "'beta-binomial' in the docstring is a label: confidence is 1-exp(-af*depth/100), not a beta-binomial test (documented, not changed)"}
json.dump(res, open(os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_liquid_biopsy.json'), 'w'), indent=1, default=str)
print(json.dumps(res, default=str)[:5000])
