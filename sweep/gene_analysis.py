"""Sweep: gene_analysis vs live NCBI Gene/UniProt/RefSeq/PubMed + Biopython + statsmodels."""
import json, math, os, sys, time, urllib.request, urllib.parse, urllib.error
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
from sugarcode.modules.gene_analysis import core as G
from sugarcode.bio import entrez

CD = os.path.join(os.path.dirname(__file__), '..', 'data', 'gene_analysis')
res = {'module': 'gene_analysis', 'sources': ['NCBI Gene/nuccore/PubMed E-utilities', 'UniProt REST', 'Reactome', 'Biopython ProtParam', 'statsmodels NormalIndPower']}
# Ground truth (NCBI Gene / UniProtKB Swiss-Prot canonical)
TRUTH = {'TP53': (7157, '17', 'P04637', 393), 'BRCA1': (672, '17', 'P38398', 1863),
         'CFTR': (1080, '7', 'P13569', 1480), 'HBB': (3043, '11', 'P68871', 147),
         'EGFR': (1956, '7', 'P00533', 1210)}
rows = []
for sym, (uid, chrom, acc, L) in TRUTH.items():
    p = G.live_gene_profile(sym)
    ng = p['sources'].get('ncbi_gene', {}); up = p['sources'].get('uniprot', {})
    seq = up.get('sequence', '')
    hyd = round(sum(1 for a in seq if a in 'AILMFWVY') / len(seq), 3) if seq else None
    rows.append({'symbol': sym, 'uid': ng.get('uid'), 'uid_ok': str(ng.get('uid')) == str(uid),
                 'chrom': ng.get('chromosome'), 'chrom_ok': ng.get('chromosome') == chrom,
                 'accession': up.get('accession'), 'acc_ok': up.get('accession') == acc,
                 'protein_length': (p.get('protein_stats') or {}).get('length'),
                 'len_ok': (p.get('protein_stats') or {}).get('length') == L,
                 'hydrophobic_fraction_recomputed_ok': hyd == (p.get('protein_stats') or {}).get('hydrophobic_fraction'),
                 'clinvar_n': p['sources'].get('clinvar', {}).get('n_variants'),
                 'reactome_n': len((p['sources'].get('reactome') or {}).get('pathways') or []),
                 'errors': {k: v['error'] for k, v in p['sources'].items() if isinstance(v, dict) and 'error' in v}})
    time.sleep(0.4)
res['live_gene_profile'] = rows

# 2. gene_profile on real RefSeq mRNAs: longest ORF must equal the annotated CDS protein
from Bio.SeqUtils.ProtParam import ProteinAnalysis
mr = []
for sym, nm, L in [('HBB', 'NM_000518.5', 147), ('TP53', 'NM_000546.6', 393), ('CFTR', 'NM_000492.4', 1480)]:
    fn = os.path.join(CD, nm + '.fa')
    if not os.path.exists(fn):
        open(fn, 'w').write(entrez.efetch_fasta('nuccore', nm)); time.sleep(0.4)
    s = ''.join(l.strip() for l in open(fn) if not l.startswith('>'))
    gp = G.gene_profile(sym, s)
    from sugarcode.bio.sequence import orfs
    o = max(orfs(s, min_aa=30), key=lambda x: x['aa_length'])
    prot = o['protein'].rstrip('*')
    bp = ProteinAnalysis(prot).molecular_weight()
    mr.append({'symbol': sym, 'refseq': nm, 'mrna_len': len(s), 'orf_aa': gp['protein']['aa_length'],
               'orf_ok': gp['protein']['aa_length'] == L, 'mw_module': gp['protein']['molecular_weight_da'],
               'mw_biopython': round(bp, 1), 'mw_rel_err': abs(gp['protein']['molecular_weight_da'] - bp) / bp,
               'strand': gp['protein']['strand'], 'n_crispr': len(gp['crispr_targets']), 'gc': gp['locus']['gc_content']})
res['refseq_profiles'] = mr

# 3. publication trend vs independent esearch
time.sleep(2)
pt = G.live_publication_trend('HBB', years=3)
time.sleep(2)
ind = {}
for y in pt['counts_by_year']:
    u = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=0&term='
         + urllib.parse.quote(f'HBB[Title/Abstract] AND {y}[dp]'))
    for k in range(5):
        try:
            ind[y] = int(json.load(urllib.request.urlopen(u))['esearchresult']['count']); break
        except urllib.error.HTTPError:
            time.sleep(2 * (k + 1))
    time.sleep(1.0)
res['pub_trend'] = {'module': pt['counts_by_year'], 'independent': ind, 'match': {str(k): pt['counts_by_year'][k] == v for k, v in ind.items()},
                    'trend': pt['trend'], 'partial': pt['partial_current_year']}

# 4. power vs statsmodels
from statsmodels.stats.power import NormalIndPower
pw = []
for d, var, pwr, a in [(0.5, .25, .8, .05), (0.3, 1, .9, .05), (1.0, 1, .5, .01), (0.2, .25, .95, .05)]:
    m = G.power_estimate(d, var, pwr, a)['replicates_per_group']
    sm = NormalIndPower().solve_power(effect_size=d / math.sqrt(var), power=pwr, alpha=a, ratio=1)
    pw.append({'delta': d, 'var': var, 'power': pwr, 'alpha': a, 'module': m, 'statsmodels_ceil': math.ceil(sm), 'ok': m == math.ceil(sm)})
res['power'] = pw

# 5. closed-form identities
fb = G.feedback_update(0.3, 7, 10); a, b = 1 + 2.4 + 7, 1 + 5.6 + 3
from scipy.stats import beta as B
res['feedback'] = {'mean_ok': abs(fb['mean'] - B(a, b).mean()) < 1e-12, 'std_ok': abs(fb['std'] - B(a, b).std()) < 1e-12}
ce = G.conformational_ensemble(1.0)
RT = 8.314462618 * 298.15 / 4184
res['ensemble'] = {'RT_module': 0.593, 'RT_298K': round(RT, 5), 'pop_sum': sum(ce['populations'].values())}

# 6. HGVS c. hotspot mapping on real RefSeq transcripts (BUG 44 regression on real data)
from sugarcode.bio.sequence import translate
HOT = [('NM_000518.5', 'NM_000518.5(HBB):c.20A>T', 'A', 'E'), ('NM_000518.5', 'c.118C>T', 'C', 'Q'),
       ('NM_000546.6', 'NM_000546.6(TP53):c.524G>A', 'G', 'R'), ('NM_000546.6', 'c.743G>A', 'G', 'R'),
       ('NM_000546.6', 'c.817C>T', 'C', 'R'), ('NM_000492.4', 'c.1652G>A', 'G', 'G')]
hs = []
for nm, v, ref, aa in HOT:
    s = ''.join(l.strip() for l in open(os.path.join(CD, nm + '.fa')) if not l.startswith('>'))
    f = G._variant_features({'variant': v, 'consequence': 'missense'}, s)
    got_aa = translate(f['codon']) if len(f['codon']) == 3 else ''
    hs.append({'refseq': nm, 'variant': v, 'locus_pos': f['position'], 'ref_base': s[f['position'] - 1],
               'ref_ok': s[f['position'] - 1] == ref, 'codon': f['codon'], 'aa': got_aa, 'aa_ok': got_aa == aa})
res['hgvs_hotspots'] = hs
res['hgvs_pre_fix'] = {'note': 'before BUG 44 fix (sugarcode 1ad2a4e): NM_000518.5(HBB):c.20A>T -> position 518, codon TTC; c.20A>T -> position 20, codon ACT (5prime UTR ignored); expected codon GAG (Glu7, sickle).'}
json.dump(res, open(os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_gene_analysis.json'), 'w'), indent=1, default=str)
print(json.dumps({k: res[k] for k in ['refseq_profiles', 'pub_trend', 'hgvs_hotspots']}, default=str)[:6000])
