"""oncocircuit TUMOR_PROMOTERS table vs Human Protein Atlas v23 pathology.tsv
(403,241 gene-cancer patient-sample rows): claimed cancer associations and
tumor_selectivity values vs observed detection and High-expression specificity."""
import sys, os, json, csv
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from collections import defaultdict
from sugarcode.modules.oncocircuit.core import TUMOR_PROMOTERS, _select_promoters, design_oncocircuit

GENE = {'hTERT': 'TERT', 'Survivin': 'BIRC5', 'AFP': 'AFP', 'PSA': 'KLK3', 'HER2_inducible': 'ERBB2'}
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'hpa', 'pathology_5genes.tsv')

rows = list(csv.DictReader(open(DATA), delimiter='\t'))
stats = defaultdict(dict)  # gene -> cancer -> dict
for r in rows:
    h, m, l, nd = int(r['High']), int(r['Medium']), int(r['Low']), int(r['Not detected'])
    stats[r['Gene name']][r['Cancer']] = {'high': h, 'med': m, 'low': l, 'nd': nd,
                                          'n': h + m + l + nd, 'det': h + m + l}

# map module cancer labels to HPA cancer names
CANCER_MAP = {'pan-cancer': None, 'lung': 'lung cancer', 'breast': 'breast cancer',
              'colorectal': 'colorectal cancer', 'hepatocellular': 'liver cancer',
              'prostate': 'prostate cancer', 'gastric': 'stomach cancer'}

res = {'promoters': {}}
for prom, ann in TUMOR_PROMOTERS.items():
    g = GENE[prom]
    per = stats[g]
    det_rates = {c: s['det'] / s['n'] for c, s in per.items() if s['n'] > 0}
    high_rates = {c: s['high'] / s['n'] for c, s in per.items() if s['n'] > 0}
    n_detected_cancers = sum(1 for v in det_rates.values() if v > 0)
    # specificity: share of all detected samples falling in the claimed cancers
    claimed = [CANCER_MAP[c] for c in ann['cancers'] if CANCER_MAP[c]]
    det_claimed = sum(s['det'] for c, s in per.items() if c in claimed)
    det_total = sum(s['det'] for s in per.values())
    high_claimed = sum(s['high'] for c, s in per.items() if c in claimed)
    high_total = sum(s['high'] for s in per.values())
    res['promoters'][prom] = {
        'gene': g,
        'claimed_cancers': ann['cancers'],
        'claimed_selectivity': ann['tumor_selectivity'],
        'n_cancers_detected': n_detected_cancers,
        'of_cancers': len(per),
        'detection_share_in_claimed': round(det_claimed / det_total, 3) if det_total else None,
        'high_share_in_claimed': round(high_claimed / high_total, 3) if high_total else None,
        'det_rate_claimed': {c: round(det_rates[c], 3) for c in claimed if c in det_rates},
        'high_rate_claimed': {c: round(high_rates[c], 3) for c in claimed if c in high_rates},
        'max_det_rate_any': round(max(det_rates.values()), 3),
        'max_high_cancer': max(high_rates, key=high_rates.get),
    }

# selector behavior on module's own terms
res['select_prostate'] = _select_promoters('prostate')
res['select_liver'] = _select_promoters('hepatocellular')
res['select_unknown_fallback'] = _select_promoters('glioblastoma')
# design sanity: two-input AND should not express payload when only one input active
d2 = design_oncocircuit('prostate', two_input=True)
d1 = design_oncocircuit('prostate', two_input=False)
res['dual_input_circuit'] = {'name': d2['circuit']['name'], 'logic': d2['circuit']['gates'][0]['logic'],
                             'promoters': [i[0] for i in d2['circuit']['gates'][0]['inputs']]}
res['dual_sim'] = d2['simulation']
res['dual_selective'] = d2['simulation']['output_tumor'] > 10 * max(d2['simulation']['output_normal'], d2['simulation']['output_single_positive'])
res['single_sim'] = d1['simulation']
res['single_input_promoters'] = [i[0] for i in d1['circuit']['gates'][0]['inputs']]
res['datasets'] = ['HPA:v23 pathology.tsv (2026-09-25 download)'] + [f'ensembl:{e}' for e in
    ['ENSG00000164362', 'ENSG00000089685', 'ENSG00000081051', 'ENSG00000142515', 'ENSG00000141736']]
res['n_datasets'] = len(res['datasets'])

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks', 'sweep_oncocircuit.json')
json.dump(res, open(out, 'w'), indent=2)
print(json.dumps(res['promoters'], indent=2))
print(res['select_prostate'], res['select_liver'], res['select_unknown_fallback'])
