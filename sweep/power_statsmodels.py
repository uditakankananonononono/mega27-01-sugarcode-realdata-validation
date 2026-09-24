"""gene_analysis.power_estimate vs statsmodels NormalIndPower on a grid."""
import json, math, os, sys, warnings, itertools
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from statsmodels.stats.power import NormalIndPower
from sugarcode.modules.gene_analysis.core import power_estimate
def old(es, v, pw):
    z = 1.96 + (.84 if pw <= .8 else 1.28); return math.ceil(2 * v * z * z / es ** 2)
rows = []
for pw, es, v in itertools.product((0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99), (0.2, 0.5, 1.0), (0.25, 1.0)):
    d = es / math.sqrt(v)
    try: sm = NormalIndPower().solve_power(effect_size=d, power=pw, alpha=0.05)
    except Exception: sm = float('nan')
    if sm != sm: continue
    rows.append({'power': pw, 'effect': es, 'variance': v, 'statsmodels_n': math.ceil(sm), 'new': power_estimate(es, v, pw)['replicates_per_group'], 'old': old(es, v, pw)})
res = {'tool': 'statsmodels NormalIndPower', 'n_cases': len(rows), 'new_equal': sum(r['new'] == r['statsmodels_n'] for r in rows),
       'old_equal': sum(r['old'] == r['statsmodels_n'] for r in rows), 'old_max_ratio': max(r['old'] / r['statsmodels_n'] for r in rows),
       'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else '', 'rows': rows}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_power_statsmodels.json'), 'w'), indent=1)
print({k: v for k, v in res.items() if k != 'rows'})
