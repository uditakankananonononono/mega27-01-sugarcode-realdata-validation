"""promoter_lib Hill occupancy / burst noise and stability_ai Wright-Fisher drift vs independent references."""
import sys, os, json, math
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import numpy as np
from scipy.special import expit
from sugarcode.modules.promoter_lib.core import thermodynamic_occupancy, expression_noise
from sugarcode.modules.stability_ai.core import stochastic_drift

# Hill occupancy vs independent expit form for h=1 and direct power form for h!=1
rows = []; worst = 0.0
for h in (0.5, 1.0, 1.5, 2.0, 3.0):
    for kd in (0.1, 1.0, 10.0):
        for c in (0.01, 0.1, 1.0, 10.0, 100.0):
            got = thermodynamic_occupancy(c, kd, h)
            ref = 1.0 / (1.0 + (kd / c) ** h)
            worst = max(worst, abs(got - ref))
            rows.append((h, kd, c, got, ref))
half = all(abs(thermodynamic_occupancy(kd, kd, h) - 0.5) < 1e-12 for h in (0.5, 1, 2, 3) for kd in (0.1, 1, 10))
# expit cross-check at h=1: occ = expit(ln(c/kd))
expit_worst = max(abs(thermodynamic_occupancy(c, kd, 1) - expit(math.log(c / kd))) for kd in (0.3, 1, 7) for c in (0.05, 0.5, 5, 50))
# burst noise vs theory: var = mean*(1+B), Fano = 1+B, cv = sqrt(mean(1+B))/mean
noise = []
for mean, B in ((0.5, 1), (2, 5), (10, 3), (40, 8)):
    r = expression_noise(mean, B)
    noise.append({'mean': mean, 'burst': B, 'fano_ok': abs(r['fano'] - (1 + B)) < 1e-9,
                  'var_ok': abs(r['variance'] - mean * (1 + B)) < 1e-9,
                  'cv_ok': abs(r['cv'] - math.sqrt(mean * (1 + B)) / mean) < 1e-9})
# Wright-Fisher drift vs independent deterministic recursion (large population => near-deterministic)
def det_traj(L, m, s, g):
    f = 1.0; out = []
    for gen in range(g + 1):
        out.append(f)
        a = f * (1 - L - m)
        f = a * (1 - s) / (1 - a * s)
    return out
constructs = [
    {'mode': 'plasmid', 'burden': 0.1, 'size_kb': 5.0, 'toxic': False},
    {'mode': 'plasmid', 'burden': 0.8, 'size_kb': 8.0, 'toxic': True},
    {'mode': 'genomic', 'burden': 0.3, 'size_kb': 3.0, 'toxic': False},
]
wf = []
for c in constructs:
    r = stochastic_drift(c, generations=200, population=2_000_000, replicates=6, seed=11)
    L, m, s = r['parameters']['loss_rate'], r['parameters']['mutation_rate'], r['parameters']['selection_against_function']
    ref = det_traj(L, m, s, 200)
    got = {t['generation']: t['mean_functional_fraction'] for t in r['trajectory']}
    diffs = [abs(got[g] - ref[g]) for g in got]
    wf.append({'construct': c, 'max_abs_diff_vs_deterministic': round(max(diffs), 5),
               'terminal_mean': round(r['terminal_mean'], 5), 'deterministic_terminal': round(ref[200], 5),
               'extinction_probability': round(r['extinction_probability'], 4)})
# seed reproducibility
a = stochastic_drift(constructs[0], generations=50, population=1000, replicates=8, seed=5)['terminal_fractions']
b = stochastic_drift(constructs[0], generations=50, population=1000, replicates=8, seed=5)['terminal_fractions']
res = {'reference': 'Hill occupancy 1/(1+(Kd/c)^h) + expit cross-check; Poisson-burst noise theory; independent deterministic Wright-Fisher recursion',
 'hill_max_abs_diff': worst, 'hill_at_kd_is_half': half, 'hill_expit_max_diff': expit_worst,
 'noise_all_ok': all(all(v for k, v in x.items() if k.endswith('_ok')) for x in noise), 'noise_rows': noise,
 'wf_rows': wf, 'wf_seed_reproducible': a == b,
 'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else ''}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_stability_promoter.json'), 'w'), indent=1)
print(json.dumps(res, default=str)[:1200])
