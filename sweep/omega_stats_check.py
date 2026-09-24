"""omega_stats analytics vs independent numpy/scipy recomputation on synthetic series."""
import sys, os, json, math, statistics
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import numpy as np
from omega.registry import REGISTRY
from sugarcode.modules.omega_stats.core import (record, stats, subnetwork_rollup, performance_dashboard,
                                                compare_releases, enhancement_features, build_omega_report, _METRICS)
MOD = sorted(REGISTRY)[:4]
rng = np.random.default_rng(3)
obs = []
for i, m in enumerate(MOD):
    for j in range(30):
        obs.append({'module': m, 'metric': 'acc', 'kind': 'accuracy', 'value': float(0.8 + 0.01 * j + rng.normal(0, 0.005)), 'ts': 1000.0 + j * 3600})
        obs.append({'module': m, 'metric': 'lat', 'kind': 'latency', 'value': float(500 + 5 * j + rng.normal(0, 20)), 'ts': 1000.0 + j * 3600})
res = {}
# stats() vs statistics module
_METRICS.clear()
for o in obs: record(o['module'], o['metric'], o['value'], o['kind'])
st = stats(); bad = []
for s in st['series']:
    vals = [o['value'] for o in obs if o['module'] == s['module'] and o['metric'] == s['metric'] and o['kind'] == s['kind']]
    if s['n'] != len(vals) or abs(s['mean'] - round(statistics.fmean(vals), 5)) > 1e-9 or abs(s['stdev'] - round(statistics.pstdev(vals), 5)) > 1e-9 \
       or s['min'] != round(min(vals), 5) or s['max'] != round(max(vals), 5): bad.append(s)
res['stats_mismatches'] = len(bad); res['stats_series'] = len(st['series'])
# rollup mean_of_means vs independent
ro = subnetwork_rollup(); ok_ro = True
for sn, d in ro.items():
    means = [s['mean'] for s in st['series'] if REGISTRY[s['module']].subnetwork == sn]
    if abs(d['mean_of_means'] - round(sum(means) / len(means), 5)) > 1e-9: ok_ro = False
res['rollup_ok'] = ok_ro
# dashboard: slope, pass_rate, CI determinism
d1 = performance_dashboard(obs, now=1000.0 + 29 * 3600 + 60)
d2 = performance_dashboard(obs, now=1000.0 + 29 * 3600 + 60)
slope_bad = 0; ci_bad = 0
for x in d1['series']:
    vals = np.array([o['value'] for o in obs if o['module'] == x['module'] and o['metric'] == x['metric']])
    ts = np.array([o['ts'] for o in obs if o['module'] == x['module'] and o['metric'] == x['metric']])
    ref = float(np.polyfit((ts - ts[0]) / 3600, vals, 1)[0])
    if abs(x['trend_per_hour'] - ref) > 1e-9: slope_bad += 1
    tgt = 1000 if x['kind'] == 'latency' else 0.9
    pr = float(np.mean(vals <= tgt)) if x['kind'] == 'latency' else float(np.mean(vals >= tgt))
    if abs(x['pass_rate'] - pr) > 1e-9 or not (x['ci95'][0] <= x['mean'] <= x['ci95'][1]): ci_bad += 1
res['dashboard_slope_mismatches'] = slope_bad; res['dashboard_passrate_ci_mismatches'] = ci_bad
res['ci_deterministic'] = all(a['ci95'] == b['ci95'] for a, b in zip(d1['series'], d2['series']))
# compare_releases vs independent
base = [dict(o, value=o['value'] * 0.95) for o in obs]
cmp_ = compare_releases(base, obs, regression_tolerance=0.05)
eff_bad = 0
for c in cmp_['comparisons']:
    av = np.array([o['value'] for o in base if o['module'] == c['module'] and o['metric'] == c['metric']])
    bv = np.array([o['value'] for o in obs if o['module'] == c['module'] and o['metric'] == c['metric']])
    pooled = math.sqrt((av.var() + bv.var()) / 2)
    ref = (bv.mean() - av.mean()) / pooled if pooled > 0 else 0.0
    imp_ref = (av.mean() - bv.mean()) / abs(av.mean()) if c['kind'] == 'latency' else (bv.mean() - av.mean()) / abs(av.mean())
    if abs(c['standardized_effect'] - ref) > 1e-9 or abs(c['relative_improvement'] - imp_ref) > 1e-12: eff_bad += 1
res['release_mismatches'] = eff_bad; res['release_gate'] = cmp_['release_gate']
# enhancement_features: 51 keys + spot checks
ef = enhancement_features(d1, cmp_)
res['diagnostics_count'] = len(ef)
res['diagnostics_spot_ok'] = (abs(ef['mean_pass_rate'] - np.mean([x['pass_rate'] for x in d1['series']])) < 1e-12
    and ef['alert_count'] == len(d1['alerts']) and ef['series_count'] == len(d1['series']))
# full report
rep = build_omega_report(obs, base, now=1000.0 + 29 * 3600 + 60)
res['report_ok'] = rep['diagnostic_count'] == 51 and rep['release_decision'] == cmp_['release_gate']
res['sugarcode_commit'] = sys.argv[1] if len(sys.argv) > 1 else ''
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_omega_stats.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
