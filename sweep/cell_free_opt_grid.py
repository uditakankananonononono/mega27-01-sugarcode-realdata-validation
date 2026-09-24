"""Validate cell_free_opt: BUG 35 grid alignment + kinetics/cost/stats vs independent recomputation."""
import json, math, os, sys
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from sugarcode.modules.cell_free_opt.core import (optimize_cfps, _yield_model, kinetics,
    cost_model, batch_normalize, resource_sensitivity, pareto_conditions, replicate_qc, REAGENTS)
res = {'module': 'cell_free_opt'}
# 1) BUG 35: pre-fix default grid missed tpl optimum; post-fix grid hits it for all steps
rows = []
for step in (1, 2, 3, 4):
    b = optimize_cfps(step)
    rows.append({'grid_step': step, 'tpl': b['template_ng_ul'], 'mg': b['mg_mm'], 'yield': b['yield_g_l'],
                 'evaluated': b['conditions_evaluated'],
                 'hits_model_max': abs(b['yield_g_l'] - 2.0) < 1e-9})
res['grid_after_fix'] = rows
res['bug35_prefix'] = {'grid_step': 2, 'tpl_returned': 10, 'yield': 1.879, 'model_max': 2.0}
# 2) grid optimum == brute-force argmax over the exact same product grid (independent)
step = 3
b = optimize_cfps(step)
mg_v = sorted(set(range(4, 21, step)) | {12}); k_v = sorted(set(range(40, 201, 20*step)) | {120})
pep_v = sorted(set(range(0, 41, 2*step)) | {20}); tpl_v = sorted(set(range(2, 31, 2*step)) | {12})
peg_v = sorted(set(range(0, 5, 1)) | {2})
import itertools
best = max(itertools.product(mg_v, k_v, pep_v, tpl_v, peg_v), key=lambda c: _yield_model(*c))
res['grid_argmax_independent'] = {'module': [b['mg_mm'], b['k_mm'], b['pep_mm'], b['template_ng_ul'], b['peg_pct']],
    'independent': list(best), 'match': [b['mg_mm'], b['k_mm'], b['pep_mm'], b['template_ng_ul'], b['peg_pct']] == list(best)}
# 3) kinetics vs analytic cumulative vmax*tau*(1-exp(-t/tau))
kin = kinetics({'yield_g_l': 2.0}, 6.0, 0.25)
vmax, tau = 1.0, 2.0
devs = [abs(p['cumulative_g_l'] - vmax*tau*(1-math.exp(-p['t_h']/tau))) for p in kin['trajectory']]
res['kinetics'] = {'max_euler_deviation_vs_analytic': round(max(devs), 4),
    'euler_expected_bias': 'left-Riemann overestimate on decreasing rate, O(vmax*dt/2) ~ 0.125 g/L',
    'plateau_time_h': kin['plateau_time_h'], 'final_g_l': kin['final_g_l'],
    'analytic_final': round(vmax*tau*(1-math.exp(-6/tau)), 3)}
# 4) cost_model independent recomputation
cond = {'mg_mm': 12, 'k_mm': 120, 'pep_mm': 20, 'template_ng_ul': 12, 'peg_pct': 2}
cm = cost_model(cond, 15.0, 96)
ind = (12*.05*15/1000 + 120*.02*15/1000 + 20*.40*15/1000 + 12*.001*15 + .35)
res['cost_model'] = {'module_per_rxn': cm['per_reaction_usd'], 'independent': round(ind, 4),
    'match': abs(cm['per_reaction_usd'] - ind) < 1e-3, 'plate': cm['plate_usd'],
    'plate_match': abs(cm['plate_usd'] - round(ind*96, 2)) < 0.01}
# 5) stats helpers
bn = batch_normalize([2.0, 4.0, 1.0], [1.0, 2.0, 0.5])
import statistics
norm = [2.0, 2.0, 2.0]
res['batch_normalize'] = {'normalized': bn['normalized_yield'], 'expected': norm,
    'match': bn['normalized_yield'] == norm, 'cv_zero': abs(bn['cv']) < 1e-12}
rq = replicate_qc([1.0, 1.1, 0.9])
m = statistics.mean([1.0, 1.1, 0.9]); cv = statistics.pstdev([1.0, 1.1, 0.9])/m
res['replicate_qc'] = {'cv_match': abs(rq['cv'] - cv) < 1e-12, 'passes': rq['passes'] == (cv <= .15)}
rs = resource_sensitivity(cond, .1)
base = _yield_model(12, 120, 20, 12, 2)
ok = True
for key, dims in [('mg_mm', (12, 5)), ('k_mm', (120, 50)), ('pep_mm', (20, 12)), ('template_ng_ul', (12, 8)), ('peg_pct', (2, 1.5))]:
    c = dict(cond); c[key] *= 1.1
    exp_el = (_yield_model(c['mg_mm'], c['k_mm'], c['pep_mm'], c['template_ng_ul'], c['peg_pct']) - base)/(base*.1)
    ok = ok and abs(rs['elasticity'][key] - exp_el) < 1e-9
res['resource_sensitivity'] = {'elasticities_match': ok, 'most_sensitive': rs['most_sensitive']}
cands = [{'yield_g_l': 2.0, 'cost_usd': 1.0}, {'yield_g_l': 1.5, 'cost_usd': 0.5},
         {'yield_g_l': 1.0, 'cost_usd': 2.0}, {'yield_g_l': 2.0, 'cost_usd': 1.0}]
pa = pareto_conditions(cands)
# exact duplicates are mutually non-dominating (no strict inequality) - both survive: correct non-dominated sort
res['pareto'] = {'kept': pa, 'expected_n': 3, 'duplicates_kept_by_design': True,
    'correct': len(pa) == 3 and pa[0]['cost_usd'] == 0.5 and all(
        not any(o['yield_g_l'] > c['yield_g_l'] and o['cost_usd'] <= c['cost_usd'] or
                o['yield_g_l'] >= c['yield_g_l'] and o['cost_usd'] < c['cost_usd'] for o in cands) for c in pa)}
json.dump(res, open('benchmarks/sweep_cell_free_opt.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
