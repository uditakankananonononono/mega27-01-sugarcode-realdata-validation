"""Sweep: phage_designer recomputation + adsorption time-dependence check."""
import json, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
from sugarcode.modules.phage_designer import core as P
res = {'module': 'phage_designer', 'sources': ['recomputation', 'first-order adsorption literature (Schlesinger 1932)']}
# 1. design_fiber recompute
for rec in ('LamB', 'OmpC', 'FhuA'):
    d = P.design_fiber(rec)
    exp_ads = round(min(0.99, 0.55 + 12.0 / d['predicted_binding_nM']), 3)
    res.setdefault('fibers', {})[rec] = {'scaffold': d['scaffold'], 'adsorption': d['predicted_adsorption_rate'],
                                         'adsorption_recomputed': exp_ads, 'ok': d['predicted_adsorption_rate'] == exp_ads,
                                         'off_target': d['off_target_risk']}
# 2. design_lysin recompute
ly = P.design_lysin('gram-', ['amidase', 'glycosidase'])
exp_pot = round((0.8 * 0.9 + 0.7 * 0.9) / 2, 3)
res['lysin'] = {'potency': ly['predicted_potency'], 'potency_recomputed': exp_pot,
                'mic_module': ly['MIC_estimate_ug_mL'], 'mic_recomputed': round(max(0.1, 8.0 * (1.1 - exp_pot)), 2),
                'om_note': ly['outer_membrane_strategy'][:60]}
# 3. adsorption_kinetics: time parameter
k1 = P.adsorption_kinetics(1e6, 1e8, 0.5, time_min=1)
k2 = P.adsorption_kinetics(1e6, 1e8, 0.5, time_min=60)
k1s = P.adsorption_kinetics(1e6, 100, 1e-3, time_min=10)
k2s = P.adsorption_kinetics(1e6, 100, 1e-3, time_min=60)
res['adsorption_time'] = {'pre_fix_bound_t1': 999999.9800000004, 'pre_fix_bound_t60': 999999.9800000004,
                          'pre_fix': 'bound identical at t=1 and t=60 (time_min ignored); sugarcode 5b4e870 parent',
                          'post_fix_small_k': {'bound_t10': k1s['bound'], 'bound_t60': k2s['bound'], 'grows': k2s['bound'] > k1s['bound']},
                          'post_fix_first_order_ok': abs(P.adsorption_kinetics(1000, 50, 0.01, time_min=20)['adsorbed_fraction'] - (1 - math.exp(-0.01*50*20))) < 1e-12,
                          'literature': 'phage adsorption is first-order: bound = P0*(1-exp(-k*B*t)) (Schlesinger 1932)'}
# 4. escape_probability recompute
e = P.escape_probability(1e-8, 2, 1e9)
res['escape'] = {'module': e, 'expected_prob': 1 - math.exp(-1e-8 * 2 * 1e9), 'expected_variants': 1e-8 * 2 * 1e9}
# 5. host_range + cocktail recompute
prof = {'E coli K12': {'LamB': 1.0, 'OmpC': 0.8}, 'E coli B': {'LamB': 0.2, 'OmpC': 0.9}}
hr = P.host_range(prof, P.design_fiber('LamB'))
res['host_range'] = hr
fn = os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_phage_designer.json')
json.dump(res, open(fn, 'w'), indent=1, default=str)
print(json.dumps(res, default=str)[:2500])
