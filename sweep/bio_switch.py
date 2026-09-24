"""Sweep: bio_switch vs PubChem descriptors, Hill identities, and internal math."""
import json, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
from sugarcode.modules.bio_switch import core as B

CACHE = os.path.join(os.path.dirname(__file__), '..', 'data', 'pubchem')
res = {'module': 'bio_switch', 'sources': ['https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<name>/property/...']}

# 1. descriptors vs PubChem
NAMES = ['tetracycline', 'lactose', 'arabinose', 'theophylline', 'benzene']
desc = []
for n in NAMES:
    pub = json.load(open(os.path.join(CACHE, n + '.json')))['PropertyTable']['Properties'][0]
    mod = B._CHEM_DEScriptors = B.__dict__.get('_CHEM_DESCRIPTORS')[n]
    desc.append({'analyte': n, 'pubchem_cid': pub['CID'],
                 'mw': {'module': mod['molecular_weight'], 'pubchem': pub['MolecularWeight'], 'match': abs(float(pub['MolecularWeight']) - mod['molecular_weight']) < 0.5},
                 'xlogp': {'module_clogp': mod['clogp'], 'pubchem_xlogp3': pub.get('XLogP'), 'note': 'clogP vs XLogP3 are different methods; |diff|>1 flagged not failed'},
                 'tpsa': {'module': mod['polar_surface_area'], 'pubchem': pub['TPSA'], 'match': abs(pub['TPSA'] - mod['polar_surface_area']) < 10,
                          'note': 'open-chain vs furanose (arabinose) or tautomer definitions differ'},
                 'hbd': {'module': mod['hbond_donors'], 'pubchem': pub['HBondDonorCount'], 'match': mod['hbond_donors'] == pub['HBondDonorCount']},
                 'hba': {'module': mod['hbond_acceptors'], 'pubchem': pub['HBondAcceptorCount'], 'match': mod['hbond_acceptors'] == pub['HBondAcceptorCount'],
                         'note': 'module uses Lipinski N+O count (tetracycline 10 = 8O+2N), PubChem uses refined acceptor definition (9)'}})
res['descriptor_check'] = desc
# ions: atomic weights (CIAAW standard) - fluoride 19.0 (F 18.998), Cu 63.546, As 74.922, Hg 200.59, Cd 112.41
res['ion_weights'] = {k: {'module': B.__dict__['_CHEM_DESCRIPTORS'][k]['molecular_weight'],
                          'standard_atomic_weight': {'fluoride': 18.998, 'copper': 63.546, 'arsenic': 74.922, 'mercury': 200.59, 'cadmium': 112.41}[k]}
                      for k in ['fluoride', 'copper', 'arsenic', 'mercury', 'cadmium']}
res['ion_weights_note'] = 'rounded to 1 decimal in module; Lipinski rules are not meaningful for ionic analytes (caveat)'

# 2. lipinski violation count (bool bug demonstration, pre/post fix)
lip = []
for k, d in B.__dict__['_CHEM_DESCRIPTORS'].items():
    true_n = int(d['molecular_weight'] > 500) + int(d['clogp'] > 5) + int(d['hbond_donors'] > 5) + int(d['hbond_acceptors'] > 10)
    out = B.cheminformatics_descriptors(k)
    expected_class = 'high' if true_n == 0 else ('moderate' if true_n == 1 else 'low')
    lip.append({'analyte': k, 'true_violations': true_n, 'module_field': out['lipinski_violations'],
                'module_field_type': type(out['lipinski_violations']).__name__,
                'expected_class': expected_class, 'module_class': out['bioavailability_class'],
                'correct': out['bioavailability_class'] == expected_class and out['lipinski_violations'] == true_n})
res['lipinski_check'] = lip

# 3. Hill identities
import numpy as np
hs = []
for kd, n in [(1.0, 1.2), (0.3, 3.0), (60.0, 1.2)]:
    for c in [0.01, 1.0, 50.0, 1000.0]:
        got = B.hill_response(c, kd, n)
        want = c ** n / (kd ** n + c ** n)
        hs.append(abs(got - want))
res['hill_response_max_abs_err'] = float(max(hs))
t = B.design_sensor_tunable('arsenic', cooperativity=2.0)
# c10/c90 are scanned on the 41-point log grid (0.1 log10 steps): recompute from the grid
grid = t['external_concentrations_uM']; rr = [B.hill_response(c, t['kd_uM'], 2.0) for c in grid]
gc10 = next(c for c, y in zip(grid, rr) if y >= 0.1); gc90 = next(c for c, y in zip(grid, rr) if y >= 0.9)
res['tunable_c90_c10'] = {'module': t['threshold_ratio'], 'grid_recompute': round(gc90 / gc10, 4),
                          'continuous_theory_81^(1/n)': 81 ** (1 / 2.0),
                          'match_grid': abs(t['threshold_ratio'] - gc90 / gc10) < 1e-3,
                          'note': 'c10/c90 quantized to the 0.1 log10 scan grid; continuous theory is 9.0 for n=2 - discretization caveat, not an error'}
res['tunable_curves_vs_hill'] = bool(np.allclose(t['response'], [c ** 2 / (t['kd_uM'] ** 2 + c ** 2) for c in t['external_concentrations_uM']], atol=1e-4))
d = B.design_biosensor('copper')
res['design_biosensor'] = {'lod': d['limit_of_detection_uM'], 'kd': d['binding_domain']['kd_uM'],
                           'lod_is_kd_over_10': abs(d['limit_of_detection_uM'] - d['binding_domain']['kd_uM'] / 10) < 1e-3,
                           'curve_vs_hill_1.2': bool(np.allclose(d['dose_response']['response'],
                                                    [1 / (1 + (d['binding_domain']['kd_uM'] / c) ** 1.2) for c in d['dose_response']['concentration_uM']], atol=2e-3)),
                           'response_time': d['response_time_min']}

# 4. reaction-diffusion: module exp profile vs finite-domain solutions
rd = B.reaction_diffusion_profile('copper', extracellular_uM=1.0)
lam = math.sqrt(10.0 / 0.1); L = 2.0
xs = rd['positions_um']
semi = [math.exp(-x / lam) for x in xs]
zero_flux = [math.cosh((L - x) / lam) / math.cosh(L / lam) for x in xs]
bathed = [math.cosh((x - L / 2) / lam) / math.cosh(L / 2 / lam) for x in xs]
res['reaction_diffusion'] = {'module_is_semi_infinite': bool(np.allclose(rd['intracellular_concentration_uM'], semi, atol=1e-6)),
                             'tip_to_base_module': rd['tip_to_base_ratio'],
                             'tip_to_base_zero_flux': round(zero_flux[-1], 6), 'tip_to_base_bathed_cell': round(bathed[-1], 6),
                             'note': 'module uses the semi-infinite slab exp(-x/lambda); a finite 2 um cell gives flatter profiles (cosh solutions). No boundary condition declared - documented approximation, not an internal contradiction'}

# 5. gillespie
g = B.gillespie_sensor_noise('arsenic', ligand_uM=1.0, sensor_copies=10)
res['gillespie'] = {'mean_bound': g['mean_bound_fraction'], 'deterministic': g['deterministic_occupancy'],
                    'expected_deterministic': 1.0 / (1.0 + 5.0), 'copies_per_fL_module': g['estimated_ligand_copies_per_fL'],
                    'copies_per_fL_check': round(1.0 * 1e-6 * 6.022e23 * 1e-15, 1),
                    'seed_reproducible': B.gillespie_sensor_noise('arsenic', ligand_uM=1.0)['mean_bound_fraction'] == g['mean_bound_fraction']}
res['gillespie']['copies_formula_ok'] = abs(g['estimated_ligand_copies_per_fL'] - res['gillespie']['copies_per_fL_check']) < 0.5

# 6. gates, filter, feedback, degradation
gate = B.multi_input_gate({'copper': 0.001, 'arsenic': 5.0}, 'AND')
va, vb = 0.001 / (0.001 + 0.001), 5.0 / (5.0 + 5.0)
res['gate_and'] = {'module': gate['gate_output'], 'expected': round(va * vb, 3)}
tf = B.temporal_filter([[0, 1.0], [10, 1.0], [20, 1.0], [30, 1.0], [40, 0], [50, 0]])
res['temporal_filter'] = {'verdict_sustained': tf['verdict'], 'peak': tf['peak_response'],
                          'decay_tail': tf['filtered_response'][-1], 'expected_tail': round(1.0 * math.exp(-10 / 5.0), 5)}
tf2 = B.temporal_filter([[0, 1.0], [1, 0], [30, 0]])
res['temporal_filter_transient'] = {'verdict': tf2['verdict'], 'peak': tf2['peak_response']}
fb = B.feedback_loop('positive', gain=2.0)
n = 3
res['feedback_positive'] = {'module': fb['output_response'], 'expected': [round(x ** n / (0.5 ** n + x ** n), 5) for x in fb['input_signal']], 'verdict': fb['verdict']}
fbn = B.feedback_loop('negative', gain=2.0)
res['feedback_negative'] = {'module': fbn['output_response'], 'expected': [round(x / (1 + 2 * x), 5) for x in fbn['input_signal']]}
dk = B.degradation_kinetics('theophylline', degradation_rate_per_h=0.2)
res['degradation'] = {'half_life': dk['half_life_h'], 'expected': round(math.log(2) / 0.2, 4),
                      'curve_max_err': float(max(abs(a - math.exp(-0.2 * t)) for a, t in zip(dk['fraction_remaining'], dk['time_points_h'])))}
ot = B.off_target_profile('copper')
res['off_target'] = {'margin': ot['specificity_margin'], 'verdict': ot['competing_interactions'],
                     'note': 'cross-reactivity proxy cross_kd = max(kd_self, kd_other) is an order-of-magnitude heuristic with no binding model; documented'}
syn = B.synthesize_biosensor('fluoride', membrane_permeability=0.1)
res['synthesize'] = {'lod_external': syn['limit_of_detection_uM_external'], 'expected': 60.0 / 10 / 0.1,
                     'ultrasensitive_uses_n3': bool(np.allclose(syn['predicted_curves']['ultrasensitive'],
                                                    [(0.1 * c) ** 3 / (60.0 ** 3 + (0.1 * c) ** 3) for c in syn['external_concentrations_uM']], atol=1e-4))}

out = sys.argv[1] if len(sys.argv) > 1 else 'benchmarks/sweep_bio_switch.json'
json.dump(res, open(out, 'w'), indent=1, default=str)
print(json.dumps({'lipinski': [(l['analyte'], l['true_violations'], l['module_field'], l['module_class'], l['correct']) for l in lip],
                  'hill_err': res['hill_response_max_abs_err'], 'c90c10': res['tunable_c90_c10'],
                  'rd': res['reaction_diffusion']['tip_to_base_module'], 'gillespie': res['gillespie'],
                  'hba_mismatches': [d['analyte'] for d in desc if not d['hba']['match']],
                  'tpsa_mismatches': [d['analyte'] for d in desc if not d['tpsa']['match']]}, indent=1))
