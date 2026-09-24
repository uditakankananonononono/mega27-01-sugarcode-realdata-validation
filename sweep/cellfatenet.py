"""Sweep: cellfatenet vs HGNC + reprogramming literature + internal dynamics math."""
import json, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
import numpy as np
from scipy.integrate import solve_ivp
from sugarcode.modules.cellfatenet import core as C

HD = os.path.join(os.path.dirname(__file__), '..', 'data', 'hgnc')
PD = os.path.join(os.path.dirname(__file__), '..', 'data', 'pmid')
res = {'module': 'cellfatenet', 'sources': ['https://rest.genenames.org/fetch/symbol/<sym>', 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed']}

# 1. HGNC gene symbol check
genes = sorted(set(C.LINEAGE_GRN) | {t for d in C.LINEAGE_GRN.values() for t in d['targets']} | {m for v in C.CELL_MARKERS.values() for m in v})
hg = []
for g in genes:
    d = json.load(open(os.path.join(HD, g + '.json')))
    docs = d['response']['docs']
    hg.append({'gene': g, 'hgnc_status': docs[0]['status'] if docs else 'not an approved symbol',
               'approved_symbol': docs[0]['symbol'] if docs else None})
oct4 = json.load(open(os.path.join(HD, 'OCT4_alias.json')))['response']['docs']
res['hgnc_check'] = hg
res['hgnc_unapproved'] = [h['gene'] for h in hg if h['hgnc_status'] != 'Approved']
res['oct4_alias'] = {'approved_symbol': oct4[0]['symbol'], 'aliases': oct4[0].get('alias_symbol'),
                     'note': 'OCT4 is an alias of POU5F1; module uses the common alias'}

# 2. reprogramming cocktails vs GRN
abs169 = open(os.path.join(PD, '16904174.txt')).read()
abs206 = open(os.path.join(PD, '20691899.txt')).read()
abs242 = open(os.path.join(PD, '24243019.txt')).read()
res['literature_anchors'] = [
 {'pmid': '16904174', 'cocktail': 'OSKM (Oct3/4, Sox2, Klf4, c-Myc)', 'in_grn': {'OCT4': True, 'SOX2': True, 'KLF4': False, 'MYC': False},
  'verified_terms': [t for t in ['Oct3/4', 'Sox2', 'Klf4', 'c-Myc'] if t in abs169],
  'note': 'KLF4 and MYC are absent from the GRN - coverage caveat'},
 {'pmid': '20691899', 'cocktail': 'GMT (Gata4, Mef2c, Tbx5) fibroblast -> cardiomyocyte, Ieda 2010 Cell',
  'in_grn': {'GATA4': True, 'MEF2C': True, 'TBX5': True}, 'verified_terms': [t for t in ['Gata4', 'Mef2c', 'Tbx5'] if t in abs206],
  'note': 'all three GMT factors are GRN nodes and drive the cardiomyocyte markers TNNT2/MYH6 through MEF2C/TBX5 edges'},
 {'pmid': '24243019', 'cocktail': 'BAM (Brn2, Ascl1, Myt1l) fibroblast -> neuron, Wapinski 2013',
  'in_grn': {'POU3F2/Brn2': False, 'ASCL1': True, 'MYT1L': True}, 'verified_terms': [t for t in ['Ascl1', 'Brn2', 'Myt1l'] if t in abs242],
  'note': 'Brn2/POU3F2 absent from the GRN - coverage caveat'}]

# 3. math recomputation
nodes, W = C.grn_matrix()
W2 = np.zeros_like(W); idx = {n: i for i, n in enumerate(nodes)}
for src, d in C.LINEAGE_GRN.items():
    for dst, sg in d['targets'].items():
        W2[idx[dst], idx[src]] = 1 if sg == '+' else -1
res['grn_matrix_match'] = bool(np.array_equal(W, W2))
cb = C.chromatin_binding([0.5, 1.0], [0.0, 0.5], 2.0)
R, T = .001987, 310.15
ka = math.exp(7.0 / (R * T))
eff = [2.0 * .5 * 1.0, 2.0 * 1.0 * .5]
occ = [ka * e / (1 + ka * e) for e in eff]
res['chromatin_binding'] = {'ka_module': cb['association_constant'], 'ka_recompute': ka,
                            'occupancy_match': bool(np.allclose(cb['occupancy'], occ, atol=1e-9))}
init = {n: 0.05 for n in nodes}; init['OCT4'] = 0.4; init['GATA4'] = 0.1
sim = C.simulate_fate(init, hours=48, controls={'ASCL1': 0.3}, methylation={'MYOD1': 0.5})
gate = np.array([0.5 if n == 'MYOD1' else 1.0 for n in nodes])
u = np.array([0.3 if n == 'ASCL1' else 0.0 for n in nodes])
x0 = np.array([init.get(x, 0.05) for x in nodes])
def rhs(_, x):
    signal = W @ (x * x / (.25 + x * x)); return gate * (1 / (1 + np.exp(-4 * (signal - .5)))) + u - .35 * x
def rhs_free(_, x):
    signal = W @ (x * x / (.25 + x * x)); return 1 / (1 + np.exp(-4 * (signal - .5))) - .35 * x
t = np.linspace(0, 48, 289)
y = solve_ivp(rhs, (0, 48), x0, t_eval=t, rtol=1e-7, atol=1e-9).y.T
res['simulate_fate_deterministic_max_abs_diff'] = float(np.max(np.abs(np.asarray(sim['mean']) - y)))
sim2 = C.simulate_fate(init, hours=48, controls={'ASCL1': 0.3}, methylation={'MYOD1': 0.5}, noise=0.02, seed=3, trajectories=48)
tr = np.asarray(sim2['trajectories'])
res['noise'] = {'nonnegative': bool((tr >= 0).all()), 'mean_close_to_deterministic_max_diff': float(np.max(np.abs(tr.mean(0) - y))),
                'trajectory_variance_positive': bool((tr.var(0)[1:] > 0).any())}
att = C.identify_attractors(starts=12, hours=150, seed=4)
stab = []
for a in att['attractors']:
    x = np.array([a['state'][n] for n in nodes])
    stab.append(float(np.max(np.abs(rhs_free(0, x)))))
res['basal_floor'] = {'tnnt2_no_regulator_steady_state': float(C.simulate_fate({n: 0.05 for n in nodes}, hours=72)['final_state']['TNNT2']),
                      'note': 'pre-fix sigmoid(2*signal) gave a 0.5 production floor -> TNNT2 steady state 1.43 with zero regulators; post-fix sigmoid(4*(signal-0.5))'}
res['attractors'] = {'count': len(att['attractors']), 'max_rhs_at_centers': max(stab),
                     'basin_fractions_sum': sum(a['basin_fraction'] for a in att['attractors']),
                     'note': 'rhs with zero control and unit gate at attractor centers (identify_attractors uses no controls/methylation)'}
gmp = C.graph_message_passing({n: 0.5 for n in nodes}, layers=3)
h = np.full(len(nodes), 0.5); deg = np.maximum(np.abs(W).sum(1), 1)
for _ in range(3): h = np.tanh(h + (W @ h) / deg)
res['message_passing_match'] = bool(np.allclose(gmp['embedding'], h, atol=1e-12))

# 4. end-to-end fibroblast -> cardiomyocyte
source = {'COL1A1': 0.9, 'VIM': 0.9}  # markers not GRN nodes; use GRN-visible proxy state
src_state = {n: 0.05 for n in nodes}; src_state.update({'OCT4': 0.02, 'NANOG': 0.02})
tgt_state = {'TNNT2': 2.0, 'MYH6': 2.0}  # above basal-driven level: requires genuine driver activation
oc = C.optimal_reprogramming(src_state, tgt_state, hours=48)
res['optimal_reprogramming_cardiac'] = {'interventions': oc['interventions'], 'objective': oc['objective'], 'converged': oc['converged']}
controls = {x['factor']: x['control'] for x in oc['interventions']}
fin = C.simulate_fate(src_state, hours=48, controls=controls)['final_state']
res['cardiac_outcome'] = {'TNNT2': fin['TNNT2'], 'MYH6': fin['MYH6'], 'MEF2C': fin['MEF2C'], 'TBX5': fin['TBX5'],
                          'OCT4': fin['OCT4']}
sv = C.stochastic_validate(src_state, tgt_state, oc['interventions'], replicates=32, noise=.04, seed=8)
res['stochastic_validate_cardiac'] = sv
rec = C.transition_recipe('fibroblast', 'cardiomyocyte')
res['transition_recipe_cardiac'] = {'steps': [(s['order'], s['action'], s['nodes']) for s in rec['recipe']],
                                    'note': 'drivers are direct regulators of target markers only (one-step lookahead)'}
out = sys.argv[1] if len(sys.argv) > 1 else 'benchmarks/sweep_cellfatenet.json'
json.dump(res, open(out, 'w'), indent=1, default=str)
print(json.dumps({'hgnc_unapproved': res['hgnc_unapproved'], 'grn_match': res['grn_matrix_match'],
                  'sim_diff': res['simulate_fate_deterministic_max_abs_diff'],
                  'attractor_stab': res['attractors']['max_rhs_at_centers'],
                  'cardiac_interventions': [(x['factor'], round(x['control'], 3)) for x in oc['interventions']],
                  'cardiac_final': {k: round(v, 3) for k, v in res['cardiac_outcome'].items()},
                  'success_p': sv['success_probability'],
                  'noise_mean_diff': res['noise']['mean_close_to_deterministic_max_diff']}, indent=1))
