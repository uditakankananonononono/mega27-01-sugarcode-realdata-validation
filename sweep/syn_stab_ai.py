"""Sweep: syn_stab_ai recomputation + BUG 49 palindrome evidence (real restriction sites)."""
import json, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
from sugarcode.modules.syn_stab_ai import core as S
res = {'module': 'syn_stab_ai', 'sources': ['recomputation', 'real inverted-repeat restriction sites (EcoRI/HindIII/SnaBI)']}
# 1. evaluate_circuit_stability exact recompute vs stability_ai formula
c = {'n_gates': 4, 'expression_level': 0.6, 'mode': 'plasmid', 'generations': 200}
ev = S.evaluate_circuit_stability(c)
burden = min(1.0, 0.25 * 4 * 0.6); loss = 0.002 * (1 + 4 * burden); size = 4.0 * 4 / 3
exp_last = (1 - loss) ** 200 * math.exp(-1e-4 * size * 200)
res['evaluate'] = {'burden_module': ev['estimated_burden'], 'burden_recomputed': round(burden, 3),
                   'failure_risk_module': ev['circuit_failure_risk'], 'failure_recomputed': round(1 - exp_last, 3),
                   'trajectory_last': ev['forecast']['trajectory'][-1]}
# 2. population_simulation deterministic limit vs closed recursion
pop = S.population_simulation(mutation_rate=.001, selection_cost=.01, drift_population=10**15, generations=100, seed=0)
f = 1.0
for g in range(100): f = f * (1 - .001) * (1 - .01) / (1 - f * .01)
res['population'] = {'final_module': pop['functional_fraction'][-1], 'final_recomputed': round(f, 8),
                     'match': abs(pop['functional_fraction'][-1] - f) < 1e-4,  # residual diffusion noise at finite N,
                     'half_life': pop['half_life_generation'],
                     'half_life_recomputed': next((i for i in range(101) if (lambda x: x)([None]) is None), None)}
# recompute half-life from the deterministic recursion
f = 1.0; hl = None
for i in range(101):
    if f < .5 and hl is None: hl = i
    f = f * (1 - .001) * (1 - .01) / (1 - f * .01)
res['population']['half_life_recomputed'] = hl
res['population']['half_life_match'] = pop['half_life_generation'] == hl
# 3. palindrome evidence: real restriction sites are inverted repeats
sites = {'EcoRI_GAATTC': 'GAATTC', 'HindIII_AAGCTT': 'AAGCTT', 'SnaBI_TACGTA': 'TACGTA', 'SfiI_GGCCNNNNNGGCC_core': 'GGCCGGCC'}
res['palindromes'] = {name: S.sequence_risk('TT' + sq + 'AA')['palindrome_count'] for name, sq in sites.items()}
res['palindromes']['pre_fix'] = 'all 0 before BUG 49 fix (plain-reverse comparison); textual GATTAG still 0 after fix'
res['palindromes']['textual_GATTAG_post_fix'] = S.sequence_risk('TTGATTAGAA')['palindrome_count']
# 4. chemical_stress / resource_burden / fitness_landscape recomputation
cs = S.chemical_stress(oxidative=.3, ph=6.9, reactive=.2)
res['chemical'] = {'module': cs['damage_rate'], 'recomputed': round(min(1, .5*.3 + .2*.5 + .3*.2), 6)}
fl = S.fitness_landscape([{'id': 'loss', 'burden_relief': .3, 'functional_loss': .8}, {'id': 'cheap_loss', 'burden_relief': .5, 'functional_loss': .1}])
res['fitness_landscape'] = {'fitness_values': [(m['id'], m['fitness']) for m in fl['mutations']],
                            'favored_first': fl['favored'][0]['id'], 'expected': 'cheap_loss'}
# 5. compare_stabilization ordering: genomic integration beats plasmid
cmp_ = S.compare_stabilization({'n_gates': 4, 'expression_level': .6, 'mode': 'plasmid'},
                               [{'name': 'integrate', 'mode': 'genomic'}, {'name': 'weaker promoter', 'expression_factor': .5}])
res['compare_stabilization'] = [{'name': x['name'], 'risk_reduction': round(x['risk_reduction'], 4)} for x in cmp_]
fn = os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_syn_stab_ai.json')
json.dump(res, open(fn, 'w'), indent=1, default=str)
print(json.dumps(res, default=str)[:3000])
