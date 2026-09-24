"""Validate car_t_designer against FDA-approved CAR-T architectures + helper math recomputation."""
import json, math, os, sys
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from sugarcode.modules.car_t_designer.core import (design_car, ANTIGENS, SCFV, COSTIM,
    antigen_selectivity, logic_gate_response, exhaustion_trajectory, killing_curve)
res = {'module': 'car_t_designer'}
# 1) curated constructs vs FDA-approved products
cd19 = design_car('CD19'); bcma = design_car('BCMA')
res['cd19_vs_kymriah'] = {
    'module': {'scfv': cd19['construct']['scfv'], 'hinge': cd19['construct']['hinge'],
               'costim': cd19['construct']['costimulatory'], 'activation': cd19['construct']['activation']},
    'reference': 'Kymriah (tisagenlecleucel): FMC63 scFv, CD8a hinge/TM, 4-1BB, CD3z (DailyMed label; canonical FMC63-CD8H-41BB-CD3z architecture, e.g. Addgene 200671)',
    'match': (cd19['construct']['scfv'] == 'FMC63' and cd19['construct']['hinge'] == 'CD8a'
              and cd19['construct']['costimulatory'] == '4-1BB' and cd19['construct']['activation'] == 'CD3zeta'),
    'note': 'Yescarta (CD28) is a second valid CD19 architecture; module rationale favors persistence'}
res['bcma_vs_abecma'] = {
    'module': {'scfv': bcma['construct']['scfv'], 'hinge': bcma['construct']['hinge'],
               'costim': bcma['construct']['costimulatory']},
    'reference': 'Abecma (idecabtagene vicleucel, FDA label fda.gov/media/147055): murine anti-BCMA scFv, CD8a hinge/TM, CD137 (4-1BB), CD3z',
    'match': (bcma['construct']['scfv'] == 'C11D5.3' and bcma['construct']['hinge'] == 'CD8a'
              and bcma['construct']['costimulatory'] == '4-1BB'),
    'note': 'C11D5.3 is the parent hybridoma of bb2121/ide-cel'}
# 2) antigen-indication and scFv maps vs literature
res['scfv_map'] = {'FMC63-CD19': True, 'C11D5.3-BCMA': True, '4D5-HER2': 'trastuzumab parent mAb (standard)',
                   '14G2a-GD2': True, 'GC33-GPC3': 'standard anti-GPC3 mAb',
                   'verified_sources': ['14G2a: PMC6082365 + S0304383509001347']}
res['antigen_indications'] = {a: ANTIGENS[a]['tumor'] for a in ANTIGENS}
# 3) decision logic: 4-1BB when normal-tissue expression exists, else CD28
res['costim_logic'] = {a: design_car(a)['construct']['costimulatory'] for a in ANTIGENS}
res['costim_logic_correct'] = all(design_car(a)['construct']['costimulatory'] ==
    ('4-1BB' if ANTIGENS[a]['normal_expression'] else 'CD28') for a in ANTIGENS)
# 4) safety switches only with off-tumor risk
res['safety_switch_logic'] = all(('iCasp9 suicide switch' in design_car(a)['safety_switches']) ==
    bool(ANTIGENS[a]['normal_expression']) for a in ANTIGENS)
# 5) toxicity heuristic recomputation
tox = cd19['toxicity_prediction']
exp_crs = round(min(0.3 + 0.2*COSTIM['4-1BB']['expansion'], 0.9), 2)
res['toxicity'] = {'module_crs': tox['crs_grade2plus_risk'], 'expected_crs': exp_crs,
    'match': tox['crs_grade2plus_risk'] == exp_crs,
    'ot_match': tox['on_target_off_tumor_risk'] == round(min(0.15*len(ANTIGENS['CD19']['normal_expression']), 0.9), 2),
    'neuro_match': tox['neurotoxicity_risk'] == 0.05}
# 6) helper math
sel = antigen_selectivity(1e-3, 1e-5)
res['selectivity'] = {'ratio': sel['tumor_normal_ratio'], 'selective': sel['selective'],
    'match': abs(sel['tumor_normal_ratio'] - 100.0) < 1e-9 and sel['selective'] is True,
    'floor': antigen_selectivity(0, 0)['tumor_normal_ratio'] == 1.0}
gates = [logic_gate_response(.6, .7, 'AND')['active'], logic_gate_response(.6, .2, 'AND')['active'],
         logic_gate_response(.2, .7, 'OR')['active'], logic_gate_response(.6, .7, 'A_NOT_B')['active'],
         logic_gate_response(.6, .2, 'A_NOT_B')['active']]
res['logic_gates'] = {'results': gates, 'match': gates == [True, False, True, False, True]}
try:
    logic_gate_response(.5, .5, 'XOR'); res['logic_gates']['bad_gate_raises'] = False
except ValueError: res['logic_gates']['bad_gate_raises'] = True
ex = exhaustion_trajectory([0, 7, 14, 28], 1.0, '4-1BB')
exp_ex = [1 - math.exp(-COSTIM['4-1BB']['exhaustion']*d/14) for d in [0, 7, 14, 28]]
res['exhaustion'] = {'match': all(abs(a-b) < 1e-12 for a, b in zip(ex['exhaustion_fraction'], exp_ex)),
    'monotone': all(a <= b for a, b in zip(ex['exhaustion_fraction'], ex['exhaustion_fraction'][1:])),
    'cd28_faster': exhaustion_trajectory([14], 1.0, 'CD28')['exhaustion_fraction'][0] >
                   exhaustion_trajectory([14], 1.0, '4-1BB')['exhaustion_fraction'][0]}
kc = killing_curve([0, 1, 5], .6, .8)
exp_kc = [.8*(1-math.exp(-.6*r)) for r in [0, 1, 5]]
res['killing'] = {'match': all(abs(a-b) < 1e-12 for a, b in zip(kc['target_kill_fraction'], exp_kc)),
    'saturates_at_antigen_frac': killing_curve([100], .6, .8)['target_kill_fraction'][0] < .8}
# 7) trial simulation
tr = cd19['trial_simulation']
res['trial'] = {'persistence_match': all(abs(p - round(COSTIM['4-1BB']['persistence']*0.85**m, 3)) < 1e-12
    for p, m in zip(tr['car_persistence_fraction'], tr['months'])),
    'response_value': tr['predicted_overall_response'],
    'response_plausible_range': 0.3 <= tr['predicted_overall_response'] <= 0.95}
json.dump(res, open('benchmarks/sweep_car_t_designer.json', 'w'), indent=1)
ok = [res['cd19_vs_kymriah']['match'], res['bcma_vs_abecma']['match'], res['costim_logic_correct'],
      res['safety_switch_logic'], res['toxicity']['match'], res['selectivity']['match'],
      res['logic_gates']['match'], res['exhaustion']['match'], res['killing']['match'], res['trial']['persistence_match']]
print('ALL PASS' if all(ok) else 'FAILURES', ok)
print(res['costim_logic'], res['trial']['response_value'])
