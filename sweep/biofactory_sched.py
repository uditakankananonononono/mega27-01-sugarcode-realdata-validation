"""Validate biofactory_1_a: BUG 36 schedule dependency fix + Bayesian/PID/Z-factor/activity math vs recomputation + protocol steps vs standard methods."""
import json, math, os, sys
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from sugarcode.modules.biofactory_1_a.core import (generate_protocol, workflow_dag, schedule_resources,
    bayesian_condition, pid_control, reagent_activity, assay_quality, reroute_on_failure, _step_minutes)
res = {'module': 'biofactory_1_a'}
# 1) BUG 36: schedule must respect linear DAG order
out = {}
for wf in ('golden_gate', 'gibson', 'pcr_screen'):
    p = generate_protocol(wf, 8)
    sch = schedule_resources(p)['schedule']
    seq = all(b['start_min'] >= a['end_min'] - 1e-9 for a, b in zip(sch, sch[1:]))
    out[wf] = {'ordered': seq, 'makespan': schedule_resources(p)['makespan_min'],
               'sequential_sum': sum(_step_minutes(x) for x in p['steps'])}
res['schedules'] = out
res['bug36_prefix'] = {'golden_gate_makespan': 120.0, 'sequential': 274,
    'symptom': 'thermocycle scheduled at t=0 in parallel with its input dispense steps'}
# 2) Bayesian condition: normal-normal conjugate update with known noise variance
bc = bayesian_condition(.7, .1, [.6, .8, .9], .05)
prec = 1/.1 + 3/.05; mean = (.7/.1 + (.6+.8+.9)/.05)/prec
res['bayesian'] = {'module_mean': bc['posterior_mean'], 'expected_mean': mean,
    'match': abs(bc['posterior_mean'] - mean) < 1e-12 and abs(bc['posterior_variance'] - 1/prec) < 1e-12,
    'shrinks_toward_data': abs(bc['posterior_mean'] - sum([.6,.8,.9])/3) < abs(.7 - sum([.6,.8,.9])/3)}
# 3) PID recomputation
pc = pid_control(10, [8, 9, 11], .5, .1, .05, 1)
integ = 0; prev = 0; exp = []
for m in [8, 9, 11]:
    err = 10 - m; integ += err; exp.append(.5*err + .1*integ + .05*(err - prev)); prev = err
res['pid'] = {'match': all(abs(a-b) < 1e-12 for a, b in zip(pc['controls'], exp)),
    'final_error': pc['final_error'], 'integral': pc['integral_error']}
# 4) reagent activity
ra = reagent_activity(1.0, 10, 5, 0.9)
exp_a = 0.9 * 2**(-10/5)
res['reagent_activity'] = {'match': abs(ra['activity_fraction'] - exp_a) < 1e-12
    and abs(ra['compensation_factor'] - 1/exp_a) < 1e-9}
# 5) Z-prime factor
aq = assay_quality([10, 11, 9.5], [1, 1.2, .8])
import statistics
ms, mb = statistics.mean([10,11,9.5]), statistics.mean([1,1.2,.8])
ss, sb = statistics.stdev([10,11,9.5]), statistics.stdev([1,1.2,.8])
exp_z = 1 - 3*(ss+sb)/abs(ms-mb)
res['z_prime'] = {'module': aq['z_prime'], 'expected': exp_z, 'match': abs(aq['z_prime'] - exp_z) < 1e-12,
    'robust_flag': aq['robust'] == (exp_z > .5)}
# 6) DAG + reroute
p = generate_protocol('gibson', 8); dag = workflow_dag(p['steps'])
rr = reroute_on_failure(dag, 2, {'op': 'incubate', 'temp_c': 50, 'min': 120})
res['dag'] = {'linear_edges': dag['edges'] == [{'source': i, 'target': i+1, 'material_flow': True} for i in range(len(p['steps'])-1)],
    'reroute_removes_failed_outgoing': all(e['source'] != 2 or e.get('condition') == 'failure' for e in rr['edges']),
    'fallback_added': rr['nodes'][-1]['min'] == 120}
# 7) protocol steps vs standard methods (text audit)
gg = generate_protocol('golden_gate', 8)
res['protocol_audit'] = {
    'golden_gate_thermocycle': gg['steps'][3]['program'],
    'reference': 'NEB Golden Gate: 30x(37C/16C), final digest + heat kill; module 50C 10min + 80C 10min - plausible variant',
    'gibson_incubation': generate_protocol('gibson', 8)['steps'][2],
    'gibson_reference': 'NEB Gibson: 50C, 15-60 min - module 60 min at upper bound',
    'transformation': '42C 45s heat shock DH5-alpha - standard',
    'plates_math': generate_protocol('golden_gate', 100)['plates'] == 2 and generate_protocol('golden_gate', 96)['plates'] == 1,
    'screening_candidates': gg['screening']['candidates'] == 64}
json.dump(res, open('benchmarks/sweep_biofactory.json', 'w'), indent=1)
checks = [all(v['ordered'] and v['makespan'] == v['sequential_sum'] for v in out.values()),
    res['bayesian']['match'], res['pid']['match'], res['reagent_activity']['match'],
    res['z_prime']['match'], res['dag']['linear_edges'], res['dag']['reroute_removes_failed_outgoing']]
print('ALL PASS' if all(checks) else 'FAIL', checks)
print(out)
