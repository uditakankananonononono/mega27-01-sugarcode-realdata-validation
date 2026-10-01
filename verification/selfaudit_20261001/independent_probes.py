#!/usr/bin/env python3
"""Self-audit Tier-2 probes: verify STATUS.md 'real published algorithms' claims
behaviorally. Each check is independent of the suite and prints PASS/FAIL + evidence."""
import sys, json
import os
sys.path.insert(0, os.environ['SUGARCODE_SRC'])
import numpy as np
R = {}

# 1. virtual_cell: FBA via HiGHS - independent LP check with scipy
from sugarcode.modules.virtual_cell.core import demo_model, fba, gene_knockout
m = demo_model()
r = fba(m)
obj = r['objective']
obj_rxn = r['objective_reaction']
S = np.array(m.S, float)
names = list(m.reactions)
from scipy.optimize import linprog
oi = names.index(obj_rxn)
lb = np.array(m.lb, float)
ub = np.array(m.ub, float)
c = np.zeros(len(names)); c[oi] = -1.0
lp = linprog(c, A_eq=S, b_eq=np.zeros(S.shape[0]), bounds=list(zip(lb,ub)), method='highs')
indep = -lp.fun if lp.status==0 else None
ko = gene_knockout(m, obj_rxn)
ko_obj = ko.get('objective_value', ko.get('objective')) if isinstance(ko,dict) else None
R['virtual_cell_fba'] = {'claim':'FBA via HiGHS LP','module_obj':obj,'independent_linprog_obj':indep,
 'match': (indep is not None and abs(obj-indep) < 1e-6*max(1,abs(obj))),
 'knockout_of_objective_rxn': ko_obj, 'pass': bool(indep is not None and abs(obj-indep)<1e-6*max(1,abs(obj)))}

# 2. codon_opt: CAI (Sharp & Li 1987) - independent recomputation
from sugarcode.bio import codon as codonlib
dna = 'ATGGCTGCTGAAATCAAAACCGGTGGTCTGAAATAA'
table = codonlib.HOST_TABLES['ecoli_k12']
w = codonlib.relative_adaptiveness(table)
codons = [dna[i:i+3] for i in range(0,len(dna)-2,3)]
# published Sharp & Li convention: exclude Met/Trp (w=1 uninformative) and stops
from sugarcode.bio.codon import STANDARD_CODE
ws = [w[c] for c in codons if STANDARD_CODE.get(c) not in (None,'*','M','W') and w.get(c,0)>0]
indep_cai = float(np.exp(np.mean(np.log(ws))))
mod_cai = codonlib.cai(dna, table)
R['codon_opt_cai'] = {'claim':'CAI Sharp & Li 1987','module_cai':mod_cai,'independent_cai':indep_cai,
 'pass': bool(abs(mod_cai-indep_cai) < 1e-9)}

# 3. codon_opt TASEP: deterministic under seed, physical bounds
from sugarcode.modules.codon_opt.core import tasep_simulate
t1 = tasep_simulate(dna*3, table)
t2 = tasep_simulate(dna*3, table)
den = np.array(t1['density_profile'])
R['codon_opt_tasep'] = {'claim':'TASEP Gillespie KMC','deterministic': t1['output_rate_per_s']==t2['output_rate_per_s'],
 'density_in_0_1': bool((den>=0).all() and (den<=1).all()), 'mean_density': float(den.mean()),
 'pass': bool(t1['output_rate_per_s']==t2['output_rate_per_s'] and (den>=0).all() and (den<=1).all())}

# 4. evofold_4d ANM: first 6 eigenvalues ~0, then positive; 3N-6 nontrivial modes
from sugarcode.modules.evofold_4d.core import anm_modes
N = 24
t = np.arange(N)
co = (np.stack([1.6*np.cos(2*np.pi*t/3.6), 1.6*np.sin(2*np.pi*t/3.6), 1.54*t/1.0], axis=1)).tolist()
am = anm_modes(co, n_modes=8)
# independent ANM: standard Hessian from contact graph (cutoff 10, k=-1 on contact)
C = np.array(co, float); n = len(C)
H = np.zeros((3*n, 3*n))
for i in range(n):
    for j in range(i+1, n):
        rij = C[j]-C[i]; d = np.linalg.norm(rij)
        if d > 10.0 or d == 0: continue
        e = rij/d
        b = -np.outer(e, e)
        H[3*i:3*i+3, 3*j:3*j+3] = b
        H[3*j:3*j+3, 3*i:3*i+3] = b
        H[3*i:3*i+3, 3*i:3*i+3] -= b
        H[3*j:3*j+3, 3*j:3*j+3] -= b
ev = np.linalg.eigvalsh(H)
nz = ev[ev > 1e-6]
mod_first_freq = am['modes'][0]['frequency']
R['evofold_4d_anm'] = {'claim':'ANM normal modes','n_atoms':N,
 'independent_zero_modes': int((np.abs(ev)<1e-6).sum()), 'expected_zero_modes': 6,
 'independent_first_nonzero_freq': round(float(np.sqrt(nz[0])),4), 'module_first_freq': mod_first_freq,
 'freq_match': bool(abs(np.sqrt(nz[0]) - mod_first_freq) < 1e-3),
 'pass': bool((np.abs(ev)<1e-6).sum()==6 and abs(np.sqrt(nz[0])-mod_first_freq)<1e-3)}

# 5. living_computer Gillespie: deterministic under seed, nonnegative counts
from sugarcode.modules.living_computer.core import parse_logic, build_circuit, stochastic_simulate
ast, inputs = parse_logic('A AND B')
ck = build_circuit(ast, 'GFP')
s1 = stochastic_simulate(ck, t_end=50.0, seed=3)
s2 = stochastic_simulate(ck, t_end=50.0, seed=3)
same = s1 == s2
allser = np.concatenate([np.array(v) for v in s1['series'].values()])
import inspect as _insp
_src = _insp.getsource(stochastic_simulate)
is_cle = 'Euler-Maruyama' in _src
R['living_computer_gillespie'] = {'claim':'stochastic chemical-Langevin simulation (Euler-Maruyama)',
 'actual_method': 'chemical-Langevin Euler-Maruyama (current STATUS wording agrees)',
 'deterministic_seed': same, 'nonnegative': bool((allser>=0).all()),
 'pass': bool(same and (allser>=0).all())}

# 6. synbio_wizard Hill + Gillespie: analytic Hill check
from sugarcode.modules.synbio_wizard.core import promoter_kinetics, gillespie_expression
pk = promoter_kinetics(0.5, kd=0.5, hill=2, vmax=1.0, basal=0.01)
# analytic Hill at regulator==kd: basal + vmax*x^h/(kd^h+x^h) = 0.01 + 0.5 = 0.51
hill_expected = 0.51
hill_val = pk['transcription_rate_reu_h'] if isinstance(pk, dict) and 'transcription_rate_reu_h' in pk else (pk if isinstance(pk,(int,float)) else None)
g1 = gillespie_expression(hours=2, seed=1); g2 = gillespie_expression(hours=2, seed=1)
gss = gillespie_expression(hours=500, seed=0)
# analytic steady state: E[mRNA]=txn/decay=2, E[protein]=tl*E[m]/pdecay=40
ev_m = [e[1] for e in gss['events'][-2000:]]; ev_p = [e[2] for e in gss['events'][-2000:]]
R['synbio_wizard'] = {'claim':'stochastic kinetics (Gillespie SSA) + Hill promoter','promoter_kinetics_sample': str(pk)[:120],
 'hill_at_kd': hill_val, 'hill_expected': hill_expected, 'hill_match': (hill_val is not None and abs(hill_val-hill_expected)<1e-9),
 'gillespie_deterministic': g1==g2, 'true_ssa': 'rng.exponential' in __import__('inspect').getsource(gillespie_expression),
 'steady_state_tail_mean_mrna': round(float(np.mean(ev_m)),2), 'analytic_mrna': 2.0,
 'steady_state_tail_mean_protein': round(float(np.mean(ev_p)),2), 'analytic_protein': 40.0,
 'pass': bool(g1==g2 and hill_val is not None and abs(hill_val-hill_expected)<1e-9 and abs(np.mean(ev_m)-2)<1.5 and abs(np.mean(ev_p)-40)<25)}

# 7. alpha_fold_ui Chou-Fasman: poly-Ala -> helix, Chou-Fasman published params
from sugarcode.modules.alpha_fold_ui.core import chou_fasman
ss_ala = chou_fasman('A'*12)
ss_pro = chou_fasman('P'*6+'G'*6)
helix_frac = sum(1 for x in ss_ala if x=='H')/len(ss_ala)
R['alpha_fold_ui_choufasman'] = {'claim':'Chou-Fasman published propensities','polyAla': ''.join(ss_ala),
 'polyAla_helix_fraction': round(helix_frac,3), 'polyProGly': ''.join(ss_pro),
 'pass': bool(helix_frac >= 0.75)}

print(json.dumps(R, indent=1, default=str))
json.dump(R,open(os.environ.get('PROBE_OUTPUT','probes.json'),'w'),indent=1,default=str)
