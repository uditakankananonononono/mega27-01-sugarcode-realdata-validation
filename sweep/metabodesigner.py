"""Sweep: metabodesigner vs KEGG ENZYME/COMPOUND records + textbook fermentation yields."""
import json, os, re, sys, time, urllib.request, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
import numpy as np
from sugarcode.modules.metabodesigner import core as M

CACHE = os.path.join(os.path.dirname(__file__), '..', 'data', 'kegg')
os.makedirs(CACHE, exist_ok=True)

def kegg(entry):
    fn = os.path.join(CACHE, entry.replace(':', '_') + '.txt')
    if not os.path.exists(fn):
        for i in range(3):
            try:
                txt = urllib.request.urlopen('https://rest.kegg.jp/get/' + entry, timeout=30).read().decode()
                break
            except Exception:
                time.sleep(2)
        else:
            raise
        open(fn, 'w').write(txt); time.sleep(0.4)
    return open(fn).read()

def field(txt, name):
    out, on = [], False
    for line in txt.splitlines():
        if line[:12].strip():
            on = line[:12].strip() == name
        if on:
            out.append(line[12:].strip())
    return out

CPD = {'glucose': 'C00031', 'g6p': 'C00092', 'f6p': 'C00085', 'f1,6bp': 'C00354', 'gap': 'C00118',
       'dhap': 'C00111', '1,3bpg': 'C00236', '3pg': 'C00197', 'pep': 'C00074', 'pyruvate': 'C00022',
       'acetyl-coa': 'C00024', 'lactate': 'C00186', 'ethanol': 'C00469', 'citrate': 'C00158',
       'succinate': 'C00042', 'glyoxylate': 'C00048', 'malate': 'C00149', '3-hydroxybutyrate': 'C01089',
       'coa': 'C00010', 'co2': 'C00011'}
ALIAS = {'glucose': ['C00738'], 'g6p': ['C00668', 'C02965'], 'f6p': ['C05345'], 'f1,6bp': ['fructofuranose 1,6-bisphosphate', 'C05378'],
         'lactate': ['C00186'], 'co2': ['C00011']}
def ids(x):
    return ([CPD[x]] if x in CPD else []) + ALIAS.get(x, [])
def inside(x, text):
    return any(i in text for i in ids(x))

def carbons(met):
    txt = kegg('cpd:' + CPD[met]); f = field(txt, 'FORMULA')[0]
    m = re.match(r'C(\d*)', f); c = int(m.group(1) or 1) if m else 0
    return c, f

res = {'module': 'metabodesigner', 'sources': ['https://rest.kegg.jp/get/ec:<EC>', 'https://rest.kegg.jp/get/cpd:<C-id>']}
# 1. EC name + substrate/product membership vs KEGG ENZYME
ec_rows = []
for r in M.REACTION_DB:
    row = {'reaction': r['id'], 'ec': r['ec'], 'enzyme': r['enzyme']}
    if r['ec'].endswith('-'):
        row['status'] = 'EC incomplete (class-level), not checkable'; ec_rows.append(row); continue
    txt = kegg('ec:' + r['ec'])
    names = [n.rstrip(';').lower() for n in field(txt, 'NAME')]
    row['kegg_names'] = names[:3]
    head = r['enzyme'].split('(')[0].split('+')[0].strip().lower()
    row['name_match'] = any(head in n or n in head for n in names) or any(w in ' '.join(names) for w in head.split() if len(w) > 6)
    subs = ' '.join(field(txt, 'SUBSTRATE')); prods = ' '.join(field(txt, 'PRODUCT'))
    rxn_all = subs + ' ' + prods
    s_ok = [inside(x, subs) for x in r['substrates']]
    p_ok = [inside(x, prods) for x in r['products']]
    s_rev = [inside(x, prods) for x in r['substrates']]
    p_rev = [inside(x, subs) for x in r['products']]
    row['kegg_substrates_match'] = all(s_ok) and all(p_ok) or (all(s_rev) and all(p_rev))
    row['unmatched'] = [x for x, ok in zip(r['substrates'] + r['products'], s_ok + p_ok) if not ok]
    row['status'] = 'exact' if row['kegg_substrates_match'] and row['name_match'] else ('name ok, lumped/abbreviated step' if row['name_match'] else 'MISMATCH')
    ec_rows.append(row)
res['ec_check'] = ec_rows
res['ec_summary'] = {k: sum(1 for x in ec_rows if x['status'] == k) for k in sorted({x['status'] for x in ec_rows})}

# 2. carbon balance per reaction (acyl-CoA counted as acyl carbons: formula C minus CoA C)
coa_c, _ = carbons('coa')
def acyl_c(m):
    c, f = carbons(m)
    return c - coa_c if m == 'acetyl-coa' else c
bal = []
for r in M.REACTION_DB:
    if any(x not in CPD for x in r['substrates'] + r['products']):
        bal.append({'reaction': r['id'], 'status': 'not checkable (polymer/unmapped)'}); continue
    cs = sum(acyl_c(x) for x in r['substrates']); cp = sum(acyl_c(x) for x in r['products'])
    bal.append({'reaction': r['id'], 'substrate_C': cs, 'product_C': cp, 'balanced': cs == cp,
                'note': '' if cs == cp else ('CO2 release not modelled' if cs > cp else 'carbon source missing from substrates (co-substrate treated as cofactor or stoichiometry missing)')})
res['carbon_balance'] = bal
res['carbon_balance_summary'] = {'balanced': sum(1 for b in bal if b.get('balanced') is True),
                                 'unbalanced': sum(1 for b in bal if b.get('balanced') is False),
                                 'not_checkable': sum(1 for b in bal if 'balanced' not in b)}

# 3. molar yields per glucose from the module's own LP vs textbook
TEXT = {'lactate': 2.0, 'ethanol': 2.0, '3-hydroxybutyrate': 1.0}
yl = []
for t, ref in TEXT.items():
    d = M.design_pathway(t); route = d['route']
    fl = M.pathway_flux(route, t, 'glucose', upper=1.0)
    S = np.asarray(fl['stoichiometry']['matrix']); mets = fl['stoichiometry']['metabolites']
    v = np.asarray(fl['fluxes'])
    net = S @ v
    glc = -net[mets.index('glucose')]; prod = net[mets.index(t)]
    yl.append({'target': t, 'steps': d['steps'], 'route': [x['reaction'] for x in route],
               'glucose_uptake': round(float(glc), 6), 'product_out': round(float(prod), 6),
               'yield_mol_per_mol': round(float(prod / glc), 6) if glc else None, 'textbook': ref,
               'match': bool(glc) and abs(prod / glc - ref) < 1e-6,
               'dead_end_byproducts': [m for i, m in enumerate(mets) if m not in ('glucose', t) and net[i] > 1e-9]})
res['yields'] = yl
res['textbook_refs'] = {'lactate': 'homolactic fermentation, glucose -> 2 lactate (Embden-Meyerhof)',
                        'ethanol': 'glucose -> 2 ethanol + 2 CO2 (Gay-Lussac)',
                        '3-hydroxybutyrate': '2 acetyl-CoA -> 1 (R)-3-hydroxybutyryl-CoA (PhaA thiolase EC 2.3.1.9 condenses two acetyl-CoA)'}

# 4. internal identities
st = M.stoichiometric_matrix(M.design_pathway('lactate')['route'])
dyn = M.dynamic_pathway([1.0, 0.5, 2.0], initial_substrate=10, hours=20)
C = np.asarray(dyn['concentrations'])
res['dynamic_mass_conservation_max_dev'] = float(np.max(np.abs(C.sum(0) - 10)))
th = M.thermodynamics(M.design_pathway('lactate')['route'])
res['thermodynamics_note'] = 'default dG = -5 + n_cofactors kJ/mol per step is a placeholder prior, not eQuilibrator values (unverified)'
res['thermo_default_total_lactate'] = th['total_dg']
bn = M.bottleneck_analysis(M.design_pathway('lactate')['route'])
res['bottleneck_top_lactate'] = bn['top_bottleneck']
res['kcat_note'] = 'class-level kcat priors (kinase 100, dehydrogenase 50 ...) are not BRENDA values (unverified)'
out = sys.argv[1] if len(sys.argv) > 1 else 'benchmarks/sweep_metabodesigner.json'
json.dump(res, open(out, 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(json.dumps({'ec_summary': res['ec_summary'], 'carbon': res['carbon_balance_summary'],
                  'yields': [(y['target'], y['yield_mol_per_mol'], y['textbook'], y['dead_end_byproducts']) for y in yl],
                  'dyn_dev': res['dynamic_mass_conservation_max_dev']}, indent=1))
for x in ec_rows: print(x['reaction'], x['ec'], x['status'], x.get('kegg_names', '')[:1] if x.get('kegg_names') else '')
for b in bal: print(b)
