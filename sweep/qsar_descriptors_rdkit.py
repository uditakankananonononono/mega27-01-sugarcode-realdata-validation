"""qsar_bench.chemistry.descriptors vs RDKit on 993 ChEMBL phase-4 molecules."""
import csv, json, os, sys, collections
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
RDLogger.DisableLog('rdApp.*')
from sugarcode.modules.qsar_bench.chemistry import descriptors
rows = list(csv.DictReader(open(os.path.expanduser('~/mega/m01/data/chembl_approved_1000.tsv')), delimiter='\t'))
def ringrank(m): return m.GetNumBonds() - m.GetNumAtoms() + len(Chem.GetMolFrags(m))
REF = {'molecular_weight': (Descriptors.MolWt, 0.05), 'heavy_atoms': (lambda m: m.GetNumHeavyAtoms(), 0),
       'hetero_atoms': (rdMolDescriptors.CalcNumHeteroatoms, 0), 'formal_charge': (Chem.GetFormalCharge, 0),
       'hbd_heuristic': (Lipinski.NumHDonors, 0), 'hba_heuristic': (Lipinski.NumHAcceptors, 0),
       'rotatable_bonds_heuristic': (lambda m: rdMolDescriptors.CalcNumRotatableBonds(m, rdMolDescriptors.NumRotatableBondsOptions.NonStrict), 0), 'ring_rank': (ringrank, 0),
       'fraction_csp3': (rdMolDescriptors.CalcFractionCSP3, 0.01)}
stat = {k: [0, 0] for k in REF}; ex = {k: [] for k in REF}; pfail = []; n = 0; strict_ok = []
for r in rows:
    m = Chem.MolFromSmiles(r['smiles'])
    if m is None: continue
    try: d = descriptors(r['smiles'])
    except Exception as e: pfail.append([r['chembl_id'], type(e).__name__]); continue
    n += 1
    strict = rdMolDescriptors.CalcNumRotatableBonds(m); strict_ok.append(strict == d['rotatable_bonds_heuristic'])
    for k, (f, tol) in REF.items():
        v = f(m); stat[k][1] += 1
        ok = abs(v - d[k]) <= tol + 1e-9
        stat[k][0] += ok
        if not ok and len(ex[k]) < 5: ex[k].append([r['chembl_id'], r['smiles'][:80], d[k], v])
res = {'reference': 'RDKit ' + __import__('rdkit').__version__, 'n_molecules_rdkit_parsed': sum(Chem.MolFromSmiles(r['smiles']) is not None for r in rows),
       'n_scored': n, 'parse_failures': len(pfail), 'parse_failure_examples': pfail[:10],
       'agreement': {k: {'match': a, 'n': b} for k, (a, b) in stat.items()}, 'mismatch_examples': ex, 'rotatable_vs_rdkit_default_strict': {'match': sum(strict_ok), 'n': len(strict_ok)}, 'note': 'rotatable_bonds reference is RDKit NonStrict (SMARTS [!$(*#*)&!D1]-&!@[!$(*#*)&!D1]); RDKit default Strict additionally excludes amide C-N and similar bonds'}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_qsar_descriptors.json'), 'w'), indent=1)
print(res['n_scored'], res['parse_failures'], res['parse_failure_examples'][:5]); print({k: f"{a}/{b}" for k, (a, b) in stat.items()})
for k, v in ex.items():
    if v: print(k, v[:2])
