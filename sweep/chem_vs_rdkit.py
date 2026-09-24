"""chem_descriptors + chem_similarity vs RDKit on 993 ChEMBL phase-4 small
molecules (held out from the module's 66-molecule ladder)."""
import csv, json, os, sys, random, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, os.environ.get("SUGARCODE_SRC", str(ROOT.parent / "sc-ai" / "src")))
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, DataStructs
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
from sugarcode.modules.chem_descriptors import compute_descriptors
from sugarcode.modules.chem_similarity import morgan as mg, similarity as sim
rows = list(csv.DictReader(open(ROOT / "data" / "chembl_approved_1000.tsv"), delimiter="\t"))
ref = {"mol_wt": Descriptors.MolWt, "exact_mol_wt": Descriptors.ExactMolWt, "hbd": Lipinski.NumHDonors,
       "hba": Lipinski.NumHAcceptors, "rotatable_bonds": Descriptors.NumRotatableBonds,
       "tpsa": rdMolDescriptors.CalcTPSA, "fraction_csp3": rdMolDescriptors.CalcFractionCSP3,
       "heavy_atom_count": lambda m: m.GetNumHeavyAtoms(), "formula": rdMolDescriptors.CalcMolFormula}
stat = {k: [0, 0] for k in ref}; parse_fail = 0; rd_fail = 0; mism = {k: [] for k in ref}
ok_mols = []
for r in rows:
    m = Chem.MolFromSmiles(r["smiles"])
    if m is None: rd_fail += 1; continue
    try: d = compute_descriptors(r["smiles"])
    except Exception: parse_fail += 1; continue
    ok_mols.append((r["smiles"], m))
    for k, f in ref.items():
        v = f(m); stat[k][1] += 1
        good = (v == d[k]) if isinstance(v, (int, str)) else abs(v - d[k]) <= 0.01 + 1e-9 * abs(v)
        stat[k][0] += good
        if not good and len(mism[k]) < 3: mism[k].append([r["chembl_id"], d[k], v])
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
random.seed(7); pairs = [random.sample(range(len(ok_mols)), 2) for _ in range(500)]
diffs = []; fp_fail = 0
for i, j in pairs:
    try: t_mod = sim.tanimoto(mg.morgan_bits(ok_mols[i][0]), mg.morgan_bits(ok_mols[j][0]))
    except Exception: fp_fail += 1; continue
    t_rd = DataStructs.TanimotoSimilarity(gen.GetFingerprint(ok_mols[i][1]), gen.GetFingerprint(ok_mols[j][1]))
    diffs.append(abs(t_mod - t_rd))
res = {"n_input": len(rows), "rdkit_parse_fail": rd_fail, "module_parse_fail": parse_fail,
       "descriptor_agreement": {k: {"agree": a, "n": n, "rate": round(a / n, 4)} for k, (a, n) in stat.items()},
       "examples_mismatch": mism,
       "tanimoto_pairs": len(diffs), "tanimoto_fp_fail": fp_fail,
       "tanimoto_exact_frac": round(sum(x < 1e-9 for x in diffs) / max(1, len(diffs)), 4),
       "tanimoto_max_abs_diff": max(diffs) if diffs else None}
json.dump(res, open(ROOT / "benchmarks" / "sweep_chem.json", "w"), indent=1)
print(json.dumps({k: v for k, v in res.items() if k != "examples_mismatch"}, indent=1))
