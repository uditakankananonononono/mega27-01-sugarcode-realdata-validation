"""Sweep: alpha_fold_ui vs real RCSB/AlphaFold structures + exact recomputation."""
import json, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
import numpy as np
from sugarcode.modules.alpha_fold_ui import core as A
from sugarcode.bio import structures as ST
res = {'module': 'alpha_fold_ui', 'sources': ['RCSB PDB (files.rcsb.org)', 'AlphaFold DB (alphafold.ebi.ac.uk)', 'Bio.PDB roundtrip', 'numpy recomputation']}

# 1. real experimental structure: 1UBQ ubiquitin
an = A.analyze_real_structure('1UBQ')
raw = ST.fetch_pdb('1ubq')
bs = [r['bfactor'] for r in raw['residues']]
res['real_1ubq'] = {'n_residues': an['n_residues'], 'expected': 76, 'chains': an['chains'],
    'b_mean_module': an['confidence_stats']['mean'], 'b_mean_recomputed': round(sum(bs)/len(bs), 2),
    'b_minmax_ok': (an['confidence_stats']['min'], an['confidence_stats']['max']) == (round(min(bs),2), round(max(bs),2)),
    'method': an.get('method'), 'resolution_A': an.get('resolution_A'), 'n_pockets': len(an['pockets'])}
# pocket density recomputation for the top pocket
coords = np.array([r['ca'] for r in raw['residues']])
dens = np.array([int(((np.linalg.norm(coords-c,axis=1)>1e-6)&(np.linalg.norm(coords-c,axis=1)<9.0)).sum()) for c in coords])
if an['pockets']:
    top = an['pockets'][0]
    res['real_1ubq']['top_pocket_density_recomputed_ok'] = abs(top['mean_density'] - float(dens[[r-1 for r in top['residues']]].mean())) < 0.15
# 2. real AlphaFold model: P68871 (HBB)
af = A.analyze_real_structure('P68871')
res['real_af_hbb'] = {'n_residues': af['n_residues'], 'expected': 147, 'mean_plddt': af['mean_plddt'],
                      'fraction_low_confidence': af['fraction_low_confidence'], 'confidence_kind': af['confidence_kind']}
# 3. Chou-Fasman on ubiquitin vs known fold (helix 23-34, 5 beta strands); CF table vs canonical 1978
seq76 = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG'
ss = A.chou_fasman(seq76)
helix_region = ''.join(ss[22:34])
CANON = {"A":(142,83,66),"R":(98,93,95),"N":(67,89,156),"D":(101,54,146),"C":(70,119,119),"Q":(111,110,98),
 "E":(151,37,74),"G":(57,75,156),"H":(100,87,95),"I":(108,160,47),"L":(121,130,59),"K":(114,74,101),
 "M":(145,105,60),"F":(113,138,60),"P":(57,55,152),"S":(77,75,143),"T":(83,119,96),"W":(108,137,96),
 "Y":(69,147,114),"V":(106,170,50)}
ref = ['C']*76
for lo,hi,ssr in [(2,7,'E'),(10,17,'E'),(23,34,'H'),(41,45,'E'),(48,49,'E'),(65,71,'E')]:
    for k in range(lo-1,hi): ref[k]=ssr
q3 = sum(1 for a,b in zip(ss,ref) if a==b)/76
res['chou_fasman'] = {'table_matches_canonical_1978': A.CF == CANON,
    'ubq_ss_post_fix': ''.join(ss), 'ubq_helix_23_34_pred': helix_region,
    'ubq_helix23_34_H_fraction': round(helix_region.count('H')/12, 2),
    'q3_vs_1ubq_literature_assignment': round(q3, 3),
    'pre_fix': 'BUG 48: runaway extension predicted HHHHHHHHHEHHH...H (70/76 helix, 92%) on ubiquitin',
    'note': 'ubiquitin has one alpha helix (residues 23-34) and a 5-strand beta sheet (experimental, 1UBQ); CF is a ~50-60% accuracy method'}
# 4. predict_structure geometry + PDB roundtrip
ps = A.predict_structure(seq76)
from Bio.PDB import PDBParser
import io
struct = PDBParser(QUIET=True).get_structure('x', io.StringIO(ps['pdb']))
cas = [a.get_vector().get_array() for a in struct.get_atoms() if a.get_id() == 'CA']
cas = np.array(cas)
consec = np.linalg.norm(np.diff(cas, axis=0), axis=1)
far = [np.linalg.norm(cas[i]-cas[j]) for i in range(len(cas)) for j in range(i+2, len(cas))]
res['geometry'] = {'pdb_ca_count': len(cas), 'expected': 76, 'consec_CA_mean': round(float(consec.mean()), 4),
    'consec_CA_minmax': [round(float(consec.min()),3), round(float(consec.max()),3)], 'target': A.CA_CA,
    'min_nonbonded_CA': round(min(far),3), 'MIN_NONBONDED_CA': A.MIN_NONBONDED_CA,
    'mean_plddt': ps['mean_plddt'], 'composition': ps['composition']}
# 5. msa_couplings exact recompute
msa = ['ACDEFG','ACDEYG','ACDEFG','ACDAFG','TCLEFG']
mc = A.msa_couplings(msa)
alpha='ACDEFGHIKLMNPQRSTVWY-'; ent=[]
for i in range(6):
    p=np.array([sum(1 for r in msa if r[i]==a)/5 for a in alpha]); p=p[p>0]; ent.append(float(-np.sum(p*np.log(p))))
mi01=0.0
from collections import Counter
j=Counter((r[0],r[1]) for r in msa); pi=Counter(r[0] for r in msa); pj=Counter(r[1] for r in msa)
mi01=sum(v/5*math.log((v/5)/((pi[a]/5)*(pj[b]/5))) for (a,b),v in j.items())
res['msa_couplings'] = {'entropy_match': all(abs(a-b)<1e-12 for a,b in zip(ent, mc['entropy'])),
    'mi01_module': mc['mutual_information'][0][1], 'mi01_recomputed': round(mi01,8), 'effective_depth': mc['effective_depth']}
# 6. refine_coordinates: energy decreases, bonds stay ~3.8
x = cas + np.random.default_rng(1).normal(0, 0.3, cas.shape)
e0 = A.physical_energy(x)['total'] + A._bond_energy(x)
ref = A.refine_coordinates(x, steps=20)
xr = np.array(ref['coordinates'])
e1 = A.physical_energy(xr)['total'] + A._bond_energy(xr)
bl = np.linalg.norm(np.diff(xr,axis=0),axis=1)
res['refine'] = {'energy_before': round(float(e0),4), 'energy_after': round(float(e1),4), 'decreased': e1 < e0,
                 'trajectory_monotone': all(a >= b for a, b in zip(ref['energy_trajectory'], ref['energy_trajectory'][1:])),
                 'max_bond_deviation_A': ref['max_bond_deviation_A'], 'bond_mean_after': round(float(bl.mean()),4)}
fn = os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_alpha_fold_ui.json')
json.dump(res, open(fn,'w'), indent=1, default=str)
print(json.dumps(res, default=str)[:4500])
