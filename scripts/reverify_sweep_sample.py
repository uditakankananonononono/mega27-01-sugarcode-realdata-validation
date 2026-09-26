#!/usr/bin/env python3
"""Live re-verification of a sample of validated-module claims (2026-09-26).
Re-executes three claims end-to-end against independent references:
  1. cfd_offtarget  vs Doench 2016 official CFD pickles (600 pairs x 16 PAMs)
  2. chem_descriptors vs RDKit 2024.09.6 (993 ChEMBL approved drugs, 6 descriptors)
  3. pwm/motif log-odds vs Biopython 1.88 (10 JASPAR motifs)
Results: verification/sweep_reverify_20260926.json"""
import sys, csv, json, pickle
sys.path.insert(0, '/home/sandbox/sugarcode-ai/src')
ROOT = '/home/sandbox/mega27-01-sugarcode-realdata-validation'
ORACLE = '/home/sandbox/wave1/oracle'

def check_cfd():
    from sugarcode.modules.cfd_offtarget import core as cfd
    mm = pickle.load(open(f'{ORACLE}/mismatch_score.pkl','rb'))
    pam = pickle.load(open(f'{ORACLE}/pam_scores.pkl','rb'))
    comp = {'A':'T','C':'G','G':'C','T':'A','U':'A'}
    def ref(wt, sg, p):
        s = 1.0
        sg = sg.upper().replace('T','U'); wt = wt.upper().replace('T','U')
        for i in range(len(sg)):
            if wt[i] != sg[i]:
                k = 'r%s:d%s,%d' % (wt[i], comp[sg[i]], i+1)
                if k not in mm: return None
                s *= mm[k]
        return s * pam[p]
    PAMS = ['GG','AG','TG','CG','GA','GC','GT','AT','AA','AC','TA','TC','CA','CT','CC','TT']
    rows = list(csv.DictReader(open(f'{ROOT}/data/crispr_cfd_oracle.tsv'), delimiter='\t'))
    n = exact = 0; mx = 0.0
    for r in rows:
        g, ot = r['guide'][:20], r['offtarget'][:20]
        for p in PAMS:
            rv = ref(g, ot, p)
            if rv is None: continue
            d = abs(cfd.score(g, ot, pam=p).score - rv)
            n += 1; exact += d < 1e-12; mx = max(mx, d)
    return {'n_scores_compared': n, 'exact_1e-12': exact, 'max_abs_diff': mx}

def check_chem():
    from sugarcode.modules.chem_descriptors import descriptors as d
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    rows = list(csv.DictReader(open(f'{ROOT}/data/chembl_approved_1000.tsv'), delimiter='\t'))
    keys = ['mol_wt','exact_mol_wt','hbd','hba','rotatable_bonds','tpsa']
    agree = {k: [0,0] for k in keys}
    for r in rows:
        rm = Chem.MolFromSmiles(r['smiles'])
        if rm is None: continue
        desc = d.compute_descriptors(d.parse_smiles(r['smiles']))
        refv = {'mol_wt': Descriptors.MolWt(rm), 'exact_mol_wt': Descriptors.ExactMolWt(rm),
                'hbd': Lipinski.NumHDonors(rm), 'hba': Lipinski.NumHAcceptors(rm),
                'rotatable_bonds': Lipinski.NumRotatableBonds(rm), 'tpsa': rdMolDescriptors.CalcTPSA(rm)}
        for k in keys:
            agree[k][0] += abs(desc[k]-refv[k]) < 1e-6; agree[k][1] += 1
    return {'n': len(rows), 'agreement': {k: {'agree': a, 'n': b, 'rate': round(a/b,4)} for k,(a,b) in agree.items()}}

def check_motif():
    from Bio import motifs
    from sugarcode.bio import pwm as sc_pwm
    claim = json.load(open(f'{ROOT}/benchmarks/sweep_motif_biopython.json'))
    names = [m['motif'].split('\t')[0] for m in claim['per_motif']]
    out = []
    for name in names:
        with open(f'{ROOT}/data/jaspar/{name}.jaspar') as fh:
            mot = motifs.read(fh, 'jaspar')
        L = mot.length
        pssm = mot.counts.normalize(0.5).log_odds()
        pwm_list = [{b: (mot.counts[b][j]+0.5)/(sum(mot.counts[x][j] for x in 'ACGT')+2.0) for b in 'ACGT'} for j in range(L)]
        sc = sc_pwm.log_odds_matrix(pwm_list, 0.25)
        mx = max(abs(pssm[b][j]-sc[j][b]) for b in 'ACGT' for j in range(L))
        out.append({'motif': name, 'len': L, 'max_lod_diff': mx})
    return out

if __name__ == '__main__':
    print(json.dumps({'cfd': check_cfd(), 'chem': check_chem(), 'motif': check_motif()}, indent=1))
