"""Validate mutdock against Platinum experimental ligand-binding ddG data.
- Builds real binding pockets from wt PDB structures (residues within 6A of ligand).
- Runs mutation_effect and compares sign + rank vs experimental RT*ln(K_mt/K_wt).
- Checks affinity_change_fold formula and mutation naming."""
import csv, json, os, sys, math, urllib.request
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import numpy as np
from sugarcode.modules.mutdock.core import mutation_effect, DDG_CLASS_CHANGE
D = 'data/platinum'
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
rows = [r for r in csv.DictReader(open(f'{D}/platinum_flat_file.csv'))]
cand = [r for r in rows if r['mut.is_single_point']=='YES' and r['mut.in_binding_site']=='YES'
        and r['affin.k_wt'] not in ('','NR') and r['affin.k_mt'] not in ('','NR')
        and r['mut.wt_pdb'] not in ('','NO') and r['affin.chain']]
def fetch(url, path):
    if not os.path.exists(path):
        with urllib.request.urlopen(url, timeout=60) as r: open(path,'wb').write(r.read())
    return path
def fetch_pdb(pdb_id, path):
    if os.path.exists(path): return path
    for u in (f'https://files.rcsb.org/download/{pdb_id}.pdb',
              f'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdb_id.lower()}.ent'):
        try:
            with urllib.request.urlopen(u, timeout=60) as r: open(path,'wb').write(r.read())
            return path
        except Exception: continue
    raise FileNotFoundError(pdb_id)
def pocket_from_pdb(pdb_path, lig_code, chain, cutoff=6.0):
    lig_atoms=[]; res={}
    for line in open(pdb_path):
        if line.startswith('HETATM') and line[17:20].strip()==lig_code:
            try: lig_atoms.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
            except ValueError: pass
        elif line.startswith('ATOM') and line[21]==chain:
            rn=line[17:20].strip()
            if rn not in AA3: continue
            try: xyz=(float(line[30:38]),float(line[38:46]),float(line[46:54])); num=int(line[22:26])
            except ValueError: continue
            res.setdefault(num, rn); 
            res.setdefault('_xyz',{}).setdefault(num,[]).append(xyz)
    if not lig_atoms or not res: return None
    keep=[]
    xyzmap=res.pop('_xyz')
    for num,rn in res.items():
        for x,y,z in xyzmap[num]:
            if any((x-a)**2+(y-b)**2+(z-c)**2 <= cutoff*cutoff for a,b,c in lig_atoms):
                keep.append((num,rn)); break
    keep.sort()
    return ''.join(AA3[r] for _,r in keep), [n for n,_ in keep]
def ligand_smiles(code):
    p=fetch(f'https://data.rcsb.org/rest/v1/core/chemcomp/{code}', f'{D}/cc_{code}.json')
    d=json.load(open(p))
    for desc in d.get('rcsb_chem_comp_descriptor',{}).keys():
        pass
    smi=(d.get('rcsb_chem_comp_descriptor') or {}).get('SMILES') or (d.get('rcsb_chem_comp_descriptor') or {}).get('smiles')
    return smi
picked=[]; per_prot={}
for r in cand:
    key=r['prot.molecule_name']+r['mut.uniprot']
    if per_prot.get(key,0)>=3: continue
    per_prot[key]=per_prot.get(key,0)+1; picked.append(r)
    if len(picked)>=40: break
res={'module':'mutdock','source':'Platinum flat file (biosig.lab.uq.edu.au, 2026-09-25)','n_candidates':len(cand),'entries':[]}
RT=0.001987*298.15
for r in picked:
    e={'mutation':r['mutation'],'pdb':r['mut.wt_pdb'],'lig':r['affin.lig_id'],'protein':r['prot.molecule_name']}
    try:
        pdb=fetch_pdb(r['mut.wt_pdb'], f'{D}/{r["mut.wt_pdb"]}.pdb')
        pk=pocket_from_pdb(pdb, r['affin.lig_id'].strip(), r['affin.chain'])
        if not pk: e['skip']='no pocket parsed'; res['entries'].append(e); continue
        seq,resnums=pk
        wt,mt=r['mutation'][0],r['mutation'][-1]; num=int(r['mutation'][1:-1])
        if num not in resnums or seq[resnums.index(num)]!=wt:
            e['skip']=f'mutated residue {r["mutation"]} not in parsed pocket (pocket resnums {resnums[0]}..{resnums[-1]})'
            res['entries'].append(e); continue
        smi=ligand_smiles(r['affin.lig_id'].strip())
        if not smi: e['skip']='no SMILES'; res['entries'].append(e); continue
        out=mutation_effect(seq, smi, resnums.index(num), mt, resnums=resnums, drug_name=r['affin.lig_name'][:30])
        kwt=float(r['affin.k_wt']); kmt=float(r['affin.k_mt'])
        exp_ddg=RT*math.log(kmt/kwt)
        e.update({'exp_ddg_kcal_mol':round(exp_ddg,3),'module_ddg':out['ddg_kcal_mol'],
                  'module_mutation_str':out['mutation'],'module_fold':out['affinity_change_fold'],
                  'exp_fold':round(kmt/kwt,2),'module_risk':out['resistance_risk'],
                  'sign_match':(exp_ddg>0)==(out['ddg_kcal_mol']>0),'pocket_len':len(seq)})
    except Exception as ex:
        e['skip']=f'{type(ex).__name__}: {ex}'
    res['entries'].append(e)
sc=[e for e in res['entries'] if 'exp_ddg_kcal_mol' in e]
res['scored']=len(sc)
res['sign_accuracy']=round(sum(e['sign_match'] for e in sc)/len(sc),3) if sc else None
if len(sc)>=8:
    from scipy.stats import spearmanr
    res['spearman_exp_vs_module']=round(float(spearmanr([e['exp_ddg_kcal_mol'] for e in sc],[e['module_ddg'] for e in sc]).correlation),3)
# formula checks (independent recomputation)
chk=[]
for e in sc:
    chk.append(abs(e['module_fold']-round(2.718**(e['module_ddg']/0.593),2))<0.011)
res['fold_formula_ok']=all(chk) if chk else None
res['mutation_str_ok']=all(e['module_mutation_str']==e['mutation'] for e in sc) if sc else None
# dataset-level audit of the class-penalty table: all single-point Platinum entries
def _cls(aa):
    for k,v in {'hydrophobic':set('AILMFWVPG'),'positive':set('KRH'),'negative':set('DE'),'polar':set('STNQCY')}.items():
        if aa in v: return k
    return 'polar'
agg={}
for r in rows:
    if r['mut.is_single_point']!='YES' or r['affin.k_wt'] in ('','NR') or r['affin.k_mt'] in ('','NR'): continue
    m=r['mutation']
    if len(m)<3 or not m[0].isalpha() or not m[-1].isalpha(): continue
    pair=(_cls(m[0]),_cls(m[-1]))
    try: dd=RT*math.log(float(r['affin.k_mt'])/float(r['affin.k_wt']))
    except Exception: continue
    agg.setdefault(pair,[]).append(dd)
table_audit={}
for pair,vals in sorted(agg.items()):
    tab=DDG_CLASS_CHANGE.get(pair)
    table_audit[f'{pair[0]}->{pair[1]}']={'n':len(vals),'mean_exp_ddg':round(sum(vals)/len(vals),3),
        'median_exp_ddg':round(sorted(vals)[len(vals)//2],3),'table_value':tab,
        'sign_agree':None if tab is None else (sum(vals)/len(vals)>0)==(tab>0)}
res['class_table_audit']={'pairs':table_audit,
    'sign_agreement_share':round(sum(1 for v in table_audit.values() if v['sign_agree'])/len(table_audit),3),
    'n_entries':sum(v['n'] for v in table_audit.values())}
json.dump(res, open('benchmarks/sweep_mutdock.json','w'), indent=1)
print(json.dumps({k:v for k,v in res.items() if k!='entries'}, indent=1))
for e in res['entries']: print(e.get('mutation'), e.get('skip') or (e['exp_ddg_kcal_mol'], e['module_ddg'], e['sign_match']))
