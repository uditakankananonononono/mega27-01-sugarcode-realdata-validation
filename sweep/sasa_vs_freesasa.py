"""sugarcode bio.structures.sasa_shrake_rupley vs FreeSASA (Shrake-Rupley, 96 pts, 1.4 A probe) on real PDB receptors."""
import freesasa, sys, glob, json, numpy as np
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.structures import parse_pdb_atoms, sasa_shrake_rupley
from scipy.stats import pearsonr
files=sorted(glob.glob('/home/sandbox/mega/item01/dock/*_rec.pdb'))[:30]; out=[]
for f in files:
    if sum(1 for l in open(f) if l.startswith('ATOM'))>6000: continue
    txt=open(f).read(); lines=[l for l in txt.splitlines() if l.startswith('ATOM') and (l[76:78].strip() or l[12:16].strip()[0])!='H']
    open('/tmp/noH.pdb','w').write('\n'.join(lines)+'\nEND\n')
    a=parse_pdb_atoms('\n'.join(lines)); r=sasa_shrake_rupley(a)
    s=freesasa.Structure('/tmp/noH.pdb'); fr=freesasa.calc(s,freesasa.Parameters({'algorithm':freesasa.ShrakeRupley,'n-points':96,'probe-radius':1.4}))
    if s.nAtoms()!=len(a): out.append({'pdb':f.split('/')[-1][:4],'error':f'atom count {len(a)} vs {s.nAtoms()}'}); continue
    fa=np.array([fr.atomArea(i) for i in range(s.nAtoms())]); sa=np.array(r['per_atom_A2'])
    out.append({'pdb':f.split('/')[-1][:4],'n_atoms':len(a),'sugarcode_total':round(r['total_A2'],1),'freesasa_total':round(fr.totalArea(),1),
      'rel_diff':round((r['total_A2']-fr.totalArea())/fr.totalArea(),4),'per_atom_pearson':round(float(pearsonr(sa,fa)[0]),4)})
    print(out[-1],flush=True)
ok=[o for o in out if 'error' not in o]
json.dump({'tool':'FreeSASA %s (Shrake-Rupley, 96 points, probe 1.4 A, default ProtOr radii)'%getattr(freesasa,'__version__','2.x'),'structures':'RCSB PDB receptors from the LP-PDBBind docking set (hydrogens removed)',
 'n':len(ok),'median_rel_diff':float(np.median([o['rel_diff'] for o in ok])),'max_abs_rel_diff':float(max(abs(o['rel_diff']) for o in ok)),
 'median_per_atom_pearson':float(np.median([o['per_atom_pearson'] for o in ok])),'per_structure':out},open('sweep_sasa_freesasa.json','w'),indent=1)
