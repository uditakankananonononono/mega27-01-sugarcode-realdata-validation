import freesasa,sys,json,numpy as np
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.structures import parse_pdb_atoms, sasa_shrake_rupley
from scipy.stats import pearsonr
out=[]
for p in ['1gi1','3po6','4b2i']:
    L=[l for l in open(f'/home/sandbox/mega/item01/dock/{p}_rec.pdb') if l.startswith('ATOM') and (l[76:78].strip() or l[12:16].strip()[0])!='H']
    open('/tmp/n.pdb','w').write(''.join(L)+'END\n'); a=parse_pdb_atoms(''.join(L)); r=sasa_shrake_rupley(a)
    s=freesasa.Structure('/tmp/n.pdb'); fr=freesasa.calc(s,freesasa.Parameters({'algorithm':freesasa.ShrakeRupley,'n-points':96}))
    o={'pdb':p,'n_sugarcode':len(a),'n_freesasa':s.nAtoms(),'sugarcode_total':round(r['total_A2'],1),'freesasa_total':round(fr.totalArea(),1)}
    if len(a)==s.nAtoms(): o['per_atom_pearson']=round(float(pearsonr(r['per_atom_A2'],[fr.atomArea(i) for i in range(s.nAtoms())])[0]),4)
    out.append(o); print(o,flush=True)
json.dump({'after_altloc_fix':out},open('sasa_altloc_recheck.json','w'),indent=1)
