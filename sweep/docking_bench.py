import csv,json,os,sys,subprocess,random,math,urllib.request,numpy as np
from rdkit import Chem
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.modules.docking_studio.vina import vina_score_pose
from vina import Vina
rows=[r for r in csv.DictReader(open('lp.csv')) if r['new_split']=='test' and r['category']=='refined' and r['covalent']=='False' and r['value']]
random.seed(7); random.shuffle(rows); rows=rows[:int(sys.argv[1])]
SKIP={'HOH','NA','CL','K','MG','CA','ZN','MN','SO4','PO4','GOL','EDO','PEG','ACT','DMS','FMT','CD','NI','CO','CU','FE','IOD','BR','NO3','TRS','MES','EPE'}
out=open('dock_results.jsonl','a'); done={json.loads(l)['pdb'] for l in open('dock_results.jsonl')} if os.path.exists('dock_results.jsonl') else set()
for r in rows:
    p=r['']; 
    if p in done: continue
    try:
        m=Chem.MolFromSmiles(r['smiles']); nh=m.GetNumHeavyAtoms()
        if not os.path.exists(p+'.pdb'): urllib.request.urlretrieve(f'https://files.rcsb.org/download/{p}.pdb',p+'.pdb')
        L=open(p+'.pdb').read().splitlines()
        res={}
        for l in L:
            if l.startswith('HETATM') and l[17:20].strip() not in SKIP and l[76:78].strip()!='H':
                res.setdefault((l[17:20],l[21],l[22:27]),[]).append(l)
        cand=sorted(res.items(),key=lambda kv:abs(len(kv[1])-nh))
        if not cand or abs(len(cand[0][1])-nh)>3: raise Exception('no ligand match')
        lig=cand[0][1]
        open(p+'_lig.pdb','w').write('\n'.join(lig)+'\nEND\n')
        open(p+'_rec.pdb','w').write('\n'.join(l for l in L if l.startswith('ATOM'))+'\nEND\n')
        subprocess.run(['obabel',p+'_rec.pdb','-xr','-h','-p','7.4','-O',p+'_rec.pdbqt'],capture_output=True,timeout=120)
        subprocess.run(['obabel',p+'_lig.pdb','-h','-p','7.4','-O',p+'_lig.pdbqt'],capture_output=True,timeout=60)
        lx=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in lig]); c=lx.mean(0)
        v=Vina(sf_name='vina',verbosity=0); v.set_receptor(p+'_rec.pdbqt'); v.set_ligand_from_file(p+'_lig.pdbqt')
        v.compute_vina_maps(center=c.tolist(),box_size=[24,24,24])
        crystal=v.score()[0]; mini=v.optimize()[0]
        # heavy-atom coords in pdbqt order
        def heavy(fn):
            return np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in open(fn) if l.startswith(('ATOM','HETATM')) and l[77:79].strip() not in ('H','HD')])
        ref=heavy(p+'_lig.pdbqt')
        v.dock(exhaustiveness=8,n_poses=9); v.write_poses(p+'_out.pdbqt',n_poses=9,overwrite=True)
        poses=open(p+'_out.pdbqt').read().split('ENDMDL'); rmsds=[]; es=v.energies(n_poses=9)[:,0].tolist()
        for blk in poses:
            xs=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in blk.splitlines() if l.startswith(('ATOM','HETATM')) and l[77:79].strip() not in ('H','HD')])
            if len(xs)==len(ref): rmsds.append(float(np.sqrt(((xs-ref)**2).sum(1).mean())))
        # sugarcode scorer on crystal pose
        el=lambda l:(l[76:78].strip() or l[12:14].strip()).capitalize()
        la=[{'element':el(l).upper() if len(el(l))==1 else el(l),'xyz':tuple(x)} for l,x in zip(lig,lx)]
        pa=[]
        for l in L:
            if l.startswith('ATOM') and l[76:78].strip()!='H':
                x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
                if np.linalg.norm(x-c)<14: pa.append({'element':l[76:78].strip(),'xyz':tuple(x)})
        t=vina_score_pose(la,pa); W={"gauss1":-0.0356,"gauss2":-0.00516,"repulsion":0.840,"hydrophobic":-0.0351,"hbond":-0.587}
        sc=sum(W[k]*t[k] for k in W)
        rec={'pdb':p,'pK':float(r['value']),'n_heavy':nh,'vina_crystal':crystal,'vina_min':mini,'vina_top_E':es[0],'rmsd_top1':rmsds[0] if rmsds else None,'rmsd_best_of_9':min(rmsds) if rmsds else None,'sugarcode_crystal':sc,'sugarcode_terms':t}
    except Exception as e:
        rec={'pdb':p,'error':str(e)[:150]}
    out.write(json.dumps(rec)+'\n'); out.flush(); print(rec.get('pdb'),rec.get('rmsd_top1'),rec.get('error',''),flush=True)
