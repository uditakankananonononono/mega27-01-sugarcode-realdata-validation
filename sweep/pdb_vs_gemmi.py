"""bio.pdb (PDB and mmCIF parsers) vs gemmi on 30 RCSB entries (the first 30 docking PDB IDs in manifests/accessions.tsv, both formats downloaded from files.rcsb.org 2026-09-25).
Every atom site, all altlocs and models: name, resname, chain (auth), resseq, icode, altloc, xyz, occupancy, B, element."""
import json, sys, glob, os, gemmi
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import pdb as P
def gem(path):
    st=gemmi.read_structure(path); st.setup_entities() if False else None
    out=[]
    for mi,m in enumerate(st):
        for ch in m:
            for r in ch:
                for a in r:
                    out.append((m.num if hasattr(m,'num') else mi+1,a.name,r.name,ch.name,r.seqid.num,(r.seqid.icode.strip() or None),(a.altloc if a.altloc!='\x00' else None),
                                round(a.pos.x,3),round(a.pos.y,3),round(a.pos.z,3),round(a.occ,2),round(a.b_iso,2),a.element.name.upper()))
    return sorted(out,key=lambda t:tuple(str(x) for x in t))
def ours(s):
    return sorted([(a['model'],a['name'],a['resname'],a['chain'],a['resseq'],a['icode'],a['altloc'],round(a['x'],3),round(a['y'],3),round(a['z'],3),
                    round(a['occupancy'] or 0,2),round(a['bfactor'] or 0,2),(a['element'] or '').upper()) for a in s['atoms']],key=lambda t:tuple(str(x) for x in t))
res=[]
for f in sorted(glob.glob('data/pdb_files/*.pdb')):
    pid=os.path.basename(f)[:-4]; row={'id':pid}
    for fmt,path,parse in (('pdb',f,P.parse_pdb),('cif',f[:-4]+'.cif',P.parse_mmcif)):
        o=ours(parse(open(path).read())); g=gem(path)
        so,sg=set(o),set(g)
        row[fmt]={'n_sugarcode':len(o),'n_gemmi':len(g),'identical':o==g,'only_sugarcode':len(so-sg),'only_gemmi':len(sg-so),
                  'ex_sc':[list(map(str,x)) for x in sorted(so-sg)[:2]],'ex_gemmi':[list(map(str,x)) for x in sorted(sg-so)[:2]]}
    res.append(row)
out={'reference':'gemmi %s'%gemmi.__version__,'n_entries':len(res),'pdb_identical':sum(r['pdb']['identical'] for r in res),'cif_identical':sum(r['cif']['identical'] for r in res),
     'atoms_total_pdb':sum(r['pdb']['n_gemmi'] for r in res),'per_entry':res}
json.dump(out,open('benchmarks/sweep_pdb_gemmi.json','w'),indent=1)
print(out['n_entries'],out['pdb_identical'],out['cif_identical'],out['atoms_total_pdb'])
for r in res:
    for fmt in ('pdb','cif'):
        if not r[fmt]['identical']: print(r['id'],fmt,{k:v for k,v in r[fmt].items() if k!='identical'})
