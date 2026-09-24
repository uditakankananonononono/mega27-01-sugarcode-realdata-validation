import sys,random; sys.path.insert(0,'sweep'); import importlib.util
spec=importlib.util.spec_from_file_location('rv','sweep/rna_vs_vienna.py')
src=open('sweep/rna_vs_vienna.py').read(); src=src[:src.index('res = {}')]
exec(src)
random.seed(7); R={}
for fam in ("RF00005","RF00001","RF00010"):
    seqs, ss = read_sto(ROOT/"data"/"rfam"/f"{fam}.sto")
    items=[project(s,ss) for s in seqs.values()]; items=[(s,r) for s,r in items if set(s)<=set("ACGU") and len(s)<=400 and len(r)>=5]
    items=random.sample(items,min(25,len(items)))
    out={k:[] for k in ('mfe','cent','mea1','mea4')}
    for s,ref in items:
        fc=RNA.fold_compound(s); mfe,_=fc.mfe(); fc.exp_params_rescale(_); fc.pf()
        out['mfe'].append(f1(db_pairs(mfe),ref)); out['cent'].append(f1(db_pairs(fc.centroid()[0]),ref))
        out['mea1'].append(f1(db_pairs(fc.MEA(1.0)[0]),ref)); out['mea4'].append(f1(db_pairs(fc.MEA(4.0)[0]),ref))
    R[fam]={k:round(sum(v)/len(v),4) for k,v in out.items()}; print(fam,R[fam],flush=True)
json.dump(R,open(str(ROOT/"benchmarks"/"sweep_rna_vienna_methods.json"),"w"),indent=1)
