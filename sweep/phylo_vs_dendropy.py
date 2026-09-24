"""sugarcode bio.phylo NJ/UPGMA vs DendroPy (nj_tree, upgma_tree) on identical K80
distance matrices from Rfam seed alignments; plus p-distance vs Biopython identity."""
import sys, json, io, random, numpy as np, dendropy
from dendropy.calculate import treecompare
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio import phylo
from sugarcode.bio.newick import write_newick
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
res={}
for fam in ['RF00001','RF00005','RF00010']:
    aln=AlignIO.read(f'/home/sandbox/mega/m01/data/rfam/{fam}.sto','stockholm')
    recs=[r for r in aln]; random.Random(1).shuffle(recs); recs=recs[:30]
    seqs={f"s{i}":str(r.seq).upper().replace('U','T').replace('.','-') for i,r in enumerate(recs)}
    names=list(seqs)
    dm=phylo.distance_matrix(seqs,model='k80'); M=dm['matrix'] if isinstance(dm,dict) and 'matrix' in dm else dm
    if isinstance(M,dict): M=[[M[a][b] for b in names] for a in names]
    M=[[float(x) for x in row] for row in M]
    model='k80'
    if any(not np.isfinite(x) for row in M for x in row):
        model='pdistance (K80 saturated: some pairs have infinite distance)'
        dm=phylo.distance_matrix(seqs,model='pdistance'); M=dm['matrix'] if isinstance(dm,dict) and 'matrix' in dm else dm
        if isinstance(M,dict): M=[[M[a][b] for b in names] for a in names]
        M=[[float(x) for x in row] for row in M]
    csv_=io.StringIO(); csv_.write(','+','.join(names)+'\n')
    for a,row in zip(names,M): csv_.write(a+','+','.join(repr(x) for x in row)+'\n')
    tns=dendropy.TaxonNamespace()
    pdm=dendropy.PhylogeneticDistanceMatrix.from_csv(io.StringIO(csv_.getvalue()),taxon_namespace=tns,delimiter=',')
    out={'distance_model':model,'n_seqs':len(names),'aln_len':aln.get_alignment_length()}
    for meth,dfun in (('nj',pdm.nj_tree),('upgma',pdm.upgma_tree)):
        ours=phylo.build_tree(names,M,method=meth); nw=write_newick(ours)
        t1=dendropy.Tree.get(data=nw.rstrip(';')+';',schema='newick',taxon_namespace=tns,rooting='force-unrooted'); t2=dfun(); t2.is_rooted=False
        t1.encode_bipartitions(); t2.encode_bipartitions()
        rf=treecompare.symmetric_difference(t1,t2)
        p1=t1.phylogenetic_distance_matrix(); p2=t2.phylogenetic_distance_matrix(); tx=list(tns)
        a=[p1.distance(tx[i],tx[j]) for i in range(len(tx)) for j in range(i+1,len(tx))]; b=[p2.distance(tx[i],tx[j]) for i in range(len(tx)) for j in range(i+1,len(tx))]
        out[meth]={'rf_symmetric_difference':int(rf),'patristic_max_abs_diff':round(float(np.max(np.abs(np.array(a)-np.array(b)))),6),'patristic_pearson':round(float(np.corrcoef(a,b)[0,1]),6)}
    # identity vs Biopython (p-distance on gap-free columns differs by gap handling; report both)
    pd_=phylo.distance_matrix(seqs,model='pdistance'); P=pd_['matrix'] if isinstance(pd_,dict) and 'matrix' in pd_ else pd_
    if isinstance(P,dict): P=[[P[a][b] for b in names] for a in names]
    from Bio.Align import MultipleSeqAlignment; from Bio.SeqRecord import SeqRecord; from Bio.Seq import Seq
    bm=DistanceCalculator('identity').get_distance(MultipleSeqAlignment([SeqRecord(Seq(seqs[n]),id=n) for n in names]))
    x=[float(P[i][j]) for i in range(len(names)) for j in range(i+1,len(names))]; y=[bm[names[i],names[j]] for i in range(len(names)) for j in range(i+1,len(names))]
    out['pdistance_note']='sugarcode uses pairwise deletion of gap columns; Biopython identity counts gaps'; out['pdistance_vs_biopython_identity_pearson']=round(float(np.corrcoef(x,y)[0,1]),4)
    res[fam]=out; print(fam,out,flush=True)
json.dump({'tools':'DendroPy %s (PhylogeneticDistanceMatrix.nj_tree/upgma_tree, treecompare.symmetric_difference); Biopython DistanceCalculator'%dendropy.__version__,'data':'Rfam seed alignments RF00001/RF00005/RF00010, 30 random sequences each, K80 distances (p-distance where K80 saturates)','results':res},open('sweep_phylo_dendropy.json','w'),indent=1)
