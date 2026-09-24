"""bio.proteinprops vs Biopython ProtParam on 500 UniProt reviewed human proteins; bio.restriction vs Biopython Restriction on pBR322 (J01749.1) and lambda (NC_001416.1)."""
import json, sys, numpy as np, Bio
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Restriction import RestrictionBatch, Analysis
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import proteinprops as pp, restriction as rs
seqs=[(r.id,str(r.seq)) for r in SeqIO.parse('data/uniprot/human_reviewed_500.fasta','fasta') if set(str(r.seq))<=set('ACDEFGHIKLMNPQRSTVWY')]
M={k:[] for k in ['mw','gravy','aromaticity','instability','pI','ext_reduced','ext_oxidized']}
for i,s in seqs:
    b=ProteinAnalysis(s)
    M['mw'].append((pp.molecular_weight(s),b.molecular_weight()))
    M['gravy'].append((pp.gravy(s),b.gravy())); M['aromaticity'].append((pp.aromaticity(s),b.aromaticity()))
    M['instability'].append((pp.instability_index(s),b.instability_index())); M['pI'].append((pp.isoelectric_point(s),b.isoelectric_point()))
    e=pp.extinction_coefficient(s); r,o=b.molar_extinction_coefficient()
    M['ext_reduced'].append((e['reduced'],r)); M['ext_oxidized'].append((e['oxidized'],o))
prot={}
for k,v in M.items():
    a=np.array(v); d=np.abs(a[:,0]-a[:,1]); rel=d/np.maximum(np.abs(a[:,1]),1e-9)
    tol=0.01 if k=='pI' else 1e-6
    prot[k]={'n':len(a),'max_abs_diff':float(d.max()),'max_rel_diff':float(rel.max()),'n_within_tol':int((d<=tol) .sum() if k=='pI' else (rel<=tol).sum()),'tol':('abs %.2f'%tol if k=='pI' else 'rel %g'%tol),'pearson':float(np.corrcoef(a[:,0],a[:,1])[0,1])}
enz=[e for e in rs.list_enzymes()]
rb=RestrictionBatch([]); common=[]
for e in enz:
    try: rb.add(e); common.append(e)
    except Exception: pass
restr={}
from Bio.Seq import Seq
for acc,circ in (('J01749.1',True),('NC_001416.1',False)):
    s=str(next(SeqIO.parse(f'data/refseq/{acc}.fa','fasta')).seq).upper()
    bio=Analysis(rb,Seq(s),linear=not circ).full()
    agree=0; dis=[]
    for e in common:
        ours=sorted({h['cut_top']+1 for h in rs.find_sites(s,e,circular=circ)})  # Biopython reports 1-based position of first base after the cut
        theirs=sorted(set(bio[rb.get(e)]))
        if [x if x>0 else x+len(s) for x in ours]==theirs or ours==theirs: agree+=1
        else: dis.append({'enzyme':e,'sugarcode':ours[:8],'biopython':theirs[:8],'n_sc':len(ours),'n_bio':len(theirs)})
    restr[acc]={'length':len(s),'circular':circ,'n_enzymes':len(common),'agree':agree,'disagree':dis[:15],'n_disagree':len(dis)}
out={'reference':'Biopython %s ProtParam + Restriction'%Bio.__version__,'proteins':'UniProtKB reviewed human, length 50-2000, first 500 of search (rest.uniprot.org, 2026-09-25)','n_proteins':len(seqs),'proteinprops':prot,'restriction':restr,'enzymes_not_in_biopython':[e for e in enz if e not in common]}
json.dump(out,open('benchmarks/sweep_proteinprops_restriction.json','w'),indent=1)
for k,v in prot.items(): print(k,v)
for k,v in restr.items(): print(k,{x:v[x] for x in ['n_enzymes','agree','n_disagree']},v['disagree'][:3])
