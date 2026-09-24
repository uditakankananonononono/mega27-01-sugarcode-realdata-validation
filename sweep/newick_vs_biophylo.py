"""bio.newick vs Biopython Bio.Phylo and DendroPy on 6 Open Tree of Life synthetic subtrees (api.opentreeoflife.org v3 tree_of_life/subtree, height_limit 7, 2026-09-25).
The height limit truncates clades, so the source trees contain empty leaf labels, e.g. "(,Conepatus_ott231614)"; both parsers must count those as unnamed leaves.
Checks leaf set, internal-node count, internal labels, and write->parse round trip."""
import json, sys, glob, io, Bio, dendropy
from Bio import Phylo
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import newick as nw
res=[]
for f in sorted(glob.glob('data/opentree/*.nwk')):
    t=open(f).read().strip(); a=nw.parse_newick(t); s=nw.stats(a)
    b=Phylo.read(io.StringIO(t),'newick')
    bl=sorted(x.name for x in b.get_terminals() if x.name); bi=len(b.get_nonterminals())
    bint=sorted(x.name for x in b.get_nonterminals() if x.name)
    d=dendropy.Tree.get(data=t,schema='newick',preserve_underscores=True,suppress_internal_node_taxa=True,suppress_leaf_node_taxa=False)
    dl=sorted(n.taxon.label for n in d.leaf_node_iter() if n.taxon is not None and n.taxon.label)
    sint=sorted(n['name'] for n in nw.preorder(a) if n['children'] and n['name'])
    rt=nw.parse_newick(nw.write_newick(a))==a
    res.append({'file':f,'root':a['name'],'leaves':s['leaves'],'unnamed_leaves':s['leaves']-len(s['leaf_names']),'unnamed_leaves_biopython':sum(1 for x in b.get_terminals() if not x.name),'internal':s['internal_nodes'],'binary':s['binary'],
      'leaves_eq_biopython':s['leaf_names']==bl,'internal_count_eq_biopython':s['internal_nodes']==bi,'internal_labels_eq_biopython':sint==bint,
      'leaves_eq_dendropy':s['leaf_names']==dl,'roundtrip':rt})
out={'reference':'Biopython %s Bio.Phylo; DendroPy %s'%(Bio.__version__,dendropy.__version__),'n_trees':len(res),'total_leaves':sum(r['leaves'] for r in res),
 'n_all_pass':sum(all(r[k] for k in r if k.endswith(('biopython','dendropy','roundtrip'))) for r in res),'per_tree':res}
json.dump(out,open('benchmarks/sweep_newick_biophylo.json','w'),indent=1)
print(out['n_trees'],out['n_all_pass'],out['total_leaves'])
for r in res: print({k:v for k,v in r.items() if k!='file'})
