"""sections_inventory.tex from manifests/module_inventory.tsv."""
import csv
R=list(csv.DictReader(open('../manifests/module_inventory.tsv'),delimiter='\t'))
esc=lambda s:s.replace('_','\\_')
v=[r for r in R if r['status']=='validated']
o=['\\section{Module inventory and validation coverage}',
f"sugarcode-ai contains {sum(r['kind']=='module' for r in R)} application modules and {sum(r['kind']=='bio' for r in R)} bio/ libraries. {len(v)} of these {len(R)} have been checked against an independent reference on real data in this work ({sum(r['kind']=='module' for r in v)} modules, {sum(r['kind']=='bio' for r in v)} libraries). The rest are not yet validated. That does not mean they are wrong, only that no claim about them is made here. Validation status is a curated list tied to committed result files (manifests/module\\_inventory.tsv, built by sweep/module\\_inventory.py), not a keyword match.",
'\\begin{center}\\scriptsize\\begin{longtable}{llrl}\\hline kind & name & py lines & status\\\\\\hline\\endhead']
for r in sorted(R,key=lambda r:(r['status']!='validated',r['kind'],r['name'])):
    st='\\textbf{validated}' if r['status']=='validated' else 'not yet'
    o.append(f"{r['kind']} & {esc(r['name'])} & {r['python_lines']} & {st}"+'\\\\')
o.append('\\hline\\end{longtable}\\end{center}')
open('sections_inventory.tex','w').write('\n'.join(o)+'\n')
