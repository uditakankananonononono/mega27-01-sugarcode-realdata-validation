"""bio.gnomad.gene_constraint on gnomAD r4 (live GraphQL API, 2026-09-25) for the 60 HGNC genes in data/hgnc.
The LOEUF cutoff is compared under the old v2-era rule (0.35) and gnomAD's v4 guidance (0.45; https://gnomad.broadinstitute.org/help/constraint).
Secondary reference: pLI > 0.9 on the same records."""
import sys, glob, json, os, time
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import gnomad as G
rows=[]; T0=time.time()
for f in sorted(glob.glob('data/hgnc/*.json')):
    g=os.path.basename(f)[:-5]
    try: c=G.gene_constraint(g, offline=True)
    except G.GnomADError:
        if time.time()-T0>75: print('partial; rerun later (gnomAD rate limit)'); sys.exit(3)
        time.sleep(6)
        try: c=G.gene_constraint(g)
        except G.GnomADError as e: print('rate-limited at',g,str(e)[:80]); sys.exit(3)
    if 'loeuf' not in c or c['loeuf'] is None: rows.append({'gene':g,'status':c.get('status')}); continue
    rows.append({'gene':g,'loeuf':c['loeuf'],'pli':c['pli'],'old_035':c['loeuf']<0.35,'new_045':c['lof_constrained'],'pli_gt_0.9':(c['pli'] or 0)>0.9})
sc=[r for r in rows if 'loeuf' in r]
out={'source':'gnomAD r4 GraphQL gene.gnomad_constraint','guidance':'https://gnomad.broadinstitute.org/help/constraint: LOEUF < 0.45 suggested for Mendelian interpretation (15th percentile of 17,063 MANE Select); v2 thresholds do not transfer to v4',
 'n_genes':len(rows),'n_scored':len(sc),'constrained_old_035':sum(r['old_035'] for r in sc),'constrained_new_045':sum(r['new_045'] for r in sc),'pli_gt_0.9':sum(r['pli_gt_0.9'] for r in sc),
 'agree_pli_old':sum(r['old_035']==r['pli_gt_0.9'] for r in sc),'agree_pli_new':sum(r['new_045']==r['pli_gt_0.9'] for r in sc),'flipped_to_constrained':[r['gene'] for r in sc if r['new_045'] and not r['old_035']],'rows':rows}
json.dump(out,open('benchmarks/sweep_gnomad_constraint.json','w'),indent=1)
print({k:v for k,v in out.items() if k not in('rows','guidance')})
