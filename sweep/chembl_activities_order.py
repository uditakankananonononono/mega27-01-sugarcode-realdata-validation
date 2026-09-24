"""bio.chembl.activities_for_target ('best first') vs ChEMBL server-side ordering (order_by=standard_value) over ALL IC50/Ki/Kd nM activities, 12 human targets, live ChEMBL API 2026-09-25."""
import sys, json, urllib.request
sys.path.insert(0, sys.argv[1] if len(sys.argv)>1 else '../sugarcode-ai/src')
from sugarcode.bio import chembl as C
T={'EGFR':'CHEMBL203','ABL1':'CHEMBL1862','BRAF':'CHEMBL5145','JAK2':'CHEMBL2971','FLT3':'CHEMBL1974','ERBB2':'CHEMBL1824','PIK3CA':'CHEMBL4005','IDH1':'CHEMBL2007625','ALK':'CHEMBL4247','KIT':'CHEMBL1936','CDK4':'CHEMBL331','DRD2':'CHEMBL217'}
B='https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={}&standard_type__in=IC50,Ki,Kd&standard_units=nM&standard_value__isnull=false&standard_relation=%3D&order_by=standard_value&limit=10'
BP='https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={}&standard_type__in=IC50,Ki,Kd&standard_units=nM&pchembl_value__isnull=false&order_by=-pchembl_value&limit=1'
rows=[]
for g,t in T.items():
    ours=C.activities_for_target(t, max_n=50)
    ref=json.loads(urllib.request.urlopen(B.format(t),timeout=60).read())
    top=[float(a['standard_value']) for a in ref['activities']]
    bp=json.loads(urllib.request.urlopen(BP.format(t),timeout=60).read())['activities'][0]
    oe=[a['value_nM'] for a in ours if a['relation']=='=']
    rows.append({'gene':g,'target':t,'total_activities':ref['page_meta']['total_count'],'ours_min_nM':min(oe) if oe else None,'global_min_nM':top[0] if top else None,
                 'best_valid_pchembl':float(bp['pchembl_value']),'best_valid_nM':float(bp['standard_value']),'ours_best_pchembl':(float(ours[0]['pchembl']) if ours and ours[0]['pchembl'] else None),'ours_top10_in_global_top10':len({a['molecule_chembl_id'] for a in ours[:10]}&{a['molecule_chembl_id'] for a in ref['activities']})})
out={'reference':'ChEMBL API order_by=standard_value (relation =)','rows':rows,'n_targets':len(rows),'n_ours_min_equals_global':sum(r['ours_min_nM']==r['global_min_nM'] for r in rows),'n_ours_best_equals_best_valid_pchembl':sum(r['ours_best_pchembl']==r['best_valid_pchembl'] for r in rows),'note':'global_min_nM includes ChEMBL data-entry artefacts (0.0 nM, negative values, 5e-9 nM) that carry no pChEMBL; best_valid_* uses only activities with a pChEMBL value'}
json.dump(out,open(sys.argv[2] if len(sys.argv)>2 else 'benchmarks/sweep_chembl_order.json','w'),indent=1)
for r in rows: print(r)
print(out['n_ours_min_equals_global'],out['n_ours_best_equals_best_valid_pchembl'])
