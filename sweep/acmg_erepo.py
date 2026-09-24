"""bayesian_acmg vs ClinGen VCEP classifications (Evidence Repository), with ACMG-2015 rule baseline."""
import csv, collections, json, os, sys
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src')); sys.path.insert(0, os.path.dirname(__file__))
from sugarcode.modules.acmg_bayesian.core import bayesian_acmg
from acmg_rules import acmg2015
csv.field_size_limit(10**8)
SRC = sys.argv[1] if len(sys.argv) > 1 else '/tmp/erepo.tsv'
r = list(csv.reader(open(SRC), delimiter='\t', quoting=csv.QUOTE_NONE))
rows = [x for x in r[1:] if len(x) == 20 and x[17] == 'false']
M = {'Pathogenic': 'Pathogenic', 'Likely Pathogenic': 'Likely pathogenic', 'Uncertain Significance': 'Uncertain significance', 'Likely Benign': 'Likely benign', 'Benign': 'Benign'}
L = {'very strong': 'VS', 'strong': 'S', 'moderate': 'M', 'supporting': 'P', 'stand alone': 'A'}
def counts(codes):
    c = dict(PVS=0, PS=0, PM=0, PP=0, BA=0, BS=0, BP=0)
    for t in codes:
        code, _, st = t.partition('_'); d = code[0]
        base = {'PVS': 'VS', 'PS': 'S', 'PM': 'M', 'PP': 'P', 'BA': 'A', 'BS': 'S', 'BP': 'P'}[code.rstrip('0123456789')]
        lv = L[st.lower()] if st else base
        c[{('P','VS'):'PVS',('P','S'):'PS',('P','M'):'PM',('P','P'):'PP',('B','A'):'BA',('B','S'):'BS',('B','P'):'BP',('B','VS'):'BA',('B','M'):'BS',('P','A'):'PVS'}[(d, lv)]] += 1
    return c
out = []; cb = collections.Counter(); cr = collections.Counter(); rej = 0; panels = collections.Counter()
for x in rows:
    codes = [t.strip() for t in x[9].split(',') if t.strip()]
    o = bayesian_acmg(codes); y = M[x[8]]
    rej += bool(o['rejected_evidence'])
    ry = acmg2015(counts(codes))
    cb[(y, o['classification'])] += 1; cr[(y, ry)] += 1; panels[x[13]] += 1
    out.append([x[1], x[4], x[13], y, x[9], o['points'], o['classification'], ry])
with open(os.path.expanduser('~/mega/m01/data/erepo/erepo_vcep_codes.tsv'), 'w') as f:
    w = csv.writer(f, delimiter='\t', lineterminator='\n')
    w.writerow(['clinvar_variation_id', 'gene', 'expert_panel', 'vcep_assertion', 'met_codes', 'points', 'bayesian_class', 'acmg2015_rules_class']); w.writerows(out)
n = len(rows)
agree_b = sum(v for k, v in cb.items() if k[0] == k[1]); agree_r = sum(v for k, v in cr.items() if k[0] == k[1])
res = {'source': 'https://erepo.clinicalgenome.org/evrepo/api/classifications/all?format=tabbed',
       'note': 'download stopped by the server connection after ~18 MB; only complete, non-retracted rows used',
       'n': n, 'n_expert_panels': len(panels), 'n_genes': len({x[4] for x in rows}),
       'bayesian_agree': agree_b, 'acmg2015_rules_agree': agree_r, 'rows_with_rejected_codes_after_fix': rej,
       'bayesian_confusion': {f'{a} -> {b}': v for (a, b), v in cb.most_common()},
       'rules_confusion': {f'{a} -> {b}': v for (a, b), v in cr.most_common()},
       'single_supporting_benign_vus_to_lb': sum(1 for q in out if q[3] == 'Uncertain significance' and q[6] == 'Likely benign' and q[5] == -1),
       'ba1_with_pathogenic_evidence': dict(collections.Counter(q[3] for q in out if 'BA1' in q[4] and any(t.strip().startswith('P') for t in q[4].split(','))))}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_acmg_erepo.json'), 'w'), indent=1)
print({k: res[k] for k in ('n', 'n_expert_panels', 'n_genes', 'bayesian_agree', 'acmg2015_rules_agree', 'rows_with_rejected_codes_after_fix', 'single_supporting_benign_vus_to_lb', 'ba1_with_pathogenic_evidence')})
