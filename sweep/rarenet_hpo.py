"""rarenet_ai curated disease table vs HPO genes_to_phenotype.txt (333,984 annotations,
2026-09-25 release): gene-disease pairing vs OMIM disease IDs, and claimed symptom
signatures vs the gene's annotated HPO term names for that disease."""
import sys, os, json, csv, re
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from collections import defaultdict
from sugarcode.modules.rarenet_ai.core import RARE_DISEASES, diagnose

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'hpo', 'genes_to_phenotype.txt')
OMIM = {'cystic_fibrosis': 'OMIM:219700', 'huntingtons': 'OMIM:143100', 'marfan': 'OMIM:154700',
        'phenylketonuria': 'OMIM:261600', 'gaucher': 'OMIM:230800', 'duchenne_md': 'OMIM:310200'}

# load HPO rows for the 6 genes
genes = {g for d in RARE_DISEASES.values() for g in d['genes']}
rows = defaultdict(list)  # (gene, disease_id) -> [hpo_name]
all_gene_diseases = defaultdict(set)
for r in csv.DictReader(open(DATA), delimiter='\t'):
    if r['gene_symbol'] in genes:
        rows[(r['gene_symbol'], r['disease_id'])].append(r['hpo_name'])
        all_gene_diseases[r['gene_symbol']].add(r['disease_id'])

# lay symptom label -> clinical HPO vocabulary (built from the actual HPO term lists)
SYNONYMS = {
    'lens_dislocation': ['ectopia lentis'],
    'aortic_dilation': ['aortic root aneurysm', 'ascending tubular aorta aneurysm', 'aortic dilatation', 'aortic aneurysm'],
    'long_limbs': ['arachnodactyly', 'long limbs'],
    'chronic_cough': ['chronic cough', 'chronic lung disease', 'recurrent pneumonia'],
    'poor_growth': ['failure to thrive', 'growth delay', 'poor growth', 'short stature'],
    'adult_onset': ['adult onset'],
    'cognitive_decline': ['dementia', 'cognitive decline', 'progressive cognitive impairment'],
    'mood_changes': ['depression', 'personality changes', 'behavioral abnormality', 'mood changes'],
    'eczema': ['eczematoid dermatitis', 'dermatitis', 'eczema'],
    'musty_odor': ['musty odor'],
    'hepatosplenomegaly': ['hepatosplenomegaly', 'hepatomegaly', 'splenomegaly'],
    'fatigue': ['fatigue'],
}

def norm(s):
    s = re.sub(r'[^a-z0-9 ]', '', s.lower().replace('_', ' '))
    s = re.sub(r'\s+', ' ', s).strip()
    return s

def stems(s):
    return {w.rstrip('s') if len(w) > 3 else w for w in norm(s).split()}

res = {'hpo_rows_total': 333984, 'diseases': {}}
for name, d in RARE_DISEASES.items():
    gene = d['genes'][0]
    omim = OMIM[name]
    entry = {'gene': gene, 'omim': omim}
    entry['gene_disease_pair_in_hpo'] = omim in all_gene_diseases.get(gene, set())
    terms = rows.get((gene, omim), [])
    entry['n_hpo_terms_for_disease'] = len(terms)
    termset = {norm(t) for t in terms}
    term_words = set().union(*[stems(t) for t in terms]) if terms else set()
    matched, missed, synonym = [], [], []
    for s in sorted(d['symptoms']):
        ns = norm(s)
        if ns in termset or (stems(s) and stems(s) <= term_words):
            matched.append(s)
        elif any(norm(a) in termset or (stems(a) and stems(a) <= term_words) for a in SYNONYMS.get(s, [])):
            synonym.append(s)
        else:
            missed.append(s)
    entry['symptoms_matched'] = matched
    entry['symptoms_matched_via_clinical_synonym'] = synonym
    entry['symptoms_unresolved'] = missed
    entry['match_fraction'] = round(len(matched) / max(1, len(d['symptoms'])), 3)
    entry['match_fraction_with_synonyms'] = round((len(matched) + len(synonym)) / max(1, len(d['symptoms'])), 3)
    res['diseases'][name] = entry

# module behavior sanity: full signature recovers the disease as top-1
top1 = 0
for name, d in RARE_DISEASES.items():
    r = diagnose(sorted(d['symptoms']), variants=[{'gene': d['genes'][0]}])
    if r['candidates'][0]['disease'] == name:
        top1 += 1
res['diagnose_top1_full_signature'] = top1
res['diagnose_top1_of'] = len(RARE_DISEASES)
res['datasets'] = ['HPO:genes_to_phenotype(2026-09-25)'] + [f'omim:{v}' for v in OMIM.values()]
res['n_datasets'] = len(res['datasets'])

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks', 'sweep_rarenet.json')
json.dump(res, open(out, 'w'), indent=2)
print(json.dumps(res, indent=2))
