"""Sweep: biogpt_lit vs live PubMed/iCite + statsmodels meta-analysis + recomputation."""
import json, math, os, sys, time, urllib.request
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
import numpy as np
from sugarcode.modules.biogpt_lit import core as L

CD = os.path.join(os.path.dirname(__file__), '..', 'data', 'icite')
os.makedirs(CD, exist_ok=True)
res = {'module': 'biogpt_lit', 'sources': ['https://eutils.ncbi.nlm.nih.gov (via bio.entrez)', 'https://icite.od.nih.gov/api/pubs?pmids=...', 'statsmodels 0.15.0']}

# 1. live ingest: PubMed abstracts + real citation counts from iCite
kg = L.KnowledgeGraph()
ing = L.ingest_pubmed(kg, 'SIRT1 resveratrol activation', retmax=8)
pmids = ing['pmids']
fn = os.path.join(CD, 'icite_' + '_'.join(pmids) + '.json')
if not os.path.exists(fn):
    urllib.request.urlretrieve('https://icite.od.nih.gov/api/pubs?pmids=' + ','.join(pmids), fn); time.sleep(0.4)
cit = {str(p['pmid']): p.get('citation_count', 0) for p in json.load(open(fn))['data']}
kg2 = L.KnowledgeGraph()
from sugarcode.bio import entrez
abs_ = entrez.pubmed_abstracts(pmids)
papers = []
for a in abs_:
    claims = L.extract_claims(a['title'] + '. ' + a['abstract'])
    papers.append({'pmid': a['pmid'], 'year': a['year'], 'citations': cit.get(a['pmid'], 0), 'n_claims': len(claims)})
    kg2.ingest({'id': f"PMID:{a['pmid']}", 'title': a['title'], 'year': int(a['year']) if a['year'].isdigit() else 2000,
                'citations': cit.get(a['pmid'], 0), 'claims': claims})
res['live_ingest'] = {'papers': papers, 'claims_no_citations': ing['claims_ingested'],
                      'claims_with_icite': len(kg2.claims),
                      'icite_citation_counts': cit}
res['extracted_claims_audit'] = [{'paper': c['paper'], 'triple': [c['subject'], c['relation'], c['object']]} for c in kg2.claims]

# 2. extractor unit checks (negation, loss-of, normalisation)
sent = {
 'SIRT1 activates PGC1A.': 1,
 'Neither drug inhibits mTOR.': 0,
 'Loss of PTEN activates AKT.': 0,
 'PTEN knockdown activates AKT.': 0,  # 'knockdown' subject is dropped
 'Drug X does not inhibit SIRT1.': 0,
 'VEGF stimulates angiogenesis.': 1,
 'The kinase phosphorylates BAD.': 0,  # 'phosphorylates' not in relation map
}
ext = {s: len(L.extract_claims(s)) for s in sent}
res['extractor_unit'] = {'expected': sent, 'got': ext, 'all_match': all(ext[s] == n for s, n in sent.items())}

# 3. evidence weight + query + consensus recomputation
kg3 = L.KnowledgeGraph()
kg3.ingest({'id': 'p1', 'year': 2003, 'citations': 1000, 'replicated': True,
            'claims': [{'subject': 'resveratrol', 'relation': 'activates', 'object': 'SIRT1', 'n': 50}]})
kg3.ingest({'id': 'p2', 'year': 2005, 'citations': 10, 'replicated': False,
            'claims': [{'subject': 'resveratrol', 'relation': 'inhibits', 'object': 'SIRT1', 'n': 5}]})
w1 = 0.3 + 0.3 * min(math.log10(1 + 1000) / 3, 1.0) + 0.25 + 0.15 * min(50 / 100, 1.0)
w2 = 0.3 + 0.3 * min(math.log10(1 + 10) / 3, 1.0) + 0.15 * min(5 / 100, 1.0)
q = kg3.query('resveratrol', 'SIRT1')
res['query_check'] = {'w1_module': kg3._weight(kg3.claims[0]), 'w1_recompute': round(w1, 3),
                      'w2_module': kg3._weight(kg3.claims[1]), 'w2_recompute': round(w2, 3),
                      'consensus': q['consensus'], 'expected_consensus': 'activates',
                      'confidence': q['confidence'], 'expected_confidence': round(w1, 3)}
contra = kg3.contradictions()
res['contradictions_check'] = {'found': len(contra), 'pair': [contra[0]['claim_a']['relation'], contra[0]['claim_b']['relation']] if contra else None}
tc = L.temporal_consensus(kg3, 'resveratrol', 'SIRT1')
res['temporal_consensus'] = {'years': [t['year'] for t in tc['trajectory']], 'status': tc['status'],
                             'agreement_last': tc['trajectory'][-1]['agreement'],
                             'expected_agreement_last': round(w1 / (w1 + w2), 3)}
cc = L.contradiction_context(kg3, 'resveratrol', 'SIRT1')
res['contradiction_context_count'] = cc['count']

# 4. meta-analysis vs statsmodels
from statsmodels.stats.meta_analysis import combine_effects
effs = [0.5, 0.2, 0.7, -0.1, 0.4]; ses = [0.1, 0.15, 0.2, 0.12, 0.25]
ma = L.meta_analysis(effs, ses)
sm = combine_effects(np.array(effs), np.array(ses) ** 2, method_re='DL')
res['meta_analysis'] = {'module_fixed': ma['fixed_effect'], 'sm_fixed': float(sm.mean_effect_fe),
                        'module_random': ma['random_effect'], 'sm_random': float(sm.mean_effect_re),
                        'module_tau2': ma['tau_squared'], 'sm_tau2': float(sm.tau2),
                        'module_ci': ma['ci95'],
                        'fixed_match': abs(ma['fixed_effect'] - sm.mean_effect_fe) < 1e-9,
                        'random_match': abs(ma['random_effect'] - sm.mean_effect_re) < 1e-9,
                        'tau2_match': abs(ma['tau_squared'] - sm.tau2) < 1e-9}

# 5. graph algorithms: multi_hop BFS vs independent enumeration, path strength, gaps, feedback
kg4 = L.KnowledgeGraph()
for pid, s, r, o, cit_n in [('a', 'ASCR', 'activates', 'BDNF', 500), ('b', 'BDNF', 'increases', 'LTP', 100),
                            ('c', 'ASCR', 'binds', 'TRKB', 50), ('d', 'TRKB', 'activates', 'LTP', 5)]:
    kg4.ingest({'id': pid, 'year': 2010, 'citations': cit_n, 'claims': [{'subject': s, 'relation': r, 'object': o}]})
mh = kg4.multi_hop('ASCR', 'LTP', max_hops=3)
paths = {tuple(p['path']) for p in mh['paths']}
res['multi_hop'] = {'paths': sorted(paths), 'expected': [('ascr', 'bdnf', 'ltp'), ('ascr', 'trkb', 'ltp')],
                    'match': paths == {('ascr', 'bdnf', 'ltp'), ('ascr', 'trkb', 'ltp')}}
ep = L.evidence_paths(kg4, 'ASCR', 'LTP')
wa = 0.3 + 0.3 * min(math.log10(501) / 3, 1) + 0.15 * 0.1
wb = 0.3 + 0.3 * min(math.log10(101) / 3, 1) + 0.15 * 0.1
res['evidence_paths_strength'] = {'top_strength': ep['paths'][0]['path_strength'],
                                  'expected_min_edge': round(min(wa, wb), 3)}
gaps = L.prioritize_gaps(kg4)
res['prioritize_gaps'] = {'top_edge': list(gaps[0]['edge']), 'leverage': gaps[0]['information_leverage'],
                          'expected_top': ['trkb', 'activates', 'ltp']}
uv = L.update_validation(kg4, 'BDNF', 'increases', 'LTP', confirmed=False, n=200)
res['update_validation'] = {'before_consensus': uv['before']['consensus'], 'after_relations': list(uv['after']['relations']),
                            'flips_to_decreases': 'decreases' in uv['after']['relations']}
out = sys.argv[1] if len(sys.argv) > 1 else 'benchmarks/sweep_biogpt_lit.json'
json.dump(res, open(out, 'w'), indent=1, default=str)
print(json.dumps({'papers': len(papers), 'claims': len(kg2.claims), 'extractor_ok': res['extractor_unit']['all_match'],
                  'meta_matches': [bool(res['meta_analysis']['fixed_match']), bool(res['meta_analysis']['random_match']), bool(res['meta_analysis']['tau2_match'])],
                  'multi_hop': res['multi_hop']['match'], 'consensus': q['consensus'],
                  'temporal': res['temporal_consensus']['status'],
                  'gaps_top': res['prioritize_gaps']['top_edge']}, indent=1))
