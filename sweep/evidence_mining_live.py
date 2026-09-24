"""evidence_mining.EvidenceMiner vs independent live calls to PubMed E-utilities and
ClinicalTrials.gov API v2 (both free, no key): ID lists, record fields, cache/offline semantics."""
import sys, os, json, time, tempfile
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
import urllib.request, urllib.parse
from sugarcode.modules.evidence_mining.core import EvidenceMiner
from sugarcode.modules.evidence_mining.client import OfflineCacheMiss

QUERY = 'ivacaftor cystic fibrosis'
LIM = 10
res = {'query': QUERY, 'limit': LIM}

def jget(url, params):
    full = url + '?' + urllib.parse.urlencode(params)
    req = urllib.request.Request(full, headers={'User-Agent': 'mega27-sweep/1.0'})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())

# independent esearch
es = jget('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',
          {'db': 'pubmed', 'term': QUERY, 'retmode': 'json', 'retmax': LIM, 'sort': 'relevance'})
ref_pmids = [str(x) for x in es['esearchresult']['idlist']]
time.sleep(0.4)
# independent CTG v2
ct = jget('https://clinicaltrials.gov/api/v2/studies',
          {'query.term': QUERY, 'pageSize': LIM, 'format': 'json'})
ref_ncts = [s['protocolSection']['identificationModule']['nctId'] for s in ct['studies']]

with tempfile.TemporaryDirectory() as cache:
    miner = EvidenceMiner(cache)
    out = miner.mine(QUERY, pubmed_limit=LIM, trial_limit=LIM)
    pm = [r for r in out['records'] if r['source'] == 'PubMed']
    tr = [r for r in out['records'] if r['source'] == 'ClinicalTrials.gov']
    got_pmids = [r['source_id'] for r in pm]
    got_ncts = [r['source_id'] for r in tr]
    res['pmid_exact_order_match'] = got_pmids == ref_pmids
    res['pmid_set_overlap'] = len(set(got_pmids) & set(ref_pmids))
    res['nct_exact_order_match'] = got_ncts == ref_ncts
    res['nct_set_overlap'] = len(set(got_ncts) & set(ref_ncts))
    res['n_pubmed'] = len(pm); res['n_trials'] = len(tr)
    res['record_count_field'] = out['record_count']
    res['record_count_matches'] = out['record_count'] == len(out['records'])
    # field integrity: URLs, IDs, no empty titles
    res['pubmed_url_format_ok'] = all(r['url'] == f"https://pubmed.ncbi.nlm.nih.gov/{r['source_id']}/" for r in pm)
    res['pmids_numeric'] = all(r['source_id'].isdigit() for r in pm)
    res['nct_format_ok'] = all(r['source_id'].startswith('NCT') for r in tr)
    res['titles_nonempty'] = sum(1 for r in out['records'] if r['title'] and r['title'] != 'N/A')
    # strict vs non-strict offline cache-miss behavior
    try:
        miner.mine(QUERY, offline=True, strict=True) if False else None
    except Exception:
        pass
    with tempfile.TemporaryDirectory() as empty:
        m2 = EvidenceMiner(empty)
        try:
            m2.mine(QUERY, offline=True, strict=True)
            res['strict_offline_raises'] = False
        except OfflineCacheMiss:
            res['strict_offline_raises'] = True
        except Exception as e:
            res['strict_offline_raises'] = f'other: {type(e).__name__}'
        m3 = EvidenceMiner(empty)
        o3 = m3.mine(QUERY, offline=True, strict=False)
        res['nonstrict_offline_errors_inband'] = sorted(o3['source_errors'].keys())
        res['nonstrict_offline_records'] = o3['record_count']
    # offline-from-cache replay equals live result
    out2 = miner.mine(QUERY, pubmed_limit=LIM, trial_limit=LIM, offline=True)
    res['offline_cache_replay_identical'] = [r['source_id'] for r in out2['records']] == [r['source_id'] for r in out['records']]

res['ref_pmids'] = ref_pmids
res['ref_ncts'] = ref_ncts
res['datasets'] = [f'pmid:{p}' for p in ref_pmids] + [f'nct:{n}' for n in ref_ncts]
res['n_datasets'] = len(res['datasets'])
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'benchmarks', 'sweep_evidence_mining.json')
json.dump(res, open(out, 'w'), indent=2)
print(json.dumps({k: v for k, v in res.items() if k not in ('ref_pmids', 'ref_ncts', 'datasets')}, indent=2))
