import os
import json, time, urllib.request, urllib.parse
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
def get(u):
    for i in range(4):
        try: return json.load(urllib.request.urlopen(u, timeout=30))
        except Exception as e: time.sleep(1.5*(i+1))
    raise SystemExit('fail '+u)
genes=['BRCA1','BRCA2','TP53','MLH1','MSH2','CFTR','APC','PTEN','LDLR','MYH7','SCN5A','RYR1']
allrec={}
for g in genes:
    ids=get(E+'esearch.fcgi?'+urllib.parse.urlencode({'db':'clinvar','term':f'{g}[gene] AND single_gene[prop]','retmax':60,'retmode':'json','sort':'relevance'}))['esearchresult']['idlist']
    time.sleep(0.4)
    s=get(E+'esummary.fcgi?'+urllib.parse.urlencode({'db':'clinvar','id':','.join(ids),'retmode':'json'}))['result']
    time.sleep(0.4)
    for i in s.get('uids',[]):
        r=s[i]; gc=r.get('germline_classification',{})
        allrec[i]={'gene':g,'title':r.get('title'),'mc':r.get('molecular_consequence_list'),'review':gc.get('review_status'),'sig':gc.get('description'),'protein_change':r.get('protein_change')}
json.dump(allrec,open(os.path.expanduser('~/mega/m01/data/clinvar_esummary/clinvar_esummary.json'),'w'),indent=0)
print(len(allrec))
