import json,os,re,gzip,sys,urllib.request,collections,numpy as np
from scipy.stats import spearmanr
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.codon import cai, STANDARD_CODE
def get(u,t=60): return urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'curl/8.0'}),timeout=t).read()
taxa=[x for x in re.findall(r'href="([0-9]+)/"',get('https://pax-db.org/downloads/latest/datasets/').decode())]
out='multi_cai.jsonl'; done={json.loads(l)['taxid'] for l in open(out)} if os.path.exists(out) else set()
def table(seqs):
    c=collections.Counter()
    for s in seqs:
        for i in range(0,len(s)-3,3): c[s[i:i+3]]+=1
    t={}
    for cod,aa in STANDARD_CODE.items():
        tot=sum(c[x] for x,a in STANDARD_CODE.items() if a==aa); t[cod]=(c[cod]+0.5)/(tot+0.5) if tot else 0
    return t
for tx in taxa:
    if tx in done: continue
    rec={'taxid':tx}
    try:
        idx=get(f'https://pax-db.org/downloads/latest/datasets/{tx}/').decode()
        f=re.search(r'href="([^"]*integrated[^"]*)"',idx)
        if not f: raise Exception('no integrated')
        rep=json.loads(get(f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{tx}/dataset_report?filters.reference_only=true&page_size=1'))
        if not rep.get('reports'): raise Exception('no reference assembly')
        r0=rep['reports'][0]; acc=r0['accession']; asm=r0['assembly_info']['assembly_name'].replace(' ','_')
        rec.update(organism=r0['organism']['organism_name'],assembly=acc,paxdb_file=f.group(1))
        n=acc.split('_')[1].split('.')[0]; base=f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{n[0:3]}/{n[3:6]}/{n[6:9]}/{acc}_{asm}/{acc}_{asm}_cds_from_genomic.fna.gz"
        cds=gzip.decompress(get(base,180)).decode()
        seqs={};names={}
        for blk in cds.split('>')[1:]:
            h,*sq=blk.split('\n'); s=''.join(sq)
            if 'pseudo=true' in h or len(s)%3 or len(s)<300 or not re.fullmatch('[ACGT]+',s): continue
            keys=re.findall(r'\[(?:locus_tag|old_locus_tag|gene)=([^\]]+)\]',h)
            g=re.search(r'\[gene=([^\]]+)\]',h)
            for k in keys:
                for kk in k.split(','): seqs.setdefault(kk.strip(),s)
            for k in keys: names[k]=g.group(1) if g else ''
        ab={}
        for l in get(f'https://pax-db.org/downloads/latest/datasets/{tx}/{f.group(1)}').decode().splitlines():
            if l.startswith('#'): continue
            p=l.split('\t')
            if len(p)<3: continue
            sid=p[1].split('.',1)[1] if '.' in p[1] else p[1]
            try: a=float(p[2])
            except: continue
            if a<=0: continue
            for k in (sid,p[0]):
                if k in seqs: ab[k]=a; break
        keys=list(ab); rec['n_matched']=len(keys)
        if len(keys)<300: raise Exception(f'only {len(keys)} matched')
        y=np.log10([ab[k] for k in keys]); S=[seqs[k] for k in keys]
        ribo=[k for k in set(seqs) if re.match(r'^(rp[lsm][A-Z]|R[Pp][LlSs][0-9])',names.get(k,'') or '')]
        rec['n_ribo']=len(ribo)
        gt=table(list(dict.fromkeys(seqs.values())))
        xg=[cai(s,gt) for s in S]
        ev=[i for i,k in enumerate(keys) if k not in set(ribo)]
        rec['rho_genome_table']=round(float(spearmanr([xg[i] for i in ev],y[ev])[0]),4)
        if len(ribo)>=15:
            rt=table([seqs[k] for k in ribo]); xr=[cai(S[i],rt) for i in ev]
            rec['rho_ribo_ref']=round(float(spearmanr(xr,y[ev])[0]),4)
        rng=np.random.default_rng(7); fo=rng.integers(0,5,len(keys)); pr=np.zeros(len(keys))
        for k5 in range(5):
            tr=np.where(fo!=k5)[0]; top=tr[np.argsort(-y[tr])[:max(20,int(.05*len(tr)))]]
            t=table([S[i] for i in top])
            for i in np.where(fo==k5)[0]: pr[i]=cai(S[i],t)
        rec['rho_top5_cv']=round(float(spearmanr(pr,y)[0]),4)
    except Exception as e: rec['error']=str(e)[:150]
    open(out,'a').write(json.dumps(rec)+'\n'); print(tx,rec.get('organism'),rec.get('rho_genome_table'),rec.get('rho_ribo_ref'),rec.get('rho_top5_cv'),rec.get('error',''),flush=True)
