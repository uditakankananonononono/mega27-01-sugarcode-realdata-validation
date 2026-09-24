import json,re,gzip,sys,urllib.request,collections,numpy as np
from scipy.stats import spearmanr
import codonbias.scores as cb
sys.path.insert(0,'/home/sandbox/mega/sc-ai/src')
from sugarcode.bio.codon import cai, STANDARD_CODE
def get(u,t=180): return urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'curl/8.0'}),timeout=t).read()
SP=[('511145','Esch_coli_K_12_MG1655','bacteria',True),('224308','Baci_subt_subtilis_168','bacteria',True),('4932','Scere3','eukaryota',False)]
raw={json.loads(l)['taxid']:json.loads(l) for l in open('multi_cai.jsonl')}
def table(seqs):
    c=collections.Counter()
    for s in seqs:
        for i in range(0,len(s)-3,3): c[s[i:i+3]]+=1
    t={}
    for cod,aa in STANDARD_CODE.items():
        tot=sum(c[x] for x,a in STANDARD_CODE.items() if a==aa); t[cod]=(c[cod]+0.5)/(tot+0.5) if tot else 0
    return t
OUT=[]
for tx,gid,dom,prok in SP:
    r0=raw[tx]; acc=r0['assembly']
    rep=json.loads(get(f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/dataset_report'))
    asm=rep['reports'][0]['assembly_info']['assembly_name'].replace(' ','_'); n=acc.split('_')[1].split('.')[0]
    cds=gzip.decompress(get(f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{acc[:3]}/{n[0:3]}/{n[3:6]}/{n[6:9]}/{acc}_{asm}/{acc}_{asm}_cds_from_genomic.fna.gz")).decode()
    seqs={};names={}
    for blk in cds.split('>')[1:]:
        h,*sq=blk.split('\n'); s=''.join(sq)
        if 'pseudo=true' in h or len(s)%3 or len(s)<300 or not re.fullmatch('[ACGT]+',s): continue
        keys=re.findall(r'\[(?:locus_tag|old_locus_tag|gene)=([^\]]+)\]',h); g=re.search(r'\[gene=([^\]]+)\]',h)
        for k in keys:
            for kk in k.split(','): seqs.setdefault(kk.strip(),s)
        for k in keys: names[k]=g.group(1) if g else ''
    ab={}
    for l in get(f"https://pax-db.org/downloads/latest/datasets/{tx}/{r0['paxdb_file']}").decode().splitlines():
        if l.startswith('#'): continue
        p=l.split('\t')
        if len(p)<3: continue
        sid=p[1].split('.',1)[1] if '.' in p[1] else p[1]
        try: a=float(p[2])
        except: continue
        if a<=0: continue
        for k in (sid,p[0]):
            if k in seqs: ab[k]=a; break
    ribo=set(k for k in seqs if re.match(r'^(rp[lsm][A-Z]|R[Pp][LlSs][0-9])',names.get(k,'') or ''))
    keys=[k for k in ab if k not in ribo]; y=np.log10([ab[k] for k in keys]); S=[seqs[k] for k in keys]
    gt=table(list(dict.fromkeys(seqs.values()))); rt=table([seqs[k] for k in ribo])
    gcn=cb.fetch_GCN_from_GtRNAdb(genome=gid,domain=dom)
    tai=cb.TrnaAdaptationIndex(tGCN=gcn,prokaryote=prok).get_score(S)
    enc=cb.EffectiveNumberOfCodons().get_score(S)
    rec={'taxid':tx,'organism':r0['organism'],'gtrnadb_genome':gid,'n_genes':len(keys),
         'rho_cai_genome':round(float(spearmanr([cai(s,gt) for s in S],y)[0]),4),
         'rho_cai_ribo':round(float(spearmanr([cai(s,rt) for s in S],y)[0]),4),
         'rho_tai':round(float(spearmanr(tai,y,nan_policy='omit')[0]),4),
         'rho_neg_enc':round(float(spearmanr(-np.array(enc),y,nan_policy='omit')[0]),4)}
    OUT.append(rec); print(rec,flush=True)
json.dump({'tools':'codon-bias (TrnaAdaptationIndex dos Reis s-values, EffectiveNumberOfCodons), GtRNAdb tRNA gene copy numbers','results':OUT},open('tai_enc.json','w'),indent=1)
