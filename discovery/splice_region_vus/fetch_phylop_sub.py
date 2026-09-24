"""phyloP100way at ~1,500 modelled variants (stratified by class x site, seed 11) via the UCSC REST API, 1 req/s, resumable."""
import csv, json, time, urllib.request, os, random
rows=[r for r in csv.DictReader(open('pwm_scores.tsv'),delimiter='\t') if r['train_gene']=='0' and r['cls'] in ('P','B')]
rng=random.Random(7); data=[]
for s in ('donor','acceptor'):
    P=[r for r in rows if r['site']==s and r['cls']=='P']; B=[r for r in rows if r['site']==s and r['cls']=='B']
    data+=P+rng.sample(B,min(len(B),4*len(P)))
r2=random.Random(11); sub=[]
for s in ('donor','acceptor'):
    for c in ('P','B'):
        g=[r for r in data if r['site']==s and r['cls']==c]; sub+=r2.sample(g,round(len(g)*1500/len(data)))
open('phylop_sub_ids.txt','w').write('\n'.join(r['vid'] for r in sub)+'\n')
snv={r[0]:r for r in csv.reader(open('splice_region_snvs.tsv'),delimiter='\t')}
done=set()
if os.path.exists('phylop_sub.tsv'): done={l.split('\t')[0] for l in open('phylop_sub.tsv')}
out=open('phylop_sub.tsv','a')
for r in sub:
    v=r['vid']
    if v in done: continue
    s=snv[v]; c='chr'+s[5]; p=int(s[6])
    u=f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=phyloP100way;chrom={c};start={p-1};end={p}"
    val='ERR'
    for a in range(3):
        try:
            j=json.load(urllib.request.urlopen(u,timeout=30)); vals=j.get('phyloP100way',[]); val=vals[0]['value'] if vals else 'NA'; break
        except Exception: time.sleep(3*(a+1))
    out.write(f"{v}\t{c}\t{p}\t{val}\n"); out.flush(); time.sleep(1.0)
print('done')
