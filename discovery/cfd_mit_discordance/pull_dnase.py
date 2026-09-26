#!/usr/bin/env python3
"""Pull ENCODE DNase-seq signal (UCSC API, hg19) at each mapped off-target/guide locus
in the cell type matched to its study. Incremental: skips seqs already in the JSON."""
import csv, json, time, urllib.request, os
hits = json.load(open('logs/seq_coords_origdata.json'))
rows = list(csv.DictReader(open('/home/sandbox/mega27-01-sugarcode-realdata-validation/data/crispor_offtargets_7studies.tsv'), delimiter='\t'))
seq2study = {}
for r in rows:
    seq2study.setdefault(r['seq'].upper(), r['name'].split('_')[0])  # 'Kim/K562' keeps cell line, others plain study
# study -> (track, confidence note)
TRACK = {
  'Cho':    ('wgEncodeOpenChromDnaseHelas3Sig',   'Cho2014 HeLa; medium confidence (no direct source pulled)'),
  'Hsu':    ('wgEncodeOpenChromDnaseHek293tSig',  'Hsu2013 HEK293T; well documented'),
  'Kim/K562':('wgEncodeOpenChromDnaseK562SigV2',  'Kim2015 K562; VERIFIED crisporPaper README'),
  'Kim/Hap1':(None, 'Hap1: no hg19 DNase track found on UCSC - excluded'),
  'Wang':   ('wgEncodeOpenChromDnaseK562SigV2',   'Wang2015 K562; medium confidence'),
  'Ran':    ('wgEncodeOpenChromDnaseHek293tSig',  'Ran2015 HEK293; 293T track as proxy'),
  'Tsai':   ('wgEncodeOpenChromDnaseHek293tSig',  'Tsai2015 U2OS+HEK293; U2OS has no track, HEK293T proxy'),
  'Frock':  (None, 'Frock2015 U2OS; no hg19 track found - excluded'),
  'Kim16':  (None, 'Kim2016 Digenome-seq cell-free assay - no chromatin by design, excluded'),
}
OUT = 'logs/dnase_signals.json'
sig = json.load(open(OUT)) if os.path.exists(OUT) else {}
done = set(sig)
todo = []
for s, hs in hits.items():
    if s in done: continue
    study = seq2study.get(s)
    track = TRACK.get(study, (None, 'unknown study'))[0]
    if not track: continue
    chrom, st, en = hs[0][0], hs[0][1], hs[0][2]
    todo.append((s, study, track, chrom, st, en))
print('todo:', len(todo), 'already:', len(done))
t0 = time.time()
n = 0
for s, study, track, chrom, st, en in todo:
    url = f'https://api.genome.ucsc.edu/getData/track?genome=hg19;track={track};chrom={chrom};start={st};end={en}'
    try:
        req = urllib.request.Request(url, headers={'User-Agent':'research/1.0'})
        d = json.load(urllib.request.urlopen(req, timeout=20))
        vals = [v['value'] for v in d.get(track, [])]
        sig[s] = {'study': study, 'track': track, 'chrom': chrom, 'start': st, 'end': en,
                  'mean': sum(vals)/len(vals) if vals else 0.0, 'max': max(vals) if vals else 0.0, 'n_bases': len(vals)}
    except Exception as e:
        sig[s] = {'study': study, 'track': track, 'error': str(e)}
    n += 1
    if n % 25 == 0:
        json.dump(sig, open(OUT,'w'))
        print(f'{n}/{len(todo)} elapsed={time.time()-t0:.0f}s', flush=True)
    if time.time() - t0 > 100:
        break
json.dump(sig, open(OUT,'w'))
print('saved', len(sig), 'total; this run', n)
