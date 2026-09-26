#!/usr/bin/env python3
"""F1 fixture extension: TP53 (NG_017013.2) guide/off-target pairs with CRISPOR-oracle
CFD (Doench 2016) and MIT (Hsu 2013) off-target scores.
Disjoint from existing BRCA1 NG_005905.2 fixtures and deepsplice training accessions.
Oracle: github.com/maximilianh/crisporWebsite CFD_Scoring pkls + calcCfdScore doctests (6/6 exact).
"""
import pickle, random, sys

SEED = 20260926
random.seed(SEED)
ORACLE = '/home/sandbox/wave1/oracle'
mm_scores = pickle.load(open(f'{ORACLE}/mismatch_score.pkl','rb'))
pam_scores = pickle.load(open(f'{ORACLE}/pam_scores.pkl','rb'))
_comp = {'A':'T','C':'G','G':'C','T':'A','U':'A'}

def calc_cfd(wt, sg, pam):
    score = 1.0
    sg = sg.upper().replace('T','U'); wt = wt.upper().replace('T','U')
    for i in range(len(sg)):
        if wt[i] != sg[i]:
            score *= mm_scores['r'+wt[i]+':d'+_comp[sg[i]]+','+str(i+1)]
    return score * pam_scores[pam.upper()]

hitScoreM = [0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]
def calc_mit(s1, s2):
    s1, s2 = s1.upper(), s2.upper()
    dists, mmCount, lastMmPos, score1 = [], 0, None, 1.0
    for pos in range(20):
        if s1[pos] != s2[pos]:
            mmCount += 1
            if lastMmPos is not None: dists.append(pos-lastMmPos)
            score1 *= 1-hitScoreM[pos]
            lastMmPos = pos
    score2 = 1.0 if mmCount < 2 else 1.0/(((19-sum(dists)/len(dists))/19.0)*4+1)
    score3 = 1.0 if mmCount == 0 else 1.0/(mmCount**2)
    return score1*score2*score3*100

def revcomp(s):
    return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

seq = ''.join(l.strip() for l in open('/home/sandbox/wave1/fixtures/tp53_NG_017013.2.fa') if not l.startswith('>')).upper()
print(f'locus length {len(seq)}', file=sys.stderr)

# mine all 20nt+NGG guides, both strands
guides = []
for strand, s in (('+', seq), ('-', revcomp(seq))):
    for i in range(len(s)-23+1):
        g, pam = s[i:i+20], s[i+20:i+23]
        if pam[1:] == 'GG' and 'N' not in g:
            guides.append((g, pam, strand, i))
print(f'{len(guides)} NGG guides mined', file=sys.stderr)
sample = random.sample(guides, 350)

BASES = 'ACGT'
def mutate(g, n):
    positions = sorted(random.sample(range(20), n))
    ot = list(g)
    for p in positions:
        ot[p] = random.choice([b for b in BASES if b != g[p]])
    return ''.join(ot), positions

rows = []
for gi, (g, pam, strand, pos) in enumerate(sample):
    for n_mm, rep in ((1,0),(1,1),(2,0),(4,0)):
        ot, mmpos = mutate(g, n_mm)
        rows.append(dict(pair_id=f'TP53_G{gi:03d}_mm{n_mm}r{rep}', guide_id=f'TP53_G{gi:03d}',
            guide_seq=g, pam=pam, strand=strand, locus_pos=pos, offtarget_seq=ot,
            n_mismatches=n_mm, mm_positions=';'.join(map(str,mmpos)),
            oracle_cfd=calc_cfd(g, ot, pam[1:]), oracle_mit=calc_mit(g, ot)))

out = '/home/sandbox/wave1/fixtures/crispr_tp53_pairs.tsv'
with open(out,'w') as f:
    hdr = ['pair_id','guide_id','guide_seq','pam','strand','locus_pos','offtarget_seq',
           'n_mismatches','mm_positions','oracle_cfd','oracle_mit']
    f.write('\t'.join(hdr)+'\n')
    for r in rows:
        f.write('\t'.join(str(r[h]) for h in hdr)+'\n')
print(f'wrote {len(rows)} pairs -> {out}', file=sys.stderr)
