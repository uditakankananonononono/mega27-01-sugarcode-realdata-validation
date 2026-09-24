"""Validate dna_to_code executable utilities against Ensembl protein + Ensembl VEP.
1) execute_translation(CDS) vs Ensembl REST protein sequence (39/40 mouse transcripts).
2) STANDARD_CODE vs NCBI translation table 1 (Biopython CodonTable).
3) execute_mutation class vs Ensembl VEP consequence for random CDS SNVs.
"""
import json, random, os, sys, time, urllib.request
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from sugarcode.modules.dna_to_code import execute_translation, execute_mutation
from sugarcode.bio.sequence import STANDARD_CODE
from Bio.Data import CodonTable
D = 'data/m6a/ensembl'; OUT = 'data/dna_to_code'
def post(url, body):
    for k in range(5):
        try:
            r = urllib.request.Request(url, data=json.dumps(body).encode(), headers={'Content-Type':'application/json','Accept':'application/json'})
            return json.load(urllib.request.urlopen(r, timeout=120))
        except Exception as e:
            print('retry', k, e); time.sleep(3 + 3*k)
    raise RuntimeError(url)
def fasta(p):
    return ''.join(l.strip() for l in open(p) if not l.startswith('>')).upper()
tx = sorted(f[:-4] for f in os.listdir(D) if f.endswith('.cds'))
cds = {}
for t in tx:
    a, b = map(int, open(f'{D}/{t}.cds').read().split()[:2])
    cds[t] = fasta(f'{D}/{t}.fa')[a:b]  # 0-based half-open (as written by rna_decoder_m6a)
# 1 protein
pf = f'{OUT}/proteins.json'
if not os.path.exists(pf):
    res = post('https://rest.ensembl.org/sequence/id?type=protein', {'ids': tx})
    json.dump({r['query']: r['seq'] for r in res}, open(pf, 'w'))
prot = json.load(open(pf))
tr = {}
for t in tx:
    p = prot.get(t)
    mine = execute_translation(cds[t])
    tr[t] = {'ensembl_len': len(p) if p else None, 'mine_len': len(mine), 'exact': mine == p,
             'cds_len_mod3': len(cds[t]) % 3, 'starts_ATG': cds[t][:3] == 'ATG'}
# 2 code table
ncbi = CodonTable.unambiguous_dna_by_id[1]
table = {c: ncbi.forward_table.get(c, '*') for c in STANDARD_CODE}
code_mism = [c for c in STANDARD_CODE if table[c] != STANDARD_CODE[c]]
# 3 VEP
random.seed(7); vf = f'{OUT}/vep.json'
var = []
for t in [t for t in tx if tr[t]['exact']][:20]:
    s = cds[t]
    pos = random.sample(range(len(s)), 8) + [0, len(s)-3, len(s)-1]
    for p in pos:
        alt = random.choice([b for b in 'ACGT' if b != s[p]])
        var.append((t, p, s[p], alt))
vep = {}
if os.path.exists(vf):
    vep = {r['input']: r for r in json.load(open(vf))}
elif os.environ.get('TRY_VEP'):
    try:
        out = []
        for i in range(0, len(var), 100):
            chunk = [f'{t}:c.{p+1}{r}>{a}' for t, p, r, a in var[i:i+100]]
            out += post('https://rest.ensembl.org/vep/mus_musculus/hgvs', {'hgvs_notations': chunk})
        json.dump(out, open(vf, 'w')); vep = {r['input']: r for r in out}
    except Exception as e:
        print('VEP unavailable', e)
def vep_cls(r, t):
    for tc in r.get('transcript_consequences', []):
        if tc['transcript_id'] == t:
            cs = tc['consequence_terms']
            for k in ['stop_gained', 'stop_lost', 'start_lost', 'stop_retained_variant', 'synonymous_variant', 'missense_variant']:
                if k in cs: return k
            return cs[0]
    return None
from Bio.Seq import Seq
def bio_cls(s, p, a):
    # independent Sequence-Ontology-style classifier on Biopython translation (codon-level, VEP term names)
    ci = p // 3; ref = s[ci*3:ci*3+3]; mut = list(ref); mut[p % 3] = a; mut = ''.join(mut)
    ra, ma = str(Seq(ref).translate()), str(Seq(mut).translate())
    if ci == 0 and ref == 'ATG' and ma != 'M': return 'start_lost'
    if ra == '*' and ma != '*': return 'stop_lost'
    if ra == '*' and ma == '*': return 'stop_retained_variant'
    if ma == '*': return 'stop_gained'
    return 'synonymous_variant' if ra == ma else 'missense_variant'
rows = []; conf = {}
for t, p, r, a in var:
    key = f'{t}:c.{p+1}{r}>{a}'
    v = vep.get(key); vc = vep_cls(v, t) if v else None
    ref_c = vc or bio_cls(cds[t], p, a)
    m = execute_mutation(cds[t], p, a)
    mc_old = 'synonymous_variant' if m['silent'] else 'stop_gained' if m['nonsense'] else 'missense_variant'
    mc = m.get('consequence', mc_old)
    rows.append({'module_legacy_flags_class': mc_old, 'hgvs': key, 'vep': vc, 'reference': ref_c, 'reference_source': 'VEP' if vc else 'biopython_SO_rules', 'module': mc,
                 'module_raw': {k: m[k] for k in ('silent', 'nonsense')}})
    conf.setdefault(ref_c, {}).setdefault(mc, 0); conf[ref_c][mc] += 1
scored = rows
agree = sum(x['reference'] == x['module'] for x in scored)
legacy_agree = sum(x['reference'] == x['module_legacy_flags_class'] for x in scored)
legacy_agree_lenient = sum(x['reference'] == x['module_legacy_flags_class'] or (x['reference'] == 'stop_retained_variant' and x['module_legacy_flags_class'] == 'synonymous_variant') for x in scored)
# template-strand check on real cDNA: mRNA == complement of reversed template
from sugarcode.modules.dna_to_code import execute_transcription
comp = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
tmpl_ok = 0; naive_ok = 0
for t in tx:
    c = fasta(f'{D}/{t}.fa'); template = c[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    m = execute_transcription(c)
    tmpl_ok += ''.join(comp[b] for b in reversed(template)) == m
    naive_ok += template.replace('T', 'U') == m
res = {'module': 'dna_to_code', 'translation': {'n': len(tx), 'exact': sum(v['exact'] for v in tr.values()), 'per_tx': tr},
       'code_table_vs_ncbi1': {'n_codons': len(STANDARD_CODE), 'mismatches': code_mism},
       'template_strand': {'n': len(tx), 'revcomp_snippet_matches_mRNA': tmpl_ok, 'old_TtoU_on_template_matches_mRNA': naive_ok},
       'vep': {'n_submitted': len(var), 'n_scored': len(scored), 'agree': agree, 'confusion_vep_to_module': conf,
               'legacy_flags_agree': legacy_agree, 'legacy_flags_agree_lenient_stop_retained': legacy_agree_lenient, 'n_vep_scored': sum(1 for x in rows if x['vep']), 'disagreements': [x for x in scored if x['reference'] != x['module']]}}
json.dump(res, open('benchmarks/sweep_dna_to_code.json', 'w'), indent=1)
print(json.dumps({k: (v if k != 'translation' else {'n': v['n'], 'exact': v['exact']}) for k, v in res.items() if k != 'vep'}, indent=0))
print(res['template_strand']); print(legacy_agree, legacy_agree_lenient); print(res['vep']['n_submitted'], res['vep']['n_scored'], agree, json.dumps(conf))
