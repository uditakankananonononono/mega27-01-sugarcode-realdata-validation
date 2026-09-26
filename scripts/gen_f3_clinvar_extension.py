#!/usr/bin/env python3
"""F3 extension: reservoir-sample 5000 ClinVar variants (2500 path-class + 2500 benign-class).
Protocol matches data/README.md of mega27-01 (GRCh38 SNVs, criteria-provided review status),
different seed; de-duplicated against the frozen 500 after repo clone.
Source: ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz, retrieved 2026-09-26.
"""
import gzip, random, csv, sys

SEED = 20260926
random.seed(SEED)
SRC = '/home/sandbox/wave1/fixtures/clinvar/variant_summary.txt.gz'
PATH_CLASS = {'Pathogenic','Likely pathogenic','Pathogenic/Likely pathogenic'}
BENIGN_CLASS = {'Benign','Likely benign','Benign/Likely benign'}
N_PER = 2500
res = {'path': [], 'ben': []}
counts = {'path': 0, 'ben': 0}
keep = lambda sig, cls: sig in cls

with gzip.open(SRC, 'rt') as f:
    rd = csv.DictReader(f, delimiter='\t')
    cols = rd.fieldnames
    for row in rd:
        if row.get('Type') != 'single nucleotide variant': continue
        if row.get('Assembly') != 'GRCh38': continue
        if 'criteria provided' not in row.get('ReviewStatus',''): continue
        sig = row.get('ClinicalSignificance','')
        if sig in PATH_CLASS: k = 'path'
        elif sig in BENIGN_CLASS: k = 'ben'
        else: continue
        counts[k] += 1
        r = res[k]
        if len(r) < N_PER: r.append(row)
        else:
            j = random.randrange(counts[k])
            if j < N_PER: r[j] = row

sel = ['VariationID','Name','GeneSymbol','Chromosome','Start','ReferenceAllele','AlternateAllele',
       'ClinicalSignificance','ReviewStatus','Assembly','Type','RS (dbSNP)'] if 'RS (dbSNP)' in cols else list(cols)[:12]
out = '/home/sandbox/wave1/fixtures/clinvar/clinvar_extension_5000.tsv'
with open(out,'w',newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(sel)
    for k in ('path','ben'):
        for row in res[k]:
            w.writerow([row.get(c,'') for c in sel])
print('eligible:', counts, '-> wrote', sum(len(v) for v in res.values()), 'to', out, file=sys.stderr)
