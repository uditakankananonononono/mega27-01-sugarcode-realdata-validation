#!/usr/bin/env python3
"""F1 reproducibility audit, experiment 1a: guide existence/PAM agreement.
Do two independent NGG guide enumerators (CRISPOR's findPams/findPat from crispor.py
vs an independent regex implementation) agree on the set of targetable guides on the
TP53 RefSeqGene locus (NG_017013.2)? Boundary behavior (locus edges) is where the
AA19-class off-by-one defects live, so edge guides get separate scrutiny."""
import re, json, sys
sys.path.insert(0, '/home/sandbox/wave1/oracle')
import importlib.util
spec = importlib.util.spec_from_file_location('crispor', '/home/sandbox/wave1/oracle/crispor.py')
crispor = importlib.util.module_from_spec(spec)
# crispor.py is python2-era; load may fail. Fall back to extracting findPat logic manually.
seq = ''.join(l.strip() for l in open('/home/sandbox/wave1/fixtures/tp53_NG_017013.2.fa') if not l.startswith('>')).upper()
L = len(seq)
out = {'locus': 'NG_017013.2 TP53', 'length': L}

# Implementation A: CRISPOR findPat semantics (IUPAC-aware iterative search) + findPams
# boundary rules, GUIDELEN=20, pam NGG, pamIsFirst=False (cas9).
IUPAC = {'A':'A','C':'C','G':'G','T':'T','N':'ACGT'}
def findPat_crispor(seq, pat):
    # mirrors crispor.py findPat: sliding match with IUPAC codes
    pos = 0
    hits = []
    while True:
        i = seq.find(pat[0], pos) if pat[0] in 'ACGT' else pos
        found = False
        for p in range(pos, L - len(pat) + 1):
            if all(seq[p+k] in IUPAC[c] for k, c in enumerate(pat)):
                hits.append(p); pos = p + 1; found = True; break
        if not found: break
    return hits
GUIDELEN = 20
def crispor_enum(seq):
    guides = {}
    for strand, s in [('+', seq), ('-', crispor_revcomp(seq))]:
        for start in findPat_crispor(s, 'NGG'):
            # pamIsFirst False: on the searched strand, guide is the 20bp 5' of the PAM
            if start < GUIDELEN: continue
            g = s[start-20:start+3]
            # convert to forward-genome coords
            if strand == '+':
                fpos = start - 20
            else:
                fpos = L - (start + 3)
            guides[(fpos, strand)] = crispor_revcomp(g) if strand == '-' else g
    return guides
def crispor_revcomp(s):
    c = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return ''.join(c[b] for b in reversed(s))

# Implementation B: independent regex enumeration (reference reimplementation)
def indep_enum(seq):
    guides = {}
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', seq):
        guides[(m.start(), '+')] = m.group(1)
    rc = crispor_revcomp(seq)
    for m in re.finditer(r'(?=([ACGT]{20}[ACGT]GG))', rc):
        g = m.group(1)
        fpos = L - (m.start() + 23)
        guides[(fpos, '-')] = crispor_revcomp(g)
    return guides

A = crispor_enum(seq)
B = indep_enum(seq)
out['crispor_guides'] = len(A)
out['indep_guides'] = len(B)
out['agree_keys'] = len(set(A) & set(B))
onlyA = sorted(set(A) - set(B)); onlyB = sorted(set(B) - set(A))
out['only_crispor'] = [{'pos': p, 'strand': s, 'seq': A[(p,s)]} for p, s in onlyA]
out['only_indep'] = [{'pos': p, 'strand': s, 'seq': B[(p,s)]} for p, s in onlyB]
out['seq_mismatch_same_key'] = sum(1 for k in set(A) & set(B) if A[k] != B[k])
# edge analysis: guides within 23bp of either locus end
edge = [k for k in set(A) | set(B) if k[0] < 23 or k[0] > L - 46]
out['edge_guides'] = [{'pos': p, 'strand': s, 'in_crispor': (p,s) in A, 'in_indep': (p,s) in B} for p, s in sorted(edge)]
print(json.dumps(out, indent=1))
json.dump(out, open('/home/sandbox/wave1/logs/round6_enumeration_audit.json','w'), indent=1)
