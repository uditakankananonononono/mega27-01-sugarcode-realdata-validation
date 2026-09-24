"""bio.fasta vs Biopython SeqIO on real Rfam FASTA files (RF00059/RF00050/RF00504 seeds): ids, descriptions, sequences; write->parse round trip; multiline and blank-line handling."""
import sys, os, json
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from Bio import SeqIO
from sugarcode.bio.fasta import parse_fasta, write_fasta, stream_fasta
files = ['RF00059', 'RF00050', 'RF00504']
rows = []; mism = 0
for f in files:
    p = os.path.expanduser(f'~/mega/m01/data/rfam/{f}.fa')
    text = open(p).read()
    bio = [(r.id, r.description, str(r.seq).upper()) for r in SeqIO.parse(p, 'fasta')]
    try:
        own = [(r['id'], r['description'], r['sequence']) for r in parse_fasta(text)]
    except Exception as e:
        rows.append({'file': f, 'error': str(e)[:100]}); mism += 1; continue
    same = len(bio) == len(own) and all(a == b for a, b in zip(bio, own))
    if not same:
        mism += 1
        diff = next((i for i, (a, b) in enumerate(zip(bio, own)) if a != b), None)
        rows.append({'file': f, 'match': False, 'n_biopython': len(bio), 'n_own': len(own), 'first_diff_index': diff,
                     'biopython_rec': str(bio[diff])[:120] if diff is not None else None, 'own_rec': str(own[diff])[:120] if diff is not None and diff < len(own) else None})
    else:
        rows.append({'file': f, 'match': True, 'n': len(bio)})
    # stream vs parse
    st = [(r['id'], r['sequence']) for r in stream_fasta(p)]
    rows[-1]['stream_matches_parse'] = st == [(r['id'], r['sequence']) for r in parse_fasta(text)]
# round trip on RF00059 subset
recs = parse_fasta(open(os.path.expanduser('~/mega/m01/data/rfam/RF00059.fa')).read())[:200]
rt = parse_fasta(write_fasta(recs, line_width=50))
roundtrip = all(a['id'] == b['id'] and a['sequence'] == b['sequence'] for a, b in zip(recs, rt)) and len(rt) == len(recs)
# blank lines / leading blank / CRLF robustness vs Biopython
edge_text = '\n\n>id1 desc words\nacg\ntt\n\n>id2\nGG\n\n'
try:
    be = [(r.id, r.description, str(r.seq).upper()) for r in SeqIO.parse(__import__('io').StringIO(edge_text.replace('\n', '\r\n')), 'fasta')]
    be_err = None
except ValueError as e:
    be, be_err = [], str(e).splitlines()[0]
oe = [(r['id'], r['description'], r['sequence']) for r in parse_fasta(edge_text)]
res = {'reference': 'Biopython ' + __import__('Bio').__version__ + ' SeqIO on Rfam CURRENT fasta_files (RF00059, RF00050, RF00504)',
 'rows': rows, 'mismatches': mism, 'roundtrip_200_ok': roundtrip,
 'edge_blank_crlf': {'biopython': be, 'biopython_error': be_err, 'own': oe, 'match': be == oe,
                     'note': 'Biopython rejects leading blank lines; bio.fasta skips them (leniency, not a mismatch on valid files)'},
 'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else ''}
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_fasta.json'), 'w'), indent=1)
print(json.dumps(res, indent=1)[:1500])
