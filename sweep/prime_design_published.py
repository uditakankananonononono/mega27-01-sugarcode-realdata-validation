"""prime_design.design_pegrna/design_edit vs published, experimentally validated pegRNAs:
HEK3 +1 CTT ins (Anzalone 2019; Addgene 132778; also Chow 2021 pegFinder, validated in HEK293T)
and RNF2 +5 G>T (Anzalone 2019). Sources: Frontiers 2025 fcell.1589034 Data Sheet 1 Table 2;
Chow et al. Nat Biomed Eng 2021 (s41551-020-00622-8) methods; Anzalone 2019 suppl (substrate AVA024
gives the HEK3 genomic context ...GGCCCAGACTGAGCACGTGATGGCAGAGGAAAGG)."""
import sys, os, json
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from sugarcode.modules.prime_design import design_pegrna, design_edit
from sugarcode.bio.sequence import reverse_complement as rc
CASES = [
 {'name': 'HEK3 +1 CTT ins', 'spacer': 'GGCCCAGACTGAGCACGTGA', 'rtt': 'TCTGCCATCAAAGC', 'pbs': 'GTGCTCAGTCTG'},
 {'name': 'RNF2 +5 G>T', 'spacer': 'GTCATCTTAGTCATTACCTG', 'rtt': 'AACGAACACATCAGG', 'pbs': 'TAATGACTAAGATG'},
]
rows = []
for c in CASES:
    d = design_pegrna(c['spacer'], rc(c['rtt']), pbs_len=len(c['pbs']), rtt_len=len(c['rtt']))
    rows.append({'name': c['name'], 'rtt_match': d['rt_template'] == c['rtt'], 'pbs_match': d['pbs'] == c['pbs'],
                 'extension_match': d['pegrna_3p_extension'] == c['rtt'] + c['pbs'],
                 'got_pbs': d['pbs'], 'got_rtt': d['rt_template'],
                 'pre_fix_pbs': c['rtt'][:len(c['pbs'])]})  # old behaviour: PBS from the RT-template region
# end-to-end: design the HEK3 +1 CTT edit from the genomic context (Anzalone substrate AVA024)
region = 'CCTGGGTCAATCCTTGGGGCCCAGACTGAGCACGTGATGGCAGAGGAAAGG'
r = design_edit(region, {'type': 'insertion', 'position': 34, 'ref': '', 'alt': 'CTT'})
top = r['pegrna_designs'][0] if r['pegrna_designs'] else {}
res = {'reference': 'Anzalone 2019 Nature suppl + Addgene 132778 + Chow 2021 Nat Biomed Eng (experimentally validated)',
 'unit_cases': rows, 'all_unit_match': all(x['rtt_match'] and x['pbs_match'] and x['extension_match'] for x in rows),
 'end_to_end_hek3': {'top_spacer': top.get('spacer'), 'spacer_match': top.get('spacer') == CASES[0]['spacer'],
   'pbs': top.get('pbs'), 'pbs_match_published': top.get('pbs') == CASES[0]['pbs'],
   'auto_pbs_len': top.get('pbs_length'), 'auto_rtt_len': top.get('rtt_length'),
   'rtt_contains_published_at_nick': (top.get('rt_template') or '').endswith(CASES[0]['rtt'])},
 'sugarcode_commit': sys.argv[1] if len(sys.argv) > 1 else ''}
res['pre_fix_pbs_matched_published'] = [x['pre_fix_pbs'] == rows[i]['got_pbs'] for i, x in enumerate(rows)]
json.dump(res, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_prime_design.json'), 'w'), indent=1)
print(json.dumps(res, indent=1))
