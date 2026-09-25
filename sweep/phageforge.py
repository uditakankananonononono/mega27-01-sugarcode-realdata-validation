"""Sweep: phageforge vs real NDM-1 plasmid FN396876 + host chromosome screen + recomputation."""
import json, math, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'sc-ai', 'src'))
from sugarcode.bio import entrez
from sugarcode.bio.sequence import reverse_complement
from sugarcode.modules.phageforge import core as PF
from sugarcode.modules.crispr_opt import design_guides, score_off_targets_cfd

res = {'module': 'phageforge', 'sources': [
    'NCBI nuccore FN396876.1 (pKpANDM-1, blaNDM-1 plasmid)',
    'NCBI nuccore NC_009648.1 (K. pneumoniae MGH 78578 chromosome, host background)',
    'NCBI nuccore NC_001604/NC_001416/NC_000866/NC_005856/NC_003287 (T7/lambda/T4/P1/M13 genomes)',
    'recomputation']}

# 1. real blaNDM-1 CDS from FN396876.1 (GBFeature CDS complement(2407..3219))
fa = entrez.efetch_fasta('nuccore', 'FN396876')
seq = ''.join(l.strip() for l in fa.split('\n')[1:] if l)
cds = reverse_complement(seq[2406:3219])
res['target'] = {'plasmid': 'FN396876.1', 'plasmid_bp': len(seq),
                 'cds_location': 'complement(2407..3219)', 'cds_bp': len(cds),
                 'starts_atg': cds[:3] == 'ATG', 'stop': cds[-3:],
                 'internal_stops': sum(1 for i in range(0, len(cds) - 3, 3) if cds[i:i+3] in ('TAA', 'TAG', 'TGA'))}

# 2. design_phage on the real CDS (default background = target gene)
out = PF.design_phage(cds)
guides_check = []
for g in out['payload']['guides']:
    ok = (cds[g['start']:g['end']] == g['guide']) if g['strand'] == '+' else (reverse_complement(cds[g['start']:g['end']]) == g['guide'])
    guides_check.append({'guide': g['guide'], 'strand': g['strand'], 'start': g['start'],
                         'in_cds': bool(ok), 'len20': len(g['guide']) == 20,
                         'on_model': g['on_target_model'], 'off_risk': g['off_target_risk']})
res['design_default'] = {'backbone': out['backbone'], 'fits': out['payload']['fits'],
                         'payload_kb': out['payload']['payload_kb'],
                         'microbiome_sparing': out['specificity']['microbiome_sparing'],
                         'background_screened': out['specificity']['background_screened'],
                         'guides': guides_check}
cands = {k: v for k, v in PF.PHAGE_BACKBONES.items() if v['capacity_kb'] >= 4.5 and v['lytic']}
pick = max(cands, key=lambda k: cands[k]['capacity_kb'])
res['backbone_logic'] = {'eligible': sorted(cands), 'expected_pick': pick,
                         'module_pick': out['backbone'], 'ok': out['backbone'] == pick}

# 3. BUG 51: perfect-match off-target was invisible pre-fix; reported post-fix
G = out['payload']['guides'][0]['guide']
bg_syn = 'T' * 50 + G + 'AGG' + 'T' * 50
hits = score_off_targets_cfd(G, bg_syn)
res['bug51_perfect_match'] = {
    'pre_fix_observed': 'score_off_targets_cfd skipped every mm==0 site: the same planted perfect-match background returned [] pre-fix (verified 2026-09-25, sugarcode 102d187)',
    'post_fix_hits': [{'position': h['position'], 'mismatches': h['mismatches'], 'cfd': h['cfd_score'], 'risk': h['risk']} for h in hits],
    'post_fix_reports_perfect': any(h['mismatches'] == 0 and h['cfd_score'] == 1.0 for h in hits)}

# 4. real host-chromosome background screen (40 kb window)
t0 = time.time()
data = entrez._get('efetch.fcgi', {'db': 'nuccore', 'id': 'NC_009648.1', 'rettype': 'fasta',
                                   'retmode': 'text', 'seq_start': '1000001', 'seq_stop': '1040000'})
host = ''.join(l.strip() for l in data.decode().split('\n') if not l.startswith('>'))
out_bg = PF.design_phage(cds, background_genome=host)
res['host_screen'] = {'background': 'NC_009648.1:1000001-1040000 (40 kb K. pneumoniae MGH 78578 chromosome)',
                      'bg_bp': len(host), 'screen_seconds': round(time.time() - t0, 1),
                      'background_screened': out_bg['specificity']['background_screened'],
                      'guides': [{'guide': g['guide'], 'off_risk': g['off_target_risk'],
                                  'n_top_off': len(g['top_off_targets'])} for g in out_bg['payload']['guides']]}

# 5. exact recomputations
kc = PF.kill_curve(3, adsorption_rate=0.7, burst_size=80, resistant_fraction=0.05)
exp_inf = 1 - math.exp(-0.7 * 3)
res['kill_curve'] = {'module': kc, 'infected_recomputed': exp_inf,
                     'surviving_recomputed': (1 - exp_inf) * 0.95 + 0.05,
                     'burst_recomputed': exp_inf * 80}
er = PF.escape_risk(3, mutation_rate=1e-8, population_size=1e10)
exp_e = 1e10 * (1e-8 ** 3)
res['escape_risk'] = {'module': er, 'expected_recomputed': exp_e,
                      'prob_recomputed': 1 - math.exp(-exp_e)}
pb = PF.payload_budget('P1', {'cas9': 4.1, 'guide_array': 0.3, 'selection': 0.1})
res['payload_budget'] = {'module': pb, 'total_recomputed': 4.5, 'headroom_recomputed': 15.0 - 4.5}
cc = PF.cocktail_coverage({'phageA': ['Kp_ST258', 'Kp_ST11'], 'phageB': ['Kp_ST11', 'Ec_O157']},
                          target_hosts=['Kp_ST258', 'Kp_ST11', 'Ec_O157', 'Sa_MRSA'])
res['cocktail'] = {'module': cc, 'coverage_recomputed': 0.75, 'uncovered_recomputed': ['Sa_MRSA'],
                   'legacy_fraction_without_targets': PF.cocktail_coverage({'phageA': ['a'], 'phageB': ['a', 'b']})['coverage_fraction']}

# 6. backbone capacities vs real genome sizes (NCBI-verified this sweep)
accs = {'T7': 'NC_001604', 'lambda': 'NC_001416', 'T4': 'NC_000866', 'P1': 'NC_005856', 'M13': 'NC_003287'}
ids = {n: entrez.esearch('nuccore', f'{a}[accn]', retmax=1)[0] for n, a in accs.items()}
summ = entrez.esummary('nuccore', list(ids.values()))
inv = {v: k for k, v in ids.items()}
genomes = {inv[u]: summ[u]['slen'] for u in ids.values()}
res['backbone_genomes'] = {n: {'capacity_kb_module': PF.PHAGE_BACKBONES[n]['capacity_kb'],
                               'genome_bp_ncbi': genomes[n],
                               'genome_kb': round(genomes[n] / 1000, 3)} for n in accs}
res['backbone_note'] = ('capacity_kb is payload headroom, not genome size: M13 capacity equals its '
                        '6407 bp genome (filamentous phages package extended genomes); lytic-phage '
                        'headrooms (T7 4.0, lambda 8.0, T4 10.0, P1 15.0 kb) are plausible packaging '
                        'headroom values but are not individually literature-sourced - documented caveat')
res['microbiome_sparing'] = {'module': out['specificity']['microbiome_sparing'],
                             'recomputed': round(0.6 + 0.3 * (1 - 0.5 * 0.7), 2)}

fn = os.path.join(os.path.dirname(__file__), '..', 'benchmarks', 'sweep_phageforge.json')
json.dump(res, open(fn, 'w'), indent=1, default=str)
print(json.dumps(res, default=str)[:3200])
