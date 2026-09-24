"""Recompute per-taxon significance of rho(CAI_ribo, abundance) at full precision.
Fisher z: z = atanh(rho) * sqrt(n-3); two-sided p; Benjamini-Hochberg (statsmodels)."""
import json, math
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
d = json.load(open('benchmarks/sweep_codon_cai_multispecies.json'))
T = [t for t in d['per_taxon'] if t.get('rho_ribo_ref') is not None]
p = [2 * norm.sf(abs(math.atanh(t['rho_ribo_ref']) * math.sqrt(t['n_matched'] - 3))) for t in T]
rej, q, _, _ = multipletests(p, alpha=0.05, method='fdr_bh')
out = {'method': 'Fisher z two-sided, BH-FDR (statsmodels multipletests fdr_bh), full precision',
       'n': len(T), 'n_significant_q05': int(rej.sum()),
       'per_taxon': [{'taxid': t['taxid'], 'organism': t['organism'], 'rho_ribo_ref': t['rho_ribo_ref'], 'n': t['n_matched'], 'p': pi, 'q': float(qi)} for t, pi, qi in zip(T, p, q)]}
json.dump(out, open('benchmarks/sweep_codon_bh_fdr.json', 'w'), indent=1)
print(out['n'], out['n_significant_q05']); [print(r['organism'][:30], '%.3g' % r['q'], d['bh_fdr_rho_ribo_gt0']['per_taxon_q'].get(r['organism'])) for r in out['per_taxon'][:6]]
