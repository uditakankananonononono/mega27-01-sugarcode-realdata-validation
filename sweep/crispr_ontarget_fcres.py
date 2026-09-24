"""crispr_opt on-target scorers vs measured activity in Doench 2016 FC+RES (Azimuth repo)."""
import csv, json, os, sys
import numpy as np
from scipy.stats import spearmanr
sys.path.insert(0, os.path.expanduser('~/mega/sc-ai/src'))
from sugarcode.modules.crispr_opt.core import doench2014_ontarget, score_on_target
from sugarcode.modules.crispr_opt.rule_set_2 import score_guides
SRC = 'https://raw.githubusercontent.com/MicrosoftResearch/Azimuth/master/azimuth/data/FC_plus_RES_withPredictions.csv'
r = list(csv.DictReader(open(os.path.expanduser('~/mega/m01/data/doench2016/FC_plus_RES_withPredictions.csv'))))
seqs = [x['30mer'] for x in r]; y = np.array([float(x['score_drug_gene_rank']) for x in r])
pp = [float(x['Percent Peptide']) for x in r]; ac = [float(x['Amino Acid Cut position']) for x in r]
genes = sorted({x['Target gene'] for x in r})
S = {'rule_set_2_with_position': score_guides(seqs, pp, ac), 'rule_set_2_sequence_only': score_guides(seqs),
     'rule_set_1_doench2014': np.array([doench2014_ontarget(s) for s in seqs]),
     'heuristic_edge_fallback': np.array([score_on_target(s[4:24]) for s in seqs]),
     'azimuth_repo_predictions_column': np.array([float(x['predictions']) for x in r])}
out = {'source': SRC, 'n_guides': len(r), 'n_genes': len(genes), 'target': 'score_drug_gene_rank (measured, rank-normalised within gene)', 'spearman': {}, 'per_gene_median_spearman': {}}
for k, v in S.items():
    out['spearman'][k] = round(float(spearmanr(v, y).correlation), 4)
    pg = [spearmanr(v[[i for i, x in enumerate(r) if x['Target gene'] == g]], y[[i for i, x in enumerate(r) if x['Target gene'] == g]]).correlation for g in genes]
    out['per_gene_median_spearman'][k] = round(float(np.nanmedian(pg)), 4)
out['caveat'] = 'Rule Set 1 and Rule Set 2 were trained on these guides (FC for RS1, FC+RES for RS2): the correlations are in-sample and optimistic; the heuristic was not fitted to them.'
json.dump(out, open(os.path.expanduser('~/mega/m01/benchmarks/sweep_crispr_ontarget_fcres.json'), 'w'), indent=1)
print(out)
