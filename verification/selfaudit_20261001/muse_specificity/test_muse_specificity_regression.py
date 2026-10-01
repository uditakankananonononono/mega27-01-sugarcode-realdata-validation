import math
from sugarcode.modules.crispr_muse import reward_decomposition,repair_outcomes
G='GAGTCCGAGCAGAAGAAGAA'
def test_exact_copy_is_observed():
 r=reward_decomposition(G,background=G+'AGG')
 assert r['off_target_hits']>=1

def test_exact_copy_penalizes_specificity():
 none=reward_decomposition(G);copy=reward_decomposition(G,background=G+'AGG')
 assert copy['components']['specificity']<none['components']['specificity']

def test_chromatin_condition_changes_repair_probabilities():
 low=repair_outcomes(G*2,chromatin=.1);high=repair_outcomes(G*2,chromatin=.9)
 assert any(not math.isclose(low[k],high[k]) for k in ('NHEJ','MMEJ','HDR'))
