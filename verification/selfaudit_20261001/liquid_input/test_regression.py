"""Independent API input validity contracts; these intentionally fail at audit HEAD."""
import pytest
from sugarcode.modules.liquid_biopsy import detect_ctdna
@pytest.mark.parametrize('signal,depth', [([-.1]*8,30000),([2]*8,30000),([float('nan')]*8,30000),([.1]*8,-100)])
def test_invalid_allele_frequencies_and_depth_rejected(signal,depth):
    with pytest.raises(ValueError): detect_ctdna(signal,depth=depth)
def test_valid_trace_positive_control():
    r=detect_ctdna([0,0,0,0,.1,0,0,0])
    assert r['ctdna_detected'] and 0<=r['estimated_sensitivity']<=1
