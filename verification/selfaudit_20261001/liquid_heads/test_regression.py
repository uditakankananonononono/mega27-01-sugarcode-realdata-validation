import pytest
from sugarcode.modules.liquid_biopsy import transformer_denoise
@pytest.mark.parametrize('heads',[3,5])
def test_positive_head_count_produces_attention(heads):
    a=transformer_denoise([[[.01,.95,.95,.95],[.1,.95,.95,.95]]],heads=heads)
    assert a['architecture']['heads']==heads
def test_default_positive_control():
    assert transformer_denoise([[[.01,.95,.95,.95],[.1,.95,.95,.95]]])['architecture']['heads']==4
