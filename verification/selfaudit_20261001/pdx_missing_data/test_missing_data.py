"""Validation evidence expected to fail at pinned product 1ccbd030.

Missing observations must not be labelled measured. This test does not
require a particular numerical fallback or impose clinical calibration.
"""
from sugarcode.modules.pdx_insight import multiomic_fidelity

def test_omitted_mutations_are_not_labelled_measured():
    result = multiomic_fidelity({}, {}, passage=0)
    assert 'mutations' not in result['layers_measured']
