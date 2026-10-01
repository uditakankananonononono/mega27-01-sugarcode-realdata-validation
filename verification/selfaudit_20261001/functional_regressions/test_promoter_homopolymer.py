"""Independent regression evidence; expected to fail at product 1ccbd030.

Run against the pinned product source with PYTHONPATH pointing to its src/.
This file is audit evidence, not a modification to the product test suite.
"""
import itertools
import pytest
from sugarcode.modules.promoter_lib import sequence_features

@pytest.mark.parametrize("sequence", ["ACGTACGT", "AAAAAAAA", "AAAACCCC", "ATATATAT", "AAAAACGT"])
def test_homopolymer_max_equals_independent_contiguous_run(sequence):
    expected = max(len(list(group)) for _, group in itertools.groupby(sequence))
    assert sequence_features(sequence)["homopolymer_max"] == expected
