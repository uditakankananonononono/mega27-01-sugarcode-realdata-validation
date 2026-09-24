from sc_validate.metrics import binary_metrics, pearson
from sc_validate.auc import auc


def test_perfect_classifier():
    m = binary_metrics([True, True, False, False], [True, True, False, False])
    assert m["accuracy"] == 1.0 and m["sensitivity"] == 1.0 and m["specificity"] == 1.0


def test_pearson_against_known_value():
    xs = [1, 2, 3, 4, 5]
    ys = [2, 4, 5, 4, 5]
    assert abs(pearson(xs, ys) - 0.7745966692) < 1e-6


def test_auc_extremes():
    assert auc([3, 4], [1, 2]) == 1.0
    assert auc([1, 2], [3, 4]) == 0.0
    assert auc([1, 3], [2, 4]) == 0.25  # (3>2) is the only win of 4 pairs
