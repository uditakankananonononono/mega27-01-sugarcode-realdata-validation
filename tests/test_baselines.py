from sc_validate.baselines import load_skempi_subset, skempi_baselines


def test_skempi_baselines_real_fixture():
    rows = load_skempi_subset()
    assert len(rows) == 150
    b = skempi_baselines(rows)
    assert b["n"] == 150
    assert b["rmse_mean_predictor"] > 0
    assert -1 <= b["pearson_abs_volume_change"] <= 1
