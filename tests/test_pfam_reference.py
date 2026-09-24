from sc_validate.pfam_reference import pfam_globin_reference


def test_pfam_globin_separates_real_fixtures():
    r = pfam_globin_reference()
    assert r["n_pos"] == 25 and r["n_neg"] == 25
    assert r["auc"] > 0.95
