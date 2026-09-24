"""The frozen oracle must reproduce the CRISPOR reference doctests exactly."""
from sc_validate.cfd_oracle import cfd_score


def test_reference_doctest_values():
    assert abs(cfd_score("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGAAAGGG") - 0.4635989007074176) < 1e-12
    assert cfd_score("GGGGGGGGGGGGGGGGGGGGGGG", "GGGGGGGGGGGGGGGGGGGGGGG") == 1.0
    assert abs(cfd_score("ATGGTCGGACTCCCTGCCAGAGG", "ATGGTGGGACTCCCTGCCAGAGG") - 0.5) < 1e-12
    assert abs(cfd_score("ATGTGGAGATTGCCACCTACCGG", "ATCTGGAGATTGCCACCTACAGG") - 0.384615385) < 1e-9


def test_frozen_fixture_self_consistent():
    import csv, pathlib
    p = pathlib.Path(__file__).resolve().parent.parent / "data" / "crispr_cfd_oracle.tsv"
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    assert len(rows) == 600
    for r in rows[:50]:
        assert abs(cfd_score(r["guide"], r["offtarget"]) - float(r["cfd"])) < 1e-9
