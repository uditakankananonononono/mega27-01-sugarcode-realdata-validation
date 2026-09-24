from sc_validate.harvest_junctions import acceptor15, intronic_decoys, junction_windows, parse_cds_exons


def test_acceptor15_puts_ag_at_12():
    w = "CTACCTTTCCCCCACCCCAGGTC"  # real NG_008021.1 acceptor window
    s = acceptor15(w)
    assert len(s) == 15 and s[12:14] == "AG" and s[14] == "G"


def test_decoys_share_dinucleotide_and_frame():
    seq = "A" * 10 + "CCCGTAAGT" + "T" * 40 + "AG" + "C" * 40 + "TTTTTTTTTTTCAG" + "GCC" + "A" * 10
    exons = [(1, 13), (seq.find("TTTTTTTTTTTCAG") + 15, len(seq))]
    acc, don = intronic_decoys(seq, exons)
    assert acc and all(a[12:14] == "AG" and len(a) == 15 for a in acc)
