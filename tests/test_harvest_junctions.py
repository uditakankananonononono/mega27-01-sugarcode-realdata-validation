from sc_validate.harvest_junctions import junction_windows


def test_windows_around_real_junction():
    seq = list("G" * 600)
    # exons 1-based inclusive: (1,100) (150,260); intron 101..149
    intron = "GT" + "C" * 45 + "AG"
    for i, ch in enumerate(intron, start=101):
        seq[i - 1] = ch
    js = junction_windows("".join(seq), [(1, 100), (150, 260)])
    assert len(js) == 2
    donor = next(j for j in js if j.kind == "donor")
    acc = next(j for j in js if j.kind == "acceptor")
    assert len(donor.window) == 9 and donor.window[3:5] == "GT"
    assert len(acc.window) == 23 and acc.window[18:20] == "AG"
    assert donor.intron_len == 49
