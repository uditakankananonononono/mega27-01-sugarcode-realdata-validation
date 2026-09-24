"""Standalone CFD off-target oracle, extracted from the public CRISPOR
reference distribution (Concordet & Haeussler 2018, Genome Biol 19:171;
github.com/maximilianh/crisporWebsite, CFD_Scoring pickles retrieved
2026-09-24). Doctest values reproduced exactly before freezing."""
import pickle, re, pathlib

_D = pathlib.Path(__file__).resolve().parents[2] / "data" / "CFD_Scoring"
_mm = pickle.load(open(_D / "mismatch_score.pkl", "rb"))
_pam = pickle.load(open(_D / "pam_scores.pkl", "rb"))


def _revcom(s):
    b = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}
    return "".join(b[x] for x in s[::-1])


def cfd_score(guide23, ot23):
    """CRISPOR calcCfdScore semantics: 23-mers (20nt guide + NGG PAM)."""
    wt = guide23.upper()
    off = ot23.upper()
    if re.search("[^ATCG]", wt) or re.search("[^ATCG]", off):
        return -1
    wtu = wt[:20].replace("T", "U")
    sgu = off[:20].replace("T", "U")
    score = 1.0
    for i, (w, s) in enumerate(zip(wtu, sgu)):
        if w != s:
            score *= _mm["r" + w + ":d" + _revcom(s) + "," + str(i + 1)]
    return score * _pam[off[-2:]]
