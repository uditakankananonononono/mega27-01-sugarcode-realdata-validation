"""Real-world reference model for profile-HMM validation: Pfam PF00042.29
(Globin) HMMER3 model, scored with pyhmmer (HMMER 3 bindings) on the frozen
globin / non-globin FASTA fixtures. Output: per-sequence bit scores + AUC."""
import pathlib
import pyhmmer
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile, Alphabet
from .auc import auc

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"


def _scores(hmm, fasta, alphabet):
    with SequenceFile(str(fasta), digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()
    best = {s.name.decode() if isinstance(s.name, bytes) else s.name: -1e9 for s in seqs}
    for hits in pyhmmer.hmmsearch([hmm], seqs, E=1e9, domE=1e9):
        for h in hits:
            n = h.name.decode() if isinstance(h.name, bytes) else h.name
            best[n] = max(best[n], h.score)
    return best


def pfam_globin_reference():
    abc = Alphabet.amino()
    with HMMFile(str(DATA / "pfam_PF00042_globin.hmm")) as hf:
        hmm = hf.read()
    pos = _scores(hmm, DATA / "globins25.fa", abc)
    neg = _scores(hmm, DATA / "nonglobins25.fa", abc)
    return {"model": "Pfam PF00042.29 (HMMER3)", "n_pos": len(pos), "n_neg": len(neg),
            "auc": auc(list(pos.values()), list(neg.values())),
            "min_pos_bits": round(min(pos.values()), 2),
            "neg_with_any_hit": sum(v > -1e8 for v in neg.values()),
            "max_neg_bits": (round(max(neg.values()), 2) if max(neg.values()) > -1e8 else None)}
