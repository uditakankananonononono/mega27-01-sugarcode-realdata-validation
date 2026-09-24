"""Harvest real GT-AG splice junction windows from RefSeqGene GenBank records.

Held-out design: the caller passes accessions disjoint from the module's
training set (its PROVENANCE lists the 29 used). Donor window: 3 exonic + 6
intronic bases (9-mer). Acceptor window: 20 intronic + 3 exonic bases (23-mer).
"""
import re
from dataclasses import dataclass


@dataclass
class Junction:
    accession: str
    kind: str          # donor | acceptor
    window: str
    intron_len: int


def parse_cds_exons(feature_line_block):
    """Parse join(...) locations from a GenBank CDS feature block -> exon spans."""
    m = re.search(r"join\(([^)]*)\)", feature_line_block.replace("\n", ""))
    if not m:
        return []
    spans = []
    for part in m.group(1).split(","):
        part = part.strip()
        if part.startswith(("complement", ">", "<")):
            part = part.replace("complement", "").strip("()<>")
        a, b = part.split("..")[-2:]
        spans.append((int(a.replace("<", "").replace(">", "")),
                      int(b.replace("<", "").replace(">", ""))))
    return sorted(spans)


def junction_windows(seq, exons, strand=1):
    """Extract donor/acceptor windows from exon boundaries (1-based inclusive)."""
    out = []
    for (a1, b1), (a2, b2) in zip(exons, exons[1:]):
        intron = seq[b1:a2 - 1] if strand == 1 else None
        if intron is None or len(intron) < 26:
            continue
        donor = seq[b1 - 3:b1 + 6]
        acceptor = seq[a2 - 21:a2 + 2]
        if intron[:2] == "GT" and intron[-2:] == "AG" and "N" not in donor + acceptor:
            out.append(Junction("", "donor", donor, len(intron)))
            out.append(Junction("", "acceptor", acceptor, len(intron)))
    return out


def acceptor15(window23):
    """deepsplice's 15-nt acceptor frame from the 23-nt MaxEntScan-style
    window (20 intron + 3 exon): 12 intron nt + AG + first exon nt, i.e.
    window[6:21]. The first validation run sliced [5:20] (AG at 13-14, no
    exon base), which scored every true acceptor off-frame (AUC 0.547)."""
    if len(window23) != 23:
        raise ValueError("acceptor window must be 23 nt")
    return window23[6:21]


def intronic_decoys(seq, exons, flank=30):
    """Real intronic AG (acceptor frame) and GT (donor frame) sites at least
    `flank` nt from either intron end: hard negatives that share the
    invariant dinucleotide, unlike shuffled windows."""
    acc, don = [], []
    for (a1, b1), (a2, b2) in zip(exons, exons[1:]):
        intron = seq[b1:a2 - 1]
        if intron[:2] != "GT" or intron[-2:] != "AG":
            continue
        for i in range(b1 + flank, a2 - 1 - flank):
            if seq[i:i + 2] == "AG":
                acc.append(seq[i - 12:i + 3])
            if seq[i:i + 2] == "GT":
                don.append(seq[i - 3:i + 6])
    return acc, don
