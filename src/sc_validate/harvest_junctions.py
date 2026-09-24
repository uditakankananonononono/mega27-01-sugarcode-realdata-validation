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
