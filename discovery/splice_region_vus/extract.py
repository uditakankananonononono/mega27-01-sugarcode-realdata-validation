"""Build transcript-oriented ref/alt windows for ClinVar splice-region SNVs
(donor +3..+6, acceptor -3..-14) from hg38; strand inferred from the canonical
GT/AG dinucleotide at the implied natural site."""
import csv, re, sys, twobitreader
RC = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s.translate(RC)[::-1]
g = twobitreader.TwoBitFile("hg38.2bit")
F = 40
out = open("windows.tsv", "w")
out.write("vid\tgene\tname\tcls\tsite\tk\tstrand\tref_ctx\talt_ctx\n")
kept = drop = 0
for r in csv.reader(open("splice_region_snvs.tsv"), delimiter="\t"):
    vid, gene, name, sig, rev, chrom, pos, ref, alt = r
    m = re.search(r"c\.[0-9]+([+-][0-9]+)", name)
    if not m or not pos.isdigit() or len(ref) != 1 or len(alt) != 1: continue
    k = int(m.group(1))
    if not ((3 <= k <= 6) or (-14 <= k <= -3)): continue
    cls = ("P" if sig in ("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic") else
           "B" if sig in ("Benign", "Likely benign", "Benign/Likely benign") else
           "VUS" if sig == "Uncertain significance" else None)
    if cls is None: continue
    ch = "chrM" if chrom == "MT" else "chr" + chrom
    try:
        p = int(pos) - 1  # 0-based
        seq = g[ch][p - F:p + F + 1].upper()
    except Exception:
        drop += 1; continue
    if len(seq) != 2 * F + 1 or seq[F] != ref:
        drop += 1; continue
    alt_seq = seq[:F] + alt + seq[F + 1:]
    ok = []
    for strand in "+-":
        s = seq if strand == "+" else rc(seq)
        # in transcript orientation the variant is always at index F
        if k > 0:   # donor: exon last base at F-k, intron starts F-k+1
            ok.append(s[F - k + 1:F - k + 3] == "GT")
        else:       # acceptor: intron ends at F-k-1 (k negative => F+|k|-1)
            e = F - k  # first exon base index
            ok.append(s[e - 2:e] == "AG")
    if ok.count(True) != 1:
        drop += 1; continue
    strand = "+" if ok[0] else "-"
    rs, as_ = (seq, alt_seq) if strand == "+" else (rc(seq), rc(alt_seq))
    out.write(f"{vid}\t{gene}\t{name}\t{cls}\t{'donor' if k>0 else 'acceptor'}\t{k}\t{strand}\t{rs}\t{as_}\n")
    kept += 1
print("kept", kept, "dropped", drop)
