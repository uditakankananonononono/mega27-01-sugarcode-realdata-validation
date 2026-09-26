# CRISPRscan AA18/AA19 resolution

CRISPOR paramsCRISPRscan places the largest negative CRISPRscan weight on AA18;
the original Moreno-Mateos 2015 supplement (MOESM640) places it on AA19.
90/91 features match at identical positions and weights, so conventions align and
this is a one-position misplacement in CRISPOR's Excel conversion. crisprScore
(Bioconductor) and sugarcode match the original supplement.

Upstream report filed 2026-09-26 with owner approval:
https://github.com/maximilianh/crisporWebsite/issues/76
