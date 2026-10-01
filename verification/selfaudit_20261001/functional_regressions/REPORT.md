# Promoter Lib wrong homopolymer-length output

Product commit: 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

sequence_features computes homopolymer_max by finding the longest substring left after splitting the sequence separately on each base. That is not a maximum contiguous repeated-base run. Four of five exact test cases disagree with an independent itertools.groupby calculation:
- ACGTACGT: module 3, expected 1.
- AAAAAAAA: module 8, expected 8 (positive control).
- AAAACCCC: module 8, expected 4.
- ATATATAT: module 8, expected 1.
- AAAAACGT: module 7, expected 5.

Repeat calls reproduce the same errors. Existing test_promoter_lib_spec.py asserts only homopolymer_max >= 1, so the 2,319-pass suite does not catch this defect. promoter_report uses sequence_features and propagates this feature.

The attached standalone regression file belongs to validation evidence and is expected to fail on the pinned source. No product mutation was performed. Product owner should implement contiguous-run length and run these exact controls plus existing suite checks.
