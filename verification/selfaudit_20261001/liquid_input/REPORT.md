# Liquid Biopsy legacy input gate
Product HEAD 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

The exported detect_ctdna API accepts negative or >1 allele fractions, NaNs, and negative depth instead of rejecting them. Negative depth yields negative confidence (-0.105) and negative estimated sensitivity (-0.016). All-NaN input is reported as not detected, conflating invalid data with a negative call. Statistical demonstration only, not a real blood/cancer test.

Independent input-rejection contracts repeated twice: 4 failed, 1 valid-input positive-control passed on each run. Failures do not modify the shipped suite's 2,319 passing count. Only the legacy exported entrypoint was checked; these results do not prove every new pipeline entrypoint shares the issue. Raw vs filtered visualization pixels remain uninspected.
