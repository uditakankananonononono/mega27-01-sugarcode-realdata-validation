# Liquid untrained-attention dimension boundary
Product 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

For a two-position/four-channel input, positive heads=3 and heads=5 pass the explicit positive-count gate but crash with matrix dimension mismatch. d_model=max(8,ceil(channels/heads)*heads) yields 8, not divisible by these heads. Default heads=4 returns. Independent contract regressions twice: 2 fail/1 pass. No product changes.

Also recorded separately in evidence: all four missing modalities produce lung tissue-of-origin and cancer_probability=0.075858 with zero reliability. This is a no-data default, not a measured cancer result. Existing research-use disclaimer does not establish diagnostic validity. No real patient data involved.
