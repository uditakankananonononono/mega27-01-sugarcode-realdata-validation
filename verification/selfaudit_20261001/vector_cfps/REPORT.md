# Vector receptor identity and CFPS initial mass
Product 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

Vector Opt target_receptor is returned as a label but never used in variant scoring. AAVR and heparan_sulfate produce identical variants/lead at seed=1. AAV2/AAV9 numeric variants also match after ID stripping: CAPSID_REGIONS is static despite comment claiming serotype-dependent contact handling. No actual capsid sequence/coordinate/docking calculation in this entrypoint.

Cell-Free Opt kinetics increments rate*dt before recording time 0. hours=0 yields 0.25 g/L produced, rather than zero; the complete zero-duration call reports final_g_l=0.25. Baseline assumes no pre-existing product.

Independent regressions twice: 2 fail/1 seeded positive control passes. Model consistency only, no wet-lab action. Product unchanged.
