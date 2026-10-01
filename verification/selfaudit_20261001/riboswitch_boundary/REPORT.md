# Riboswitch allowed input crashes; design mode has no effect

Product commit: 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

response_curve validates 0 <= basal < maximum, then computes dynamic_range=maximum/basal without zero handling. basal=0, Kd=1, concentrations=[0,1,10] throws ZeroDivisionError. Default basal=.05 succeeds, reports dynamic range 20 and half-occupancy signal .525. Independent boundary regression fails on the unexpected arithmetic exception; it does not mandate an infinity/null fallback versus a documented validation rejection.

optimize_switch accepts and echoes on/off mode but uses the same sequence search and objective for either. Theophylline with four candidates produces identical rankings and selected sequence for both modes. This is a confirmed absence of design-mode effect, not evidence that either generated sensor works. The downstream response curve does use mode for signal direction.

No product changes or biochemical action. Product owner should define the zero-basal contract, handle it safely and either implement mode-specific design or label search mode-independent. No biological efficacy established.
