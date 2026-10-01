# CellPainter similarity depends on feature dictionary order

Product commit: 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

nearest_mechanism constructs numpy arrays from dictionary .values() for both query/reference without matching keys. A dna_damage reference gives cosine 1.0 against the identical query. Reversing only the reference dictionary insertion order, while dictionary equality remains true, gives 0.0. Repeat calls reproduce that result.

This is a named-feature alignment defect. The feature mapping is unchanged; a serialization/import order should not alter the numerical similarity. Existing tests construct references in the same ordering, so the passing product suite does not exercise this case. cellpainting_report also converts dictionary values into a matrix and should be reviewed for the same issue; only nearest_mechanism is independently proven here.

Expected correction direction: align vectors by explicit shared feature keys and define missing-feature handling, then run the regression and existing suite. No product edits performed. The standalone regression is audit evidence expected to fail at the pinned source.
