# PDX Insight absent data reported as perfect measured fidelity

Product commit: 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai

multiomic_fidelity({}, {}, passage=0) returns Model Fidelity Index 1.0 and layers_measured=['mutations'] with no mutation or other omics observations supplied. Explicitly empty mutation lists return the same score. A nonempty identical-pair positive control also returns 1.0, so the score does not distinguish these evidence states. All three cases reproduce exactly.

Source defaults missing mutations to empty lists and assigns retention 1 when the patient set is empty. It always includes the mutation layer's weight and calls it measured. Distinguish omitted/unassessed mutations from a verified zero-mutation assay; neither should silently establish clinically meaningful perfect fidelity. The narrow independent regression checks only the incorrect measured label on omitted data.

Legacy fidelity_assessment({}, {},0) also returns MFI 1.0, mutation retention basis 'measured', and verdict 'high fidelity - results translate'. The newer multiomic route explicitly disclaims translational/clinical validation, but the measured-label defect remains.

No clinical harm, calibration or real specimen outcome established. No product changes. Product owner should define missing-data semantics and block unsupported assurance; the exact fallback requires a documented contract rather than an improvised audit fix.
