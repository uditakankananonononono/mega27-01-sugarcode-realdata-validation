# Deeper exported-route checks

Product HEAD: 1ccbd0305cd9b639b06324fa136654d686483b7a.
Source: https://github.com/uditakankananonononono/sugarcode-ai
These findings are independent of the 2,319 passing repository tests. No product edits were made.

## NeuroPlan imaging-claim provenance
The exported plan_surgery(tumor, image_size) route takes geometry, not MRI or diffusion inputs. It returns segmentation.confidence=0.9, named tissue labels, and tractography.method="deterministic streamline, FA threshold 0.2". Inspection of core.py lines 16-42 and 74-80 shows spherical volume calculation, constant confidence and static tract names, not those imaging operations. Repeat inputs produce identical output. analyze_neurosurgical_case also inherits these fields and marks segmentation/tractography capabilities complete; its mechanism/non-autonomous disclaimer does not identify these fields as geometry-only estimates.

Separate segment_mri and analyze_tractography functions genuinely process supplied intensity volumes and streamline geometry. This finding does not say those functions are absent. It says the exported planning routes' output labels overstate the computations actually performed on their supplied inputs.

Proposed review direction: distinguish geometry estimates from image-derived segmentation, remove unsupported confidence/method labels or mark them unavailable, and test the provenance of every output field. Product owner decides the change.

## RareNet accepted but inert inheritance input
explainable_rank validates inheritance in {dominant, recessive, x_linked, None} and echoes it in inheritance_assumption. Source inspection shows no use of it in evidence contributions or posterior/ranking calculations. Identical full candidate lists are returned for all three valid values on the same phenotype input; repeat runs confirm this.

This is an API clarity/behavior gap, not a demonstrated clinical misdiagnosis. The original aspirational spec does not explicitly require inheritance-conditioned scoring. A caller should nevertheless not have to discover silently that an accepted assumption does nothing. Proposed review direction: document it as metadata only, remove the inert parameter, or implement evidence-based handling with an explicit verification plan. Do not improvise inheritance rules from memory.

## Limits
Neither check establishes clinical validity, clinical harm, or benchmark superiority. Neither modifies the product. Existing suite success remains true; it does not cover these provenance/API questions. Per-sentence verification of all 77 aspirational spec files remains open.
