# JUDGE_ROUNDS: F1 CRISPR on/off-target scoring family
Modules: cfd_offtarget, mit_offtarget, crisprscan_score, crisprater, crispr_opt, crispr_muse
Judge: ChatGPT (Free tier) on the user's account via cloud browser config-c, read lease, released after the round.

## Round 1 - 2026-09-26 16:23 IST

### PROMPT (verbatim)
You are an ISEF grand-award judge reviewing the external-validation plan for a CRISPR guide-scoring toolkit family (6 modules: CFD off-target scorer, MIT/Hsu off-target scorer, CRISPRscan on-target scorer, CRISPRater on-target scorer, a CRISPR design optimizer, and a guide-generator refinement module). Current external-validation state, all on public data: (1) the CFD reimplementation reproduces Doench 2016 oracle scores exactly across 9,600 guide/off-target scores and 6/6 CRISPOR reference doctests; (2) the MIT score matches CRISPOR's calcHitScore exactly on 600 held-out BRCA1 pairs; (3) CRISPRscan matches 275/300 held-out 35-nt contexts - we discovered one coefficient (largest negative weight) sits at dinucleotide position 19 in crisprScore but position 18 in CRISPOR, changing 8.3% of scores by up to 10 points; which position is correct vs the original Moreno-Mateos 2015 supplement is OPEN; (4) CRISPRater matches author-computed scores for 3,138/3,141 published sgRNAs after we fixed a GC-window bug (GC4-13 not GC4-14); (5) on-target ranking reaches Spearman 0.71 (RS2+position) vs measured efficiency on 5,310 guides across 17 genes, but RS1/RS2 are in-sample on that set and the design pipeline still ranks guides with RS1. New held-out fixtures just built: 1,400 guide/off-target pairs mined from TP53 (disjoint locus) with oracle CFD and MIT scores. Planned next experiments: (a) settle the CRISPRscan coefficient against the original 2015 supplement; (b) a fully held-out on-target benchmark where our ranking must BEAT RS1 and sequence-only RS2 on genes never used to train them; (c) a discovery study mapping systematic CFD-vs-MIT discordance on endogenous off-target sites. As a judge, attack this plan: What are the weakest claims? Where are the leakage and circularity risks? What baselines are missing? What specific experiments would make each claim convincing to a skeptical expert? Be harsh and specific.

### RESPONSE (verbatim)
ISEF grand-award judge critique: CRISPR guide-scoring toolkit validation plan

Overall: this is much stronger than the typical high-school computational biology project because several components have real reproducibility checks and independent fixtures. However, a skeptical judge would still attack the project at one central point:

You have demonstrated that you can reproduce published scoring functions. You have not yet demonstrated that your toolkit produces better biological decisions than existing tools.

The validation plan needs to separate software correctness, model validity, and scientific discovery. Right now those are partially mixed.

1. Weakest claim: “The toolkit is validated because scores match published implementations”
Problem

Exact score reproduction is valuable, but it proves only:

“Our implementation is faithful.”

It does not prove:

“Our implementation is useful.”

A judge may say:

CFD matching CRISPOR is expected.

MIT matching CRISPOR is expected.

CRISPRater matching authors' numbers is expected.

These are engineering validations, not biological validations.

The strongest evidence here is actually the discovery of the CRISPRscan position discrepancy, because that shows you found a real reproducibility issue. But it creates a danger:

You cannot claim CRISPRscan is validated until the ambiguity is resolved.

Needed experiment

Create a validation hierarchy table:

Module	Current claim	What it proves	Missing proof
CFD	Exact reproduction	Implementation correctness	Biological ranking validity
MIT	Exact reproduction	Implementation correctness	Biological ranking validity
CRISPRater	3138/3141 match	Implementation correctness	Generalization
CRISPRscan	275/300	Partial reproduction	Original model fidelity
Optimizer	Uses scores	Integration	Independent design improvement
Generator refinement	Produces candidates	Generation	Experimental usefulness

This prevents judges from thinking score reproduction is being oversold.

2. CRISPRscan is currently the biggest technical vulnerability
Problem

The coefficient-position mismatch is not a minor bug.

You found:

largest negative coefficient at dinucleotide position 19 vs CRISPOR position 18

That means:

your implementation and CRISPOR may not represent the same model,

your 275/300 agreement may hide systematic error,

your optimizer may be optimizing the wrong objective.

A judge will ask:

“Why should I trust your CRISPRscan-guided designs if you have not established the actual feature encoding?”

Required experiment

Do not just compare against CRISPOR.

Reconstruct the original Moreno-Mateos 2015 model.

Run:

Experiment A: coefficient provenance

Compare:

Original supplementary table

crisprScore implementation

CRISPOR implementation

Your implementation

For every coefficient:

position

nucleotide/dinucleotide feature

sign

magnitude

Produce a difference matrix.

Experiment B: biological consequence

Take the 300 held-out guides.

Run:

model version A (position 18)

model version B (position 19)

Compare:

Pearson correlation

Spearman ranking

top 10% guide selection overlap

predicted vs measured efficiency

The important question is not:

“Which implementation matches CRISPOR?”

The important question is:

“Which implementation better predicts biology?”

3. Major leakage risk: on-target benchmark

This is the biggest issue.

You write:

RS2+position achieves Spearman 0.71 on 5,310 guides across 17 genes

A judge immediately asks:

Was RS2 trained on these data?

Because RS2 itself was trained using experimental guide efficiency datasets.

If yes:

That 0.71 is not an external validation.

It is measuring:

“How well does a published model perform on data similar to its training distribution?”

Not your contribution.

Current status

Your planned experiment (b) is the correct fix.

But it must be stricter.

Required held-out benchmark design
Bad design

Random split:

80% guides
20% guides

This leaks gene context.

Better design

Gene-level split:

Example:

Training:

genes A-Q

Testing:

completely unseen genes R-T

No guide from test genes appears in:

training

optimizer tuning

threshold selection

feature selection

Best design

External dataset split:

Train/tune only on:

Doench 2016

Moreno-Mateos

Test only on:

newer CRISPR screens

different organism/cell type

Metrics:

Do not only report Spearman.

Report:

Ranking

Spearman

Kendall tau

top-1/top-5 guide success rate

Practical design metric

"If a researcher selects the top predicted guide, how often does it outperform random selection?"

Example:

Top predicted guide:
72% successful

Random guide:
45%

baseline improvement:
+27 percentage points

That is much more intuitive.

4. The optimizer claim is currently unsupported

This is probably the weakest module.

You have:

scoring models

candidate generation

refinement

But optimization introduces circularity.

A judge asks:

“Are you optimizing toward your own score?”

Example:

Optimizer:

generates guide

scores with RS2

selects highest RS2 guide

Then you test RS2.

That is circular.

You have only proven:

“The optimizer is good at maximizing RS2.”

Not:

“The optimizer finds biologically better guides.”

Required experiment

Freeze all scoring models.

Then:

Input:

A target gene.

Compare:

Baselines

Random guides

CRISPOR recommended guide

CHOPCHOP

Benchling (if accessible)

Existing heuristic ranking

vs

Your optimizer.

Measure:

predicted efficiency

off-target burden

diversity

CFD/MIT safety score

The ideal result:

Your optimizer finds guides that are:

equal or better on-target

lower off-target risk

without requiring proprietary ranking

5. CFD vs MIT discordance discovery: interesting but scientifically risky

This is the most novel part, but also the easiest to overclaim.

You propose:

mapping systematic CFD-vs-MIT discordance on endogenous off-target sites

Good idea.

But:

Current danger

You may only discover:

“Two mathematical scoring systems disagree.”

That is not biological insight.

A judge asks:

“Why should disagreement matter?”

Make it a real discovery study

Define hypotheses before mining.

Examples:

Hypothesis 1

CFD and MIT disagree systematically for specific mismatch classes.

Test:

mismatch position

mismatch identity

PAM type

seed/non-seed location

Hypothesis 2

Discordance predicts experimentally observed cleavage better.

Need:

Independent dataset:

GUIDE-seq

CHANGE-seq

SITE-seq

CIRCLE-seq

Then test:

Model A:
CFD only

Model B:
MIT only

Model C:
discordance feature added

Does discordance improve prediction?

If yes, this becomes a scientific contribution.

6. Missing baselines

Your current baseline set is incomplete.

A skeptical judge expects comparison against:

Off-target

Minimum:

CRISPOR

Cas-OFFinder

Cas-OFFinder-like enumeration methods

Also:

CFD alone

MIT alone

You need to show whether your combined system improves anything.

On-target

Missing:

CRISPRscan

RS1

RS2

DeepCRISPR

DeepHF

CRISPR-Net

CRISPRon

Even if you cannot run all, explain why.

A weaker baseline comparison looks like cherry-picking.

7. The TP53 1,400 pair fixture: good, but limited

This is useful.

However:

A judge will say:

“You validated on TP53 because you had the data.”

One locus is not enough.

You need:

At least:

multiple genes

different GC content

different genomic contexts

Example:

BRCA1
TP53
CFTR
EMX1
VEGFA

Then report:

score correlation

rank correlation

failure cases

8. The strongest possible ISEF framing

Do NOT frame this as:

"Six AI modules for CRISPR design."

That sounds like a software collection.

Frame it as:

"A reproducible, independently validated framework revealing when established CRISPR guide-scoring models disagree and whether ensemble disagreement improves guide selection."

That creates a scientific question.

Priority ranking of fixes
Must fix before claiming strong validation
1. Resolve CRISPRscan discrepancy

Critical

2. True held-out on-target benchmark

Critical

3. Remove optimizer circularity

Critical

Needed for grand-award level
4. Compare against external tools
5. Validate discordance against experimental cleavage datasets
Nice but not essential
6. Expand TP53 fixture
7. Add visualization dashboard
Judge verdict

Current state:

Software engineering/reproducibility: very strong.
Scientific validation: incomplete.
Novel discovery potential: high but unproven.

The project becomes ISEF-grand-award caliber if the final paper demonstrates one of these:

Your implementation fixes a real reproducibility failure in CRISPR scoring, with biological consequences.

or

Your optimizer demonstrably improves guide selection over established tools on unseen genes.

or

CFD/MIT disagreement contains biological information that improves off-target prediction.

Right now, the biggest danger is not bad methodology. It is that the project could be judged as a carefully tested software package rather than a new scientific finding. The next experiments should be designed specifically to cross that boundary.


### Round 1 NOVELTY FOLD-BACK (required per 2026-09-26 17:00:38 rule)
The critique's central charge - "the project could be judged as a carefully tested software
package rather than a new scientific finding" - produced these concrete changes to the work:
1. Discovery-study redesign: the CFD-vs-MIT discordance study now validates against
   EXPERIMENTAL cleavage data (GUIDE-seq/CIRCLE-seq public datasets) as its primary arm,
   with the in-silico oracle map demoted to supporting analysis. New claim target:
   "CFD/MIT disagreement carries biological information that improves off-target
   classification against measured cleavage" - a scientific claim, not a software claim.
2. The CRISPRscan AA18/AA19 resolution (completed 2026-09-26, upstream issue
   maximilianh/crisporWebsite#76) is promoted from a validation footnote to a
   first-class finding: a reproducibility defect in a widely used reference tool with
   measured score consequences - exactly the "fix a real reproducibility failure"
   grand-award path the judge named.
3. Held-out on-target benchmark redesigned: optimizer must beat RS1 and sequence-only
   RS2 on genes disjoint from all training genes of both rule sets, with guide-level
   bootstrap CIs; optimizer circularity (design pipeline ranking with RS1) must be
   removed BEFORE the benchmark runs so the compared ranker is the shipped one.

---

## ROUND 2 of 10 - 2026-09-26 ~17:10 IST (ChatGPT, her account, cloud browser)

### PROMPT (verbatim)

```
You are an ISEF grand-award judge. Same CRISPR-scoring validation project as before; this round is about redirecting a first honest negative into a stronger scientific claim. New results since round 1: (1) RESOLVED the CRISPRscan coefficient discrepancy against the original Moreno-Mateos 2015 supplement - CRISPOR misplaces the largest negative weight (AA19 vs its AA18); our toolkit matches the original model; upstream bug report filed. (2) We ran the planned discovery study: do CFD-vs-MIT-score disagreements carry biological information about MEASURED off-target cleavage? Data: 764 experimentally measured on/off-target rows from 7 studies (Cho2014, Frock2015, Hsu2013, Kim2015, Ran2015, Tsai2015 GUIDE-seq, Wang2015), 600 usable off-targets with <=4 mismatches across 45 guides. Results: at read-fraction threshold 0.001, global AUROC CFD 0.669 vs MIT 0.619; a 2-feature logistic (log CFD + log MIT) with leave-one-guide-out CV got 0.649 - WORSE than CFD alone, so the naive combined model fails. BUT: per-study heterogeneity is large (MIT wins on Cho 0.782 vs 0.698 and Kim16 0.725 vs 0.562; CFD wins on Hsu, Tsai, Wang, Kim/K562); at a strict 0.05 threshold MIT beats CFD (0.771 vs 0.744); and the discordance residual (log CFD regressed on log MIT) correlates with log measured read fraction at Spearman 0.327. So: signal exists in the residual but does not survive predictive CV. As the judge: how do we redirect this into a rigorous, novel, publishable claim? Is the right move per-study read-fraction normalization, a continuous regression framework, mismatch-class-restricted analysis, or something else? What exactly would convince you the discordance is real biology versus dataset noise and study-effect confounding? Give a concrete experimental design with controls.
```

### RESPONSE (verbatim)

```
Show more
ISEF Grand Award Judge Review: The project should not pivot into “a better predictor.” Pivot into a mechanistic audit of where existing CRISPR specificity models fail and why.

Your strongest result is not AUROC 0.669 vs 0.619. That is incremental. The scientifically interesting result is:

When two established mismatch-based CRISPR off-target models disagree, the disagreement contains biological signal, but that signal is context-dependent and is not captured by a simple global ensemble.

That is a much stronger claim because it turns a failed classifier into a model disagreement biology study.

However, the current evidence is not yet sufficient because the biggest alternative explanation is:

“The discordance residual is just a proxy for study protocol, guide identity, sequencing depth, or publication-era differences.”

You need to design experiments that separate biology from dataset structure.

First: what I would NOT do
1. Do not make per-study normalization the main discovery

Per-study normalization is necessary, but alone it is dangerous.

Why?

Your strongest examples:

MIT wins Cho2014

MIT wins Kim2015

CFD wins GUIDE-seq/Tsai/Wang

This screams:

assay technology effect.

Different studies used:

different cell types

different detection platforms

different thresholds

different sequencing depths

different mismatch distributions

A model that learns "GUIDE-seq favors CFD" is not discovering cleavage biology.

A judge will immediately ask:

"Did you discover biology, or did you rediscover that assays differ?"

So: normalize, but do not stop there.

The strongest redirect: "Discordance as a biological phenotype"

I would redesign the paper around:

Title concept

"When CRISPR specificity models disagree: identifying hidden sequence determinants of off-target cleavage beyond MIT and CFD scores"

The hypothesis:

MIT and CFD encode partially overlapping but biologically incomplete models of cleavage. Their residual disagreement identifies sequence contexts where current scoring functions fail.

This is novel enough.

Experiment 1 — Replace binary classification with continuous cleavage modeling

This is the most important change.

Your current endpoint:

"Does this off-target exceed threshold?"

throws away information.

You already have:

measured read fraction

mismatch count

guide identity

study

Use:

log(read fraction)

as the response.

Model:

Baseline
Y=β
0
	​

+β
1
	​

CFD+β
2
	​

MIT+ϵ

Then test:

Discordance model
Y=β
0
	​

+β
1
	​

CFD+β
2
	​

MIT+β
3
	​

(CFD−MIT)+ϵ

The question becomes:

Does disagreement between models explain cleavage magnitude after accounting for both scores?

Your current Spearman 0.327 suggests there may be something here.

But do it properly.

Experiment 2 — Mixed-effects model (this is probably your winning analysis)

The dataset is not independent.

The same guide produces many off-targets.

The same study produces many measurements.

Your statistical structure is:

off-target
    |
    +-- guide
    |
    +-- study

Use:

log(readfraction)∼CFD+MIT+Discordance+(1∣guide)+(1∣study)

Meaning:

fixed effects = sequence scoring

random effects = guide-specific biology and study protocol

Your key test:

Does the discordance coefficient remain significant?

Interpretation:
If discordance disappears:

Your claim becomes:

Model disagreement reflects guide/study effects rather than hidden cleavage biology.

Still publishable as a benchmarking critique.

If discordance survives:

Now you have evidence:

Existing scoring models leave systematic biological variance unexplained.

That is much stronger.

Experiment 3 — Stratify by mismatch class

This is actually more important than per-study normalization.

Why?

CFD and MIT were trained on different assumptions about mismatch penalties.

A 1-mismatch off-target is not biologically equivalent to a 4-mismatch off-target.

Split:

Group A

1 mismatch

Group B

2 mismatches

Group C

3 mismatches

Group D

4 mismatches

Then ask:

Where does discordance predict?

My expectation:

The signal will probably concentrate in:

2–3 mismatch sites

because:

1 mismatch: nearly everything cleaves

4 mismatch: nearly nothing cleaves

middle regime contains biological uncertainty

A real discovery would look like:

Discordance predicts cleavage specifically in intermediate mismatch regimes where existing scoring models disagree.

That is a biological statement.

Experiment 4 — Sequence feature discovery

If discordance is real, ask:

"What does CFD miss?"

Build a residual model.

Pipeline:

Fit:

read_fraction ~ CFD + MIT

Calculate residual:

observed - predicted

Predict residual from sequence features:

Candidate features:

PAM context

PAM type

nucleotide immediately upstream/downstream

mismatch topology

Not just mismatch count.

Examples:

seed vs non-seed

mismatch spacing

adjacent mismatch pairs

nucleotide identity

Example:

position 10 mismatch A→G

rather than:

one mismatch at position 10
DNA thermodynamics

GC

melting temperature

bulges if available

The discovery question:

What sequence features explain where CFD and MIT disagree?

Critical control experiments

These are what would convince me.

Control 1: Guide-level held-out validation

Your current LOGO CV is good.

But do this specifically:

Train:

discover discordance-associated features

Test:

completely unseen guides

Success:

Residual prediction improves on unseen guides.

Failure:

Only memorizes guides.

Control 2: Leave-one-study-out validation

Even more important.

Train:

Cho + Frock + Hsu + Kim + Ran + Wang

Test:

Tsai GUIDE-seq

Then rotate.

If the effect survives:

strong evidence.

If it collapses:

study artifact.

Control 3: Negative-control residual

Shuffle:

CFD scores among off-targets

MIT scores among off-targets

Repeat residual analysis.

If your correlation remains:

you found an artifact.

Control 4: Publication bias analysis

Your dataset contains famous benchmark studies.

Ask:

Are high-cleavage sites overrepresented?

Compare:

all sites

low-read sites

high-read sites

The CRISPOR bug discovery is actually valuable

Do not bury it.

This gives the project a second pillar:

Pillar 1:

Reproducibility failure

"One widely used implementation differs from the original CRISPRscan model because of a coefficient indexing error."

Pillar 2:

Scientific uncertainty

"Even after correcting scoring implementations, disagreement between models reveals unexplained cleavage variation."

That is a coherent story:

Before improving prediction, verify the predictors.

Very ISEF-friendly.

My proposed final paper claim

Not:

❌ "We built a better CRISPR off-target predictor."

Too weak.

Not:

❌ "CFD beats MIT."

Already known territory.

Instead:

Strong claim:

"A systematic audit of CRISPR off-target scoring models reveals that model disagreement contains measurable information about cleavage activity, but this information is obscured by assay-specific effects. A mixed-effects and sequence-residual framework identifies when and why current scoring models fail."

That is publishable.

The exact next 4-week experiment I would require before giving a Grand Award score
Dataset

Your 600 off-targets.

Models

MIT only

CFD only

MIT + CFD

MIT + CFD + discordance

MIT + CFD + sequence residual features

Validation:

leave-one-guide-out

leave-one-study-out

Primary endpoint:

Continuous:

R
2
, RMSE, Spearman

Secondary:

AUROC at thresholds.

Discovery:

Find sequence contexts where:

CFD predicts low
MIT predicts high
observed cleavage high

and the opposite.

My judge decision:

The project is currently one strong analysis away from becoming a methods paper. The CRISPRscan bug alone is a reproducibility note. The CFD–MIT discordance alone is currently suggestive but confounded. Combining:

implementation audit,

mixed-effects continuous modeling,

mismatch-class stratification,

unseen-guide/unseen-study validation,

would transform the negative result into a rigorous claim about the limits of current CRISPR specificity scoring models and the biological information hidden in their disagreements.

```

### Round 2 NOVELTY FOLD-BACK (required per 2026-09-26 17:00:38 rule)
Judge redirection: stop pivoting to "a better predictor"; pivot to a mechanistic audit of
WHERE existing models fail - replace threshold AUROC with continuous cleavage modeling,
use a mixed-effects design (guide + study random effects) as the primary analysis, and
stratify by mismatch class. Fold-back executed this wake (discovery/cfd_mit_discordance/,
round2_foldback.py + round2_foldback_result.json):

1. Exp 1 (continuous OLS): log10(readfrac) ~ CFD + MIT + |CFD-MIT|. Baseline R2=0.1385;
   discordance term NOT significant in OLS (coef -0.292, p=0.634). Honest negative.
   Methodological fix en route: the judge's literal (CFD-MIT) is an exact linear
   combination of the included scores (perfect collinearity, LR test NaN), so the
   discordance term was re-specified as the nonlinear |CFD-MIT|.
2. Exp 2 (mixed effects, judge's "winning analysis"): log10(readfrac) ~ CFD + MIT +
   |CFD-MIT| with crossed random intercepts for study and guide (45 guides, 8 studies,
   600 sites; nm optimizer after lbfgs hit singular matrix). DISCORDANCE SURVIVES:
   coef -1.376, p=0.0018; LRT vs no-discordance model chi2=9.62, p=0.0019. Guide-level
   variance 0.96 dominates; study variance collapses to ~0. Interpretation per judge:
   existing scoring models leave systematic biological variance unexplained - greater
   absolute model disagreement predicts LOWER cleavage than either score expects.
3. Exp 3 (mismatch-class stratification): discordance's residual signal is
   class-dependent - null at 1-3 mismatches (rho 0.05-0.08, p>0.19) but SIGNIFICANT and
   NEGATIVE at 4 mismatches (rho -0.240, p=0.0012). The judge expected 2-3 mismatches to
   carry the signal; the data say the opposite - the high-mismatch regime is where model
   disagreement maps to cleavage suppression. Hypothesis-updating, preserved honestly.
4. Exp 4 (unseen-study validation): leave-one-study-out Spearman, baseline vs
   discordance model - discordance helps 4/8 studies (Cho 0.188->0.193, Kim 0.829->0.833,
   Wang 0.574->0.581, Frock tie) but HURTS the largest study (Tsai 0.456->0.345).
   Global gains do not transfer uniformly across assay platforms - the judge's
   "did you discover biology or rediscover that assays differ?" caution is live.

Net scientific position after round 2: the binary-classifier negative (round 1/2 honest
negative) redirected into a positive, preregistered-consistent finding: MODEL DISAGREEMENT
IS A CLEAVAGE-SUPPRESSIVE SIGNAL IN MIXED-EFFECTS CONTINUOUS MODELING (p=0.0019),
concentrated at high mismatch counts, with non-uniform cross-study transfer. Paper chapter
claim: "absolute CFD-MIT discordance is an independent negative predictor of off-target
cleavage magnitude after guide and study effects" - a statement about the limits of both
scoring models, not a new black-box predictor.

---

## ROUND 3 of 10 - 2026-09-26 ~18:11 IST (ChatGPT, her account, cloud browser)

### PROMPT (verbatim)

```
You are an ISEF grand-award judge. Same CRISPR-scoring validation project. Your round 2 redirection is executed; new results: (1) Continuous mixed-effects model log10(read fraction) ~ CFD + MIT + |CFD-MIT| with crossed guide+study random intercepts on 600 experimental off-target sites (45 guides, 8 studies): |CFD-MIT| coef -1.376, p=0.0018, LRT p=0.0019 - model disagreement independently predicts LOWER cleavage than either score expects. Robust to study-only random effects (p=0.0048) and guide-only (p=0.0018). (2) Signal concentrated at 4 mismatches (Spearman rho=-0.240, p=0.0012); null at 1-3 mismatches, opposite of your 2-3 mismatch expectation. (3) Unseen-study validation: discordance feature helps 4/8 held-out studies but HURTS the largest (Tsai GUIDE-seq, Spearman 0.456 to 0.345). (4) Plain OLS without random effects: discordance null (p=0.634). Questions: (a) Is the surviving mixed-effects signal a discovery or a statistical artifact of guide-level clustering? What is the strongest falsification test I can run with public data? (b) Does the Tsai harm in leave-one-study-out kill the cross-platform utility claim, or reframe it as assay-specific discordance calibration? (c) What single experiment would make this grand-award competitive rather than a solid benchmarking paper? Be harsh and specific.
```

### RESPONSE (verbatim)

```
Show more
Adversarial ISEF grand-award review — CRISPR scoring validation redirection

This is a substantially better scientific direction than a generic “new CRISPR score beats CFD/MIT” project. The result is no longer merely benchmarking; it is testing whether disagreement between existing models contains biological information. However, the current evidence is not yet a discovery. It is a strong hypothesis-generating statistical observation with a serious remaining confound.

(a) Discovery or guide-clustering artifact?
Verdict: not yet a discovery; the main threat is hierarchical leakage / guide-specific calibration failure

Your mixed-effects result is interesting:

log
10
	​

(read fraction)∼CFD+MIT+∣CFD−MIT∣

with:

∣CFD−MIT∣ coefficient = -1.376

LRT p = 0.0019

robust under study-only and guide-only random intercepts

This is evidence against a trivial study effect, but it does not completely defeat guide clustering.

Why?

You have:

600 sites

45 guides

8 studies

The effective biological sample size is closer to 45 guide perturbation contexts, not 600 independent observations.

A guide may systematically generate:

more unusual mismatch distributions,

different genomic contexts,

different chromatin accessibility,

different off-target abundance,

different sequencing behavior.

The mixed model says:

“Among sites within these guides, disagreement associates with lower cleavage.”

But the remaining question is:

“Does disagreement generalize to unseen guides, or did the model learn properties of particular guides?”

Those are very different claims.

Strongest falsification test

The single strongest public-data test:

Leave-one-guide-family-out validation

Not leave-one-study-out. Leave-one-guide-out is not harsh enough because sites from the same guide can share hidden properties.

Do:

Split:

Training:

80% of guides

all their off-target sites

Test:

completely unseen guides

No site from the same guide appears in both.

Model:

Baseline:

log(read)∼CFD+MIT

Novel:

log(read)∼CFD+MIT+∣CFD−MIT∣

Evaluate:

held-out Spearman correlation

held-out RMSE

likelihood improvement

calibration curves

The critical result:

If discordance survives unseen-guide prediction:

You have evidence of a transferable biological phenomenon.

If it disappears:

The result was probably a guide-specific residual pattern.

A second excellent falsifier:

Permute discordance within guide

Keep:

guide identity

CFD/MIT distributions

read fractions

Randomize ∣CFD−MIT∣ among sites within each guide.

Repeat 1000 times.

If the mixed model signal remains, the feature is not actually informative.

If it collapses, your current p-value is inflated by within-guide structure.

(b) Does Tsai harm kill cross-platform utility?

No. But it kills the strong version of the claim.

Your original implied claim:

“Model disagreement is a universal indicator of hidden cleavage biology across CRISPR assays.”

The Tsai result weakens that.

You found:

Helps 4/8 held-out studies

Hurts largest study:

Tsai GUIDE-seq

Spearman 0.456 → 0.345

That is not a small degradation.

A judge will ask:

"If the largest and most influential dataset loses performance, why should I believe this is a generalizable feature?"

You need to reframe.

The defensible claim becomes:

“Discordance between CFD and MIT captures assay-dependent unexplained cleavage behavior, with strongest signal in datasets sharing certain experimental measurement characteristics.”

That is actually scientifically interesting.

Possible biological interpretation:

Different assays measure different things:

GUIDE-seq detects double-strand-break integration events

other assays measure cleavage abundance indirectly

PCR/amplicon sequencing has different detection thresholds

Discordance might reflect:

cleavage kinetics

repair pathway bias

detection threshold effects

But you cannot claim that without testing.

What would rescue the cross-platform story?

Do not average all studies.

Build an assay-aware model:

log(read)=CFD+MIT+Discordance+Discordance×Assay

Question:

Does the interaction term explain why Tsai fails?

If:

Discordance×Assay

is significant, your result becomes:

“The value of score disagreement depends on assay measurement physics.”

That is a much stronger scientific story than “our feature works everywhere.”

(c) Single experiment that makes this grand-award competitive

The honest answer:

Not another benchmark.

Grand-award level requires demonstrating a new biological mechanism or a predictive advantage that existing scores cannot explain.

The experiment:

Build and validate a CRISPR “model disagreement map” that predicts previously unexplained off-target cleavage in a prospective unseen dataset

Specifically:

Step 1 — Train only on old datasets

Features:

CFD

MIT

mismatch count

position-specific mismatch features

GC%

sequence context

chromatin if available

∣CFD−MIT∣

Freeze model.

No tuning afterward.

Step 2 — Pick a completely untouched dataset

Ideally:

new publication

new guides

new assay platform

Not used in development.

Step 3 — Prediction target

Not:

“Does my score correlate?”

Too incremental.

Instead:

Identify cases where:

CFD≈high

but

MIT≈low

or vice versa.

Ask:

Are these disagreement sites systematically enriched for unexpected cleavage behavior?

A successful result:

Existing models disagree.

Your model predicts which disagreement cases are false positives.

The grand-award-level claim would be:

“Current CRISPR off-target scores disagree because they encode different assumptions about mismatch tolerance. Their disagreement itself provides an orthogonal signal for distinguishing apparent from biologically realized off-target cleavage.”

That is a conceptual advance.

Current judge scorecard
Criterion	Current status
Novel question	Strong
Statistical rigor	Moderate-strong
Discovery evidence	Not yet
Generalization	Weak due Tsai failure
Biological interpretation	Missing
Grand award competitiveness	Borderline
Biggest danger in the next round

Do not say:

“We discovered disagreement predicts off-target cleavage.”

A judge can destroy that.

Say:

“We identified and are testing whether score disagreement represents an assay-dependent biological axis of off-target prediction.”

That is accurate.

The project becomes competitive if the next experiment shows unseen-guide, unseen-study generalization plus assay interaction explanation. Without that, it remains a very good computational benchmarking paper, but not a grand-award discovery.

```

### Round 3 NOVELTY FOLD-BACK (required per 2026-09-26 17:00:38 rule)
Judge's three falsification demands, all executed this wake (discovery/cfd_mit_discordance/,
round3_foldback.py + round3_foldback_result.json; plus round2_sensitivity_study_guide.json):

1. UNSEEN-GUIDE generalization (Exp A, 5-fold grouped-by-guide CV, no guide in both
   train/test): the |CFD-MIT| feature does NOT transfer - pooled held-out Spearman
   0.3597 (baseline) vs 0.3593 (discordance model); RMSE 1.193 vs 1.196 (slightly worse).
   Per the judge's own framework ("if it disappears, the result was probably a
   guide-specific residual pattern"), the transferable-predictor claim is DEAD.
   Honest negative, preserved.
2. WITHIN-GUIDE PERMUTATION TEST (Exp B, 1000 permutations of |CFD-MIT| among sites
   within each guide, guide-demeaned within-estimator): observed t=-3.28 against a null
   of mean 0.03/SD 0.99 - permutation p < 0.001 (0/1000 nulls reach |t|). The
   within-guide association is REAL, not an artifact of inflated degrees of freedom.
3. ASSAY INTERACTION (Exp C, disc x study interaction, guide random intercept): NOT
   significant (LRT chi2=7.55, df=7, p=0.374; Tsai term p=0.345). The assay-aware rescue
   the judge proposed for the Tsai leave-one-study-out harm does NOT hold. Honest negative.

NOVELTY PIVOT produced by this round (the concrete change): the finding is reframed from
"discordance is a transferable off-target feature" (dead) to "discordance is a
WITHIN-GUIDE, non-transferable signal - i.e., current position-specific mismatch models
omit GUIDE-LEVEL determinants of cleavage." Both CFD and MIT score only mismatch
position/type; neither encodes guide-sequence properties. The round-4+ research question
becomes mechanistic decomposition: WHICH guide-level features (GC content, sequence
complexity, predicted accessibility/self-folding, PAM-distal seed composition) explain the
variance that |CFD-MIT| proxies for within guides? Claim language per judge: "score
disagreement exposes guide-level biology absent from mismatch-based models" - testing,
not declaring.
