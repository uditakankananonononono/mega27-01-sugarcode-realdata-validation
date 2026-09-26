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
