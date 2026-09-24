[1] Review status in ClinVar — ncbi.nlm.nih.gov
https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/

Review status in ClinVar

![U.S. flag](https://www.ncbi.nlm.nih.gov/coreutils/uswds/img/favicons/favicon-57.png)

An official website of the United States government

Here's how you know

# ClinVar

## Relating variation to medicine

## Review status in ClinVar

- Review status
- - Definition
- Review status on submitted records (SCV)
- Review status on aggregate records (VCV and RCV)
- Web display
- FTP
- Assertion criteria
- - Definition
- How to provide new assertion criteria
- How to specify assertion criteria for your submission
- Web display
- FTP

# Review status

## Definition

ClinVar reports the level of review supporting the variant classification for the variation as review status. Gold stars provide a graphical representation of the aggregate review status on web pages. Review status is reported only in text format in ClinVar's products available by FTP.

## Review status on submitted records (SCV)

Review status on submitted records (SCVs) is assigned as follows:

| Number of gold stars | Review status | Description |
| --- | --- | --- |
| four | practice guideline | There is a submitted record with a classification from a practice guideline |
| three | reviewed by expert panel | There is a submitted record with a classification from an expert panel |
| one | criteria provided, single submitter | There is a single submitted record with a classification, where assertion criteria and evidence for the classification (or a public contact) were provided. |
| none | no assertion criteria provided | There are one or more submitted records with a classification but without assertion criteria and evidence for the classification (or a public contact). |
| none | no classification provided | There are one or more submitted records without a classification. |

## Review status on aggregate records (VCV and RCV)

ClinVar calculates a review status for each type of classification on aggregate records, both VCV (variation) and RCV (variation-condition) records.

### Germline classification and oncogenicity

Review status for germline classification and oncogenicity on aggregate records (VCV and RCV) is assigned based on data provided in the contributing SCVs:

1. if there is a submitted record from a practice guideline or an expert panel, the review status from the submitted record (SCV) is used for the VCV and RCV records
2. otherwise, the review status is based on a. transparency: whether any submitted record provided assertion criteria and evidence for the classification (or a public contact) and b. consensus: whether the submitted records agree on the classification

The review statuses for germline classification and oncogenicity on aggregate records are as follows:

| Number of gold stars | Review status | Description |
| --- | --- | --- |
| four | practice guideline | There is a submitted record with a classification from a practice guideline |
| three | reviewed by expert panel | There is a submitted record with a classification from an expert panel |
| two | criteria provided, multiple submitters, no conflicts | There are multiple submitted records with a classification. Assertion criteria and evidence for the classification (or a public contact) were provided and the classifications agree. |
| one | criteria provided, conflicting classifications | There are multiple submitted records with a classification, where assertion criteria and evidence for the classification (or a public contact) were provided. However there are conflicting classifications. The conflicting values for the classification are enumerated. |
| one | criteria provided, single submitter | There is a single submitted record with a classification, where assertion criteria and evidence for the classification (or a public contact) were provided. |
| none | no assertion criteria provided | There are one or more submitted records with a classification but without assertion criteria and evidence for the classification (or a public contact). |
| none | no classification provided | There are one or more submitted records without a classification. |
| none | no classification for the individual variant | The variant was not classified directly in any submitted record; it was submitted to ClinVar only as part of a haplotype or a genotype. |

### Somatic classification of clinical impact

Review status for somatic classification of clinical impact on aggregate records (VCV and RCV) is assigned based on data provided in the contributing SCVs.

It differs from the calculation for germline and oncogenicity in that consensus in the classification is not considered. The rationale is that for a given variant, it will be typical for different tiers to be reported for different types of assertions (therapeutic, diagnostic, or prognostic) and different tumor types, so conflicts are not informative. Thus the review status for somatic clinical impact is based on the following data provided in the contributing SCVs:

1. if there is a submitted record from a practice guideline or an expert panel, the review status from the submitted record (SCV) is used for the VCV and RCV records
2. otherwise, the review status is based on a. transparency: whether any submitted record provided assertion criteria and evidence for the classification (or a public contact) b. multiplicity: whether there are multiple submitted records (consensus is not required)

The review statuses for somatic classification of clinical impact on aggregate records are as follows:

| Number of gold stars | Review status | Description |
| --- | --- | --- |
| four | practice guideline | There is a submitted record with a classification from a [practice guideline](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_guidelines/) |
| three | reviewed by expert panel | There is a submitted record with a classification from an [expert panel](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_guidelines/) |
| two | criteria provided, multiple submitters | There are multiple submitted records with a somatic classification of clinical impact. Assertion criteria and evidence for the classification (or a public contact) were provided. |
| one | criteria provided, single submitter | There is a single submitted record with a classification, where assertion criteria and evidence for the classification (or a public contact) were provided. |
| none | no assertion criteria provided | There are one or more submitted records with a classification but without assertion criteria and evidence for the classification (or a public contact). |
| none | no classification provided | There are one or more submitted records without a classification. |
| none | no classification for the individual variant | The variant was not classified directly in any submitted record; it was submitted to ClinVar only as part of a haplotype or a genotype. |

## Web display

The aggregate review status for each type of classification for a VCV record is displayed at the top the VCV page in the Classification summary, with words and with gold stars.

Note that there are multiple review statuses represented by zero and one stars; see the correlation between review status and stars above.

## FTP

The review status is included in FTP files as shown below.

| File | How review status is reported |
| --- | --- |
| ClinVarVCVRelease XML This is the new format of VCV XML; we encourage you to switch to this file in early 2024. See the README file on the FTP site for more information. | Using germline classification as an example, review status for an aggregate classification in the VCV is reported as: ... criteria provided, single submitter ... Review status for each submitted record (SCV) is reported as: ... ... criteria provided, single submitter ... ... |
| ClinVarRCVRelease XML This is the new format of RCV XML; we encourage you to switch to this file in early 2024. See the README file on the FTP site for more information. | Using germline classification as an example, review status for an aggregate classification in the RCV is reported as: ... criteria provided, single submitter ... Elements for SomaticClinicalImpact and OncogenicityClassification are at the same level as GermlineClassification. Review status for each submitted record (SCV) is reported as: ...; criteria provided, single submitter ... |
| ClinVarFullRelease XML This is the old format of VCV XML which will be phased out in 2024. See the README file on the FTP site for more information. | Review status for the RCV is reported as: ... criteria provided, single submitter ... Review status for the SCV is reported as: ...; criteria provided, single submitter ... |
| ClinVarVariationRelease XML This is the old format of VCV XML which will be phased out in 2024. See the README file on the FTP site for more information. | Review status for the VCV is reported as: ... criteria provided, single submitter ... Review status for the SCV is reported as: ... ... criteria provided, single submitter ... ... |
| VCF | Review statuses for the variant are included in the VCF file as three INFO tags: CLNREVSTAT: review status for the aggregate germline classification ONCREVSTAT: review status for the aggregate oncogenicity classification SCIREVSTAT: review status for the aggregate somatic classification of clinical impact |
| variant_summary.txt | The file reports the aggregate review status per classification type for each variant. ReviewStatus: review status for the aggregate germline classification ReviewStatusOncogenicity: review status for the aggregate oncogenicity classification ReviewStatusClinicalImpact: review status for the aggregate somatic classification of clinical impact |
| summary_of_conflicting_interpretations.txt | The review status for each of a pair of conflicting submissions is reported as Submitter1_ReviewStatus and Submitter2_ReviewStatus. |
| submission_summary.txt | The review status for each submission is reported as ReviewStatus. |

# Assertion criteria

## Definition

Assertion criteria refers to a citation (PubMed ID, PubMedCentral ID, or DOI) or an electronic document (a Word document or PDF) that describes:

- the variant classification terms used by the submitter (e.g. pathogenic, uncertain significance, benign, or appropriate terms for other types of variation) and
- the criteria required to assign a variant to each category

This document describes the categories and criteria that are used by the submitter to guide their classification process, not the specific evidence for an individual variant classification. This document is made available to ClinVar users on Clinvar web displays.

Note that when assertion criteria are provided, it is presumed that the submitter attests that the variant was classified according to a comprehensive review of evidence consistent with, or more thorough than, current practice guidelines. Read more details about meeting the requirements for the "criteria provided, single submitter" review status.

### Representative documentation of assertion criteria

The following table provides examples to help submitters understand what is meant by documentation of assertion criteria and to guide you in creating your own documentation for your specific assertion criteria. These example documents are not intended to be submitted to ClinVar by other groups; please do not include either of these documents in your submission (unless your group provided the sample documentation).

Examples of submitted assertion criteria.

| Type of documentation | Example |
| --- | --- |
| Assertion criteria | GeneDx Variant Classification Process June 2021 |
| Expert panel | 2013-08_InSiGHT_VIC_v1.9.pdf |

## How to provide new assertion criteria

For a citation:

- Wizard or file submission: in the ClinVar Submission Portal, on the submission's "Assertion method" tab, select the option to "use a new assertion criteria".
- API submission: a citation is provided in the JSON, as described in the API documentation.

For a file:

Files for assertion criteria are uploaded independent of submissions. On your organization's landing page in the ClinVar Submission Portal use the button "View/add assertion criteria files" below your organization's name to upload a new file. ClinVar generates a URL to make the file publicly accessible on ClinVar web displays.

Some recommendations to name your assertion criteria file:

- use an informative name, including the name of your organization
- include a version number or a date; assertion criteria are likely to be updated over time and ClinVar requires distinct names for distinct files
- if you have multiple documents for different scenarios, for example one document for autosomal recessive variants and one for autosomal dominant variants, include that information in the name
- if your methods are published but not open access, please excerpt the variant interpretation methods from the publication and submit with reference to the full publication

## How to specify assertion criteria for your submission

Wizard or file submission: in the ClinVar Submission Portal, on the submission's "Assertion method" tab:

- select the option to "use one of my existing assertion criteria"
- choose the desired citation or file

API submission: a citation or URL for a file already uploaded in the Submission Portal is provided in the JSON, as described in the API documentation.

A submission must reference one set of assertion criteria only, and these criteria must apply to every variant classification in the submission.

Remember that the review status of "criteria provided, single submitter" requires this documentation AND either the evidence for each interpretation OR a public contact for your organization.

## Web display

Links to assertion criteria for each submission are displayed on the Variation ID pages in the "Submitted interpretations and evidence" section of the "Variant details" tab as shown in Figure 2.

Figure 2. Links to assertion criteria on a VCV record.

![Display of review status in the evidence section](https://www.ncbi.nlm.nih.gov/core/assets/clinvar/images/revstat_Fig2.png)

Assertion criteria can also be accessed per organization on each submitter page as shown in Figure 3.

Figure 3. An example of a submitter page with links to assertion criteria for that organization.

![Display of an organization page with assertion criteria](https://www.ncbi.nlm.nih.gov/core/assets/clinvar/images/revstat_Fig3.png)

## FTP

Assertion criteria are availabe in the RCV and VCV XML files as shown in the following table.

Assertion criteria in ClinVar's FTP files.

| File | How assertion criteria are reported |
| --- | --- |
| ClinVarVCVRelease XML ClinVarVariationRelease XML | Assertion criteria for each submitted (SCV) is reported as: ... Ambry Autosomal Dominant and X-Linked criteria (10/2015) https://submit.ncbi.nlm.nih.gov/ft/byid/aH4uL7vL/mid-7377_ambry_classification_scheme_oct_2015.pdf ... |
| ClinVarRCVRelease XML ClinVarFullRelease XML | Assertion criteria for each submitted record (SCV) is reported as: ... Ambry Autosomal Dominant and X-Linked criteria (10/2015) https://submit.ncbi.nlm.nih.gov/ft/byid/aH4uL7vL/mid-7377_ambry_classification_scheme_oct_2015.pdf ... |

Follow NCBI

National Library of Medicine 8600 Rockville Pike Bethesda, MD 20894

Last updated: 2024-04-18T13:38:49Z
