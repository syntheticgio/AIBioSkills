---
name: statistical-methods-agent
description: Fetch a biomedical research paper and produce a structured analysis of every statistical method used — including assumptions, plain-language interpretation, and critique
---

Analyze the statistical methods used in this paper: **$ARGUMENTS**

Use the MCP tools available to you to fetch the paper, gather citation context, and then analyze every statistical method in depth. Follow the steps below in order. If a step fails or returns no data, note the gap and continue.

## Data Gathering Steps

### 1. Fetch Paper Content

Call `paper_fetch` with:
- `identifier`: the user's input (`$ARGUMENTS`)
- `sections`: `["methods", "results", "discussion", "all"]`

Extract the metadata (title, authors, journal, date, PMID, DOI, PMC ID) and full-text sections. If only the abstract is available, note this prominently — the analysis will be limited.

### 2. Field Context and Citation Impact

Call `semantic_scholar_search` with the DOI (prefixed with `DOI:`) if available, or the paper title.

Extract: citation count, influential citation count, fields of study, TLDR summary. This provides field context for the significance assessment.

### 3. Related Work for Field Context

Call `pubmed_search` with key terms from the paper title plus `"review[pt]"`, limited to 5 results. Focus on reviews published within 5 years of the paper.

This helps assess the paper's significance within its field.

### 4. Analyze Statistical Methods

Read through the Methods and Results sections carefully. For **every** statistical test, model, or analytical approach mentioned, identify:

- **Test/Method name** — the exact name (e.g., "two-tailed Student's t-test", "Cox proportional hazards model", "Kaplan-Meier survival analysis", "Fisher's exact test", "linear mixed-effects model", "Benjamini-Hochberg correction")
- **What data it was applied to** — the specific comparison or dataset
- **Explicit assumptions** — any assumptions stated by the authors
- **Implicit assumptions** — assumptions required by the test but not stated (e.g., normality for t-test, proportional hazards for Cox, independence of observations)
- **Result** — the test statistic, p-value, confidence interval, effect size, or other reported outcome
- **Authors' interpretation** — what the authors concluded from this test
- **Plain-language significance** — translate the statistical result into plain language (e.g., "There is less than a 5% probability that this difference occurred by random chance alone")

### 5. Statistical Critique

Apply statistical knowledge to assess:

- **Appropriateness of each test** — Was this the right test for this data type and study design?
- **Assumption validity** — Are the assumptions of each test likely met given the study design?
- **Multiple comparisons** — Were corrections applied when needed?
- **Effect sizes** — Were effect sizes reported alongside p-values?
- **Sample sizes** — Were power analyses or sample size justifications provided?
- **Interpretation fairness** — Are the conclusions well-supported, or do they overreach?
- **Alternative interpretations** — Could the data support different conclusions?
- **Missing analyses** — What analyses should they have done but didn't?

### 6. Save Report

After generating the complete report below, save it as a markdown file in the current working directory using the Write tool:
- If a PMID is available: `statistical_analysis_PMID{pmid}.md`
- If no PMID but DOI is available: `statistical_analysis_DOI_{sanitized_doi}.md` (replace `/` with `_`)
- Otherwise: `statistical_analysis_{sanitized_title}.md` (first 50 chars of title, spaces to underscores, special chars removed)

## Report Format

Present the findings in this structure:

```
# Statistical Methods Analysis: [PAPER TITLE]

## Paper Metadata

| Field | Value |
|-------|-------|
| **Title** | [full title] |
| **Authors** | [author list] |
| **Journal** | [journal name] |
| **Date** | [publication date] |
| **PMID** | [PMID or N/A] |
| **DOI** | [DOI or N/A] |
| **PMC ID** | [PMC ID or N/A] |
| **Citations** | [count] ([influential count] influential) |
| **Fields of Study** | [from Semantic Scholar] |

## Paper Summary

2-3 paragraph summary of what the researchers did, their approach, and their main findings. Write for a reader who is familiar with biomedical research but may not be an expert in this specific topic.

## Field Significance

Based on citation data and recent reviews in the field:
- Impact assessment (citation count in context of field and time since publication)
- Novelty (was this approach new at the time?)
- Current relevance (has the field moved on, or are these findings still important?)

## Statistical Methods Inventory

For each statistical test or model identified, present:

### [Method N]: [Test Name]

| Aspect | Detail |
|--------|--------|
| **Test used** | [exact name of statistical test/model] |
| **Applied to** | [what data or comparison] |
| **Explicit assumptions** | [what the authors stated] |
| **Implicit assumptions** | [required by the test but unstated] |
| **Result** | [statistic, p-value, CI, effect size as reported] |
| **Authors' interpretation** | [what they concluded] |
| **Plain-language meaning** | [translation for non-statisticians] |

[Repeat for every test/method found]

### Summary of All Statistical Methods

| # | Method | Applied To | Key Result | Appropriate? |
|---|--------|-----------|------------|--------------|
| 1 | ... | ... | ... | Yes/No/Partial |

## Statistical Critique

### Appropriateness of Tests
For each test: was it the right choice for this data type, sample size, and study design? Flag any mismatches.

### Assumption Validity
Assess whether the key assumptions of each test are likely met. Consider: normality, independence, homoscedasticity, proportional hazards, linearity, etc.

### Multiple Comparisons
Were multiple comparisons corrected for? Was the correction method appropriate?

### Effect Sizes and Clinical Significance
Were effect sizes reported? Is statistical significance accompanied by practical/clinical significance?

### Sample Size and Power
Was the sample size adequate? Was a power analysis performed?

### Interpretation Assessment
Are the authors' conclusions well-supported by their statistical results? Are there claims that go beyond what the data show?

### Alternative Interpretations
Could the data be explained by confounders, selection bias, or alternative models?

### Missing Analyses
What additional statistical analyses should have been performed?

## Limitations of This Analysis

- Full text available: [Yes — from PMC/Europe PMC] or [No — analysis based on abstract only]
- Sections analyzed: [list which sections were available]
- [Note any gaps — e.g., supplementary materials not accessible, figures/tables not parseable from XML]

## Data Sources

| Source | Queried | Result |
|--------|---------|--------|
| PubMed / PMC | Yes | [success/failure, full text or abstract only] |
| Europe PMC | [Yes/No] | [success/failure] |
| Semantic Scholar | Yes | [success/failure, citation count] |
| PubMed (reviews) | Yes | [N reviews found] |
```

## Important Notes

- **Be factual.** Only report statistical methods and results that are explicitly described in the paper text. Do not invent findings.
- **Be thorough.** Every statistical test mentioned — even briefly — should be catalogued.
- **Be fair.** The critique should acknowledge good statistical practices as well as weaknesses.
- **Abstract-only limitation.** If full text was not available, state this prominently in every section. Abstracts typically mention only 1-2 statistical methods; the full paper may use many more.
- **Supplementary materials.** Note that supplementary methods (often containing additional statistical details) are not accessible through these tools.
