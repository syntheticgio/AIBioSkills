---
name: literature-agent
description: Search PubMed for relevant biomedical literature on a gene, variant, disease, or topic and produce a structured literature summary
---

Search PubMed for relevant literature on: **$ARGUMENTS**

Use the MCP tools available to you to find, contextualize, and summarize biomedical publications. Follow the steps below in order. If a step fails or returns no data, note the gap and continue.

## Data Gathering Steps

### 1. Identify the Topic Type

Determine whether the input is:
- **A gene symbol** (e.g., `BRCA1`, `TP53`) — if so, also gather gene context in step 2
- **A variant** (e.g., `BRAF V600E`, `NM_007294.4:c.5266dupC`) — if so, also gather variant context
- **A disease or general topic** (e.g., `Lynch syndrome`, `CRISPR gene therapy`) — skip to step 3

### 2. Gene/Variant Context (if applicable)

If the topic includes a gene symbol:
- Call `datasets_summary_gene` with the gene symbol (taxon: human) to get the gene's full name and summary.
- Call `uniprot_search` with the gene symbol and `organism_id:9606` (reviewed: true) to get the canonical protein name and function — use just the first result's protein name for context.

If the topic includes a specific variant:
- Call `clinvar_search` with the variant identifier to get clinical significance and associated conditions.

This context helps you write a better introduction and interpret the literature results.

### 3. Literature Search — Primary Query

- Call `pubmed_search` with the user's topic as the query, `max_results: 15`.
- Review the returned articles: titles, journals, dates, and abstracts.

### 4. Literature Search — Refined Queries (if needed)

If the primary search returns fewer than 5 results, or the topic is broad, run 1-2 additional targeted searches:
- For a gene: try `"[GENE] AND review[pt]"` to find review articles
- For a variant: try `"[GENE] AND [VARIANT] AND pathogenic"`
- For a disease: try `"[DISEASE] AND molecular mechanism"`
- For a therapy: try `"[THERAPY] AND clinical trial"`

Call `pubmed_search` for each refined query with `max_results: 10`.

### 5. Deduplicate and Rank

Remove duplicate articles and rank by relevance:
- Recency (prioritize last 5 years unless older landmark papers are relevant)
- Citation count (if available through Semantic Scholar integration)
- Journal impact
- Directly addresses the topic

Present the top 8-12 articles in the report.

## Report Format

```
# Literature Summary: [TOPIC]

## Topic Overview

### Query
User's search term: [topic]

### Context
If the topic is a gene or variant, provide:
- Gene: [symbol], Full name: [name], UniProt: [accession]
- Function: [brief description from UniProt or NCBI]
- Clinical relevance: [summary of disease associations if applicable]

Or if it's a disease/topic:
- Condition: [name]
- Key characteristics: [brief description]

## Literature Search Results

### Search Strategy
- Primary query: [what was searched]
- Refined queries (if any): [additional searches]
- Total articles found: [N] (showing top [M] most relevant)

## Key Articles

| # | Title | Authors | Journal | Year | PMID | Key Finding |
|---|-------|---------|---------|------|------|-------------|
| 1 | ... | ... | ... | ... | ... | ... |
| 2 | ... | ... | ... | ... | ... | ... |

[For each of top 5-8 articles, include:]

### Article [N]: [Short Title]
- **PMID:** [ID]
- **Authors:** [first author et al.]
- **Journal:** [journal name], Year: [year]
- **Abstract Summary:** [1-2 sentences describing the research]
- **Key Methods:** [study design, sample size if applicable]
- **Key Findings:** [main results relevant to your topic]
- **Relevance:** [why this is important for understanding your topic]

## Thematic Analysis

### Main Research Themes
- [Theme 1]: [N articles] — [brief description of what these papers focus on]
- [Theme 2]: [N articles] — [brief description]
- [Theme 3]: [N articles] — [brief description]

### Consensus Findings
[1-2 paragraphs summarizing what the literature agrees on]

### Areas of Debate or Uncertainty
[Controversies, conflicting findings, or open questions in the literature]

### Recent Advances
[What's new in the field in the last 1-2 years, if applicable]

## Recommended Reading Path

For someone new to this topic, read in this order:
1. [Review article or overview paper] — comprehensive background
2. [Landmark or high-impact paper] — foundational understanding
3. [Recent paper] — current state of the field

## Summary & Next Steps

### Key Takeaways
- [3-5 bullet points summarizing the literature]

### Questions for Further Investigation
- [Open questions or gaps in current knowledge]

### Resources
- PubMed search strategy for deeper exploration: [query to use]
- Related topics to explore: [suggestions]

## Data Sources
- PubMed: [number of results, success/failure]
- Semantic Scholar (if integrated): [citation data availability]
```

## Important Notes

- **Be accurate.** Report only what is stated in article titles and abstracts. Do not invent findings.
- **Be comprehensive.** Capture the breadth of research on this topic — show themes, not just individual papers.
- **Be fair.** Acknowledge both supportive and critical perspectives in the literature.
- **Be concise.** Summarize in plain language suitable for a researcher who may not specialize in this area.
