---
name: literature-agent
description: Search PubMed for relevant biomedical literature on a gene, variant, disease, or topic and produce a structured literature summary
user-invocable: true
args:
  - name: topic
    description: Gene symbol, variant, disease, or research topic (e.g., "BRCA1 PARP inhibitor resistance", "TP53 gain-of-function", "CRISPR gene therapy")
    required: true
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

Combine results from all searches. Remove duplicates by PMID. Rank by relevance:
1. Articles whose title/abstract directly addresses the query topic
2. Review articles (broad overviews)
3. Recent articles (newer publications)
4. Highly specific mechanistic or clinical studies

Select the top 10-15 most relevant articles for the report.

## Report Format

Present the findings as a structured literature summary:

```
# Literature Report: [TOPIC]

## Background
Brief context paragraph about the gene/variant/topic (drawn from gene/variant context gathered in step 2, or from your general knowledge if it's a broad topic). 1-3 sentences.

## Key Findings from the Literature

Synthesize the main themes and findings across the retrieved articles. Organize by theme, not by individual paper. For example:
- **Mechanism of action** — What do the papers say about how this works?
- **Clinical relevance** — What clinical implications are discussed?
- **Therapeutic approaches** — What treatments or interventions are studied?
- **Open questions** — What gaps or controversies are noted?

Use 3-5 thematic subsections as appropriate for the topic. Cite papers inline as (Author et al., Year; PMID: XXXXX).

## Selected References

For each of the top 10-15 articles, list:
| # | PMID | Title | Authors | Journal | Year |
|---|------|-------|---------|---------|------|

Include the first author + "et al." for papers with multiple authors.

## Abstract Highlights

For the 3-5 most relevant papers, include a 2-3 sentence summary of the abstract, highlighting the key finding or conclusion. Format as:

### [Title] (PMID: XXXXX)
**Authors:** First Author et al. | **Journal:** Journal Name | **Year:** YYYY
Summary of key findings from the abstract.

## Data Sources
- PubMed: [number] articles retrieved
- Gene context: [NCBI / UniProt / ClinVar — note which were queried]
```

Keep the report factual — only include information from the retrieved articles and tool results. Do not hallucinate citations or findings. If abstracts were not available for some articles, note this.
