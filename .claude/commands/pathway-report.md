---
name: pathway-report
description: Generate a pathway deep-dive report combining Reactome pathway data, member gene summaries, ClinVar variants, and PubMed literature
user-invocable: true
args:
  - name: pathway
    description: Reactome pathway ID (e.g., R-HSA-1640170), pathway name, or gene name to find related pathways
    required: true
---

Generate a comprehensive pathway deep-dive report for: **$ARGUMENTS**

Use the MCP tools available to you to gather data from all relevant sources, then synthesize a single structured report. Follow the steps below in order. If a step fails or returns no data, note the gap and continue — do not stop the report.

## Data Gathering Steps

### 1. Pathway Identification
- If the input looks like a Reactome stable ID (starts with `R-`), call `reactome_get_pathway` with `pathway_id` set to that ID, with `include_participants: true` and `include_hierarchy: true`.
- Otherwise, call `reactome_get_pathway` with `query` set to the input to search for matching pathways. Present the top hits and pick the most relevant one, then do a direct lookup with `include_participants: true` and `include_hierarchy: true`.

### 2. Pathway Overview
From the Reactome data, extract:
- Pathway name, stable ID, and species
- Summation/description
- Position in the pathway hierarchy (parent pathways)
- Sub-events (child pathways and reactions) — list the first 10-15
- Compartments (cellular locations)
- Whether it is disease-associated

### 3. Key Member Genes
- From the participants list, identify the gene/protein participants (filter for UniProt/ReferenceGeneProduct entries).
- Select up to 10 key genes from the participant list.
- Call `batch_gene_summary` with those gene symbols (taxon: human) to get brief summaries.
- For the top 3-5 most important genes, call `uniprot_get_protein` to get function descriptions and GO terms.

### 4. Protein Interactions Within the Pathway
- Pick the 2-3 most central genes from the pathway.
- Call `string_get_interactions` for each to see how pathway members interact. Focus on interactions *between* pathway members rather than external interactions.

### 5. Clinical Relevance
- For the top 3-5 key genes, call `clinvar_search` to find pathogenic/likely pathogenic variants.
- Summarize the most clinically significant variants (up to 3 per gene, 10 total max).
- Note which diseases are associated with variants in pathway members.

### 6. KEGG Cross-Reference
- Call `kegg_get_pathway` with the pathway name or key gene to find the equivalent KEGG pathway (if any).
- Note any additional member genes or connections found in KEGG but not Reactome.

### 7. Literature
- Call `pubmed_search` with the pathway name to find recent relevant publications (limit: 10).
- If few results, also search with the names of 2-3 key member genes combined with the pathway topic.

## Report Format

Present the report in this structure:

### Pathway: [Name] ([Reactome ID])

**Summary:** [1-2 sentence pathway description]

**Hierarchy:** [Top-level pathway] > [Parent pathway] > **[This pathway]**

**Compartments:** [cellular locations]

**Disease-associated:** Yes/No

---

#### Sub-pathways & Reactions
[Numbered list of child events, noting their type (Pathway vs Reaction)]

#### Key Member Genes ([N] total participants)
| Gene | Full Name | Function Summary |
|------|-----------|-----------------|
[Table of top 10 genes with brief descriptions]

#### Protein Interaction Network
[Summary of how key pathway members interact, noting interaction scores and key hubs]

#### Clinical Variants in Pathway Members
| Gene | Variant | Clinical Significance | Condition |
|------|---------|----------------------|-----------|
[Table of significant ClinVar variants]

**Disease associations:** [Summary of diseases linked to this pathway through variant data]

#### KEGG Cross-Reference
[KEGG pathway ID and name if found, plus any additional insights]

#### Recent Literature
| # | Title | Year | Journal | PMID |
|---|-------|------|---------|------|
[Table of top 5-8 relevant publications]

#### Key Takeaways
[3-5 bullet points summarizing the most important findings: biological role, clinical significance, key genes, and any notable patterns]
