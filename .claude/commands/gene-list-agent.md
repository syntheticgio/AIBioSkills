---
name: gene-list-agent
description: Functional analysis of a gene list — batch summaries, pathway mapping, protein interactions, tissue expression, and phenotype associations
user-invocable: true
args:
  - name: genes
    description: "Comma-separated gene symbols (e.g., 'TP53,BRCA1,EGFR,PTEN,RB1') or a description of the gene set"
    required: true
---

Perform a functional analysis of the gene list: **$ARGUMENTS**

Use the MCP tools to characterize this gene set — summarize each gene, find shared pathways, map interactions, assess tissue expression patterns, and identify phenotype associations. Follow the steps below. If a step fails for a specific gene, note the gap and continue.

## Input Parsing

Extract gene symbols from the input. The user may provide:
- A comma-separated list: `TP53, BRCA1, EGFR, PTEN`
- A space-separated list: `TP53 BRCA1 EGFR PTEN`
- A description with embedded genes — extract the gene symbols

Normalize all symbols to uppercase. If more than 20 genes are provided, work with the first 20 and note the truncation.

## Data Gathering Steps

### 1. Batch Gene Summary

- Call `batch_gene_summary` with all gene symbols (comma-separated) and taxon `human`.
- Extract: full name, description, chromosome location, gene type for each gene.
- If batch_gene_summary fails, fall back to individual `datasets_summary_gene` calls for each gene.

### 2. Pathway Analysis

For each gene (up to 10):
- Call `kegg_get_pathway` with the gene symbol to find associated pathways.
- Collect all pathway IDs and names.
- After all genes are queried, identify **shared pathways** — pathways that appear for 2+ genes in the list.
- Rank shared pathways by the number of input genes they contain.

### 3. Protein Interaction Network

- Call `string_get_interactions` with all gene symbols joined by commas (species: 9606, required_score: 700, limit: 10).
- This shows interactions **between** the input genes and with other proteins.
- Identify which input genes interact directly with each other.

### 4. Tissue Expression Patterns

For each gene (up to 10):
- Call `gtex_get_expression` with the gene symbol.
- Record the top 5 tissues by expression level for each gene.
- After all genes are queried, identify **shared high-expression tissues** — tissues that appear in the top 5 for multiple genes.

### 5. Functional Annotation

For each gene (up to 5 — focus on the most interesting or least familiar):
- Call `uniprot_search` with the gene symbol and `organism_id:9606` (reviewed: true).
- Call `uniprot_get_protein` with the accession to get GO terms and function descriptions.
- Collect all GO Biological Process terms across the gene set.
- Identify recurring GO terms (appearing for 2+ genes) as **enriched biological themes**.

### 6. Clinical Relevance

- Call `clinvar_search` with each gene symbol to check for pathogenic variants.
- Count pathogenic/likely pathogenic variants per gene.
- Note which genes have the most clinical variant data.

### 7. Phenotype Associations

For each gene (up to 10):
- Call `hpo_search` with the gene symbol (category: search) to find associated phenotypes.
- Collect phenotype terms across the set and identify **recurring phenotypes**.

## Report Format

```
# Gene List Functional Analysis

## Input
- Number of genes analyzed: N
- Gene symbols: [list]

## Gene Summary Table

| Gene | Full Name | Chr | Type | ClinVar Pathogenic Variants |
|------|-----------|-----|------|-----------------------------|
| TP53 | Tumor protein p53 | 17 | protein-coding | 1,247 |
| ... | ... | ... | ... | ... |

## Shared Pathway Analysis

### Top Shared Pathways
Pathways containing 2+ genes from the input list:

| Pathway | KEGG ID | Genes in List | Total Genes in Pathway |
|---------|---------|---------------|----------------------|
| p53 signaling | hsa04115 | TP53, BRCA1, PTEN | — |
| ... | ... | ... | ... |

### Per-Gene Pathway Summary
Brief note on unique pathways for each gene not shared with others.

## Protein Interaction Network
- Direct interactions between input genes (from STRING)
- Key hub genes (genes with most connections)
- Notable interaction partners outside the input list

## Tissue Expression Patterns
- Shared high-expression tissues across the gene set
- Tissue expression summary table:
  | Tissue | Genes Highly Expressed | Notable |
  |--------|----------------------|---------|

## Functional Themes
- Recurring GO Biological Process terms (enriched themes)
- Recurring GO Molecular Function terms
- Common subcellular localizations

## Phenotype Associations
- Recurring HPO phenotype terms across the gene set
- Disease associations shared by multiple genes

## Clinical Summary
- Genes with highest ClinVar pathogenic variant burden
- Shared disease associations

## Interpretation
Brief synthesis:
- What biological theme(s) unite this gene set?
- Are they co-expressed in specific tissues?
- Do they share pathway memberships or interact directly?
- What clinical significance does this gene set have collectively?

## Data Sources
List which databases were queried and their success/failure status.
```

Keep the report factual — only include data returned by the tools. For pathway and enrichment analysis, clearly state that you are reporting observed overlaps, not formal statistical enrichment (no p-values). If a gene returns no data from a source, note it in the relevant section.
