---
name: variant-report
description: Generate a variant annotation report combining Ensembl, ClinVar, and UniProt data for a gene or genomic region
user-invocable: true
args:
  - name: target
    description: "Gene symbol (e.g., BRCA1), genomic region (e.g., 7:140453136-140453236), or specific variant (e.g., NM_007294.4:c.5266dupC)"
    required: true
---

Generate a comprehensive variant annotation report for: **$ARGUMENTS**

Use the MCP tools available to you to gather variant data from all relevant sources, then synthesize a single structured report. Follow the steps below. If a step fails, note the gap and continue.

## Determine Input Type

First, classify the input:
- **Gene symbol** (e.g., TP53, BRCA1): Look up genomic coordinates first, then find variants in that region.
- **Genomic region** (e.g., 7:140453136-140453236): Use directly for variant lookup.
- **Specific variant** (e.g., NM_007294.4:c.5266dupC or rs ID): Search ClinVar directly.

## Data Gathering Steps

### For Gene Symbol Input:

#### 1. Gene Coordinates
- Call `ensembl_lookup_gene` with the gene symbol to get chromosomal coordinates.
- Note the chr:start-end for variant queries.

#### 2. Known Variants in the Region
- Call `ensembl_get_variants` with the gene's region (species: homo_sapiens, limit: 100).
- Summarize: total variants found, breakdown by consequence type, any with clinical significance.

#### 3. ClinVar Annotations
- Call `clinvar_search` with the gene symbol (limit: 20) to get clinical interpretations.
- Call `clinvar_search` with "[gene] AND pathogenic" to specifically find pathogenic variants.
- Summarize: total entries, breakdown by clinical significance (pathogenic, likely pathogenic, VUS, benign).

#### 4. Protein Feature Context
- Call `uniprot_search` with the gene symbol and `organism_id:9606` (reviewed: true) to find the UniProt accession.
- Call `uniprot_get_features` with the accession to get domain boundaries and functional sites.
- This gives context for interpreting which variants fall in functionally important regions.

#### 5. Protein Sequence (for reference)
- Call `ensembl_get_sequence` with the Ensembl gene ID (seq_type: protein) to get the protein sequence length for position mapping context.

### For Genomic Region Input:

Follow steps 2-4 above, using the provided region directly. For step 3, search ClinVar with the region or nearby gene name.

### For Specific Variant Input:

#### 1. ClinVar Lookup
- Call `clinvar_search` with the exact variant description.
- Extract clinical significance, review status, associated conditions, and supporting evidence count.

#### 2. Gene Context
- From the ClinVar result, identify the gene. Then call `ensembl_lookup_gene` to get full gene context.
- Call `uniprot_get_features` to see if the variant falls in a known functional domain.

## Report Format

Present the report with these sections:

```
# Variant Report: [TARGET]

## Summary
One-paragraph overview: what gene/region this covers, key clinical findings, and notable variants.

## Gene / Region Context
- Gene symbol, Ensembl ID, chromosomal location
- Gene function (brief, from prior lookup)
- Protein length and key domains

## Variant Landscape
- Total known variants in the region (from Ensembl)
- Breakdown by consequence type (missense, synonymous, frameshift, etc.)
- Variants with clinical significance annotations

## Clinical Significance (ClinVar)
- Total ClinVar entries for this gene/region
- Breakdown: Pathogenic | Likely Pathogenic | VUS | Likely Benign | Benign
- Top pathogenic variants (up to 10):
  | Variant | Clinical Significance | Condition(s) | Review Status |
  |---------|----------------------|---------------|---------------|
  | ...     | ...                  | ...           | ...           |

## Protein Domain Context
- Domain map showing where key variants fall relative to functional domains
- Highlight variants in active sites, binding domains, or known hotspots

## Interpretation Notes
- Which regions of this gene are most clinically relevant
- Any variant hotspots (clusters of pathogenic variants)
- Note if certain domains are enriched for pathogenic variants

## Data Sources
List which databases were queried and whether each returned data successfully.
```

Keep the report factual — only include data returned by the tools. Do not hallucinate variant interpretations. If clinical significance is uncertain (VUS), state that clearly. If a section has no data, write "No data available from [source]."
