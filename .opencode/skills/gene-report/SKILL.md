---
name: gene-report
description: Generate a comprehensive gene report by combining data from NCBI, Ensembl, UniProt, ClinVar, PDB, InterPro, STRING, and KEGG
---

Generate a comprehensive, multi-database gene report for: **$ARGUMENTS**

Use the MCP tools available to you to gather data from all relevant sources, then synthesize a single structured report. Follow the steps below in order. If a step fails or returns no data, note the gap and continue — do not stop the report.

## Data Gathering Steps

### 1. Gene Identity & Genomic Context
- Call `datasets_summary_gene` with the gene symbol or ID (taxon: human unless otherwise specified) to get NCBI gene metadata.
- Call `ensembl_lookup_gene` with the symbol (or Ensembl ID if provided) to get Ensembl coordinates, biotype, and transcript count.

### 2. Protein Function & Annotation
- Call `uniprot_search` with the gene symbol and `organism_id:9606` (reviewed: true) to find the canonical Swiss-Prot entry.
- Using the accession from step above, call `uniprot_get_protein` to get function descriptions, GO terms, and subcellular localization.
- Call `interpro_get_domains` with the same UniProt accession to get domain architecture.

### 3. Protein Features
- Call `uniprot_get_features` with the UniProt accession (no type filter) to get domains, active sites, binding sites, and modified residues. Summarize the key features — do not dump the full list.

### 4. Clinical Variants
- Call `clinvar_search` with the gene symbol to find clinical variant interpretations. Summarize the top pathogenic/likely pathogenic variants (up to 5).

### 5. 3D Structures
- Call `pdb_search` with the gene symbol (limit: 5) to find available crystal/cryo-EM structures. Report PDB IDs, titles, methods, and resolutions.

### 6. Protein Interactions
- Call `string_get_interactions` with the gene symbol (species: 9606, limit: 10, required_score: 700) to find high-confidence interaction partners.

### 7. Pathways
- Call `kegg_get_pathway` with the gene symbol to find associated KEGG pathways. List the top pathway names and IDs.

## Report Format

Present the gathered data as a structured report with these sections:

```
# Gene Report: [GENE SYMBOL]

## Summary
One-paragraph overview: what this gene encodes, its primary function, and clinical relevance.

## Gene Identity
- NCBI Gene ID, Ensembl ID, UniProt accession
- Chromosomal location, strand, coordinates
- Biotype, transcript count

## Protein Function
- Full name and alternative names
- Functional description (from UniProt)
- Subcellular localization
- Key GO terms (top 5 Biological Process, top 5 Molecular Function)

## Domain Architecture
- List of InterPro/Pfam domains with positions
- Key functional sites (active sites, binding sites)

## Clinical Significance
- Number of ClinVar entries
- Notable pathogenic variants (up to 5) with conditions
- Associated diseases/phenotypes

## 3D Structures
- Available PDB structures with method and resolution
- Best resolution structure highlighted

## Protein Interactions
- Top interaction partners from STRING with confidence scores
- Brief note on key interactions

## Pathways
- KEGG pathways this gene participates in

## Data Sources
List which databases were queried and whether each returned data successfully.
```

Keep the report factual — only include data returned by the tools. Do not hallucinate annotations. If a section has no data, write "No data available from [source]."
