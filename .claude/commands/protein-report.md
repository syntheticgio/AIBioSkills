---
name: protein-report
description: Generate a protein-centric report combining UniProt, InterPro, AlphaFold, PDB, STRING, and post-translational modification data
user-invocable: true
args:
  - name: protein
    description: UniProt accession (e.g., P04637), gene symbol (e.g., TP53), or protein name
    required: true
---

Generate a comprehensive protein-centric report for: **$ARGUMENTS**

Use the MCP tools available to you to gather data from all relevant sources, then synthesize a single structured report. Follow the steps below in order. If a step fails or returns no data, note the gap and continue — do not stop the report.

## Input Parsing

Determine the input type:
1. **UniProt accession** (e.g., `P04637`, `Q9Y6K9`) — use directly
2. **Gene symbol** (e.g., `TP53`, `BRCA1`) — resolve to UniProt accession in step 1
3. **Protein name** (e.g., `tumor protein p53`) — search UniProt in step 1

## Data Gathering Steps

### 1. Protein Identity

- If the input is a gene symbol or protein name, call `uniprot_search` with the query and `organism_id:9606` (reviewed: true) to find the canonical Swiss-Prot accession.
- Call `uniprot_get_protein` with the accession to get full protein details: function description, gene name, organism, sequence length, GO terms, subcellular localization.
- Note the UniProt accession, gene symbol, and protein name for subsequent queries.

### 2. Domain Architecture (InterPro)

- Call `interpro_get_domains` with the UniProt accession.
- Extract all domain entries: InterPro accession, name, type (domain/family/repeat), source database (Pfam, PROSITE, CDD, etc.), and amino acid positions.

### 3. Protein Features & Post-Translational Modifications

- Call `uniprot_get_features` with the UniProt accession and **no type filter** to get all features.
- Categorize features into:
  - **Domains**: Domain, Region
  - **Active/binding sites**: Active site, Binding site, Metal binding, Nucleotide binding
  - **Secondary structure**: Helix, Beta strand, Turn
  - **PTMs**: Modified residue, Cross-link, Lipidation, Glycosylation site, Disulfide bond
  - **Processing**: Signal peptide, Transit peptide, Chain, Propeptide
  - **Variants**: Natural variant, Mutagenesis
- Summarize each category — count entries and highlight the most important ones.

### 4. 3D Structures (PDB)

- Call `pdb_search` with the gene symbol (limit: 10) to find all available experimental structures.
- For the top 3 structures by resolution, call `pdb_get_structure` with the PDB ID to get detailed metadata (method, resolution, entity details).
- Note the best-resolution structure and any ligand-bound structures.

### 5. AlphaFold Predicted Structure

- Call `alphafold_get_prediction` with the UniProt accession.
- Extract: overall pLDDT confidence score, fraction of residues in each confidence bin (very high/confident/low/very low), model URLs.
- Note regions of high vs low confidence — low-confidence regions are often disordered or flexible.

### 6. Protein Interactions (STRING)

- Call `string_get_interactions` with the gene symbol (species: 9606, limit: 15, required_score: 700).
- Extract high-confidence interaction partners with scores and evidence types.

### 7. Tissue Expression (GTEx)

- Call `gtex_get_expression` with the gene symbol.
- Extract the top 10 tissues by expression level and the bottom 5 (tissues where the protein is minimally expressed).

### 8. Clinical Variants

- Call `clinvar_search` with the gene symbol to find pathogenic/likely pathogenic variants.
- Focus on variants that affect protein function (missense, nonsense, frameshift in the coding region).

### 9. Literature

- Call `pubmed_search` with the protein name or gene symbol + "protein structure function" (max_results: 5).
- Extract the most relevant recent publications.

## Report Format

Present the gathered data as a structured report:

```
# Protein Report: [PROTEIN NAME] ([UNIPROT ACCESSION])

## Summary
One-paragraph overview: what this protein does, its key domains, clinical relevance, and structural knowledge.

## Protein Identity
- UniProt accession, entry name, gene symbol(s)
- Protein name (recommended name + alternative names)
- Organism, sequence length (amino acids)
- UniProt entry status (reviewed/unreviewed)

## Function
- Functional description (from UniProt)
- Catalytic activity (if enzyme)
- Subcellular localization(s)
- Key GO terms:
  - Biological Process (top 5)
  - Molecular Function (top 5)
  - Cellular Component (top 3)

## Domain Architecture
Visual-style representation of domains along the sequence:
```
[1]---[Signal peptide: 1-22]---[Domain A: 50-150]---[Domain B: 200-350]---[393]
```
Table of all InterPro domains with positions, source database, and description.

## Post-Translational Modifications
- Phosphorylation sites (count and key residues)
- Ubiquitination, acetylation, methylation sites
- Disulfide bonds
- Glycosylation sites
- Other modifications
- Signal/transit peptides and processing events

## Active Sites & Binding
- Active site residues with positions
- Metal binding sites
- Nucleotide/cofactor binding sites
- Key binding partners and interaction interfaces

## 3D Structures
### Experimental Structures (PDB)
| PDB ID | Title | Method | Resolution | Ligands |
|--------|-------|--------|------------|---------|
Best structure highlighted.

### AlphaFold Prediction
- Overall confidence (mean pLDDT)
- Confidence breakdown by region
- Disordered regions (pLDDT < 50)
- Model URLs (PDB, CIF)

## Protein Interactions
- Top interaction partners from STRING (sorted by confidence)
  | Partner | Score | Evidence Types |
  |---------|-------|---------------|
- Key functional interaction clusters

## Tissue Expression
- Top 10 tissues by expression level (median TPM from GTEx)
- Expression pattern interpretation (ubiquitous vs tissue-specific)

## Clinical Variants
- Pathogenic/likely pathogenic variants affecting this protein
- Variant hotspots (positions with multiple pathogenic variants)
- Domains most affected by pathogenic variation

## Key Publications
Top 3-5 recent papers on this protein's structure and function.

## Data Sources
List which databases were queried and whether each returned data successfully.
```

Keep the report factual — only include data returned by the tools. Do not hallucinate annotations, structures, or modifications. If a section has no data, write "No data available from [source]."
