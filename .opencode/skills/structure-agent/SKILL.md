---
name: structure-agent
description: Structural biology deep-dive — all PDB structures, AlphaFold prediction, variant mapping onto structure, domain architecture, and literature
---

Perform a comprehensive structural biology analysis for: **$ARGUMENTS**

Use the MCP tools to gather structural data from all relevant sources, then synthesize a detailed structural report. Follow the steps below in order. If a step fails or returns no data, note the gap and continue.

## Input Parsing

Determine the input type:
1. **PDB ID** (4 characters, e.g., `1TUP`, `6W9C`) — start with PDB lookup, then identify the protein
2. **UniProt accession** (e.g., `P04637`) — resolve gene symbol, then search PDB
3. **Gene symbol** (e.g., `TP53`) — resolve to UniProt, then search PDB

## Data Gathering Steps

### 1. Protein Identity

- If starting from a gene symbol, call `uniprot_search` with the symbol and `organism_id:9606` (reviewed: true) to get the canonical UniProt accession.
- Call `uniprot_get_protein` with the accession to get protein length, function, and gene name.
- If starting from a PDB ID, call `pdb_get_structure` first, extract the UniProt mapping from entity details, then call `uniprot_get_protein`.

### 2. Domain Architecture

- Call `interpro_get_domains` with the UniProt accession.
- Extract all domains with their positions — these will be mapped onto structures later.
- Call `uniprot_get_features` with the UniProt accession to get detailed feature annotations (active sites, binding sites, disulfide bonds, metal binding).

### 3. All Experimental Structures (PDB)

- Call `pdb_search` with the gene symbol or protein name (limit: 20) to find all available structures.
- For the **top 5 structures by resolution**, call `pdb_get_structure` to get detailed metadata:
  - Experimental method (X-ray, cryo-EM, NMR)
  - Resolution
  - Entity details (chains, sequence ranges covered)
  - Ligands and small molecules bound
  - Deposition and release dates
  - Source organism
- Group structures by:
  - **Apo structures** (no ligands)
  - **Ligand-bound** (with drug/inhibitor/cofactor)
  - **Complex structures** (with DNA, RNA, or other proteins)
  - **Mutant structures** (engineered mutations)

### 4. Structure Coverage Analysis

- For each detailed PDB structure, note which residue ranges of the full-length protein are resolved.
- Identify gaps — regions of the protein with NO experimental structure coverage.
- Map domain boundaries onto PDB coverage to identify which domains have been structurally characterized.

### 5. AlphaFold Predicted Structure

- Call `alphafold_get_prediction` with the UniProt accession.
- Extract:
  - Overall mean pLDDT confidence
  - Per-region confidence breakdown (very high >90, confident 70-90, low 50-70, very low <50)
  - Identify disordered regions (pLDDT < 50) — these are likely intrinsically disordered
- Compare AlphaFold coverage with PDB coverage — AlphaFold covers the full sequence, so note regions only predicted (not experimentally determined).

### 6. Pathogenic Variant Mapping

- Call `clinvar_search` with the gene symbol to get pathogenic/likely pathogenic variants.
- Map variant positions onto:
  - **Domains**: Which domains harbor the most pathogenic variants?
  - **Active/binding sites**: Do any variants hit catalytic or functional residues?
  - **Structural context**: Are variants in structured (high pLDDT) or disordered regions?
- Identify mutation hotspots — positions or domains with multiple pathogenic variants.

### 7. Protein Interactions & Interfaces

- Call `string_get_interactions` with the gene symbol (species: 9606, limit: 10, required_score: 700).
- Cross-reference interaction partners with PDB complex structures — are any STRING partners seen in co-crystal structures?
- Note which structural regions are involved in protein-protein interfaces (from PDB entity details).

### 8. Literature

- Call `pubmed_search` with the gene symbol + "crystal structure" OR "cryo-EM structure" (max_results: 8).
- Focus on recent structural papers, especially those describing new structures, drug binding, or mechanistic insights.

## Report Format

```
# Structural Biology Report: [PROTEIN NAME] ([GENE SYMBOL])

## Summary
One-paragraph overview: structural knowledge for this protein, key structural features, druggable pockets, and gaps in structural coverage.

## Protein Overview
- Gene: [symbol], Protein: [name], UniProt: [accession]
- Length: [N] amino acids
- Function: [brief description]
- Key domains: [list with positions]

## Domain Architecture
[1]---[Domain A: 50-150]---[Domain B: 200-350]---[N]
Table of all domains with InterPro/Pfam IDs, positions, and descriptions.

## Experimental Structures

### Structure Inventory
| PDB ID | Method | Resolution | Chains | Coverage | Ligands | Year |
|--------|--------|-----------|--------|----------|---------|------|
[All structures found, sorted by resolution]

**Total structures:** [N] | **Best resolution:** [X.X Å] | **Methods:** X-ray ([N]), Cryo-EM ([N]), NMR ([N])

### Key Structures Highlighted

#### Best Overall: [PDB ID] — [Title]
- Method: [X-ray/Cryo-EM], Resolution: [X.X Å]
- Coverage: residues [start]-[end] ([%] of protein)
- Key features: [what makes this structure notable]

#### Best Ligand-Bound: [PDB ID] — [Title]
- Ligand: [name/identifier]
- Binding site residues: [list]
- Drug relevance: [if applicable]

#### Best Complex: [PDB ID] — [Title]
- Complex partners: [proteins/DNA/RNA]
- Interface residues: [key contacts]

### Structure Coverage Map
Visual representation of which protein regions are covered by experimental structures:
```
Protein:    [1]=========================[N]
Domain A:        [====50-150====]
Domain B:                    [====200-350====]
PDB 1XXX:     [--30---------180--]
PDB 2YYY:              [--120--------300--]
PDB 3ZZZ:                          [--280---350--]
Gap:         [1-29]              [181-199]       [351-N]
```

## AlphaFold Prediction
- Mean pLDDT: [score]
- High-confidence regions (pLDDT > 90): [ranges]
- Confident regions (70-90): [ranges]
- Low-confidence / disordered (< 50): [ranges]
- Comparison with PDB: [which gaps in PDB are predicted well by AlphaFold?]
- Model URLs: [PDB file], [CIF file]

## Pathogenic Variants on Structure

### Variant Hotspots
| Domain/Region | Variants | Significance |
|--------------|----------|-------------|
[Domains ranked by number of pathogenic variants]

### Key Pathogenic Variants
| Variant | Position | Domain | ClinVar Significance | Structural Context |
|---------|----------|--------|---------------------|-------------------|
[Top 10 pathogenic variants with structural interpretation]

**Interpretation:** [Where do pathogenic variants cluster? Do they affect active sites, stability, or interfaces?]

## Protein Interaction Interfaces
- Top interaction partners (STRING):
  | Partner | Score | Co-crystal Structure? |
  |---------|-------|--------------------|
- Structurally characterized interfaces: [which interactions have PDB structures?]
- Key interface residues: [from complex structures]

## Structural Literature
| # | Title | Year | Journal | PMID |
|---|-------|------|---------|------|
[Top 5-8 structural papers]

## Druggability Assessment
- Known ligand-binding sites from PDB structures
- Druggable pockets: [based on ligand-bound structures]
- Allosteric sites: [if any structures reveal allosteric mechanisms]
- Drug candidates: [if any PDB ligands are known drugs or clinical candidates]

## Key Takeaways
- [3-5 bullet points summarizing structural knowledge, gaps, variant hotspots, and druggability]

## Data Sources
List which databases were queried and whether each returned data successfully.
```

Keep the report factual — only include data returned by the tools. Do not fabricate PDB IDs, resolutions, or structural details. If a section has no data, write "No data available from [source]."
