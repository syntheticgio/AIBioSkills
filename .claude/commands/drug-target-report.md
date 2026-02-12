---
name: drug-target-report
description: Druggability assessment — protein function, 3D structures, ligand-bound PDBs, binding sites, interaction network, disease associations, and literature
user-invocable: true
args:
  - name: target
    description: "Gene symbol (EGFR), UniProt accession (P00533), or protein name"
    required: true
---

Generate a druggability and drug-target assessment report for: **$ARGUMENTS**

Use the MCP tools to gather data from all relevant sources, then synthesize a structured drug-target report. Follow the steps below in order. If a step fails or returns no data, note the gap and continue.

## Input Parsing

Determine the input type:
1. **UniProt accession** (e.g., `P00533`) — use directly
2. **Gene symbol** (e.g., `EGFR`) — resolve to UniProt accession
3. **Protein name** — search UniProt

## Data Gathering Steps

### 1. Target Identity & Function

- If input is a gene symbol or name, call `uniprot_search` with the query and `organism_id:9606` (reviewed: true).
- Call `uniprot_get_protein` with the accession to get:
  - Function description, catalytic activity (enzyme classification)
  - GO terms (molecular function — look for kinase, receptor, enzyme, transporter, channel, protease activities)
  - Subcellular localization (membrane, extracellular, cytoplasmic — affects drug accessibility)
  - Tissue specificity
- Call `datasets_summary_gene` with the gene symbol for NCBI summary and aliases.

### 2. Domain Architecture & Functional Sites

- Call `interpro_get_domains` with the UniProt accession.
- Call `uniprot_get_features` with the UniProt accession to get all features.
- Identify druggable features:
  - **Active sites** — enzymatic targets
  - **Binding sites** — substrate/cofactor pockets
  - **Transmembrane domains** — receptor targets
  - **Signal peptides / extracellular domains** — antibody-accessible regions

### 3. 3D Structures & Ligand-Bound Structures

- Call `pdb_search` with the gene symbol (limit: 15).
- For the **top 5 most relevant structures** (prioritize ligand-bound over apo), call `pdb_get_structure`.
- Categorize:
  - **Ligand-bound structures** — identify the ligand name, type (drug, inhibitor, substrate analog, cofactor), and binding site residues
  - **Apo structures** — useful for pocket identification
  - **Antibody/nanobody complexes** — relevant for biologics
- Note the best-resolution ligand-bound structure — this is the most informative for drug design.

### 4. AlphaFold Prediction

- Call `alphafold_get_prediction` with the UniProt accession.
- Focus on:
  - Confidence at known binding site residues
  - Disordered regions (poor drug targets)
  - Regions with no PDB coverage but high AlphaFold confidence (potential novel binding sites)

### 5. Disease Associations

- Call `omim_search` with the gene symbol to find Mendelian disease associations.
- Call `clinvar_search` with the gene symbol to find pathogenic variants.
- Call `hpo_search` with the gene symbol (category: search) for phenotype associations.
- Note: Genes causing disease when mutated are often good drug targets. Gain-of-function mutations suggest inhibitor targets; loss-of-function suggests activator/replacement targets.

### 6. Protein Interaction Network

- Call `string_get_interactions` with the gene symbol (species: 9606, limit: 15, required_score: 700).
- Identify:
  - **Pathway context** — is this target in a druggable pathway?
  - **Upstream regulators** — alternative targets that control this protein
  - **Downstream effectors** — alternative targets if this protein is undruggable
  - **Co-targeted proteins** — for combination therapy rationale

### 7. Pathway Context

- Call `kegg_get_pathway` with the gene symbol to find pathway memberships.
- Call `reactome_get_pathway` with the gene symbol (limit: 5) for Reactome pathway context.
- Note which pathways this target participates in — drugging a target in a critical signaling pathway has broader therapeutic implications.

### 8. Somatic Mutations (Cancer Targets)

- Call `cosmic_search` with the gene symbol (limit: 20).
- Identify mutation hotspots — recurrently mutated positions suggest oncogenic drivers.
- Map hotspots to domains and binding sites.

### 9. Literature — Drug Development

- Call `pubmed_search` with the gene symbol + "drug target" OR "inhibitor" OR "therapeutic" (max_results: 10).
- Focus on papers about:
  - Known drugs/inhibitors for this target
  - Clinical trials
  - Drug resistance mechanisms
  - Structure-based drug design studies

## Report Format

```
# Drug Target Report: [PROTEIN NAME] ([GENE SYMBOL])

## Executive Summary
One-paragraph druggability assessment: target class, existing drugs (if any), structural knowledge, disease rationale, and overall druggability score (High/Medium/Low/Undruggable).

## Target Identity
- Gene: [symbol], Protein: [name], UniProt: [accession]
- Length: [N] amino acids
- Target class: [kinase / GPCR / nuclear receptor / protease / ion channel / enzyme / other]
- Enzyme classification: [EC number if applicable]
- Localization: [membrane/extracellular/cytoplasmic/nuclear]

## Target Function
- Functional description
- Catalytic activity (if enzyme)
- Key molecular functions (GO terms)
- Biological processes regulated
- Tissue expression pattern (from UniProt)

## Disease Rationale

### Mendelian Diseases (OMIM)
| MIM Number | Disease | Type |
|-----------|---------|------|

### Somatic Cancer Associations (COSMIC)
- Total COSMIC mutations: [N]
- Mutation hotspots: [list top positions with counts]
- Dominant mutation type: [missense/nonsense/frameshift]
| Hotspot | AA Change | Count | Histology |
|---------|-----------|-------|-----------|

### Clinical Variants (ClinVar)
- Pathogenic variants: [N]
- Top conditions: [list]

### Phenotype Associations (HPO)
- Key phenotypes associated with this gene

**Therapeutic hypothesis:** [Based on disease data, what is the rationale for targeting this protein? Inhibitor for gain-of-function, activator for loss-of-function, etc.]

## Structural Druggability

### Structure Inventory
| PDB ID | Method | Resolution | Ligand | Ligand Type | Year |
|--------|--------|-----------|--------|-------------|------|

### Ligand-Bound Structures (Key for Drug Design)
For each ligand-bound structure:
- **[PDB ID]**: [Ligand name] ([drug/inhibitor/substrate analog])
  - Binding site residues: [list]
  - Binding mode: [competitive/allosteric/covalent]
  - Drug relevance: [approved drug / clinical candidate / tool compound]

### Binding Site Analysis
- **Primary binding site**: [residues, domain, pocket characteristics]
- **Allosteric sites**: [if any]
- **Protein-protein interaction interfaces**: [potentially druggable?]

### AlphaFold Confidence at Binding Sites
- pLDDT at active site: [score] — [high confidence = reliable pocket / low = flexible]
- Disordered regions: [list — typically undruggable]

### Druggability Assessment
| Feature | Assessment | Evidence |
|---------|-----------|---------|
| Target class | [e.g., Kinase — highly druggable] | [GO terms, UniProt] |
| Active site accessibility | [exposed / buried / membrane] | [PDB structures] |
| Ligand-bound structures | [N structures with ligands] | [PDB search] |
| Known drugs | [Yes: drug names / No] | [Literature] |
| Structural coverage | [% of protein with PDB structures] | [PDB coverage] |
| AlphaFold confidence | [high / medium / low at binding sites] | [AlphaFold] |

**Overall druggability:** [High / Medium / Low / Challenging]

## Interaction Network & Alternative Targets
- Top interaction partners:
  | Partner | Score | Role | Also Druggable? |
  |---------|-------|------|----------------|
- Pathway context: [KEGG/Reactome pathways]
- Alternative/combination targets: [upstream regulators or downstream effectors]

## Literature — Drug Development Status
| # | Title | Year | Journal | PMID |
|---|-------|------|---------|------|
[Key drug development papers]

**Current status:** [Pre-clinical / Phase I / Phase II / Phase III / Approved drugs exist / No known drug programs]

**Known drugs/inhibitors:** [List any known drugs, their mechanism, and approval status]

## Key Takeaways
- [5-7 bullet points: druggability verdict, best structural starting point, disease rationale, key risks, suggested approach (small molecule vs biologic), and any resistance concerns]

## Data Sources
List which databases were queried and whether each returned data successfully.
```

Keep the report factual — only include data returned by the tools. Do not fabricate drug names, clinical trial status, or binding site details. If a section has no data, write "No data available from [source]." For drug/clinical trial information, rely solely on what PubMed search returns — do not speculate about approval status.
