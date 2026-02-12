---
name: comparative-genomics-agent
description: Compare a gene across multiple species — find orthologs, retrieve sequences, compute alignments, and summarize conservation
---

Perform a multi-species comparative genomics analysis for: **$ARGUMENTS**

Use the MCP tools available to you to find orthologs, retrieve sequences, align them, and report on conservation. Follow the steps below in order. If a step fails for a particular species, note the gap and continue with the remaining species.

## Input Parsing

Parse the user input to identify:
1. **Gene identifier** — a gene symbol (e.g., `TP53`), Ensembl gene ID (e.g., `ENSG00000141510`), or NCBI Gene ID
2. **Species list** — extract species names. Convert common names to Ensembl species names:
   - human → homo_sapiens
   - mouse → mus_musculus
   - rat → rattus_norvegicus
   - zebrafish → danio_rerio
   - chicken → gallus_gallus
   - frog → xenopus_tropicalis
   - fly / fruit fly → drosophila_melanogaster
   - worm → caenorhabditis_elegans
   - dog → canis_lupus_familiaris
   - cat → felis_catus
   - pig → sus_scrofa
   - cow → bos_taurus

If the user says "across vertebrates" or similar, use: human, mouse, zebrafish, chicken (4 representative species).
If no species are specified, default to: human, mouse, zebrafish.

## Data Gathering Steps

### 1. Reference Gene Information

- If the input is a gene symbol, call `ensembl_lookup_gene` with the symbol and species `homo_sapiens` (or the first species listed) to get the Ensembl gene ID.
- Call `datasets_summary_gene` with the gene symbol (taxon: human) for NCBI gene metadata (full name, summary).
- Note the reference Ensembl gene ID for the next step.

### 2. Find Orthologs

- Call `ensembl_get_homologs` with the reference Ensembl gene ID, `homology_type: "orthologues"`.
- From the results, extract the ortholog Ensembl gene IDs for each of the requested target species.
- If a requested species has no ortholog in the results, note it as "No ortholog found."
- Record the percent identity values reported by Ensembl for each ortholog pair.

### 3. Retrieve Protein Sequences

For the reference gene and each ortholog found:
- Call `ensembl_get_sequence` with the Ensembl gene ID, `seq_type: "protein"`, `format: "json"`.
- Store the protein sequence and its length.

If a protein sequence is not available for a gene ID, try the canonical transcript ID instead.

### 4. Pairwise Alignments

For each ortholog protein sequence, align it against the reference (human) protein:
- Call `sequence_align` with the two protein sequences, `sequence_type: "protein"`, `mode: "global"`.
- Record: alignment score, percent identity, gap count, alignment length.

If there are 3+ species, also consider one key pairwise comparison between distant species (e.g., mouse vs zebrafish) to show the range of divergence.

### 5. Sequence Statistics (Optional)

For the reference protein:
- Call `sequence_stats` with the protein sequence to get molecular weight, amino acid composition.
- Note any unusual composition differences across species if evident from the alignments.

### 6. Domain Conservation Check (Optional)

- Call `interpro_get_domains` with the UniProt accession (if known) or look up via `uniprot_search` for the reference gene.
- Note which key functional domains exist — these regions are expected to be highly conserved.

## Report Format

Present the analysis as a structured comparative genomics report:

```
# Comparative Genomics Report: [GENE SYMBOL]

## Gene Overview
- Full name, function summary (from NCBI/Ensembl)
- Reference species and Ensembl gene ID
- Number of species analyzed

## Ortholog Summary

| Species | Ensembl Gene ID | Protein Length | % Identity to [Reference] | % Positives |
|---------|-----------------|----------------|---------------------------|-------------|
| Human (reference) | ENSG... | 393 aa | — | — |
| Mouse | ENSMUSG... | 390 aa | 77.8% | 85.2% |
| Zebrafish | ENSDARG... | 373 aa | 52.1% | 66.3% |

## Pairwise Alignments

For each species pair aligned:
- **[Reference] vs [Species]**: X% identity, Y gaps, alignment length Z
- Key observations: conserved regions, notable insertions/deletions

## Conservation Analysis

Summarize the overall conservation pattern:
- Which regions are most conserved (relate to known domains if domain data was retrieved)
- Which regions show the most divergence
- Overall trend: is this gene highly conserved, moderately conserved, or rapidly evolving?
- Note any species-specific insertions or deletions

## Functional Domain Context
If domain data was retrieved:
- List key domains with positions
- Note whether these domains span the conserved regions

## Evolutionary Insights
Brief interpretation:
- What does the conservation pattern suggest about functional constraints?
- Are there species-specific adaptations visible in the sequence differences?
- How does the conservation level compare to expectations for this gene family?

## Data Sources
List which databases were queried and whether each returned data successfully for each species.
```

Keep the report factual — only include data returned by the tools. Do not hallucinate sequences, identity scores, or ortholog relationships. If alignment data is unavailable for a species, note "Alignment not performed — sequence unavailable."
