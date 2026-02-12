---
name: project-summary
description: Summarize all MCP tools, capabilities, and usage examples for this project
user-invocable: true
---

Read the file PROJECT_INVENTORY.md from the project root and present a formatted summary including:

1. **Overview** — What this MCP server does (1-2 sentences)
2. **Available Tools** — List all MCP tools grouped by category (NCBI Datasets CLI, BioPython Analysis, Ensembl REST API, UniProt REST API, ClinVar, PDB/RCSB, InterPro, STRING, KEGG, NCBI BLAST, PubMed, Ensembl Regulation, AlphaFold, gnomAD, GTEx, HPO, Genome Coordinates, Reactome, Gene Ontology Enrichment, COSMIC, OMIM, Expression Atlas, NCBI Gene Links, ENCODE, Primer Design, PharmGKB, Disease Ontology, Paper Retrieval / Citation), with a one-line description of each
3. **Quick-Start Examples** — Show 7-8 practical usage examples such as:
   - Looking up a gene: `datasets_summary_gene` with symbol BRCA1 and taxon human
   - Aligning two sequences: `sequence_align` with two protein sequences
   - Ensembl gene lookup: `ensembl_lookup_gene` with symbol TP53
   - Searching UniProt: `uniprot_search` for BRCA1 in human
   - ClinVar search: `clinvar_search` for pathogenic BRCA1 variants
   - PDB search: `pdb_search` for TP53 structures
   - STRING interactions: `string_get_interactions` for TP53
   - BLAST search: `blast_search` with a protein sequence
   - PubMed search: `pubmed_search` for CRISPR gene therapy papers
   - AlphaFold prediction: `alphafold_get_prediction` for P04637
4. **Custom Skills** — Mention all seven skills: `/project-summary`, `/gene-report <gene>` (comprehensive multi-database gene report), `/variant-report <target>` (variant annotation report), `/protein-report <protein>` (protein-centric report with structure, PTMs, interactions, expression), `/pathway-report <pathway>` (pathway deep-dive with Reactome data, member genes, clinical variants, literature), `/drug-target-report <target>` (druggability assessment with structures, binding sites, disease rationale, COSMIC mutations), and `/lab-notebook <subcommand>` (research session lab notebook — start, annotate, update, report, status)
   - Also mention the six agents: `/literature-agent <topic>` (PubMed literature search and summary), `/comparative-genomics-agent <gene> <species...>` (multi-species gene comparison with alignments), `/clinical-variant-agent <variant>` (full variant workup), `/gene-list-agent <genes>` (functional analysis of a gene list), `/structure-agent <target>` (structural biology deep-dive with all PDB structures, AlphaFold, variant mapping), `/statistical-methods-agent <paper>` (paper statistical analysis — fetches full text, inventories every statistical method, critiques assumptions and appropriateness)
5. **Setup Notes** — Brief setup requirements (datasets CLI for NCBI tools, httpx + biopython installed via pip)

Keep the output concise and well-formatted with markdown.
