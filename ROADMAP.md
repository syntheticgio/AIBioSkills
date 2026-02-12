# Datasets MCP Server — Roadmap

Tracking document for all implemented tools, skills, and agents.

## Implemented Tools

| Tool | Source | Description |
|------|--------|-------------|
| `datasets_summary_genome` | NCBI CLI | Get genome assembly summaries by accession/taxon |
| `datasets_summary_gene` | NCBI CLI | Get gene summaries by ID/symbol |
| `datasets_download_genome` | NCBI CLI | Download genome assembly data |
| `datasets_download_gene` | NCBI CLI | Download gene data |
| `datasets_version` | NCBI CLI | Get datasets CLI version |
| `batch_gene_summary` | NCBI CLI | Query multiple genes in one call |
| `sequence_align` | BioPython | Pairwise sequence alignment |
| `sequence_stats` | BioPython | Sequence statistics (GC%, MW, etc.) |
| `parse_fasta` | BioPython | Parse and filter FASTA files |
| `sequence_translate` | BioPython | 6-frame nucleotide-to-protein translation |
| `primer_design` | BioPython / Primer3 | PCR primer design for a target sequence |
| `ensembl_lookup_gene` | Ensembl REST | Look up gene by ID or symbol |
| `ensembl_get_sequence` | Ensembl REST | Retrieve nucleotide/protein sequence |
| `ensembl_search` | Ensembl REST | Search cross-references by symbol |
| `ensembl_get_variants` | Ensembl REST | Get variants in a genomic region |
| `ensembl_get_homologs` | Ensembl REST | Find orthologs/paralogs |
| `ensembl_get_regulation` | Ensembl REST | Regulatory features in a genomic region |
| `liftover_coordinates` | Ensembl REST | Convert coordinates between genome assemblies (hg19 ↔ hg38) |
| `uniprot_search` | UniProt REST | Search proteins (Lucene syntax) |
| `uniprot_get_protein` | UniProt REST | Get protein details, GO terms, function |
| `uniprot_get_features` | UniProt REST | Get protein domains/sites/variants |
| `clinvar_search` | ClinVar/NCBI E-utils | Search clinical variant interpretations |
| `pdb_get_structure` | PDB/RCSB REST | Get 3D structure metadata |
| `pdb_search` | PDB/RCSB REST | Search PDB by text/keyword |
| `interpro_get_domains` | InterPro REST | Protein domain/family annotations |
| `string_get_interactions` | STRING REST | Protein-protein interaction networks |
| `kegg_get_pathway` | KEGG REST | Metabolic/signaling pathway info |
| `blast_search` | NCBI BLAST REST | Remote sequence similarity search |
| `pubmed_search` | NCBI E-utils/PubMed | Search PubMed for biomedical literature |
| `alphafold_get_prediction` | AlphaFold REST | Predicted 3D structure and confidence scores |
| `gnomad_get_variant` | gnomAD GraphQL | Population allele frequencies from gnomAD |
| `gtex_get_expression` | GTEx REST | Tissue-specific gene expression levels |
| `hpo_search` | HPO REST (Jax) | Human Phenotype Ontology — gene-phenotype mapping |
| `reactome_get_pathway` | Reactome REST | Detailed reaction-level pathway data from Reactome |
| `gene_ontology_enrich` | g:Profiler REST | GO enrichment analysis for a gene list |
| `cosmic_search` | COSMIC/NLM REST | Somatic mutation data for cancer research |
| `omim_search` | NCBI E-utils/OMIM | Genetic disease associations from OMIM |
| `expression_atlas_search` | EBI Expression Atlas | Baseline/differential expression experiments across species |
| `ncbi_gene_links` | NCBI E-utils | Gene-gene relationships (neighbors, co-expression) |
| `encode_get_experiments` | ENCODE REST | Search ENCODE for ChIP-seq, ATAC-seq, RNA-seq datasets |
| `pharmgkb_search` | PharmGKB REST | Pharmacogenomics gene-drug-variant clinical annotations |
| `disease_ontology_search` | Disease Ontology REST | Standardized disease terms, definitions, and cross-references |
| `paper_fetch` | PubMed/PMC/Europe PMC | Retrieve paper metadata and full text with section extraction |
| `semantic_scholar_search` | Semantic Scholar REST | Paper metadata, citations, impact, and open-access links |
| `lab_notebook_annotate` | Local (stdlib) | Add manual annotations to the lab notebook session log |

## Implemented Skills

| Skill | Description |
|-------|-------------|
| `/project-summary` | Summarize all MCP tools, capabilities, and usage examples |
| `/gene-report` | Comprehensive gene report combining NCBI, Ensembl, UniProt, ClinVar, PDB, InterPro, STRING, and KEGG |
| `/variant-report` | Variant annotation report combining Ensembl variants, ClinVar clinical data, and UniProt protein features |
| `/protein-report` | Protein-centric report: UniProt function, InterPro domains, AlphaFold prediction, PDB structures, STRING interactions, PTMs, GTEx expression |
| `/pathway-report` | Pathway deep-dive: Reactome pathway data, member genes, ClinVar variants, KEGG cross-reference, PubMed literature |
| `/drug-target-report` | Druggability assessment: function, structures, ligand-bound PDBs, disease associations, interaction network, COSMIC mutations, literature |
| `/lab-notebook` | Research session lab notebook — start sessions, annotate, update, report, and status via PostToolUse hook logging |

## Implemented Agents

| Agent | Description |
|-------|-------------|
| `/literature-agent` | Search PubMed for relevant papers given a gene, variant, disease, or topic — produces a structured literature summary |
| `/comparative-genomics-agent` | Multi-species gene comparison using Ensembl homologs, sequence retrieval, and pairwise alignment |
| `/clinical-variant-agent` | Full variant workup: gnomAD frequency, ClinVar significance, domain impact, AlphaFold structure, HPO phenotypes, literature, ACMG considerations |
| `/gene-list-agent` | Functional analysis of a gene list: batch summaries, KEGG pathways, STRING network, GTEx expression, HPO phenotypes |
| `/structure-agent` | Structural biology deep-dive: all PDB structures, AlphaFold, variant mapping, domain architecture, interaction interfaces, literature |
| `/statistical-methods-agent` | Paper statistical analysis: fetches paper full text, inventories every statistical method, analyzes assumptions and appropriateness, produces structured critique |
