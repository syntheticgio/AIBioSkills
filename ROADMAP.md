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

---

## Planned Features

### NCBI Datasets Enhancements

#### Tier 1 — High Priority (Core Missing Features)

| Feature | Type | Description | Rationale |
|---------|------|-------------|-----------|
| **Virus genome support** | Tool | `datasets_summary_virus` and `datasets_download_virus` for viral genome metadata and downloads | Critical for infectious disease research, COVID-19 variants, viral evolution studies |
| **Taxonomy tools** | Tool | `datasets_summary_taxonomy` for detailed taxonomic metadata and classifications | Fundamental for taxonomic queries, phylogenetic analysis, species identification |
| **Chromosome filtering** | Parameter | Add `--chromosomes` parameter to genome downloads (e.g., `chromosomes: ["X", "Y"]`) | Commonly needed for sex chromosome research, reduces download size |

#### Tier 2 — Medium Priority (Enhanced Functionality)

| Feature | Type | Description | Rationale |
|---------|------|-------------|-----------|
| **Gene accession queries** | Parameter | Support gene queries by protein accession (e.g., `NP_000483.3`) in addition to symbol/ID | Needed when starting from protein data |
| **Batch genome queries** | Tool | Support multiple genome accessions in single query (like `batch_gene_summary`) | Efficiency for comparative genomics |
| **Assembly filtering** | Parameter | Add assembly level (`complete`, `chromosome`, `scaffold`, `contig`) and category filters (`reference`, `representative`) | Better control over data quality |
| **Exclusion flags** | Parameter | Add `--exclude-gff3`, `--exclude-rna`, `--exclude-protein` to download commands | Bandwidth/storage optimization |
| **Virus host filtering** | Parameter | Add `--host` parameter to virus queries for host-specific viral genomes (e.g., dog, human) | Important for zoonotic disease research |

#### Tier 3 — Low Priority (Advanced Features)

| Feature | Type | Description | Rationale |
|---------|------|-------------|-----------|
| **Rehydrate command** | Tool | `datasets_rehydrate` to expand dehydrated datasets | For working with previously downloaded dehydrated packages |
| **Dehydrated downloads** | Parameter | Add `--dehydrated` flag to download commands | Faster initial downloads, can rehydrate later |
| **API key support** | Configuration | Support `--api-key` parameter for higher NCBI rate limits | Power users with NCBI API keys |
| **Output format control** | Parameter | TSV/CSV output options for summary commands | Integration with other analysis pipelines |
| **Assembly source filtering** | Parameter | Filter by GenBank vs RefSeq assemblies | Data source preference control |

### Future Tool Ideas

| Tool Idea | Category | Description |
|-----------|----------|-------------|
| `datasets_download_taxonomy` | NCBI Datasets | Download complete taxonomy data packages |
| `datasets_summary_virus_protein` | NCBI Datasets | Get viral protein summaries |
| `datasets_download_virus_protein` | NCBI Datasets | Download specific viral proteins |
| Chromosome-level genome browser | Visualization | Integration with genome browsers for viewing regions |
| Batch variant annotation | Clinical | Annotate multiple variants at once combining ClinVar, gnomAD, UniProt |

### Known Limitations

- **NCBI CLI dependency**: Most `datasets_*` tools require the NCBI Datasets CLI to be installed. Consider adding REST API fallbacks for core functionality.
- **Download size constraints**: Large genome downloads may timeout or exceed storage limits in some environments.
- **Rate limiting**: Some APIs (Ensembl, UniProt) have rate limits; current implementation has basic retry logic but could be improved.
- **Virus data gap**: Complete absence of viral genome support is a significant limitation for infectious disease researchers.

### References

- [NCBI Datasets Documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/)
- [datasets command-line reference](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/)
- [Exploring and retrieving sequence and metadata with NCBI Datasets (Nature, 2024)](https://www.nature.com/articles/s41597-024-03571-y)
- [NCBI Taxonomy enhanced access (NAR, 2025)](https://academic.oup.com/nar/article/53/D1/D1711/7848842)
