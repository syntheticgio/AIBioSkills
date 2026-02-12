# Datasets MCP Server — Project Inventory

An MCP (Model Context Protocol) server for bioinformatics that wraps the NCBI Datasets CLI, BioPython, Ensembl REST API, and UniProt REST API. Gives AI assistants unified access to genomic data retrieval, protein annotation, and computational sequence analysis.

## MCP Tools

### NCBI Datasets CLI Tools

These tools require the [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) to be installed.

#### `datasets_summary_genome`
Get summary information for genome assemblies.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `accession` | string | one of accession/taxon | Genome assembly accession (e.g., `GCF_000001405.40`) |
| `taxon` | string | one of accession/taxon | Taxon name or NCBI Taxonomy ID |
| `limit` | integer | no | Max results (default: 10) |

#### `datasets_summary_gene`
Get summary information for genes.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gene_id` | string | one of gene_id/symbol | NCBI Gene ID (e.g., `672`) |
| `symbol` | string | one of gene_id/symbol | Gene symbol (e.g., `BRCA1`) |
| `taxon` | string | no | Taxon filter |
| `limit` | integer | no | Max results (default: 10) |

#### `datasets_download_genome`
Download genome assembly data (sequences, annotations, metadata).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `accession` | string | one of accession/taxon | Assembly accession |
| `taxon` | string | one of accession/taxon | Taxon name or ID |
| `include` | array[string] | no | Data types: genome, rna, protein, cds, gff3, gtf, gbff, seq-report |
| `filename` | string | no | Output filename (default: `ncbi_dataset.zip`) |
| `dry_run` | boolean | no | Preview without downloading (default: false) |

#### `datasets_download_gene`
Download gene data (sequences and annotations).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gene_id` | string | one of gene_id/symbol | Gene ID or comma-separated list |
| `symbol` | string | one of gene_id/symbol | Gene symbol or comma-separated list |
| `taxon` | string | no | Taxon filter |
| `include` | array[string] | no | Data types: gene, rna, protein, cds |
| `filename` | string | no | Output filename (default: `ncbi_dataset.zip`) |
| `dry_run` | boolean | no | Preview without downloading (default: false) |

#### `datasets_version`
Get the version of the datasets CLI tool. No parameters.

#### `batch_gene_summary`
Query multiple genes in a single call (requires datasets CLI).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gene_ids` | string | one of gene_ids/symbols | Comma-separated Gene IDs (e.g., `672,675,7157`) |
| `symbols` | string | one of gene_ids/symbols | Comma-separated symbols (e.g., `BRCA1,BRCA2,TP53`) |
| `taxon` | string | no | Taxon filter |
| `limit` | integer | no | Results per gene (default: 10) |

### BioPython Analysis Tools

These tools use BioPython for local computation and do **not** require the datasets CLI.

#### `sequence_align`
Pairwise sequence alignment using BioPython's PairwiseAligner.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence1` | string | yes | First sequence (raw or FASTA format) |
| `sequence2` | string | yes | Second sequence (raw or FASTA format) |
| `sequence_type` | string | no | `protein` or `nucleotide` (default: auto-detect) |
| `mode` | string | no | `global` or `local` (default: `global`) |

**Returns:** aligned sequences, score, identity %, gap count, alignment length.

#### `sequence_stats`
Compute statistics for a nucleotide or protein sequence.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | yes | Raw sequence or path to FASTA file |
| `sequence_type` | string | no | `protein` or `nucleotide` (default: auto-detect) |

**Returns (nucleotide):** length, GC content, base counts, codon usage table.
**Returns (protein):** length, molecular weight, amino acid composition, isoelectric point, aromaticity.

#### `parse_fasta`
Parse a FASTA file and extract sequence metadata.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `file_path` | string | yes | Path to FASTA file |
| `filter_pattern` | string | no | Regex to filter IDs/descriptions |
| `include_sequences` | boolean | no | Include full sequences (default: false) |
| `limit` | integer | no | Max sequences to return (default: 100) |

**Returns:** list of sequence IDs, descriptions, lengths, and optionally sequences.

#### `sequence_translate`
Translate a nucleotide sequence to protein in all 6 reading frames with ORF detection.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | yes | Nucleotide sequence (raw or FASTA format, or file path) |
| `frames` | string | no | `all` (default, 6 frames), `forward` (+1,+2,+3), `reverse` (-1,-2,-3), or specific frame (`+1`, `-2`) |
| `table` | integer | no | NCBI genetic code table (default: 1 = Standard; 2 = Vertebrate Mitochondrial; 11 = Bacterial) |

**Returns:** per-frame translation with protein sequence, stop codon count, and ORFs (>=10 aa) with start position, length, and sequence preview. Also reports the longest ORF across all frames.

### Ensembl REST API Tools

These tools query the [Ensembl REST API](https://rest.ensembl.org/) and do **not** require the datasets CLI. Responses are cached for 5 minutes; rate limits are handled automatically.

#### `ensembl_lookup_gene`
Look up a gene by Ensembl stable ID or by symbol + species.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `id` | string | one of id/symbol | Ensembl stable ID (e.g., `ENSG00000141510`) |
| `symbol` | string | one of id/symbol | Gene symbol (e.g., `TP53`). Requires `species`. |
| `species` | string | no | Species name (default: `homo_sapiens`) |

**Returns:** Ensembl ID, display name, description, species, biotype, genomic coordinates (chr, start, end, strand), assembly name, transcript count, canonical transcript ID.

#### `ensembl_get_sequence`
Retrieve a nucleotide or protein sequence by Ensembl stable ID.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `id` | string | yes | Ensembl stable ID (gene, transcript, or protein) |
| `seq_type` | string | no | `genomic`, `cdna`, `cds`, or `protein` (default: `cdna`) |
| `format` | string | no | `json` or `fasta` (default: `json`) |

**Returns (json):** ID, molecule type, length, sequence string.
**Returns (fasta):** raw FASTA text.

#### `ensembl_search`
Search Ensembl cross-references for a gene symbol in a given species.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `symbol` | string | yes | Gene symbol (e.g., `BRCA1`) |
| `species` | string | no | Species name (default: `homo_sapiens`) |

**Returns:** list of matching IDs with type and database display name.

#### `ensembl_get_variants`
Get known genetic variants in a genomic region.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `region` | string | yes | Genomic region as `chr:start-end` (e.g., `7:140453136-140624564`) |
| `species` | string | no | Species name (default: `homo_sapiens`) |
| `limit` | integer | no | Max variants to return (default: 100) |

**Returns:** list of variants with ID, position, alleles, consequence type, clinical significance, source.

#### `ensembl_get_homologs`
Find homologous genes (orthologs/paralogs) for an Ensembl gene ID.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `id` | string | yes | Ensembl gene ID (e.g., `ENSG00000141510`) |
| `target_species` | string | no | Filter to a specific target species |
| `homology_type` | string | no | `orthologues`, `paralogues`, or `all` (default: `all`) |

**Returns:** list of homologs with type, target ID, target species, target protein ID, percent identity, percent positives.

### UniProt REST API Tools

These tools query the [UniProt REST API](https://rest.uniprot.org/) and do **not** require the datasets CLI. Responses are cached for 5 minutes; rate limits are handled automatically.

#### `uniprot_search`
Search UniProt for proteins using Lucene query syntax.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | UniProt search query (e.g., `BRCA1 AND organism_id:9606`) |
| `limit` | integer | no | Max results (default: 10) |
| `reviewed` | boolean | no | If true, only return Swiss-Prot (reviewed) entries |

**Returns:** list of entries with accession, gene name, protein name, organism, reviewed status, sequence length.

#### `uniprot_get_protein`
Get detailed protein information by UniProt accession.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `accession` | string | yes | UniProt accession (e.g., `P04637`) |

**Returns:** accession, gene name, protein name, organism, length, function descriptions, subcellular locations, GO terms (up to 30), entry type.

#### `uniprot_get_features`
Get protein sequence features (domains, active sites, modifications, variants).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `accession` | string | yes | UniProt accession (e.g., `P04637`) |
| `feature_types` | array[string] | no | Filter by feature types (e.g., `["Domain", "Active site"]`). Omit for all. |

**Returns:** list of features with type, description, start/end positions, evidence count.

### ClinVar / NCBI E-utilities Tools

These tools query ClinVar via [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/) and do **not** require the datasets CLI.

#### `clinvar_search`
Search NCBI ClinVar for clinical variant interpretations.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | ClinVar search query (e.g., `BRCA1 AND pathogenic`) |
| `limit` | integer | no | Max results (default: 20) |

**Returns:** list of variants with UID, title, accession, clinical significance, review status, gene name, associated traits/conditions, variation details.

### PDB / RCSB REST API Tools

These tools query the [RCSB PDB APIs](https://data.rcsb.org/) and do **not** require the datasets CLI.

#### `pdb_get_structure`
Get detailed metadata for a PDB structure by its 4-character PDB ID.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pdb_id` | string | yes | 4-character PDB ID (e.g., `1TUP`, `4HHB`) |

**Returns:** PDB ID, title, experimental method, resolution, deposition/release dates, polymer entity count, molecular weight, and first entity details (description, type, organism, gene names, sequence length).

#### `pdb_search`
Search the RCSB Protein Data Bank for 3D structures by free-text query.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Search text (e.g., `TP53`, `insulin receptor`) |
| `limit` | integer | no | Max results (default: 10) |

**Returns:** total structure count, list of hits with PDB ID, title, experimental method, resolution, search relevance score.

### InterPro REST API Tools

These tools query the [InterPro API](https://www.ebi.ac.uk/interpro/api/) and do **not** require the datasets CLI.

#### `interpro_get_domains`
Get protein domain and family annotations from InterPro by UniProt accession.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `accession` | string | yes | UniProt accession (e.g., `P04637`) |

**Returns:** protein name, length, organism, list of domains with InterPro accession, name, type (family/domain/repeat/etc.), source database, description, and amino acid locations.

### STRING REST API Tools

These tools query the [STRING database API](https://string-db.org/help/api/) and do **not** require the datasets CLI.

#### `string_get_interactions`
Get protein-protein interaction partners from the STRING database.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `identifiers` | string | yes | Protein name(s), comma-separated (e.g., `TP53` or `TP53,MDM2`) |
| `species` | integer | no | NCBI taxonomy ID (default: 9606 = human) |
| `limit` | integer | no | Max partners per query protein (default: 10) |
| `required_score` | integer | no | Min combined score 0-1000 (default: 400 = medium confidence) |
| `network_type` | string | no | `functional` (default) or `physical` |

**Returns:** list of interactions with protein names, combined score, and individual evidence channel scores (experimental, database, textmining, coexpression).

### KEGG REST API Tools

These tools query the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html) and do **not** require the datasets CLI. Note: KEGG API is for academic use only.

#### `kegg_get_pathway`
Search for KEGG pathways by keyword or retrieve details for a specific pathway ID.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Keyword (e.g., `apoptosis`) or KEGG pathway ID (e.g., `hsa04210`) |
| `organism` | string | no | KEGG organism code (default: `hsa` = human) |

**Returns (keyword search):** list of matching pathway IDs and names.
**Returns (pathway ID):** pathway name, description, class, organism, gene list with IDs and descriptions.

### NCBI BLAST REST API Tools

These tools use the [NCBI BLAST Common URL API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html). Searches are asynchronous and may take 30s to several minutes.

#### `blast_search`
Submit a sequence similarity search to NCBI BLAST and retrieve results.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Query sequence (FASTA or raw) |
| `program` | string | no | `blastn`, `blastp`, `blastx`, `tblastn`, `tblastx` (default: `blastp`) |
| `database` | string | no | Target database: `nr`, `swissprot`, `nt`, `refseq_protein`, `pdb` (default: `nr`) |
| `evalue` | number | no | E-value threshold (default: 0.01) |
| `max_hits` | integer | no | Max hits to return (default: 10) |

**Returns:** RID, program, database, query info, total hits, list of hits with accession, title, organism, E-value, bit score, identity %, alignment length, query/subject coverage ranges.

### gnomAD GraphQL API Tools

These tools query the [gnomAD API](https://gnomad.broadinstitute.org/) via GraphQL and do **not** require the datasets CLI.

#### `gnomad_get_variant`
Get population allele frequencies for a genetic variant from gnomAD.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `variant_id` | string | yes | Variant ID in chrom-pos-ref-alt format (e.g., `7-140753336-A-T`) or rsID (e.g., `rs11549407`) |
| `dataset` | string | no | gnomAD dataset version (default: `gnomad_r4`) |

**Returns:** variant ID, rsIDs, chromosome/position/ref/alt, exome and genome allele counts (AC, AN, AF), per-population breakdowns, homozygote/hemizygote counts, filter flags, and top transcript consequence (gene symbol, HGVS, LoF annotation).

### GTEx REST API Tools

These tools query the [GTEx Portal API v2](https://gtexportal.org/api/v2/redoc) and do **not** require the datasets CLI.

#### `gtex_get_expression`
Get tissue-specific gene expression data from GTEx (median TPM across ~54 human tissues).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gene` | string | yes | Gene symbol (e.g., `TP53`) or Ensembl gene ID (e.g., `ENSG00000141510`) |
| `dataset` | string | no | GTEx dataset: `gtex_v8` (default), `gtex_v10`, or `gtex_v7` |

**Returns:** gene symbol, GENCODE ID, description, dataset version, total tissue count, and list of tissues sorted by expression level (highest first) with tissue ID, median TPM, and unit.

### HPO REST API Tools

These tools query the [Human Phenotype Ontology API](https://hpo.jax.org/) and do **not** require the datasets CLI.

#### `hpo_search`
Search the Human Phenotype Ontology for phenotype terms, or get genes/diseases associated with a phenotype.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Search keyword (e.g., `seizure`), HPO term ID (e.g., `HP:0001250`), or gene symbol |
| `category` | string | no | `search` (default), `term` (details for an HPO ID), `genes` (genes for an HPO term), `diseases` (diseases for an HPO term) |
| `limit` | integer | no | Max results (default: 20) |

**Returns (search):** matching HPO terms (ID, name, definition), genes (symbol, Entrez ID), and diseases (ID, name, database).
**Returns (term):** HPO ID, name, definition, synonyms, parent terms, children count.
**Returns (genes):** list of genes associated with the HPO term, with gene symbol, Entrez ID, and disease count.
**Returns (diseases):** list of diseases associated with the HPO term, with disease ID, name, and source database.

### PubMed / NCBI E-utilities Tools

These tools query [PubMed via NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/) and do **not** require the datasets CLI.

#### `pubmed_search`
Search PubMed for biomedical literature with abstracts.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | PubMed search query (e.g., `BRCA1 PARP inhibitor resistance`) |
| `max_results` | integer | no | Max results (default: 10, max: 50) |

**Returns:** list of articles with PMID, title, authors, journal, publication date, DOI, and abstract text.

### Ensembl Regulation Tools

These tools query the [Ensembl REST API](https://rest.ensembl.org/) regulatory feature endpoints.

#### `ensembl_get_regulation`
Get regulatory features (promoters, enhancers, CTCF binding sites, etc.) overlapping a genomic region.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `region` | string | yes | Genomic region as `chr:start-end` (e.g., `17:7661779-7687550`) |
| `species` | string | no | Species name (default: `homo_sapiens`) |

**Returns:** list of regulatory features with ID, feature type, description, coordinates (start, end, strand), activity evidence, and bound motifs (if any).

### AlphaFold REST API Tools

These tools query the [AlphaFold Protein Structure Database API](https://alphafold.ebi.ac.uk/).

#### `alphafold_get_prediction`
Get AlphaFold predicted 3D structure and per-residue confidence scores by UniProt accession.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `uniprot_accession` | string | yes | UniProt accession (e.g., `P04637`) |

**Returns:** UniProt accession, gene name, organism, model creation date, PDB/CIF/BCIF URLs, PAE image URL, mean pLDDT confidence score, per-residue pLDDT breakdown (very high/confident/low/very low), and sequence length.

### Genome Coordinate Tools

These tools use the [Ensembl REST API](https://rest.ensembl.org/) coordinate mapping endpoint.

#### `liftover_coordinates`
Convert genomic coordinates between genome assemblies (e.g., GRCh37/hg19 to GRCh38/hg38).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `region` | string | yes | Genomic region as `chr:start-end` (e.g., `17:41196312-41277500`) |
| `source_assembly` | string | no | Source assembly: `GRCh37` or `GRCh38` (default: `GRCh37`) |
| `target_assembly` | string | no | Target assembly: `GRCh37` or `GRCh38` (default: `GRCh38`) |
| `species` | string | no | Species name (default: `human`) |

**Returns:** list of coordinate mappings with original and mapped regions (chromosome, start, end, strand), source and target assembly names. Multiple mappings may be returned if the region maps to separate locations.

#### `reactome_get_pathway`
Look up or search Reactome pathways. Supports direct pathway lookup by stable ID, keyword/gene search, and retrieval of participating molecules and pathway hierarchy.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pathway_id` | string | no | Reactome stable ID (e.g., `R-HSA-1640170`) for direct lookup |
| `query` | string | no | Search term (gene name, keyword, pathway name). Used when `pathway_id` not provided |
| `species` | string | no | Species name (default: `Homo sapiens`) |
| `include_participants` | boolean | no | Also retrieve participating molecules (default: `false`) |
| `include_hierarchy` | boolean | no | Also retrieve ancestor pathways (default: `false`) |
| `limit` | integer | no | Max search results (default: `10`) |

**Returns:** For lookup: pathway name, description, species, sub-events, compartments, disease association, and optionally participants (genes/proteins with identifiers) and hierarchy (ancestor pathways). For search: list of matching pathways with names and summaries.

#### `gene_ontology_enrich`
Perform GO enrichment analysis on a gene list using g:Profiler. Returns statistically enriched terms from GO (Biological Process, Molecular Function, Cellular Component) and optionally KEGG, Reactome, WikiPathways, and other sources.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `genes` | string | yes | Comma-separated gene symbols (e.g., `TP53,BRCA1,EGFR`) |
| `organism` | string | no | Organism (default: `hsapiens`). Others: `mmusculus`, `dmelanogaster`, `scerevisiae` |
| `sources` | string | no | Comma-separated sources (default: `GO:BP,GO:MF,GO:CC`). Options: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, HP, CORUM, TF, MIRNA |
| `threshold` | number | no | Significance threshold (default: `0.05`) |
| `correction_method` | string | no | Multiple testing correction: `g_SCS` (default), `fdr`, `bonferroni` |
| `no_iea` | boolean | no | Exclude electronically inferred annotations (default: `false`) |
| `ordered` | boolean | no | Treat gene list as ordered/ranked (default: `false`) |
| `background` | string | no | Comma-separated background gene list for custom statistical domain |
| `limit` | integer | no | Max enriched terms to return (default: `25`) |

**Returns:** List of enriched terms sorted by p-value. Each term includes: source (e.g., GO:BP), term ID (e.g., GO:0006915), term name, p-value, term size, query size, intersection size, precision, recall, and intersecting genes from the query.

#### `cosmic_search`
Search the COSMIC (Catalogue Of Somatic Mutations In Cancer) database for somatic mutations. Data sourced via the NLM Clinical Tables API (COSMIC V89, GRCh37).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Search term — gene name, mutation, or keyword |
| `gene` | string | no | Filter results to a specific gene |
| `limit` | integer | no | Max mutations to return (default: `20`, max: `500`) |

**Returns:** List of mutations with: mutation ID, gene name, CDS mutation, amino acid change, description (e.g., substitution - missense), genome position, strand, primary histology, primary tissue site, and PubMed ID.

#### `omim_search`
Search OMIM (Online Mendelian Inheritance in Man) for genetic disease associations. Uses NCBI E-utilities to find OMIM entries and gene→OMIM cross-references.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | no | Gene symbol, disease name, or MIM number |
| `gene_id` | string | no | NCBI Gene ID for gene→OMIM linkage |
| `limit` | integer | no | Max entries to return (default: `20`) |

**Returns:** List of OMIM entries with: MIM number, title, entry type (gene, phenotype, gene_and_phenotype), and alternative titles. Entry types are derived from OMIM prefix conventions (* = gene, # = phenotype, + = gene+phenotype, % = phenotype only).

### Expression Atlas REST API Tools

These tools query the [EBI Expression Atlas](https://www.ebi.ac.uk/gxa/) and do **not** require the datasets CLI.

#### `expression_atlas_search`
Search the Expression Atlas for baseline and differential gene expression experiments across species.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `keyword` | string | no | Search keyword (e.g., `brain`, `cancer`) |
| `species` | string | no | Species filter (e.g., `Homo sapiens`, `Mus musculus`) |
| `experiment_type` | string | no | `baseline`, `differential`, or `rnaseq` |
| `limit` | integer | no | Max results (default: `20`) |

**Returns:** List of experiments with accession, description, species, experiment type, technology platform, number of assays, and experimental factors.

### NCBI Gene Links Tools

These tools query [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/) gene-gene relationships and do **not** require the datasets CLI.

#### `ncbi_gene_links`
Find gene-gene relationships (genomic neighbors, co-expression) via NCBI E-utilities elink.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gene_id` | string | one of gene_id/gene_symbol | NCBI Gene ID (e.g., `7157`) |
| `gene_symbol` | string | one of gene_id/gene_symbol | Gene symbol (e.g., `TP53`). Auto-resolved to Gene ID. |
| `species` | string | no | Species for symbol resolution (default: `human`) |
| `limit` | integer | no | Max neighbor genes to return (default: `10`) |

**Returns:** Query gene ID and symbol, list of neighbor genes with Gene ID, symbol, name, chromosome, and map location.

### ENCODE REST API Tools

These tools query the [ENCODE Project](https://www.encodeproject.org/) search API and do **not** require the datasets CLI.

#### `encode_get_experiments`
Search the ENCODE portal for ChIP-seq, ATAC-seq, RNA-seq, and other functional genomics experiments.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Search term (e.g., `TP53`, `CTCF`) |
| `assay` | string | no | Assay type filter (e.g., `ChIP-seq`, `ATAC-seq`, `RNA-seq`) |
| `organism` | string | no | Organism filter (e.g., `Homo sapiens`, `Mus musculus`) |
| `limit` | integer | no | Max results (default: `10`) |

**Returns:** List of experiments with accession, assay title, description, target (e.g., antibody target), biosample summary, lab, and file count.

### Primer Design Tools (BioPython / Primer3)

These tools use the [primer3-py](https://pypi.org/project/primer3-py/) library for local primer design and do **not** require the datasets CLI. Requires `pip install primer3-py`.

#### `primer_design`
Design PCR primers for a target region within a template sequence.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sequence` | string | yes | Template nucleotide sequence |
| `target_start` | integer | yes | Start position of the target region (0-based) |
| `target_length` | integer | yes | Length of the target region |
| `primer_min_size` | integer | no | Minimum primer length (default: `18`) |
| `primer_opt_size` | integer | no | Optimal primer length (default: `20`) |
| `primer_max_size` | integer | no | Maximum primer length (default: `25`) |
| `primer_min_tm` | number | no | Minimum melting temperature (default: `57.0`) |
| `primer_opt_tm` | number | no | Optimal melting temperature (default: `60.0`) |
| `primer_max_tm` | number | no | Maximum melting temperature (default: `63.0`) |
| `num_primers` | integer | no | Number of primer pairs to return (default: `3`) |

**Returns:** List of primer pairs with left/right primer sequence, start position, length, Tm, GC%, product size, and pair penalty score.

### PharmGKB REST API Tools

These tools query the [PharmGKB API](https://api.pharmgkb.org/) and do **not** require the datasets CLI.

#### `pharmgkb_search`
Search PharmGKB for pharmacogenomics clinical annotations (gene-drug-variant associations).

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `gene` | string | no | Gene symbol to search (e.g., `CYP2D6`) |
| `variant` | string | no | Variant name (e.g., `rs1065852`) |
| `drug` | string | no | Drug name filter (client-side) |
| `limit` | integer | no | Max results (default: `20`) |

**Returns:** List of clinical annotations with annotation ID, name, genes, variant, chemicals (drugs), annotation types, and level of evidence.

### Disease Ontology REST API Tools

These tools query the [Disease Ontology API](https://api.disease-ontology.org/) and do **not** require the datasets CLI.

#### `disease_ontology_search`
Search for standardized disease terms, definitions, synonyms, and cross-references.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Disease name (e.g., `breast cancer`) or DOID (e.g., `DOID:1612`) |
| `limit` | integer | no | Max results (default: `20`) |

**Returns:** For search: list of diseases with DOID, name, definition, synonyms, and cross-references. For direct DOID lookup: full disease entry with definition, synonyms, xrefs, and parent terms.

### Paper Retrieval / Citation Tools

These tools retrieve paper metadata, full text, and citation data from PubMed Central, Europe PMC, and Semantic Scholar. They do **not** require the datasets CLI.

#### `paper_fetch`
Retrieve a biomedical research paper's metadata and full text (when available via PubMed Central or Europe PMC). Extracts specific sections from open-access full-text XML. Falls back to abstract if full text is unavailable.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `identifier` | string | yes | PMID (e.g., `33057194`), DOI (e.g., `10.1038/s41586-020-2649-2`), PMC ID (e.g., `PMC7505768`), or article title |
| `sections` | array[string] | no | Sections to extract: `methods`, `results`, `discussion`, `introduction`, `conclusions`, `all` (default: `["methods", "results"]`) |

**Returns:** PMID, PMC ID, DOI, title, authors, journal, publication date, abstract, full_text_available (boolean), full_text_source (`pmc`, `europe_pmc`, or `abstract_only`), sections (dict mapping section name to extracted text), and section_titles_found (list of all section titles in the paper).

#### `semantic_scholar_search`
Search Semantic Scholar for paper metadata, citation counts, influential citations, TLDRs, and open-access PDF links.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `query` | string | yes | Paper title, `DOI:10.1038/...` for DOI lookup, `PMID:33057194` for PMID lookup, or keyword search |
| `limit` | integer | no | Max results for keyword search (default: 5, ignored for DOI/PMID lookup) |
| `fields_of_study` | string | no | Filter by field: `Medicine`, `Biology`, `Computer Science`, etc. |

**Returns:** List of papers with Semantic Scholar ID, title, authors, year, venue, publication date, abstract (truncated to 2000 chars), TLDR summary, citation count, influential citation count, open-access PDF URL, fields of study, and external IDs (DOI, PMID, etc.).

## Custom Skills

### `/project-summary`
Reads this inventory file and presents a formatted summary of all available tools, capabilities, and usage examples.

### `/gene-report <gene>`
Generates a comprehensive multi-database gene report. Accepts a gene symbol (e.g., `TP53`), Ensembl ID (e.g., `ENSG00000141510`), or NCBI Gene ID (e.g., `7157`).

Orchestrates these tools in sequence:
1. `datasets_summary_gene` — NCBI gene metadata
2. `ensembl_lookup_gene` — Genomic coordinates, transcripts
3. `uniprot_search` + `uniprot_get_protein` — Protein function, GO terms, localization
4. `interpro_get_domains` — Domain architecture
5. `uniprot_get_features` — Active sites, binding sites, modifications
6. `clinvar_search` — Clinical variant interpretations
7. `pdb_search` — Available 3D structures
8. `string_get_interactions` — Protein interaction partners
9. `kegg_get_pathway` — Associated metabolic/signaling pathways

**Output:** Structured markdown report with sections for Gene Identity, Protein Function, Domain Architecture, Clinical Significance, 3D Structures, Protein Interactions, and Pathways.

### `/variant-report <target>`
Generates a variant annotation report. Accepts a gene symbol (e.g., `BRCA1`), genomic region (e.g., `7:140453136-140453236`), or specific variant (e.g., `NM_007294.4:c.5266dupC`).

Orchestrates these tools based on input type:
1. `ensembl_lookup_gene` — Gene coordinates (for gene symbol input)
2. `ensembl_get_variants` — Known variants in the genomic region
3. `clinvar_search` — Clinical significance annotations
4. `uniprot_search` + `uniprot_get_features` — Protein domain context for variant interpretation
5. `ensembl_get_sequence` — Protein sequence for position mapping

**Output:** Structured markdown report with sections for Gene/Region Context, Variant Landscape, Clinical Significance table, Protein Domain Context, and Interpretation Notes.

### `/protein-report <protein>`
Generates a protein-centric report. Accepts a UniProt accession (e.g., `P04637`), gene symbol (e.g., `TP53`), or protein name.

Orchestrates these tools:
1. `uniprot_search` + `uniprot_get_protein` — Protein identity, function, GO terms, localization
2. `interpro_get_domains` — Domain architecture (InterPro/Pfam)
3. `uniprot_get_features` — PTMs, active sites, binding sites, processing, variants
4. `pdb_search` + `pdb_get_structure` — Experimental 3D structures with metadata
5. `alphafold_get_prediction` — Predicted structure and per-residue confidence
6. `string_get_interactions` — Protein interaction network
7. `gtex_get_expression` — Tissue expression profile
8. `clinvar_search` — Pathogenic variants in the gene
9. `pubmed_search` — Recent structure/function publications

**Output:** Structured markdown report with sections for Protein Identity, Function, Domain Architecture, Post-Translational Modifications, Active Sites & Binding, 3D Structures (PDB + AlphaFold), Protein Interactions, Tissue Expression, Clinical Variants, and Key Publications.

### `/pathway-report <pathway>`
Generates a pathway deep-dive report. Accepts a Reactome pathway ID (e.g., `R-HSA-1640170`), pathway name, or gene name.

Orchestrates these tools:
1. `reactome_get_pathway` — Pathway details, hierarchy, and participants
2. `batch_gene_summary` — Summaries for key member genes
3. `uniprot_get_protein` (×3-5) — Protein function for top genes
4. `string_get_interactions` (×2-3) — Interactions between pathway members
5. `clinvar_search` (×3-5) — Pathogenic variants in pathway genes
6. `kegg_get_pathway` — KEGG cross-reference
7. `pubmed_search` — Recent pathway literature

**Output:** Structured markdown report with sections for Pathway Overview, Sub-pathways & Reactions, Key Member Genes, Protein Interaction Network, Clinical Variants in Pathway Members, KEGG Cross-Reference, Recent Literature, and Key Takeaways.

### `/drug-target-report <target>`
Generates a druggability and drug-target assessment. Accepts a gene symbol (e.g., `EGFR`), UniProt accession, or protein name.

Orchestrates these tools:
1. `uniprot_search` + `uniprot_get_protein` + `datasets_summary_gene` — Target identity, function, target class
2. `interpro_get_domains` + `uniprot_get_features` — Domain architecture, druggable functional sites
3. `pdb_search` + `pdb_get_structure` (×5) — Structural inventory, ligand-bound structures and binding sites
4. `alphafold_get_prediction` — Structural confidence at binding sites
5. `omim_search` + `clinvar_search` + `hpo_search` — Disease rationale
6. `cosmic_search` — Somatic mutation hotspots (cancer targets)
7. `string_get_interactions` — Interaction network, alternative targets
8. `kegg_get_pathway` + `reactome_get_pathway` — Pathway context
9. `pubmed_search` — Drug development literature

**Output:** Structured markdown report with sections for Executive Summary, Target Identity, Disease Rationale (OMIM, COSMIC, ClinVar, HPO), Structural Druggability (ligand-bound structures, binding sites, AlphaFold confidence), Interaction Network & Alternative Targets, Literature — Drug Development Status, and Key Takeaways.

## Agents

### `/literature-agent <topic>`
Searches PubMed for relevant biomedical literature and produces a structured summary with thematic synthesis and citations. Accepts a gene symbol, variant, disease, or research topic.

Orchestrates these tools:
1. `datasets_summary_gene` — Gene context (if topic is a gene)
2. `uniprot_search` — Protein context (if topic is a gene)
3. `clinvar_search` — Variant context (if topic is a variant)
4. `pubmed_search` (×1-3) — Primary + refined literature searches

**Output:** Structured markdown report with sections for Background, Key Findings from the Literature (thematic synthesis), Selected References table, and Abstract Highlights.

### `/comparative-genomics-agent <gene> <species...>`
Compares a gene across multiple species by finding orthologs, retrieving protein sequences, computing pairwise alignments, and summarizing conservation patterns.

Orchestrates these tools:
1. `ensembl_lookup_gene` — Reference gene information
2. `datasets_summary_gene` — NCBI gene metadata
3. `ensembl_get_homologs` — Find orthologs in target species
4. `ensembl_get_sequence` (×N) — Retrieve protein sequences per species
5. `sequence_align` (×N) — Pairwise global alignments against reference
6. `sequence_stats` — Reference protein statistics
7. `interpro_get_domains` — Domain context for conservation analysis

**Output:** Structured markdown report with sections for Gene Overview, Ortholog Summary table, Pairwise Alignments, Conservation Analysis, Functional Domain Context, and Evolutionary Insights.

### `/clinical-variant-agent <variant>`
Full clinical variant workup. Accepts rsID, chrom-pos-ref-alt, HGVS, or gene+variant shorthand.

Orchestrates these tools:
1. `gnomad_get_variant` — Population allele frequencies
2. `clinvar_search` — Clinical significance, conditions
3. `datasets_summary_gene` + `ensembl_lookup_gene` — Gene context
4. `uniprot_search` + `uniprot_get_features` + `interpro_get_domains` — Protein domain impact
5. `alphafold_get_prediction` — Structural context and confidence
6. `hpo_search` — Phenotype associations
7. `pubmed_search` — Relevant literature

**Output:** Structured markdown report with Variant Identity, Population Frequency table, Clinical Significance, Protein Impact Assessment, Phenotype Associations, Literature, and ACMG Classification Considerations.

### `/gene-list-agent <genes>`
Functional analysis of a gene list (comma-separated symbols, up to 20 genes).

Orchestrates these tools:
1. `batch_gene_summary` — Gene metadata
2. `kegg_get_pathway` (×N) — Shared pathway analysis
3. `string_get_interactions` — Protein interaction network
4. `gtex_get_expression` (×N) — Tissue expression patterns
5. `uniprot_search` + `uniprot_get_protein` (×N) — GO terms, functional themes
6. `clinvar_search` (×N) — Clinical variant burden
7. `hpo_search` (×N) — Phenotype associations

**Output:** Structured markdown report with Gene Summary table, Shared Pathway Analysis, Protein Interaction Network, Tissue Expression Patterns, Functional Themes, Phenotype Associations, and Interpretation.

### `/structure-agent <target>`
Structural biology deep-dive. Accepts a gene symbol (e.g., `TP53`), UniProt accession (e.g., `P04637`), or PDB ID (e.g., `1TUP`).

Orchestrates these tools:
1. `uniprot_search` + `uniprot_get_protein` — Protein identity and function
2. `interpro_get_domains` + `uniprot_get_features` — Domain architecture, functional sites
3. `pdb_search` + `pdb_get_structure` (×5) — All experimental structures with detailed metadata
4. `alphafold_get_prediction` — Full-length predicted structure and confidence
5. `clinvar_search` — Pathogenic variant mapping onto structure
6. `string_get_interactions` — Interaction partners, interface analysis
7. `pubmed_search` — Structural biology literature

**Output:** Structured markdown report with Domain Architecture, Structure Inventory (all PDB entries), Key Structures Highlighted, Structure Coverage Map, AlphaFold Prediction, Pathogenic Variants on Structure, Protein Interaction Interfaces, Druggability Assessment, and Structural Literature.

### `/statistical-methods-agent <paper>`
Fetches a biomedical research paper's full text, inventories every statistical method used, analyzes assumptions, and produces a structured critique. Accepts a PMID, DOI, PMC ID, or article title. Saves the report as a markdown file in the working directory.

Orchestrates these tools:
1. `paper_fetch` — Retrieve paper metadata and full text (Methods, Results, Discussion)
2. `semantic_scholar_search` — Citation count, fields of study, TLDR
3. `pubmed_search` — Recent reviews for field context
4. Claude reads the text to identify and critique all statistical methods

**Output:** Structured markdown report with Paper Metadata, Paper Summary, Field Significance, Statistical Methods Inventory (per-test table with assumptions and plain-language meaning), Statistical Critique (appropriateness, assumption validity, multiple comparisons, effect sizes, missing analyses), and Limitations. Report also saved as `statistical_analysis_PMID{pmid}.md`.

### Lab Notebook Tools

#### `lab_notebook_annotate`
Add a manual annotation to the lab notebook session log. Used to record research context, decisions, or observations that automated tool call logging can't capture.

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `note` | string | yes | The annotation text to add to the log |
| `session_id` | string | no | Session ID to annotate. If omitted, uses the most recent session. |

**Returns:** Confirmation with session ID and the annotation text.

### Custom Skills (continued)

#### `/lab-notebook <subcommand>`
Research session lab notebook — passively records tool calls via a PostToolUse hook, then synthesizes a structured narrative report. Subcommands: `start <title>`, `annotate <note>`, `update [note]`, `report`, `status`.

## Setup

```bash
# Install dependencies
pip install -e .

# Ensure NCBI datasets CLI is installed (only needed for NCBI tools)
# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

# Run the server
datasets-mcp-server
```

## Architecture

- **`src/datasets_mcp/server.py`** — MCP server with all tool definitions and handlers, organized in sections: NCBI CLI, BioPython, Ensembl REST, UniProt REST, ClinVar, PDB/RCSB, InterPro, STRING, KEGG, BLAST, PubMed, Ensembl Regulation, AlphaFold, gnomAD, GTEx, HPO, Genome Coordinates, Reactome, GO Enrichment, COSMIC, OMIM, Expression Atlas, NCBI Gene Links, ENCODE, Primer Design, PharmGKB, Disease Ontology, Paper Retrieval/Citation
- **`pyproject.toml`** — Project metadata, dependencies (mcp, biopython, httpx)
- **`.claude/commands/project-summary.md`** — Slash command skill definition
- **`CLAUDE.md`** — Project instructions and documentation update checklists
