# Datasets MCP Server — Usage Examples

Example queries you can use with an AI assistant (Claude Desktop, Claude Code, etc.) once the MCP server is configured. The assistant selects the right tool automatically.

## Setup

After installing and configuring the server (see README.md), restart your client. The datasets MCP server tools should appear in your available tools.

---

## NCBI Datasets CLI Tools

### Genome Assembly Information

**Query:** "Get summary information for the human genome"

Uses `datasets_summary_genome`:
```json
{"taxon": "human", "limit": 5}
```

**Query:** "Show me information about genome assembly GCF_000001405.40"
```json
{"accession": "GCF_000001405.40"}
```

**Query:** "What genome assemblies are available for *Drosophila melanogaster*?"
```json
{"taxon": "drosophila melanogaster", "limit": 5}
```

### Gene Information

**Query:** "Tell me about the BRCA1 gene in humans"

Uses `datasets_summary_gene`:
```json
{"symbol": "BRCA1", "taxon": "human"}
```

**Query:** "Get information for gene ID 672"
```json
{"gene_id": "672"}
```

### Batch Gene Queries

**Query:** "Give me a summary of TP53, BRCA1, and EGFR all at once"

Uses `batch_gene_summary`:
```json
{"symbols": "TP53,BRCA1,EGFR", "taxon": "human"}
```

**Query:** "Look up NCBI gene IDs 672, 675, and 7157 together"
```json
{"gene_ids": "672,675,7157"}
```

### Downloading Data

**Query:** "Download the E. coli genome assembly with GFF3 annotations (dry run)"

Uses `datasets_download_genome`:
```json
{"taxon": "escherichia coli", "include": ["genome", "gff3"], "dry_run": true}
```

**Query:** "Download BRCA1 and BRCA2 gene sequences for human"

Uses `datasets_download_gene`:
```json
{"symbol": "BRCA1,BRCA2", "taxon": "human", "include": ["gene", "rna", "protein"]}
```

### Version Check

**Query:** "What version of the NCBI datasets CLI is installed?"

Uses `datasets_version`:
```json
{}
```

---

## BioPython Sequence Analysis Tools

### Pairwise Alignment

**Query:** "Align these two protein sequences and tell me how similar they are:
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH"

Uses `sequence_align`:
```json
{
  "sequence1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
  "sequence2": "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH",
  "sequence_type": "protein",
  "mode": "global"
}
```

**Query:** "Do a local alignment of ATCGATCGATCG and ATCAATCGATCG"
```json
{
  "sequence1": "ATCGATCGATCG",
  "sequence2": "ATCAATCGATCG",
  "mode": "local"
}
```

### Sequence Statistics

**Query:** "What is the GC content of ATCGATCGATCGATCG?"

Uses `sequence_stats`:
```json
{"sequence": "ATCGATCGATCGATCG", "sequence_type": "nucleotide"}
```

**Query:** "Compute the molecular weight and amino acid composition of this protein: MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
```json
{"sequence": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"}
```

**Query:** "Analyze the first sequence in `/data/my_sequences.fasta`"
```json
{"sequence": "/data/my_sequences.fasta"}
```

### FASTA Parsing

**Query:** "List all sequences in `/data/proteins.fasta` that have 'kinase' in the description"

Uses `parse_fasta`:
```json
{"file_path": "/data/proteins.fasta", "filter_pattern": "kinase", "limit": 50}
```

**Query:** "How many sequences are in `/data/genome.fasta`? Show me the first 10."
```json
{"file_path": "/data/genome.fasta", "limit": 10}
```

### 6-Frame Translation

**Query:** "Translate this nucleotide sequence in all 6 reading frames: ATGGATCCAGTGGTCCAGGAGTTTGAGGCCATGCTG"

Uses `sequence_translate`:
```json
{"sequence": "ATGGATCCAGTGGTCCAGGAGTTTGAGGCCATGCTG", "frames": "all"}
```

**Query:** "What are the open reading frames in the forward frames of this sequence?"
```json
{"sequence": "ATGGATCCAGTGGTCCAGGAGTTTGAGGCCATGCTGAAATAA", "frames": "forward"}
```

**Query:** "Translate this sequence using the mitochondrial genetic code"
```json
{"sequence": "ATGGATCCAGTGGTCCAGGAG", "table": 2}
```

---

## Ensembl REST API Tools

### Gene Lookup

**Query:** "Look up the TP53 gene in Ensembl"

Uses `ensembl_lookup_gene`:
```json
{"symbol": "TP53", "species": "homo_sapiens"}
```

**Query:** "What does Ensembl say about ENSG00000141510?"
```json
{"id": "ENSG00000141510"}
```

**Query:** "Look up the insulin gene in mouse via Ensembl"
```json
{"symbol": "Ins2", "species": "mus_musculus"}
```

### Sequence Retrieval

**Query:** "Get the protein sequence for Ensembl gene ENSG00000141510"

Uses `ensembl_get_sequence`:
```json
{"id": "ENSG00000141510", "seq_type": "protein"}
```

**Query:** "Retrieve the cDNA sequence for transcript ENST00000269305 in FASTA format"
```json
{"id": "ENST00000269305", "seq_type": "cdna", "format": "fasta"}
```

**Query:** "Get the CDS sequence for ENST00000275493"
```json
{"id": "ENST00000275493", "seq_type": "cds"}
```

### Cross-Reference Search

**Query:** "Search Ensembl for cross-references to BRAF in human"

Uses `ensembl_search`:
```json
{"symbol": "BRAF", "species": "homo_sapiens"}
```

**Query:** "What Ensembl IDs correspond to the gene symbol EGFR?"
```json
{"symbol": "EGFR"}
```

### Variant Lookup

**Query:** "What known variants are in the BRAF V600 region (chromosome 7:140453136-140453236)?"

Uses `ensembl_get_variants`:
```json
{"region": "7:140453136-140453236", "species": "homo_sapiens", "limit": 50}
```

**Query:** "Show me variants on chromosome 17 from position 7661779 to 7687550"
```json
{"region": "17:7661779-7687550"}
```

### Homolog Discovery

**Query:** "Find mouse orthologs for human TP53 (ENSG00000141510)"

Uses `ensembl_get_homologs`:
```json
{"id": "ENSG00000141510", "target_species": "mus_musculus", "homology_type": "orthologues"}
```

**Query:** "What are the paralogs of ENSG00000141510?"
```json
{"id": "ENSG00000141510", "homology_type": "paralogues"}
```

**Query:** "Show me all homologs of ENSG00000157764 (BRAF) across all species"
```json
{"id": "ENSG00000157764"}
```

---

## UniProt REST API Tools

### Protein Search

**Query:** "Search UniProt for reviewed BRCA1 entries in human"

Uses `uniprot_search`:
```json
{"query": "BRCA1 AND organism_id:9606", "reviewed": true, "limit": 5}
```

**Query:** "Find proteins annotated with 'tumor suppressor' in UniProt for organism 9606"
```json
{"query": "tumor suppressor AND organism_id:9606", "limit": 10}
```

**Query:** "Search UniProt for human kinases reviewed in Swiss-Prot"
```json
{"query": "kinase AND organism_id:9606", "reviewed": true, "limit": 20}
```

### Protein Details

**Query:** "Get the full UniProt entry for human P53 (accession P04637) — function, GO terms, localization"

Uses `uniprot_get_protein`:
```json
{"accession": "P04637"}
```

**Query:** "What does UniProt say about Q9Y6K9?"
```json
{"accession": "Q9Y6K9"}
```

### Protein Features

**Query:** "Show me the protein domains and active sites for UniProt entry P04637"

Uses `uniprot_get_features`:
```json
{"accession": "P04637", "feature_types": ["Domain", "Active site"]}
```

**Query:** "What post-translational modifications are annotated for P38398?"
```json
{"accession": "P38398", "feature_types": ["Modified residue", "Cross-link"]}
```

**Query:** "List all annotated features for human insulin (P01308)"
```json
{"accession": "P01308"}
```

---

## ClinVar Tools

### Clinical Variant Search

**Query:** "Search ClinVar for pathogenic variants in BRCA1"

Uses `clinvar_search`:
```json
{"query": "BRCA1 AND pathogenic", "limit": 10}
```

**Query:** "What is the clinical significance of NM_007294.4:c.5266dupC?"
```json
{"query": "NM_007294.4:c.5266dupC"}
```

**Query:** "Find ClinVar entries related to Lynch syndrome"
```json
{"query": "Lynch syndrome", "limit": 20}
```

**Query:** "Are there any pathogenic variants in TP53 associated with Li-Fraumeni syndrome?"
```json
{"query": "TP53 AND Li-Fraumeni AND pathogenic"}
```

---

## PDB / RCSB Tools

### Structure Details

**Query:** "Get the structure details for PDB entry 1TUP (TP53 bound to DNA)"

Uses `pdb_get_structure`:
```json
{"pdb_id": "1TUP"}
```

**Query:** "What is the resolution and experimental method for hemoglobin structure 4HHB?"
```json
{"pdb_id": "4HHB"}
```

**Query:** "Show me details of the SARS-CoV-2 main protease structure 6LU7"
```json
{"pdb_id": "6LU7"}
```

### Structure Search

**Query:** "Search PDB for structures of the SARS-CoV-2 spike protein"

Uses `pdb_search`:
```json
{"query": "SARS-CoV-2 spike protein", "limit": 5}
```

**Query:** "Find crystal structures of insulin receptor"
```json
{"query": "insulin receptor", "limit": 10}
```

**Query:** "Search for PDB structures of BRCA1 BRCT domain"
```json
{"query": "BRCA1 BRCT domain", "limit": 5}
```

---

## InterPro Tools

### Domain Architecture

**Query:** "What protein domains does TP53 (P04637) have according to InterPro?"

Uses `interpro_get_domains`:
```json
{"accession": "P04637"}
```

**Query:** "Show me the domain architecture for BRCA1 (P38398)"
```json
{"accession": "P38398"}
```

**Query:** "What Pfam domains are in human insulin (P01308)?"
```json
{"accession": "P01308"}
```

---

## STRING Tools

### Protein-Protein Interactions

**Query:** "What proteins interact with TP53 in human?"

Uses `string_get_interactions`:
```json
{"identifiers": "TP53", "species": 9606, "limit": 10}
```

**Query:** "Show me the interaction network for TP53, MDM2, and CDKN2A with high confidence"
```json
{"identifiers": "TP53,MDM2,CDKN2A", "species": 9606, "required_score": 700}
```

**Query:** "What are the physical interaction partners of EGFR in human?"
```json
{"identifiers": "EGFR", "species": 9606, "network_type": "physical", "limit": 15}
```

**Query:** "Find interaction partners for mouse Trp53"
```json
{"identifiers": "Trp53", "species": 10090, "limit": 10}
```

---

## KEGG Tools

### Pathway Search and Retrieval

**Query:** "Search KEGG for pathways related to apoptosis"

Uses `kegg_get_pathway`:
```json
{"query": "apoptosis"}
```

**Query:** "Get details for the human p53 signaling pathway"
```json
{"query": "hsa04115"}
```

**Query:** "Search KEGG for insulin signaling pathways"
```json
{"query": "insulin signaling"}
```

**Query:** "Get details for the mouse MAPK signaling pathway"
```json
{"query": "mmu04010"}
```

---

## NCBI BLAST Tools

### Sequence Similarity Search

**Query:** "BLAST this protein sequence against SwissProt: MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS"

Uses `blast_search`:
```json
{"query": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS", "program": "blastp", "database": "swissprot", "max_hits": 5}
```

**Query:** "Search for nucleotide sequences similar to ATCGATCGATCGATCGATCGATCG in the nt database"
```json
{"query": "ATCGATCGATCGATCGATCGATCG", "program": "blastn", "database": "nt", "max_hits": 10}
```

**Query:** "Find proteins similar to this sequence in PDB: MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"
```json
{"query": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH", "program": "blastp", "database": "pdb", "max_hits": 5}
```

---

## PubMed Tools

### Literature Search

**Query:** "Find recent papers about CRISPR-Cas9 gene therapy"

Uses `pubmed_search`:
```json
{"query": "CRISPR Cas9 gene therapy", "max_results": 10}
```

**Query:** "Search PubMed for reviews on BRCA1 and PARP inhibitor resistance"
```json
{"query": "BRCA1 PARP inhibitor resistance review", "max_results": 5}
```

**Query:** "Find papers about TP53 gain-of-function mutations published recently"
```json
{"query": "TP53 gain-of-function mutations", "max_results": 10}
```

---

## Ensembl Regulation Tools

### Regulatory Features

**Query:** "What regulatory features are in the region around TP53 on chromosome 17?"

Uses `ensembl_get_regulation`:
```json
{"region": "17:7661779-7687550", "species": "homo_sapiens"}
```

**Query:** "Show me promoters and enhancers near BRAF on chromosome 7"
```json
{"region": "7:140719327-140924929", "species": "homo_sapiens"}
```

**Query:** "Are there any regulatory elements in the BRCA1 region?"
```json
{"region": "17:43044295-43170245", "species": "homo_sapiens"}
```

---

## AlphaFold Tools

### Predicted 3D Structures

**Query:** "Get the AlphaFold predicted structure for human TP53 (P04637)"

Uses `alphafold_get_prediction`:
```json
{"uniprot_accession": "P04637"}
```

**Query:** "What is the confidence score for the AlphaFold model of BRCA1 (P38398)?"
```json
{"uniprot_accession": "P38398"}
```

**Query:** "Show me the AlphaFold prediction for human insulin (P01308)"
```json
{"uniprot_accession": "P01308"}
```

---

## Genome Coordinate Tools

### Coordinate Liftover

**Query:** "Convert the BRCA1 region 17:41196312-41277500 from hg19 to hg38"

Uses `liftover_coordinates`:
```json
{"region": "17:41196312-41277500", "source_assembly": "GRCh37", "target_assembly": "GRCh38"}
```

**Query:** "Lift over coordinates 7:140453136-140624564 from GRCh37 to GRCh38"
```json
{"region": "7:140453136-140624564"}
```

**Query:** "Convert this hg38 position back to hg19: 17:43044295-43170245"
```json
{"region": "17:43044295-43170245", "source_assembly": "GRCh38", "target_assembly": "GRCh37"}
```

---

## Reactome Pathway Tools

### Pathway Lookup

**Query:** "Look up Reactome pathway R-HSA-1640170 with its participants"

Uses `reactome_get_pathway`:
```json
{"pathway_id": "R-HSA-1640170", "include_participants": true, "include_hierarchy": true}
```

**Query:** "Search Reactome for pathways involving TP53"
```json
{"query": "TP53", "species": "Homo sapiens", "limit": 5}
```

**Query:** "Find apoptosis-related pathways in Reactome"
```json
{"query": "apoptosis", "limit": 10}
```

**Query:** "Get the cell cycle pathway with its sub-pathways and gene participants"
```json
{"pathway_id": "R-HSA-1640170", "include_participants": true}
```

---

## Gene Ontology Enrichment

### GO Enrichment Analysis

**Query:** "Run GO enrichment on these cancer genes: TP53, BRCA1, EGFR, KRAS, MYC, RB1, APC, PTEN"

Uses `gene_ontology_enrich`:
```json
{"genes": "TP53,BRCA1,EGFR,KRAS,MYC,RB1,APC,PTEN"}
```

**Query:** "What biological processes are enriched in CDK2, CDK4, CDK6, CCND1, CCNE1, RB1?"
```json
{"genes": "CDK2,CDK4,CDK6,CCND1,CCNE1,RB1", "sources": "GO:BP"}
```

**Query:** "Do enrichment analysis including KEGG and Reactome pathways for these inflammatory genes"
```json
{"genes": "IL6,TNF,IL1B,CXCL8,CCL2,IFNG", "sources": "GO:BP,GO:MF,GO:CC,KEGG,REAC", "threshold": 0.01}
```

**Query:** "Run GO enrichment for these mouse genes"
```json
{"genes": "Trp53,Brca1,Egfr,Kras,Myc", "organism": "mmusculus"}
```

---

## COSMIC Somatic Mutations

### Cancer Mutation Search

**Query:** "Search COSMIC for BRAF V600E mutations"

Uses `cosmic_search`:
```json
{"query": "V600E", "gene": "BRAF"}
```

**Query:** "Find somatic mutations in TP53 from lung cancer"
```json
{"query": "TP53", "limit": 30}
```

**Query:** "What KRAS mutations are in COSMIC?"
```json
{"query": "KRAS", "limit": 20}
```

**Query:** "Search for EGFR T790M resistance mutation"
```json
{"query": "T790M", "gene": "EGFR"}
```

---

## OMIM Gene-Disease Associations

### Mendelian Disease Search

**Query:** "What genetic diseases are associated with BRCA1?"

Uses `omim_search`:
```json
{"query": "BRCA1"}
```

**Query:** "Search OMIM for Marfan syndrome"
```json
{"query": "Marfan syndrome"}
```

**Query:** "Get all OMIM entries linked to TP53 (Gene ID 7157)"
```json
{"gene_id": "7157"}
```

**Query:** "Look up OMIM entry 191170"
```json
{"query": "191170"}
```

---

## gnomAD Tools

### Population Allele Frequencies

**Query:** "What is the population frequency of the BRAF V600E variant in gnomAD?"

Uses `gnomad_get_variant`:
```json
{"variant_id": "7-140753336-A-T"}
```

**Query:** "Look up rs11549407 in gnomAD for population frequencies"
```json
{"variant_id": "rs11549407"}
```

**Query:** "Check if variant 17-43093449-A-G is common in gnomAD"
```json
{"variant_id": "17-43093449-A-G"}
```

---

## GTEx Tools

### Tissue Expression

**Query:** "In which tissues is TP53 most highly expressed?"

Uses `gtex_get_expression`:
```json
{"gene": "TP53"}
```

**Query:** "Show me the tissue expression profile for BRCA1"
```json
{"gene": "BRCA1", "dataset": "gtex_v8"}
```

**Query:** "What is the expression of ENSG00000141510 across human tissues?"
```json
{"gene": "ENSG00000141510"}
```

---

## HPO Tools

### Phenotype Search

**Query:** "Search HPO for phenotypes related to seizures"

Uses `hpo_search`:
```json
{"query": "seizure"}
```

**Query:** "What genes are associated with the HPO term for seizures (HP:0001250)?"
```json
{"query": "HP:0001250", "category": "genes", "limit": 20}
```

**Query:** "What diseases are associated with intellectual disability?"
```json
{"query": "HP:0001249", "category": "diseases", "limit": 20}
```

**Query:** "Get details for HPO term HP:0001250"
```json
{"query": "HP:0001250", "category": "term"}
```

---

## Expression Atlas Tools

### Expression Experiment Search

**Query:** "Search Expression Atlas for RNA-seq experiments in human"

Uses `expression_atlas_search`:
```json
{"species": "Homo sapiens", "experiment_type": "rnaseq", "limit": 10}
```

**Query:** "Find differential expression experiments for mouse"
```json
{"species": "Mus musculus", "experiment_type": "differential", "limit": 10}
```

**Query:** "Search for expression experiments related to brain"
```json
{"keyword": "brain", "species": "Homo sapiens"}
```

---

## NCBI Gene Links Tools

### Gene-Gene Relationships

**Query:** "What genes are neighbors of TP53?"

Uses `ncbi_gene_links`:
```json
{"gene_symbol": "TP53", "species": "human"}
```

**Query:** "Find genes related to BRCA1 via NCBI gene links"
```json
{"gene_symbol": "BRCA1", "limit": 15}
```

**Query:** "Get gene neighbors for NCBI Gene ID 7157"
```json
{"gene_id": "7157", "limit": 10}
```

---

## ENCODE Tools

### ENCODE Experiment Search

**Query:** "Search ENCODE for ChIP-seq experiments targeting TP53"

Uses `encode_get_experiments`:
```json
{"query": "TP53", "assay": "ChIP-seq", "organism": "Homo sapiens", "limit": 10}
```

**Query:** "Find ATAC-seq datasets in ENCODE"
```json
{"query": "ATAC-seq", "organism": "Homo sapiens", "limit": 5}
```

**Query:** "Search ENCODE for RNA-seq experiments in mouse"
```json
{"query": "RNA-seq", "organism": "Mus musculus", "limit": 10}
```

---

## Primer Design Tools

### PCR Primer Design

**Query:** "Design PCR primers for this sequence targeting positions 100-200"

Uses `primer_design`:
```json
{"sequence": "ATGGATCCAGTGGTCCAGGAGTTTGAGGCCATGCTG...", "target_start": 100, "target_length": 100}
```

**Query:** "Design 5 primer pairs for this template"
```json
{"sequence": "ATGGATCCAGTGGTCCAGGAG...", "target_start": 50, "target_length": 80, "num_primers": 5}
```

---

## PharmGKB Tools

### Pharmacogenomics Annotations

**Query:** "What pharmacogenomics annotations exist for CYP2D6?"

Uses `pharmgkb_search`:
```json
{"gene": "CYP2D6"}
```

**Query:** "Search PharmGKB for warfarin-related clinical annotations"
```json
{"drug": "warfarin"}
```

**Query:** "Find clinical annotations for the VKORC1 gene"
```json
{"gene": "VKORC1"}
```

---

## Disease Ontology Tools

### Disease Term Search

**Query:** "Search the Disease Ontology for breast cancer"

Uses `disease_ontology_search`:
```json
{"query": "breast cancer"}
```

**Query:** "Look up DOID:1612 in the Disease Ontology"
```json
{"query": "DOID:1612"}
```

**Query:** "Search for Parkinson's disease in the Disease Ontology"
```json
{"query": "Parkinson disease", "limit": 10}
```

---

## Paper Retrieval / Citation Tools

### Full-Text Paper Retrieval

**Query:** "Fetch the full text of PMID 33057194 — I need the Methods and Results sections"

Uses `paper_fetch`:
```json
{"identifier": "33057194", "sections": ["methods", "results"]}
```

**Query:** "Get the complete paper for DOI 10.1038/s41586-020-2649-2 with all sections"
```json
{"identifier": "10.1038/s41586-020-2649-2", "sections": ["all"]}
```

**Query:** "Retrieve the paper titled 'Highly accurate protein structure prediction with AlphaFold'"
```json
{"identifier": "Highly accurate protein structure prediction with AlphaFold", "sections": ["methods", "results"]}
```

**Query:** "Get the Discussion section from PMC7505768"
```json
{"identifier": "PMC7505768", "sections": ["discussion"]}
```

### Citation and Impact Search

**Query:** "Look up citation data for DOI 10.1038/s41586-020-2649-2 on Semantic Scholar"

Uses `semantic_scholar_search`:
```json
{"query": "DOI:10.1038/s41586-020-2649-2"}
```

**Query:** "Search Semantic Scholar for papers on CRISPR base editing"
```json
{"query": "CRISPR base editing", "limit": 5}
```

**Query:** "Find highly cited papers about transformer models in protein structure"
```json
{"query": "transformer protein structure prediction", "limit": 10, "fields_of_study": "Biology"}
```

**Query:** "Get Semantic Scholar metadata for PMID 33057194"
```json
{"query": "PMID:33057194"}
```

---

## Multi-Tool Workflow Examples

These examples show how the assistant combines multiple tools to answer complex questions.

### Cross-Database Gene Investigation

**Query:** "I'm researching TP53. Can you get me its NCBI summary, Ensembl coordinates and transcripts, and UniProt protein function?"

The assistant will call:
1. `datasets_summary_gene` — NCBI gene metadata
2. `ensembl_lookup_gene` — genomic coordinates, biotype, transcripts
3. `uniprot_search` + `uniprot_get_protein` — protein function, GO terms, localization

### Comparative Sequence Analysis

**Query:** "Find the mouse ortholog of human TP53, get both protein sequences from Ensembl, and align them"

The assistant will call:
1. `ensembl_get_homologs` — find the mouse ortholog ID
2. `ensembl_get_sequence` (×2) — get human and mouse protein sequences
3. `sequence_align` — align and compute identity

### Variant Region Analysis

**Query:** "Show me the known variants in the TP53 gene region, then get the UniProt features for the same protein to see which domains they overlap"

The assistant will call:
1. `ensembl_lookup_gene` — get genomic coordinates for TP53
2. `ensembl_get_variants` — get variants in that region
3. `uniprot_get_features` — get domain boundaries to cross-reference

### Protein Family Exploration

**Query:** "Search UniProt for human deubiquitinases, pick the top 3, and show me their domain architecture"

The assistant will call:
1. `uniprot_search` — find deubiquitinase entries
2. `uniprot_get_features` (×3) — get domain features for each

### Batch Gene Download with Preview

**Query:** "I need sequences for BRCA1, BRCA2, and TP53. First show me what would be downloaded, then download if it looks reasonable"

The assistant will call:
1. `datasets_download_gene` with `dry_run: true` — preview
2. `datasets_download_gene` — actual download (after user confirms)

---

## Skill Examples

Skills are slash commands that orchestrate multiple tools into comprehensive reports. Use them in Claude Code by typing the command directly.

### `/gene-report` — Comprehensive Gene Report

**Usage:** `/gene-report TP53`

This generates a full report by automatically calling 9+ tools:
1. `datasets_summary_gene` — NCBI gene metadata
2. `ensembl_lookup_gene` — Genomic coordinates, transcripts
3. `uniprot_search` + `uniprot_get_protein` — Protein function, GO terms
4. `interpro_get_domains` — Domain architecture
5. `uniprot_get_features` — Active sites, binding sites
6. `clinvar_search` — Pathogenic variants
7. `pdb_search` — 3D structures
8. `string_get_interactions` — Protein interaction network
9. `kegg_get_pathway` — Metabolic/signaling pathways

**More examples:**
- `/gene-report BRCA1` — Breast cancer susceptibility gene
- `/gene-report EGFR` — Epidermal growth factor receptor
- `/gene-report ENSG00000141510` — Using Ensembl ID directly
- `/gene-report 7157` — Using NCBI Gene ID

**When to use:** When you want a one-stop summary of everything known about a gene across all major databases. Great for starting research on a new gene or preparing a gene-focused presentation.

### `/variant-report` — Variant Annotation Report

**Usage:** `/variant-report BRCA1`

This generates a variant-focused report by calling:
1. `ensembl_lookup_gene` — Gene coordinates
2. `ensembl_get_variants` — All known variants in the region
3. `clinvar_search` — Clinical significance annotations
4. `uniprot_search` + `uniprot_get_features` — Protein domain context
5. `ensembl_get_sequence` — Protein sequence for position mapping

**More examples:**
- `/variant-report BRCA1` — All variants in BRCA1 with clinical significance
- `/variant-report 7:140453136-140453236` — Variants in a specific BRAF region
- `/variant-report NM_007294.4:c.5266dupC` — Look up a specific variant
- `/variant-report TP53` — Variant landscape for TP53

**When to use:** When investigating the clinical relevance of variants in a gene or region. Useful for understanding which variants are pathogenic, what conditions they cause, and which protein domains are affected.

### `/protein-report` — Protein-Centric Report

**Usage:** `/protein-report P04637`

This generates a protein-focused report by automatically calling 9+ tools:
1. `uniprot_search` + `uniprot_get_protein` — Protein identity, function, GO terms
2. `interpro_get_domains` — Domain architecture
3. `uniprot_get_features` — PTMs, active sites, binding sites, variants
4. `pdb_search` + `pdb_get_structure` — Experimental 3D structures
5. `alphafold_get_prediction` — Predicted structure and confidence
6. `string_get_interactions` — Protein interaction network
7. `gtex_get_expression` — Tissue expression profile
8. `clinvar_search` — Pathogenic variants
9. `pubmed_search` — Recent publications

**More examples:**
- `/protein-report P04637` — Human p53 protein by UniProt accession
- `/protein-report TP53` — Same protein by gene symbol
- `/protein-report P38398` — BRCA1 protein
- `/protein-report Q9Y6K9` — A less well-known protein

**When to use:** When you want a deep dive into a specific protein — its structure, modifications, interactions, and expression. Ideal when starting from a UniProt accession rather than a gene symbol, or when you need detailed post-translational modification and structural data.

---

### `/pathway-report` — Pathway Deep-Dive Report

**Usage:** `/pathway-report R-HSA-1640170`

This skill orchestrates multiple tools for a comprehensive pathway analysis:
1. `reactome_get_pathway` — Pathway details, hierarchy, and participants
2. `batch_gene_summary` — Summaries for key member genes
3. `uniprot_get_protein` (×3-5) — Protein function for top genes
4. `string_get_interactions` (×2-3) — Interactions between pathway members
5. `clinvar_search` (×3-5) — Pathogenic variants in pathway genes
6. `kegg_get_pathway` — KEGG cross-reference
7. `pubmed_search` — Recent literature on the pathway

**More examples:**
- `/pathway-report R-HSA-1640170` — Cell Cycle pathway by Reactome ID
- `/pathway-report apoptosis` — Search for apoptosis pathways
- `/pathway-report TP53` — Find pathways involving TP53
- `/pathway-report R-HSA-109582` — Hemostasis pathway

**When to use:** When you want a deep dive into a biological pathway — its member genes, clinical relevance, protein interactions, and recent literature. Great for understanding pathway biology, identifying druggable nodes, or contextualizing hits from pathway enrichment analysis.

---

### `/drug-target-report` — Drug Target Assessment

**Usage:** `/drug-target-report EGFR`

This skill orchestrates multiple tools for a comprehensive druggability assessment:
1. `uniprot_search` + `uniprot_get_protein` — Target identity, function, target class
2. `interpro_get_domains` + `uniprot_get_features` — Domain architecture, functional sites
3. `pdb_search` + `pdb_get_structure` (×5) — Structural inventory, ligand-bound structures
4. `alphafold_get_prediction` — Structural confidence at binding sites
5. `omim_search` — Mendelian disease associations
6. `clinvar_search` — Pathogenic variants
7. `hpo_search` — Phenotype associations
8. `cosmic_search` — Somatic mutation hotspots
9. `string_get_interactions` — Interaction network, alternative targets
10. `kegg_get_pathway` + `reactome_get_pathway` — Pathway context
11. `pubmed_search` — Drug development literature

**More examples:**
- `/drug-target-report EGFR` — Classic kinase drug target
- `/drug-target-report BRAF` — Oncogenic kinase with known inhibitors
- `/drug-target-report P00533` — EGFR by UniProt accession
- `/drug-target-report CDK4` — Cell cycle kinase target

**When to use:** When evaluating a protein as a potential drug target. Provides disease rationale, structural druggability (ligand-bound structures, binding sites), somatic mutation hotspots, pathway context, and drug development literature.

---

## Agent Examples

Agents are slash commands that perform autonomous, multi-step research by orchestrating several tools.

### `/literature-agent` — Literature Search & Summary

**Usage:** `/literature-agent BRCA1 PARP inhibitor resistance`

This agent searches PubMed for relevant literature and produces a structured summary:
1. `datasets_summary_gene` — Gene context (full name, summary)
2. `uniprot_search` — Protein context (canonical entry)
3. `pubmed_search` (×1-3) — Primary + refined literature searches
4. Deduplication, ranking, and thematic synthesis

**More examples:**
- `/literature-agent TP53 gain-of-function mutations` — Oncology research topic
- `/literature-agent CRISPR Cas9 gene therapy` — Broad topic search
- `/literature-agent BRAF V600E melanoma` — Variant-specific literature
- `/literature-agent Lynch syndrome` — Disease-focused search

**When to use:** When you need a quick literature review on a gene, variant, disease, or research topic. Great for getting up to speed on recent findings, preparing a literature review section, or identifying key papers.

### `/comparative-genomics-agent` — Multi-Species Gene Comparison

**Usage:** `/comparative-genomics-agent TP53 human mouse zebrafish`

This agent compares a gene across species by:
1. `ensembl_lookup_gene` — Reference gene info
2. `datasets_summary_gene` — NCBI gene metadata
3. `ensembl_get_homologs` — Find orthologs in target species
4. `ensembl_get_sequence` (×N) — Retrieve protein sequences for each species
5. `sequence_align` (×N) — Pairwise alignments against reference
6. `interpro_get_domains` — Domain context for conservation analysis

**More examples:**
- `/comparative-genomics-agent BRCA1 human mouse rat chicken` — Breast cancer gene across 4 species
- `/comparative-genomics-agent ENSG00000141510 across vertebrates` — TP53 by Ensembl ID, auto-selects representative vertebrates
- `/comparative-genomics-agent INS human mouse zebrafish` — Insulin gene conservation

**When to use:** When you want to understand how conserved a gene is across species, identify conserved functional domains, or explore evolutionary relationships. Useful for assessing whether model organisms are appropriate for studying a particular gene.

### `/clinical-variant-agent` — Full Clinical Variant Workup

**Usage:** `/clinical-variant-agent BRAF V600E`

This agent performs a comprehensive variant assessment:
1. `gnomad_get_variant` — Population allele frequency
2. `clinvar_search` — Clinical significance and conditions
3. `datasets_summary_gene` + `ensembl_lookup_gene` — Gene context
4. `uniprot_get_features` + `interpro_get_domains` — Domain impact
5. `alphafold_get_prediction` — Structural context
6. `hpo_search` — Phenotype associations
7. `pubmed_search` — Relevant literature

**More examples:**
- `/clinical-variant-agent rs11549407` — By rsID
- `/clinical-variant-agent 7-140753336-A-T` — By genomic coordinates
- `/clinical-variant-agent NM_007294.4:c.5266dupC` — By HGVS notation

**When to use:** When you need a complete workup of a specific variant — population frequency, clinical classification, protein impact, and literature context. Includes ACMG evidence considerations.

### `/gene-list-agent` — Gene List Functional Analysis

**Usage:** `/gene-list-agent TP53,BRCA1,EGFR,PTEN,RB1`

This agent analyzes a gene set across multiple databases:
1. `batch_gene_summary` — Gene metadata for all genes
2. `kegg_get_pathway` (×N) — Shared pathway analysis
3. `string_get_interactions` — Protein interaction network
4. `gtex_get_expression` (×N) — Tissue expression patterns
5. `uniprot_get_protein` (×N) — GO terms and functional themes
6. `clinvar_search` (×N) — Clinical variant burden
7. `hpo_search` (×N) — Phenotype associations

**More examples:**
- `/gene-list-agent MDM2,CDK4,CCND1,RB1,E2F1` — Cell cycle regulators
- `/gene-list-agent KRAS,BRAF,MEK1,ERK1,ERK2` — MAPK pathway members
- `/gene-list-agent IL6,TNF,IL1B,CXCL8,CCL2` — Inflammatory cytokines

**When to use:** When you have a gene list from an experiment (DE analysis, screen hits, pathway members) and want to understand shared pathways, interactions, expression patterns, and clinical relevance.

---

### `/structure-agent` — Structural Biology Deep-Dive

**Usage:** `/structure-agent TP53`

This agent orchestrates multiple tools for a comprehensive structural analysis:
1. `uniprot_search` + `uniprot_get_protein` — Protein identity and function
2. `interpro_get_domains` + `uniprot_get_features` — Domain architecture, active/binding sites
3. `pdb_search` + `pdb_get_structure` (×5) — All experimental structures with detailed metadata
4. `alphafold_get_prediction` — Full-length predicted structure and confidence
5. `clinvar_search` — Pathogenic variant mapping onto structure
6. `string_get_interactions` — Interaction partners, co-crystal structures
7. `pubmed_search` — Structural biology literature

**More examples:**
- `/structure-agent TP53` — Tumor suppressor with many structures
- `/structure-agent P04637` — Same protein by UniProt accession
- `/structure-agent 1TUP` — Start from a specific PDB structure
- `/structure-agent EGFR` — Receptor tyrosine kinase with drug-bound structures

**When to use:** When you need a deep structural analysis — all available PDB structures, structure coverage mapping, AlphaFold predictions, variant hotspots on structure, and druggability. Ideal for structural biology projects, drug design, and understanding structure-function relationships.

---

### `/statistical-methods-agent` — Paper Statistical Analysis

**Usage:** `/statistical-methods-agent 33057194`

This agent fetches a paper's full text and produces a structured analysis of every statistical method:
1. `paper_fetch` — Retrieve paper metadata and full text (Methods, Results, Discussion)
2. `semantic_scholar_search` — Citation count, fields of study, TLDR
3. `pubmed_search` — Recent reviews for field context
4. Claude analyzes the text to inventory and critique all statistical methods
5. Report saved as a markdown file (e.g., `statistical_analysis_PMID33057194.md`)

**More examples:**
- `/statistical-methods-agent 33057194` — By PMID
- `/statistical-methods-agent 10.1038/s41586-020-2649-2` — By DOI
- `/statistical-methods-agent PMC7505768` — By PMC ID
- `/statistical-methods-agent Highly accurate protein structure prediction with AlphaFold` — By title

**When to use:** When you need to critically evaluate the statistical methodology of a paper — for journal club presentations, peer review, systematic reviews, or learning which statistical tests are appropriate for different study designs.

---

## Lab Notebook

The lab notebook feature passively records your research session and synthesizes it into a structured report.

### Starting a Session

**Query:** "Start a lab notebook for my BRAF V600E investigation"

This creates a session log and reminds you to enable the PostToolUse hook.

### Adding Annotations

**Query:** "Note in my lab notebook that the ClinVar results confirmed pathogenicity"

Uses `lab_notebook_annotate`:
```json
{"note": "ClinVar confirmed BRAF V600E is pathogenic with strong evidence (847 submissions)"}
```

### Generating a Report

**Query:** "Generate my lab notebook report"

This reads the full session JSONL log, groups tool calls into research steps, identifies abandoned branches, and produces a polished markdown document with timeline, key findings, and session metadata.

### Checking Status

**Query:** "What's the status of my lab notebook session?"

Shows tool call count, unique tools used, session duration, and start time.

### Full Workflow Example

```
/lab-notebook start Investigating BRAF V600E in melanoma
```
> ... research using gene lookup, variant search, structure analysis ...
```
/lab-notebook annotate The gnomAD frequency is near-zero in healthy populations, confirming somatic origin
```
> ... more research ...
```
/lab-notebook report
```
> Produces `lab-notebook-2026-02-12_001.md` with the full narrative.

---

## Tips for Best Results

1. **Be specific about the organism** — Use common names (human, mouse) or taxonomic IDs (9606, 10090) or Ensembl species names (homo_sapiens, mus_musculus)

2. **Use dry_run for downloads** — Always preview downloads first with `dry_run: true`

3. **Limit results** — For broad queries, set a limit to avoid overwhelming output

4. **Specify data types** — When downloading, be explicit: "with GFF3 and protein sequences"

5. **Chain queries** — Ask multi-part questions; the assistant will orchestrate multiple tools

6. **Use accession IDs when you have them** — Direct lookups by ID (Ensembl, UniProt, NCBI) are faster and more precise than symbol searches

---

## Troubleshooting

**If the assistant doesn't use the tools:**
- Restart your client after configuration changes
- Check that the MCP server is listed in available tools
- For NCBI tools, verify `datasets` CLI is installed: `datasets version`

**If REST API queries fail:**
- Ensembl and UniProt require internet access
- Rate limits are handled automatically, but persistent errors mean you should wait a minute
- Check that species names use Ensembl format: `homo_sapiens`, not `human`

**If downloads fail:**
- Use `dry_run: true` first to preview
- Check available disk space
- Large genomes may take time — try smaller datasets first
