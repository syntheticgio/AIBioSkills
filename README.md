# Datasets MCP Server

An MCP (Model Context Protocol) server for bioinformatics that provides unified access to NCBI Datasets, Ensembl, and UniProt — plus local sequence analysis via BioPython. AI assistants can search, summarize, and download genomic and protein data without the user needing to know which database to query.

**Compatible with:** Claude Desktop, [OpenCode](https://opencode.ai), and any MCP-compatible client.

## 📖 Database Reference Documentation

Comprehensive documentation for all 24+ integrated bioinformatics databases is available in the **[docs/](docs/)** directory:

- **NCBI Databases:** [Datasets](docs/databases/ncbi-datasets.html), [BLAST](docs/databases/ncbi-blast.html), [ClinVar](docs/databases/ncbi-clinvar.html), [PubMed](docs/databases/ncbi-pubmed.html), [COSMIC](docs/databases/ncbi-cosmic.html), [OMIM](docs/databases/ncbi-omim.html), [Gene Links](docs/databases/ncbi-gene-links.html)
- **Genomic & Sequence:** [Ensembl](docs/databases/ensembl.html), [UniProt](docs/databases/uniprot.html), [RCSB PDB](docs/databases/rcsb-pdb.html), [InterPro](docs/databases/interpro.html), [AlphaFold](docs/databases/alphafold.html)
- **Pathways & Interactions:** [STRING](docs/databases/string.html), [KEGG](docs/databases/kegg.html), [Reactome](docs/databases/reactome.html)
- **Expression & Population:** [gnomAD](docs/databases/gnomad.html), [GTEx](docs/databases/gtex.html), [Expression Atlas](docs/databases/expression-atlas.html), [ENCODE](docs/databases/encode.html)
- **Annotations & Literature:** [HPO](docs/databases/hpo.html), [Disease Ontology](docs/databases/disease-ontology.html), [PubMed Central](docs/databases/pmc.html), [Semantic Scholar](docs/databases/semantic-scholar.html)
- **Pharmacogenomics:** [PharmGKB](docs/databases/pharmgkb.html)

👉 **[View all database docs →](docs/)**

## Features

### NCBI Datasets CLI Tools
- **datasets_summary_genome** — Genome assembly summaries by accession or taxon
- **datasets_summary_gene** — Gene summaries by ID or symbol
- **datasets_download_genome** — Download genome assembly data
- **datasets_download_gene** — Download gene data
- **datasets_version** — Check the datasets CLI version
- **batch_gene_summary** — Query multiple genes in one call

### BioPython Sequence Analysis Tools
- **sequence_align** — Pairwise sequence alignment (global/local, protein/nucleotide)
- **sequence_stats** — Sequence statistics (GC%, molecular weight, codon usage, etc.)
- **parse_fasta** — Parse and filter FASTA files
- **sequence_translate** — 6-frame nucleotide-to-protein translation with ORF detection

### Ensembl REST API Tools
- **ensembl_lookup_gene** — Look up a gene by Ensembl ID or symbol
- **ensembl_get_sequence** — Retrieve genomic, cDNA, CDS, or protein sequence
- **ensembl_search** — Search cross-references for a gene symbol
- **ensembl_get_variants** — Get known variants in a genomic region
- **ensembl_get_homologs** — Find orthologs and paralogs

### UniProt REST API Tools
- **uniprot_search** — Search proteins with Lucene query syntax
- **uniprot_get_protein** — Get protein function, GO terms, subcellular location
- **uniprot_get_features** — Get protein domains, active sites, modifications

### ClinVar Tools
- **clinvar_search** — Search NCBI ClinVar for clinical variant interpretations

### PDB / RCSB Tools
- **pdb_get_structure** — Get 3D structure metadata by PDB ID
- **pdb_search** — Search the Protein Data Bank by text/keyword

### InterPro Tools
- **interpro_get_domains** — Protein domain/family annotations (Pfam, PROSITE, CDD, etc.)

### STRING Tools
- **string_get_interactions** — Protein-protein interaction partners and scores

### KEGG Tools
- **kegg_get_pathway** — Search/retrieve metabolic and signaling pathway info

### NCBI BLAST Tools
- **blast_search** — Remote sequence similarity search (async, supports all BLAST programs)

### PubMed Tools
- **pubmed_search** — Search PubMed for biomedical literature with abstracts

### Ensembl Regulation Tools
- **ensembl_get_regulation** — Regulatory features (promoters, enhancers, CTCF sites) in a region

### AlphaFold Tools
- **alphafold_get_prediction** — Predicted 3D structure and confidence scores by UniProt accession

### gnomAD Tools
- **gnomad_get_variant** — Population allele frequencies from gnomAD (exome + genome, per-population breakdowns)

### GTEx Tools
- **gtex_get_expression** — Tissue-specific gene expression (median TPM across ~54 human tissues)

### HPO Tools
- **hpo_search** — Human Phenotype Ontology search — gene-phenotype mapping, disease associations

### Genome Coordinate Tools
- **liftover_coordinates** — Convert coordinates between genome assemblies (GRCh37/hg19 ↔ GRCh38/hg38)

### Reactome Pathway Tools
- **reactome_get_pathway** — Look up or search Reactome pathways; retrieve pathway details, participants, and hierarchy

### Gene Ontology Enrichment
- **gene_ontology_enrich** — GO enrichment analysis on a gene list via g:Profiler (GO:BP, GO:MF, GO:CC, plus optional KEGG/Reactome/WikiPathways)

### COSMIC Somatic Mutations
- **cosmic_search** — Search COSMIC for somatic mutations in cancer (gene, mutation, histology, tissue site)

### OMIM Gene-Disease Associations
- **omim_search** — Search OMIM for Mendelian disease associations linked to genes

### Expression Atlas Tools
- **expression_atlas_search** — Search EBI Expression Atlas for baseline/differential expression experiments

### NCBI Gene Links Tools
- **ncbi_gene_links** — Gene-gene relationships (neighbors, co-expression) via NCBI E-utilities

### ENCODE Tools
- **encode_get_experiments** — Search ENCODE for ChIP-seq, ATAC-seq, RNA-seq datasets

### Primer Design Tools
- **primer_design** — PCR primer design for a target sequence using Primer3

### PharmGKB Tools
- **pharmgkb_search** — Pharmacogenomics gene-drug-variant clinical annotations

### Disease Ontology Tools
- **disease_ontology_search** — Standardized disease terms, definitions, and cross-references

### Paper Retrieval / Citation Tools
- **paper_fetch** — Retrieve paper metadata and full text (Methods, Results, Discussion) from PubMed Central / Europe PMC
- **semantic_scholar_search** — Paper metadata, citation counts, influential citations, TLDRs, and open-access PDF links

### Lab Notebook Tools
- **lab_notebook_annotate** — Add manual annotations to the lab notebook session log for recording research context and decisions

## Prerequisites

1. **Python 3.10+**

2. **NCBI Datasets CLI** (required only for `datasets_*` and `batch_gene_summary` tools):

   ```bash
   # Conda (recommended)
   conda install -c conda-forge ncbi-datasets-cli

   # Direct download — Linux (amd64)
   curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
   chmod +x datasets && sudo mv datasets /usr/local/bin/

   # Direct download — macOS
   curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
   chmod +x datasets && sudo mv datasets /usr/local/bin/

   # Direct download — Windows (PowerShell)
   curl -o datasets.exe https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe
   ```

   See the [NCBI install guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) for ARM64, Linux 32-bit, and other platforms.

   The Ensembl, UniProt, and BioPython tools work without the NCBI CLI installed.

## Installation

### 1. Install the MCP Server

**From PyPI** (when published):

```bash
pip install datasets-mcp-server
```

**From GitHub:**

```bash
pip install git+https://github.com/youruser/datasets-mcp.git
```

**For development (editable install):**

```bash
git clone https://github.com/youruser/datasets-mcp.git
cd datasets-mcp
pip install -e .
```

### 2. Configure the MCP Server

#### Claude Desktop

Add to your config file:

**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
**Linux**: `~/.config/Claude/claude_desktop_config.json`

```json
{
  "mcpServers": {
    "datasets": {
      "command": "datasets-mcp-server"
    }
  }
}
```

#### Claude Code

Add to your project's `.mcp.json` or global `~/.claude/settings.json`:

```json
{
  "mcpServers": {
    "datasets": {
      "command": "datasets-mcp-server"
    }
  }
}
```

### 3. Install Skills and Agents (Claude Code only)

The skills (`/gene-report`, `/variant-report`, etc.) and agents (`/literature-agent`, `/clinical-variant-agent`, etc.) are Claude Code slash commands defined as `.md` files in `.claude/commands/`.

**Option A — Work inside the cloned repo** (simplest):

Skills are automatically available when your working directory is this repository.

**Option B — Install globally** (available in any project):

```bash
# Copy all skills/agents to your user-level commands directory
mkdir -p ~/.claude/commands
cp .claude/commands/*.md ~/.claude/commands/
```

**Option C — Symlink** (stays up-to-date with git pulls):

```bash
mkdir -p ~/.claude/commands
for f in /path/to/datasets-mcp/.claude/commands/*.md; do
  ln -sf "$f" ~/.claude/commands/
done
```

> **Note:** The MCP server (tools) works with any MCP-compatible client. The skills and agents are specific to Claude Code.

## Usage Examples

The examples below show the kind of natural-language prompts you can give an AI assistant (Claude Desktop, Claude Code, etc.) that has this MCP server configured. The assistant picks the right tool automatically.

---

### NCBI Datasets CLI Tools

#### datasets_summary_genome

> "Get me the genome assembly summary for the human reference genome."

> "What genome assemblies are available for *Drosophila melanogaster*?"

> "Look up genome assembly GCF_000001405.40."

```json
{"taxon": "homo_sapiens", "limit": 3}
{"accession": "GCF_000001405.40"}
```

#### datasets_summary_gene

> "What does NCBI know about the BRCA1 gene?"

> "Get gene info for NCBI gene ID 7157."

```json
{"symbol": "BRCA1", "taxon": "human"}
{"gene_id": "7157"}
```

#### datasets_download_genome

> "Download the GRCh38 human reference genome with GFF3 annotations."

> "Show me what would be downloaded for the mouse reference genome without actually downloading."

```json
{"accession": "GCF_000001405.40", "include": ["genome", "gff3"]}
{"taxon": "mus_musculus", "dry_run": true}
```

#### datasets_download_gene

> "Download the BRCA1 and BRCA2 gene sequences for human, including RNA and protein."

```json
{"symbol": "BRCA1,BRCA2", "taxon": "human", "include": ["gene", "rna", "protein"]}
```

#### datasets_version

> "What version of the NCBI datasets CLI is installed?"

```json
{}
```

#### batch_gene_summary

> "Give me a summary of TP53, BRCA1, and EGFR all at once."

> "Look up NCBI gene IDs 672, 675, and 7157 together."

```json
{"symbols": "TP53,BRCA1,EGFR", "taxon": "human"}
{"gene_ids": "672,675,7157"}
```

---

### BioPython Sequence Analysis Tools

#### sequence_align

> "Align these two protein sequences and tell me how similar they are:
> MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
> MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH"

> "Do a local alignment of these two DNA sequences: ATCGATCGATCG and ATCAATCGATCG."

```json
{
  "sequence1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH",
  "sequence2": "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH",
  "sequence_type": "protein",
  "mode": "global"
}
```

#### sequence_stats

> "What is the GC content of this sequence: ATCGATCGATCGATCG?"

> "Compute the molecular weight and amino acid composition of this protein: MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH."

> "Analyze the first sequence in `/data/my_sequences.fasta`."

```json
{"sequence": "ATCGATCGATCGATCG", "sequence_type": "nucleotide"}
{"sequence": "/data/my_sequences.fasta"}
```

#### parse_fasta

> "List all sequences in `/data/proteins.fasta` that have 'kinase' in the description."

> "How many sequences are in `/data/genome.fasta`? Show me the first 10."

```json
{"file_path": "/data/proteins.fasta", "filter_pattern": "kinase", "limit": 50}
{"file_path": "/data/genome.fasta", "limit": 10}
```

#### sequence_translate

> "Translate this nucleotide sequence in all 6 frames: ATGGATCCAGTGGTCCAGGAGTTT..."

> "What are the open reading frames in this DNA sequence?"

```json
{"sequence": "ATGGATCCAGTGGTCCAGGAGTTTGAGGCCATGCTG", "frames": "all"}
{"sequence": "ATGGATCCAGTGGTCCAGGAGTTT", "frames": "forward", "table": 1}
```

---

### Ensembl REST API Tools

#### ensembl_lookup_gene

> "Look up the TP53 gene in Ensembl."

> "What does Ensembl say about ENSG00000141510?"

```json
{"symbol": "TP53", "species": "homo_sapiens"}
{"id": "ENSG00000141510"}
```

#### ensembl_get_sequence

> "Get the protein sequence for Ensembl gene ENSG00000141510."

> "Retrieve the cDNA sequence for ENST00000269305."

```json
{"id": "ENSG00000141510", "seq_type": "protein"}
{"id": "ENST00000269305", "seq_type": "cdna", "format": "fasta"}
```

#### ensembl_search

> "Search Ensembl for cross-references to BRAF in human."

> "What Ensembl IDs correspond to the gene symbol EGFR?"

```json
{"symbol": "BRAF", "species": "homo_sapiens"}
{"symbol": "EGFR"}
```

#### ensembl_get_variants

> "What known variants are in the BRAF V600 region (chromosome 7:140453136-140453236) in human?"

> "Show me variants on chromosome 17 from position 7661779 to 7687550."

```json
{"region": "7:140453136-140453236", "species": "homo_sapiens", "limit": 50}
{"region": "17:7661779-7687550"}
```

#### ensembl_get_homologs

> "Find mouse orthologs for human TP53 (ENSG00000141510)."

> "What are the paralogs of ENSG00000141510?"

```json
{"id": "ENSG00000141510", "target_species": "mus_musculus", "homology_type": "orthologues"}
{"id": "ENSG00000141510", "homology_type": "paralogues"}
```

---

### UniProt REST API Tools

#### uniprot_search

> "Search UniProt for reviewed BRCA1 entries in human."

> "Find proteins annotated with 'tumor suppressor' in UniProt for organism 9606."

```json
{"query": "BRCA1 AND organism_id:9606", "reviewed": true, "limit": 5}
{"query": "tumor suppressor AND organism_id:9606", "limit": 10}
```

#### uniprot_get_protein

> "Get the full UniProt entry for human P53 (accession P04637) — function, GO terms, localization."

> "What does UniProt say about Q9Y6K9?"

```json
{"accession": "P04637"}
{"accession": "Q9Y6K9"}
```

#### uniprot_get_features

> "Show me the protein domains and active sites for UniProt entry P04637."

> "What post-translational modifications are annotated for P38398?"

```json
{"accession": "P04637", "feature_types": ["Domain", "Active site"]}
{"accession": "P38398", "feature_types": ["Modified residue", "Cross-link"]}
```

---

### ClinVar Tools

#### clinvar_search

> "Search ClinVar for pathogenic variants in BRCA1."

> "What is the clinical significance of NM_007294.4:c.5266dupC?"

> "Find ClinVar entries related to Lynch syndrome."

```json
{"query": "BRCA1 AND pathogenic", "limit": 10}
{"query": "NM_007294.4:c.5266dupC"}
{"query": "Lynch syndrome", "limit": 20}
```

---

### PDB / RCSB Tools

#### pdb_get_structure

> "Get the structure details for PDB entry 1TUP (TP53 bound to DNA)."

> "What is the resolution and experimental method for 4HHB?"

```json
{"pdb_id": "1TUP"}
{"pdb_id": "4HHB"}
```

#### pdb_search

> "Search PDB for structures of the SARS-CoV-2 spike protein."

> "Find crystal structures of insulin receptor."

```json
{"query": "SARS-CoV-2 spike protein", "limit": 5}
{"query": "insulin receptor", "limit": 10}
```

---

### InterPro Tools

#### interpro_get_domains

> "What protein domains does TP53 (P04637) have according to InterPro?"

> "Show me the domain architecture for BRCA1 (P38398)."

```json
{"accession": "P04637"}
{"accession": "P38398"}
```

---

### STRING Tools

#### string_get_interactions

> "What proteins interact with TP53 in human?"

> "Show me the protein-protein interaction network for TP53, MDM2, and CDKN2A."

```json
{"identifiers": "TP53", "species": 9606, "limit": 10}
{"identifiers": "TP53,MDM2,CDKN2A", "species": 9606, "required_score": 700}
```

---

### KEGG Tools

#### kegg_get_pathway

> "Search KEGG for pathways related to apoptosis."

> "Get details for the human p53 signaling pathway (hsa04115)."

```json
{"query": "apoptosis"}
{"query": "hsa04115"}
```

---

### NCBI BLAST Tools

#### blast_search

> "BLAST this protein sequence against SwissProt: MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS..."

> "Search for nucleotide sequences similar to ATCGATCGATCGATCG in the nt database."

```json
{"query": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPS", "program": "blastp", "database": "swissprot", "max_hits": 5}
{"query": "ATCGATCGATCGATCG", "program": "blastn", "database": "nt", "max_hits": 10}
```

---

### PubMed Tools

#### pubmed_search

> "Find recent papers about CRISPR-Cas9 gene therapy."

> "Search PubMed for reviews on BRCA1 and PARP inhibitor resistance."

```json
{"query": "CRISPR Cas9 gene therapy", "max_results": 10}
{"query": "BRCA1 PARP inhibitor resistance review", "max_results": 5}
```

---

### Ensembl Regulation Tools

#### ensembl_get_regulation

> "What regulatory features are in the region around TP53 on chromosome 17?"

> "Show me promoters and enhancers near BRAF (chromosome 7:140719327-140924929)."

```json
{"region": "17:7661779-7687550", "species": "homo_sapiens"}
{"region": "7:140719327-140924929", "species": "homo_sapiens"}
```

---

### AlphaFold Tools

#### alphafold_get_prediction

> "Get the AlphaFold predicted structure for human TP53 (P04637)."

> "What is the confidence score for the AlphaFold model of BRCA1 (P38398)?"

```json
{"uniprot_accession": "P04637"}
{"uniprot_accession": "P38398"}
```

---

### gnomAD Tools

#### gnomad_get_variant

> "What is the population frequency of the BRAF V600E variant (7-140753336-A-T) in gnomAD?"

> "Look up rs11549407 in gnomAD."

```json
{"variant_id": "7-140753336-A-T"}
{"variant_id": "rs11549407"}
```

---

### GTEx Tools

#### gtex_get_expression

> "In which tissues is TP53 most highly expressed?"

> "Show me the tissue expression profile for BRCA1."

```json
{"gene": "TP53"}
{"gene": "BRCA1", "dataset": "gtex_v8"}
```

---

### HPO Tools

#### hpo_search

> "Search HPO for phenotypes related to seizures."

> "What genes are associated with the HPO term HP:0001250 (Seizure)?"

> "What diseases are associated with HP:0001250?"

```json
{"query": "seizure"}
{"query": "HP:0001250", "category": "genes", "limit": 20}
{"query": "HP:0001250", "category": "diseases"}
```

---

### Genome Coordinate Tools

#### liftover_coordinates

> "Convert the BRCA1 region 17:41196312-41277500 from hg19 to hg38."

> "Lift over coordinates 7:140453136-140624564 from GRCh37 to GRCh38."

```json
{"region": "17:41196312-41277500", "source_assembly": "GRCh37", "target_assembly": "GRCh38"}
{"region": "7:140453136-140624564"}
```

### Reactome Pathway Tools

#### reactome_get_pathway

> "Look up the Reactome pathway R-HSA-1640170 and show me the participating genes."

> "Search Reactome for pathways involving TP53."

> "What Reactome pathways are related to apoptosis in humans?"

```json
{"pathway_id": "R-HSA-1640170", "include_participants": true, "include_hierarchy": true}
{"query": "TP53", "species": "Homo sapiens", "limit": 5}
{"query": "apoptosis"}
```

### Gene Ontology Enrichment

#### gene_ontology_enrich

> "Run GO enrichment analysis on this gene list: TP53, BRCA1, EGFR, KRAS, MYC, RB1, APC, PTEN"

> "What biological processes are enriched in these genes: CDK2, CDK4, CDK6, CCND1, CCNE1, RB1?"

> "Do GO enrichment including KEGG and Reactome pathways for: IL6, TNF, IL1B, CXCL8, CCL2"

```json
{"genes": "TP53,BRCA1,EGFR,KRAS,MYC,RB1,APC,PTEN"}
{"genes": "CDK2,CDK4,CDK6,CCND1,CCNE1,RB1", "sources": "GO:BP,GO:CC,KEGG,REAC"}
{"genes": "IL6,TNF,IL1B,CXCL8,CCL2", "sources": "GO:BP,GO:MF,GO:CC,KEGG,REAC,WP", "threshold": 0.01}
```

### COSMIC Somatic Mutations

#### cosmic_search

> "Search COSMIC for BRAF V600E mutations in cancer."

> "What somatic mutations in TP53 are found in lung cancer?"

> "Find COSMIC mutations in KRAS."

```json
{"query": "V600E", "gene": "BRAF"}
{"query": "TP53", "limit": 30}
{"query": "KRAS"}
```

### OMIM Gene-Disease Associations

#### omim_search

> "What genetic diseases are associated with BRCA1?"

> "Search OMIM for Marfan syndrome."

> "Look up OMIM entries linked to NCBI Gene ID 7157 (TP53)."

```json
{"query": "BRCA1"}
{"query": "Marfan syndrome"}
{"gene_id": "7157"}
```

### Expression Atlas Tools

#### expression_atlas_search

> "Search Expression Atlas for RNA-seq experiments in human brain tissue."

> "What expression experiments are available for *Mus musculus*?"

```json
{"keyword": "brain", "species": "Homo sapiens", "experiment_type": "rnaseq"}
{"species": "Mus musculus", "limit": 10}
```

### NCBI Gene Links Tools

#### ncbi_gene_links

> "What genes are neighbors of TP53 in the genome?"

> "Find genes related to BRCA1 via NCBI gene links."

```json
{"gene_symbol": "TP53", "species": "human"}
{"gene_id": "672", "limit": 10}
```

### ENCODE Tools

#### encode_get_experiments

> "Search ENCODE for ChIP-seq experiments targeting TP53."

> "Find ATAC-seq datasets in human cell lines from ENCODE."

```json
{"query": "TP53", "assay": "ChIP-seq", "organism": "Homo sapiens", "limit": 10}
{"query": "ATAC-seq", "organism": "Homo sapiens", "limit": 5}
```

### Primer Design Tools

#### primer_design

> "Design PCR primers for this sequence: ATGGATCCAGTGGTCCAGGAGTTTGAGGCCATGCTG..."

> "Design primers targeting positions 100-200 of this template."

```json
{"sequence": "ATGGATCCAGTGGTCCAGGAGTTTGAGGCCATGCTG...", "target_start": 100, "target_length": 100}
{"sequence": "ATGGATCCAGTGGTCCAGGAG...", "target_start": 50, "target_length": 80, "num_primers": 5}
```

### PharmGKB Tools

#### pharmgkb_search

> "What pharmacogenomics annotations exist for CYP2D6?"

> "Search PharmGKB for warfarin-related clinical annotations."

```json
{"gene": "CYP2D6"}
{"drug": "warfarin"}
```

### Disease Ontology Tools

#### disease_ontology_search

> "Search the Disease Ontology for breast cancer."

> "Look up DOID:1612 in the Disease Ontology."

```json
{"query": "breast cancer"}
{"query": "DOID:1612"}
```

---

### Paper Retrieval / Citation Tools

#### paper_fetch

> "Fetch the full text of this paper: 33057194"

> "Get the Methods and Results sections from DOI 10.1038/s41586-020-2649-2."

```json
{"identifier": "33057194", "sections": ["methods", "results"]}
{"identifier": "10.1038/s41586-020-2649-2", "sections": ["methods", "results", "discussion"]}
```

#### semantic_scholar_search

> "Look up citation data for DOI 10.1038/s41586-020-2649-2 on Semantic Scholar."

> "Search Semantic Scholar for papers about transformer models in protein structure prediction."

```json
{"query": "DOI:10.1038/s41586-020-2649-2"}
{"query": "transformer protein structure prediction", "limit": 5, "fields_of_study": "Biology"}
```

---

### Lab Notebook

**Query:** "Add a note to my lab notebook that I'm switching from TP53 to BRCA1 analysis"

Uses `lab_notebook_annotate`:
```json
{"note": "Switching focus from TP53 to BRCA1 — TP53 variant was benign"}
```

**Query:** "Record that the ClinVar results were inconclusive"
```json
{"note": "ClinVar search returned VUS only — no clear pathogenic variants for this region"}
```

---

## Skills

Skills are slash commands that orchestrate multiple tools to produce comprehensive reports. Use them in Claude Code by typing the command.

### `/gene-report <gene>`

Generates a comprehensive multi-database gene report. Queries NCBI, Ensembl, UniProt, InterPro, ClinVar, PDB, STRING, and KEGG — then synthesizes a structured markdown report covering identity, function, domains, clinical variants, structures, interactions, and pathways.

> `/gene-report TP53`

> `/gene-report BRCA1`

> `/gene-report ENSG00000141510`

**Example output sections:**
- Gene Identity (NCBI + Ensembl coordinates, transcripts)
- Protein Function (UniProt function, GO terms, localization)
- Domain Architecture (InterPro/Pfam domains with positions)
- Clinical Significance (ClinVar pathogenic variants, conditions)
- 3D Structures (PDB entries with resolution)
- Protein Interactions (STRING high-confidence partners)
- Pathways (KEGG pathway memberships)

### `/variant-report <target>`

Generates a variant annotation report. Accepts a gene symbol, genomic region, or specific variant identifier. Queries Ensembl for known variants, ClinVar for clinical significance, and UniProt for protein domain context.

> `/variant-report BRCA1`

> `/variant-report 7:140453136-140453236`

> `/variant-report NM_007294.4:c.5266dupC`

**Example output sections:**
- Gene/Region Context (coordinates, protein length)
- Variant Landscape (total variants, consequence type breakdown)
- Clinical Significance table (pathogenic variants with conditions and review status)
- Protein Domain Context (where variants fall relative to functional domains)
- Interpretation Notes (hotspots, enriched domains)

### `/protein-report <protein>`

Generates a protein-centric report. Accepts a UniProt accession (e.g., `P04637`), gene symbol, or protein name. Queries UniProt for function and features, InterPro for domains, AlphaFold for predicted structure, PDB for experimental structures, STRING for interactions, GTEx for tissue expression, ClinVar for pathogenic variants, and PubMed for literature.

> `/protein-report P04637`

> `/protein-report TP53`

> `/protein-report BRCA1`

**Example output sections:**
- Protein Identity (accession, name, organism, length)
- Function (description, GO terms, localization)
- Domain Architecture (InterPro/Pfam domains with positions)
- Post-Translational Modifications (phosphorylation, ubiquitination, disulfide bonds, etc.)
- Active Sites & Binding (catalytic residues, metal/cofactor binding)
- 3D Structures (PDB experimental + AlphaFold prediction with confidence)
- Protein Interactions (STRING high-confidence partners)
- Tissue Expression (GTEx top tissues)
- Clinical Variants (pathogenic variants, domain hotspots)

### `/pathway-report <pathway>`

Generates a pathway deep-dive report. Accepts a Reactome pathway ID (e.g., `R-HSA-1640170`), pathway name, or gene name. Queries Reactome for pathway details and participants, batch gene summaries for member genes, STRING for protein interactions, ClinVar for clinical variants in pathway members, KEGG for cross-reference, and PubMed for literature.

> `/pathway-report R-HSA-1640170`

> `/pathway-report apoptosis`

> `/pathway-report TP53`

**Example output sections:**
- Pathway Overview (name, description, hierarchy, compartments)
- Sub-pathways & Reactions
- Key Member Genes (table with function summaries)
- Protein Interaction Network (interactions between pathway members)
- Clinical Variants in Pathway Members (pathogenic variants, conditions)
- KEGG Cross-Reference
- Recent Literature (PubMed publications)
- Key Takeaways

### `/drug-target-report <target>`

Generates a druggability and drug-target assessment. Accepts a gene symbol (e.g., `EGFR`), UniProt accession, or protein name. Queries UniProt for function and target class, PDB for ligand-bound structures and binding sites, AlphaFold for structural confidence, OMIM/ClinVar/HPO for disease rationale, COSMIC for somatic mutation hotspots, STRING/KEGG/Reactome for pathway context, and PubMed for drug development literature.

> `/drug-target-report EGFR`

> `/drug-target-report BRAF`

> `/drug-target-report P00533`

**Example output sections:**
- Executive Summary (druggability verdict)
- Target Identity (class, localization, enzyme activity)
- Disease Rationale (OMIM, COSMIC hotspots, ClinVar, HPO)
- Structural Druggability (ligand-bound PDBs, binding sites, AlphaFold confidence)
- Interaction Network & Alternative Targets
- Literature — Drug Development Status
- Key Takeaways

### `/lab-notebook <subcommand>`

Research session lab notebook that passively records tool calls via a PostToolUse hook, then synthesizes them into a structured narrative report.

**Subcommands:**
- `start <title>` — Start a new notebook session with a title/objective
- `annotate <note>` — Add a manual context annotation
- `update [note]` — Generate a narrative summary of recent activity
- `report` — Produce a polished final lab notebook as markdown
- `status` — Show session stats (tool count, duration, tools used)

**Setup:** Add the PostToolUse hook to your Claude Code settings:
```json
{
  "hooks": {
    "PostToolUse": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "python3 scripts/lab_notebook_logger.py"
          }
        ]
      }
    ]
  }
}
```

**Example workflow:**
> "Start a lab notebook for investigating BRAF V600E in melanoma"
> ... do research using various tools ...
> "Add a note that the gnomAD frequency was surprisingly high"
> "Generate my lab notebook report"

---

## Agents

Agents are slash commands that perform autonomous, multi-step research tasks by orchestrating several tools.

### `/literature-agent <topic>`

Searches PubMed for relevant literature on a gene, variant, disease, or research topic. Gathers gene/variant context, searches PubMed (with refined queries if needed), deduplicates, and produces a structured literature summary with thematic synthesis and citations.

> `/literature-agent BRCA1 PARP inhibitor resistance`

> `/literature-agent TP53 gain-of-function mutations`

> `/literature-agent CRISPR gene therapy`

**Example output sections:**
- Background (gene/variant context from NCBI, UniProt, ClinVar)
- Key Findings from the Literature (thematic synthesis with inline citations)
- Selected References table (PMID, title, authors, journal, year)
- Abstract Highlights (detailed summaries of top 3-5 papers)

### `/comparative-genomics-agent <gene> <species...>`

Compares a gene across multiple species by finding orthologs, retrieving protein sequences, and computing pairwise alignments. Reports on conservation, divergence, and functional domain context.

> `/comparative-genomics-agent TP53 human mouse zebrafish`

> `/comparative-genomics-agent BRCA1 human mouse rat chicken`

> `/comparative-genomics-agent ENSG00000141510 across vertebrates`

**Example output sections:**
- Gene Overview (name, function, reference species)
- Ortholog Summary table (species, Ensembl ID, protein length, % identity)
- Pairwise Alignments (identity, gaps, key observations)
- Conservation Analysis (conserved vs divergent regions)
- Functional Domain Context (InterPro domains vs conservation)
- Evolutionary Insights (functional constraints, species-specific adaptations)

### `/clinical-variant-agent <variant>`

Full clinical variant workup. Accepts an rsID, chrom-pos-ref-alt, HGVS notation, or gene+variant shorthand. Queries gnomAD for population frequency, ClinVar for clinical significance, UniProt/InterPro for domain impact, AlphaFold for structural context, HPO for phenotype associations, and PubMed for literature.

> `/clinical-variant-agent BRAF V600E`

> `/clinical-variant-agent rs11549407`

> `/clinical-variant-agent 7-140753336-A-T`

**Example output sections:**
- Variant Identity (coordinates, rsID, HGVS, consequence type)
- Population Frequency (gnomAD overall + per-population table)
- Clinical Significance (ClinVar classification, conditions, review status)
- Protein Impact Assessment (domain context, structural confidence)
- Phenotype Associations (HPO terms, related diseases)
- Relevant Literature (top publications)
- ACMG Classification Considerations (evidence summary)

### `/gene-list-agent <genes>`

Functional analysis of a gene list. Accepts comma-separated gene symbols. Runs batch summaries, shared pathway analysis (KEGG), protein interaction network (STRING), tissue expression patterns (GTEx), functional annotation (UniProt GO terms), clinical relevance (ClinVar), and phenotype associations (HPO).

> `/gene-list-agent TP53,BRCA1,EGFR,PTEN,RB1`

> `/gene-list-agent MDM2,CDK4,CCND1,RB1,E2F1`

**Example output sections:**
- Gene Summary table (name, chromosome, ClinVar variants)
- Shared Pathway Analysis (KEGG pathways containing 2+ input genes)
- Protein Interaction Network (STRING direct interactions)
- Tissue Expression Patterns (shared high-expression tissues)
- Functional Themes (recurring GO terms)
- Phenotype Associations (shared HPO terms)
- Interpretation (biological themes uniting the gene set)

### `/structure-agent <target>`

Structural biology deep-dive. Accepts a gene symbol (e.g., `TP53`), UniProt accession (e.g., `P04637`), or PDB ID (e.g., `1TUP`). Retrieves all PDB structures, AlphaFold prediction, domain architecture, pathogenic variant mapping onto structure, protein interaction interfaces, and structural literature.

> `/structure-agent TP53`

> `/structure-agent P04637`

> `/structure-agent 1TUP`

**Example output sections:**
- Domain Architecture (InterPro domains with positions)
- Structure Inventory (all PDB entries: method, resolution, ligands)
- Key Structures Highlighted (best overall, ligand-bound, complex)
- Structure Coverage Map (which protein regions are experimentally resolved)
- AlphaFold Prediction (confidence, disordered regions)
- Pathogenic Variants on Structure (hotspot mapping)
- Protein Interaction Interfaces
- Druggability Assessment (ligand-binding sites, pockets)
- Structural Literature

### `/statistical-methods-agent <paper>`

Analyzes every statistical method used in a biomedical research paper. Fetches the paper full text (via PMC/Europe PMC), gathers citation context from Semantic Scholar, identifies all statistical tests and models, and produces a structured critique covering assumptions, appropriateness, and interpretation. Saves the report as a markdown file.

> `/statistical-methods-agent 33057194`

> `/statistical-methods-agent 10.1038/s41586-020-2649-2`

> `/statistical-methods-agent PMC7505768`

**Example output sections:**
- Paper Metadata (title, authors, journal, citations, fields of study)
- Paper Summary (what the researchers did and found)
- Field Significance (citation context, novelty, current relevance)
- Statistical Methods Inventory (per-test table: test name, data applied to, assumptions, result, plain-language meaning)
- Statistical Critique (appropriateness, assumption validity, multiple comparisons, effect sizes, missing analyses)
- Limitations of This Analysis

---

## Development

### Project Structure

```
datasets-mcp/
├── pyproject.toml
├── README.md
├── ROADMAP.md
└── src/
    └── datasets_mcp/
        ├── __init__.py
        └── server.py
```

### Testing

Test with the MCP Inspector:

```bash
npx @modelcontextprotocol/inspector datasets-mcp-server
```

Or test CLI tools directly:

```bash
datasets version
datasets summary genome taxon human --limit 1
```

## Troubleshooting

### "datasets CLI tool is not installed"

This only affects NCBI tools. Ensembl, UniProt, and BioPython tools work without it.

```bash
which datasets
# Should output the path to the datasets executable
```

If not found, install via conda or download from NCBI.

### Command timeouts

- NCBI CLI commands have a 60-second timeout
- REST API calls (Ensembl, UniProt) have a 30-second timeout
- Use `dry_run: true` to preview downloads
- REST responses are cached for 5 minutes to avoid redundant calls

### Rate limiting

Ensembl and UniProt enforce rate limits. The server handles HTTP 429 responses automatically by waiting and retrying once. If you see persistent rate-limit errors, wait a minute before retrying.

## License

This MCP server is provided as-is. NCBI Datasets, Ensembl, and UniProt each have their own terms of use.

## Contributing

Contributions welcome! See [ROADMAP.md](ROADMAP.md) for planned features. Feel free to open issues or submit pull requests.
