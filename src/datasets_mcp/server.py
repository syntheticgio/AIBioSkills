#!/usr/bin/env python3
"""
MCP Server for NCBI Datasets CLI Tool

This server exposes the NCBI datasets command-line tool functionality
through the Model Context Protocol, allowing AI assistants to search,
summarize, and download genomic and gene data.
"""

import asyncio
import io
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import httpx

from Bio.Align import PairwiseAligner
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data.CodonTable import standard_dna_table

from mcp.server import Server
from mcp.types import Tool, TextContent, ImageContent, EmbeddedResource

# Initialize MCP server
app = Server("datasets-mcp-server")

# REST API constants
ENSEMBL_BASE_URL = "https://rest.ensembl.org"
UNIPROT_BASE_URL = "https://rest.uniprot.org"
EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PDB_DATA_URL = "https://data.rcsb.org/rest/v1/core"
PDB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
INTERPRO_BASE_URL = "https://www.ebi.ac.uk/interpro/api"
STRING_BASE_URL = "https://string-db.org/api"
KEGG_BASE_URL = "https://rest.kegg.jp"
BLAST_BASE_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
ALPHAFOLD_BASE_URL = "https://alphafold.ebi.ac.uk/api"
GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
GTEX_BASE_URL = "https://gtexportal.org/api/v2"
HPO_BASE_URL = "https://hpo.jax.org/api/hpo"
REACTOME_BASE_URL = "https://reactome.org/ContentService"
GPROFILER_BASE_URL = "https://biit.cs.ut.ee/gprofiler/api"
COSMIC_API_URL = "https://clinicaltables.nlm.nih.gov/api/cosmic/v3/search"
EXPRESSION_ATLAS_URL = "https://www.ebi.ac.uk/gxa/json/experiments"
ENCODE_BASE_URL = "https://www.encodeproject.org"
PHARMGKB_BASE_URL = "https://api.pharmgkb.org/v1"
DISEASE_ONTOLOGY_BASE_URL = "https://api.disease-ontology.org/v1"
SEMANTIC_SCHOLAR_BASE_URL = "https://api.semanticscholar.org/graph/v1"
EUROPE_PMC_BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest"
UNPAYWALL_BASE_URL = "https://api.unpaywall.org/v2"

# Simple in-memory cache: {url: (timestamp, data)}
_rest_cache: dict[str, tuple[float, Any]] = {}
_REST_CACHE_TTL = 300  # 5 minutes


async def _rest_get(
    url: str,
    params: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    timeout: float = 30,
    use_cache: bool = True,
) -> dict[str, Any]:
    """Perform an HTTP GET with caching, rate-limit handling, and error wrapping."""
    cache_key = url + (json.dumps(params, sort_keys=True) if params else "")

    if use_cache and cache_key in _rest_cache:
        ts, cached_data = _rest_cache[cache_key]
        if time.time() - ts < _REST_CACHE_TTL:
            return {"success": True, "data": cached_data}
        del _rest_cache[cache_key]

    default_headers = {"Accept": "application/json"}
    if headers:
        default_headers.update(headers)

    try:
        async with httpx.AsyncClient(timeout=timeout) as client:
            resp = await client.get(url, params=params, headers=default_headers)

            if resp.status_code == 429:
                retry_after = float(resp.headers.get("Retry-After", "2"))
                await asyncio.sleep(retry_after)
                resp = await client.get(url, params=params, headers=default_headers)

            if resp.status_code >= 400:
                return {
                    "success": False,
                    "data": None,
                    "error": f"HTTP {resp.status_code}: {resp.text[:500]}",
                }

            data = resp.json()
            if use_cache:
                _rest_cache[cache_key] = (time.time(), data)
            return {"success": True, "data": data}

    except httpx.TimeoutException:
        return {"success": False, "data": None, "error": "Request timed out"}
    except Exception as e:
        return {"success": False, "data": None, "error": str(e)}


def check_datasets_installed() -> bool:
    """Check if the datasets CLI is installed and accessible."""
    return shutil.which("datasets") is not None


async def run_datasets_command(args: list[str]) -> dict[str, Any]:
    """
    Execute a datasets CLI command and return the result.

    Args:
        args: Command arguments to pass to datasets CLI

    Returns:
        Dictionary with success status, output, and any errors
    """
    try:
        cmd = ["datasets"] + args
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )

        return {
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "returncode": result.returncode
        }
    except subprocess.TimeoutExpired:
        return {
            "success": False,
            "stdout": "",
            "stderr": "Command timed out after 60 seconds",
            "returncode": -1
        }
    except Exception as e:
        return {
            "success": False,
            "stdout": "",
            "stderr": str(e),
            "returncode": -1
        }


@app.list_tools()
async def list_tools() -> list[Tool]:
    """List available tools."""
    return [
        Tool(
            name="datasets_summary_genome",
            description="Get summary information for genome assemblies by accession, taxon, or bioproject. Returns metadata including assembly stats, organism info, and annotation details.",
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "Genome assembly accession (e.g., GCF_000001405.40)"
                    },
                    "taxon": {
                        "type": "string",
                        "description": "Taxon name or NCBI Taxonomy ID (e.g., 'human' or '9606')"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results to return",
                        "default": 10
                    }
                },
            }
        ),
        Tool(
            name="datasets_summary_gene",
            description="Get summary information for genes by gene ID, symbol, or taxon. Returns gene metadata, genomic coordinates, and annotation details.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_id": {
                        "type": "string",
                        "description": "NCBI Gene ID (e.g., '672')"
                    },
                    "symbol": {
                        "type": "string",
                        "description": "Gene symbol (e.g., 'BRCA1')"
                    },
                    "taxon": {
                        "type": "string",
                        "description": "Taxon name or ID to filter results (e.g., 'human' or '9606')"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results to return",
                        "default": 10
                    }
                },
            }
        ),
        Tool(
            name="datasets_download_genome",
            description="Download genome assembly data including sequences, annotations, and metadata. Use --dry-run to see what would be downloaded without actually downloading.",
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "Genome assembly accession (e.g., GCF_000001405.40)"
                    },
                    "taxon": {
                        "type": "string",
                        "description": "Taxon name or NCBI Taxonomy ID"
                    },
                    "include": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Data types to include: genome, rna, protein, cds, gff3, gtf, gbff, seq-report",
                        "default": ["genome"]
                    },
                    "filename": {
                        "type": "string",
                        "description": "Output filename (default: ncbi_dataset.zip)",
                        "default": "ncbi_dataset.zip"
                    },
                    "dry_run": {
                        "type": "boolean",
                        "description": "Show what would be downloaded without downloading",
                        "default": False
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="datasets_download_gene",
            description="Download gene data including sequences and annotations. Use --dry-run to preview the download.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_id": {
                        "type": "string",
                        "description": "NCBI Gene ID or comma-separated list of IDs"
                    },
                    "symbol": {
                        "type": "string",
                        "description": "Gene symbol or comma-separated list of symbols"
                    },
                    "taxon": {
                        "type": "string",
                        "description": "Taxon name or ID to filter results"
                    },
                    "include": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Data types to include: gene, rna, protein, cds",
                        "default": ["gene"]
                    },
                    "filename": {
                        "type": "string",
                        "description": "Output filename (default: ncbi_dataset.zip)",
                        "default": "ncbi_dataset.zip"
                    },
                    "dry_run": {
                        "type": "boolean",
                        "description": "Show what would be downloaded without downloading",
                        "default": False
                    }
                },
                "required": []
            }
        ),
        Tool(
            name="datasets_version",
            description="Get the version of the datasets CLI tool",
            inputSchema={
                "type": "object",
                "properties": {},
            }
        ),
        Tool(
            name="sequence_align",
            description="Perform pairwise sequence alignment using BioPython. Returns aligned sequences, identity percentage, gap count, and alignment score. Useful for comparing protein or nucleotide sequences accurately.",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence1": {
                        "type": "string",
                        "description": "First sequence (raw sequence or FASTA format)"
                    },
                    "sequence2": {
                        "type": "string",
                        "description": "Second sequence (raw sequence or FASTA format)"
                    },
                    "sequence_type": {
                        "type": "string",
                        "description": "Type of sequence: 'protein' or 'nucleotide' (default: auto-detect)",
                        "enum": ["protein", "nucleotide"]
                    },
                    "mode": {
                        "type": "string",
                        "description": "Alignment mode: 'global' or 'local' (default: 'global')",
                        "enum": ["global", "local"],
                        "default": "global"
                    }
                },
                "required": ["sequence1", "sequence2"]
            }
        ),
        Tool(
            name="sequence_stats",
            description="Compute statistics for a nucleotide or protein sequence. Returns length, GC content (nucleotide), codon usage (nucleotide), amino acid composition (protein), and molecular weight (protein). Accepts raw sequence or file path to a FASTA file.",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "Raw sequence string OR path to a FASTA file"
                    },
                    "sequence_type": {
                        "type": "string",
                        "description": "Type of sequence: 'protein' or 'nucleotide' (default: auto-detect)",
                        "enum": ["protein", "nucleotide"]
                    }
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="parse_fasta",
            description="Parse a FASTA file and return sequence IDs, descriptions, lengths, and optionally full sequences. Supports filtering by regex pattern on IDs/descriptions.",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the FASTA file"
                    },
                    "filter_pattern": {
                        "type": "string",
                        "description": "Regex pattern to filter sequence IDs or descriptions"
                    },
                    "include_sequences": {
                        "type": "boolean",
                        "description": "Include full sequences in output (default: false)",
                        "default": False
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of sequences to return (default: 100)",
                        "default": 100
                    }
                },
                "required": ["file_path"]
            }
        ),
        Tool(
            name="sequence_translate",
            description="Translate a nucleotide sequence to protein in all 6 reading frames (3 forward + 3 reverse complement). Returns amino acid sequences with stop codons marked as '*'. Useful for identifying open reading frames.",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "Nucleotide sequence (raw or FASTA format)"
                    },
                    "frames": {
                        "type": "string",
                        "description": "Which frames to translate: 'all' (default, all 6), 'forward' (frames +1,+2,+3), 'reverse' (frames -1,-2,-3), or a specific frame like '+1', '-2'",
                        "default": "all"
                    },
                    "table": {
                        "type": "integer",
                        "description": "NCBI genetic code table number (default: 1 = Standard). Use 2 for Vertebrate Mitochondrial, 11 for Bacterial/Plant Plastid, etc.",
                        "default": 1
                    },
                },
                "required": ["sequence"]
            }
        ),
        Tool(
            name="batch_gene_summary",
            description="Query multiple genes at once by gene IDs or symbols. Returns consolidated summary results for all genes in a single call, avoiding multiple round-trips.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_ids": {
                        "type": "string",
                        "description": "Comma-separated NCBI Gene IDs (e.g., '672,675,7157')"
                    },
                    "symbols": {
                        "type": "string",
                        "description": "Comma-separated gene symbols (e.g., 'BRCA1,BRCA2,TP53')"
                    },
                    "taxon": {
                        "type": "string",
                        "description": "Taxon name or ID to filter results"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum results per gene (default: 10)",
                        "default": 10
                    }
                },
            }
        ),

        # ── Section 3: Ensembl REST API Tools ──────────────────────────
        Tool(
            name="ensembl_lookup_gene",
            description="Look up a gene via Ensembl by stable ID or by symbol+species. Returns genomic coordinates, biotype, description, and cross-references.",
            inputSchema={
                "type": "object",
                "properties": {
                    "id": {
                        "type": "string",
                        "description": "Ensembl stable ID (e.g., 'ENSG00000141510')"
                    },
                    "symbol": {
                        "type": "string",
                        "description": "Gene symbol (e.g., 'TP53'). Requires 'species'."
                    },
                    "species": {
                        "type": "string",
                        "description": "Species name (e.g., 'homo_sapiens'). Required when using 'symbol'.",
                        "default": "homo_sapiens"
                    },
                },
            }
        ),
        Tool(
            name="ensembl_get_sequence",
            description="Retrieve a nucleotide or protein sequence from Ensembl by stable ID. Supports genomic, cDNA, CDS, and protein sequence types.",
            inputSchema={
                "type": "object",
                "properties": {
                    "id": {
                        "type": "string",
                        "description": "Ensembl stable ID (gene, transcript, or protein)"
                    },
                    "seq_type": {
                        "type": "string",
                        "description": "Sequence type to retrieve",
                        "enum": ["genomic", "cdna", "cds", "protein"],
                        "default": "cdna"
                    },
                    "format": {
                        "type": "string",
                        "description": "Output format",
                        "enum": ["json", "fasta"],
                        "default": "json"
                    },
                },
                "required": ["id"]
            }
        ),
        Tool(
            name="ensembl_search",
            description="Search Ensembl cross-references for a gene symbol in a given species. Returns matching Ensembl IDs and external database links.",
            inputSchema={
                "type": "object",
                "properties": {
                    "symbol": {
                        "type": "string",
                        "description": "Gene symbol to search (e.g., 'BRCA1')"
                    },
                    "species": {
                        "type": "string",
                        "description": "Species name (e.g., 'homo_sapiens')",
                        "default": "homo_sapiens"
                    },
                },
                "required": ["symbol"]
            }
        ),
        Tool(
            name="ensembl_get_variants",
            description="Get known genetic variants in a genomic region from Ensembl. Returns variant IDs, alleles, consequences, and clinical significance.",
            inputSchema={
                "type": "object",
                "properties": {
                    "region": {
                        "type": "string",
                        "description": "Genomic region as chr:start-end (e.g., '7:140453136-140624564')"
                    },
                    "species": {
                        "type": "string",
                        "description": "Species name (e.g., 'homo_sapiens')",
                        "default": "homo_sapiens"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of variants to return (default: 100)",
                        "default": 100
                    },
                },
                "required": ["region"]
            }
        ),
        Tool(
            name="ensembl_get_homologs",
            description="Find homologous genes (orthologs/paralogs) for an Ensembl gene ID. Useful for cross-species gene comparisons and evolutionary analysis.",
            inputSchema={
                "type": "object",
                "properties": {
                    "id": {
                        "type": "string",
                        "description": "Ensembl gene ID (e.g., 'ENSG00000141510')"
                    },
                    "target_species": {
                        "type": "string",
                        "description": "Target species for orthologs (e.g., 'mus_musculus'). Omit for all species."
                    },
                    "homology_type": {
                        "type": "string",
                        "description": "Filter by homology type",
                        "enum": ["orthologues", "paralogues", "all"],
                        "default": "all"
                    },
                },
                "required": ["id"]
            }
        ),

        # ── Section 4: UniProt REST API Tools ──────────────────────────
        Tool(
            name="uniprot_search",
            description="Search UniProt for proteins using Lucene query syntax. Supports filtering by organism, gene name, protein name, GO terms, and more.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "UniProt search query in Lucene syntax (e.g., 'BRCA1 AND organism_id:9606')"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results (default: 10)",
                        "default": 10
                    },
                    "reviewed": {
                        "type": "boolean",
                        "description": "If true, only return Swiss-Prot (reviewed) entries"
                    },
                },
                "required": ["query"]
            }
        ),
        Tool(
            name="uniprot_get_protein",
            description="Get detailed protein information from UniProt by accession. Returns function, subcellular location, GO terms, and key annotations.",
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "UniProt accession (e.g., 'P04637' for human TP53)"
                    },
                },
                "required": ["accession"]
            }
        ),
        Tool(
            name="uniprot_get_features",
            description="Get protein sequence features (domains, active sites, modifications, variants) from UniProt by accession.",
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "UniProt accession (e.g., 'P04637')"
                    },
                    "feature_types": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Filter by feature types (e.g., ['Domain', 'Active site', 'Modified residue']). Omit for all features."
                    },
                },
                "required": ["accession"]
            }
        ),

        # ── Section 5: ClinVar / NCBI E-utilities Tools ───────────────
        Tool(
            name="clinvar_search",
            description="Search NCBI ClinVar for clinical variant interpretations. Query by gene, variant, disease, or condition. Returns variant accessions, clinical significance, conditions, and review status.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "ClinVar search query (e.g., 'BRCA1', 'BRCA1 AND pathogenic', 'NM_007294.4:c.5266dupC', 'breast cancer')"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results (default: 20)",
                        "default": 20
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 6: PDB / RCSB REST API Tools ──────────────────────
        Tool(
            name="pdb_get_structure",
            description="Get detailed metadata for a PDB structure by its 4-character PDB ID. Returns title, experimental method, resolution, organism, deposition date, and polymer entity details.",
            inputSchema={
                "type": "object",
                "properties": {
                    "pdb_id": {
                        "type": "string",
                        "description": "4-character PDB ID (e.g., '1TUP', '4HHB', '6LU7')"
                    },
                },
                "required": ["pdb_id"]
            }
        ),
        Tool(
            name="pdb_search",
            description="Search the RCSB Protein Data Bank for 3D structures by free-text query. Supports searching by protein name, gene name, organism, ligand, or any keyword.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search text (e.g., 'TP53', 'insulin receptor', 'SARS-CoV-2 spike protein')"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results (default: 10)",
                        "default": 10
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 7: InterPro REST API Tools ────────────────────────
        Tool(
            name="interpro_get_domains",
            description="Get protein domain and family annotations from InterPro by UniProt accession. Returns domain architecture including Pfam, PROSITE, CDD, and other member database matches.",
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "UniProt accession (e.g., 'P04637' for human TP53)"
                    },
                },
                "required": ["accession"]
            }
        ),

        # ── Section 8: STRING REST API Tools ──────────────────────────
        Tool(
            name="string_get_interactions",
            description="Get protein-protein interaction partners from STRING database. Returns interaction scores, evidence channels, and network neighbors for one or more proteins.",
            inputSchema={
                "type": "object",
                "properties": {
                    "identifiers": {
                        "type": "string",
                        "description": "Protein name(s), separated by newline or comma (e.g., 'TP53' or 'TP53,MDM2,CDKN2A')"
                    },
                    "species": {
                        "type": "integer",
                        "description": "NCBI taxonomy ID (e.g., 9606 for human, 10090 for mouse)",
                        "default": 9606
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max interaction partners per query protein (default: 10)",
                        "default": 10
                    },
                    "required_score": {
                        "type": "integer",
                        "description": "Minimum combined score (0-1000, default: 400 = medium confidence)",
                        "default": 400
                    },
                    "network_type": {
                        "type": "string",
                        "description": "Type of interactions: 'functional' (default) or 'physical'",
                        "enum": ["functional", "physical"],
                        "default": "functional"
                    },
                },
                "required": ["identifiers"]
            }
        ),

        # ── Section 9: KEGG REST API Tools ────────────────────────────
        Tool(
            name="kegg_get_pathway",
            description="Search and retrieve KEGG pathway information. Can search for pathways by keyword or retrieve details for a specific pathway ID. Returns pathway names, descriptions, genes, and linked entries.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Keyword search (e.g., 'apoptosis', 'insulin signaling') or a KEGG pathway ID (e.g., 'hsa04210', 'map00010')"
                    },
                    "organism": {
                        "type": "string",
                        "description": "KEGG organism code (e.g., 'hsa' for human, 'mmu' for mouse). Used for keyword search.",
                        "default": "hsa"
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 10: NCBI BLAST REST API Tools ─────────────────────
        Tool(
            name="blast_search",
            description="Submit a sequence similarity search to NCBI BLAST and retrieve results. Supports blastn (DNA), blastp (protein), blastx, tblastn, tblastx. Note: BLAST searches are asynchronous and may take 30 seconds to several minutes.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Query sequence in FASTA format or as a raw sequence string"
                    },
                    "program": {
                        "type": "string",
                        "description": "BLAST program to use",
                        "enum": ["blastn", "blastp", "blastx", "tblastn", "tblastx"],
                        "default": "blastp"
                    },
                    "database": {
                        "type": "string",
                        "description": "Target database (e.g., 'nr', 'swissprot', 'nt', 'refseq_protein', 'pdb')",
                        "default": "nr"
                    },
                    "evalue": {
                        "type": "number",
                        "description": "E-value threshold (default: 0.01)",
                        "default": 0.01
                    },
                    "max_hits": {
                        "type": "integer",
                        "description": "Maximum number of hits to return (default: 10)",
                        "default": 10
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 11: PubMed / NCBI E-utilities Tools ───────────────
        Tool(
            name="pubmed_search",
            description="Search PubMed for biomedical literature. Returns article titles, authors, journals, dates, abstracts, DOIs, and PMIDs. Supports full PubMed query syntax.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "PubMed search query (e.g., 'BRCA1 AND PARP inhibitor', 'TP53 cancer therapy 2024')"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results (default: 10)",
                        "default": 10
                    },
                    "sort": {
                        "type": "string",
                        "description": "Sort order: 'relevance' (default), 'date', or 'author'",
                        "enum": ["relevance", "date", "author"],
                        "default": "relevance"
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 12: Ensembl Regulation Tools ──────────────────────
        Tool(
            name="ensembl_get_regulation",
            description="Get regulatory features (promoters, enhancers, CTCF binding sites, open chromatin, TF binding sites) in a genomic region from Ensembl. Useful for understanding non-coding regulatory elements near a gene.",
            inputSchema={
                "type": "object",
                "properties": {
                    "region": {
                        "type": "string",
                        "description": "Genomic region as chr:start-end (e.g., '17:7661779-7687550')"
                    },
                    "species": {
                        "type": "string",
                        "description": "Species name (e.g., 'homo_sapiens')",
                        "default": "homo_sapiens"
                    },
                },
                "required": ["region"]
            }
        ),

        # ── Section 13: AlphaFold REST API Tools ──────────────────────
        Tool(
            name="alphafold_get_prediction",
            description="Get AlphaFold predicted 3D structure information for a protein by UniProt accession. Returns confidence scores (pLDDT), model URLs, gene info, and sequence coverage.",
            inputSchema={
                "type": "object",
                "properties": {
                    "accession": {
                        "type": "string",
                        "description": "UniProt accession (e.g., 'P04637' for human TP53)"
                    },
                },
                "required": ["accession"]
            }
        ),

        # ── Section 14: gnomAD GraphQL API Tools ─────────────────────────
        Tool(
            name="gnomad_get_variant",
            description="Get population allele frequencies for a genetic variant from gnomAD. Accepts a variant ID (chrom-pos-ref-alt format like '1-55516888-G-A') or rsID (like 'rs11549407'). Returns allele counts, frequencies, and population breakdowns from exome and genome data.",
            inputSchema={
                "type": "object",
                "properties": {
                    "variant_id": {
                        "type": "string",
                        "description": "Variant ID in chrom-pos-ref-alt format (e.g., '1-55516888-G-A') or rsID (e.g., 'rs11549407')"
                    },
                    "dataset": {
                        "type": "string",
                        "description": "gnomAD dataset version (default: 'gnomad_r4')",
                        "default": "gnomad_r4"
                    },
                },
                "required": ["variant_id"]
            }
        ),

        # ── Section 15: GTEx REST API Tools ──────────────────────────────
        Tool(
            name="gtex_get_expression",
            description="Get tissue-specific gene expression data from GTEx. Returns median TPM expression values across ~54 human tissues. Accepts a gene symbol or Ensembl gene ID.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol (e.g., 'TP53') or Ensembl gene ID (e.g., 'ENSG00000141510')"
                    },
                    "dataset": {
                        "type": "string",
                        "description": "GTEx dataset version: 'gtex_v8' (default) or 'gtex_v10'",
                        "default": "gtex_v8"
                    },
                },
                "required": ["gene"]
            }
        ),

        # ── Section 16: HPO REST API Tools ───────────────────────────────
        Tool(
            name="hpo_search",
            description="Search the Human Phenotype Ontology (HPO) for phenotype terms, or get genes/diseases associated with a phenotype. Use for gene-phenotype mapping, clinical phenotyping, and rare disease research.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query — a keyword (e.g., 'seizure'), HPO term ID (e.g., 'HP:0001250'), or gene symbol (e.g., 'BRCA1')"
                    },
                    "category": {
                        "type": "string",
                        "description": "What to retrieve: 'search' (search terms/genes/diseases), 'genes' (genes for an HPO term), 'diseases' (diseases for an HPO term), 'term' (details of an HPO term). Default: 'search'",
                        "default": "search"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum results to return (default: 20)",
                        "default": 20
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 17: Genome Coordinate Tools ─────────────────────────
        Tool(
            name="liftover_coordinates",
            description="Convert genomic coordinates between genome assemblies (e.g., GRCh37/hg19 to GRCh38/hg38 or vice versa) using the Ensembl coordinate mapping API. Essential when working with data from different assembly versions.",
            inputSchema={
                "type": "object",
                "properties": {
                    "region": {
                        "type": "string",
                        "description": "Genomic region as chr:start-end or chr:start..end (e.g., '17:43044295-43170245' or 'X:1000000-1001000')"
                    },
                    "source_assembly": {
                        "type": "string",
                        "description": "Source assembly: 'GRCh37' (hg19) or 'GRCh38' (hg38). Default: 'GRCh37'",
                        "default": "GRCh37"
                    },
                    "target_assembly": {
                        "type": "string",
                        "description": "Target assembly: 'GRCh37' (hg19) or 'GRCh38' (hg38). Default: 'GRCh38'",
                        "default": "GRCh38"
                    },
                    "species": {
                        "type": "string",
                        "description": "Species name (default: 'human')",
                        "default": "human"
                    },
                },
                "required": ["region"]
            }
        ),

        # ── Section 18: Reactome Pathway Tools ───────────────────────────
        Tool(
            name="reactome_get_pathway",
            description="Get detailed pathway information from Reactome — the curated knowledgebase of biological pathways and reactions. Supports lookup by Reactome stable ID (e.g., R-HSA-1640170), keyword/gene search, and retrieval of pathway participants (genes/proteins). Returns pathway details, hierarchy, and participating molecules.",
            inputSchema={
                "type": "object",
                "properties": {
                    "pathway_id": {
                        "type": "string",
                        "description": "Reactome stable ID (e.g., 'R-HSA-1640170'). Use for direct pathway lookup."
                    },
                    "query": {
                        "type": "string",
                        "description": "Search term (gene name, keyword, or pathway name) to search for pathways. Used when pathway_id is not provided."
                    },
                    "species": {
                        "type": "string",
                        "description": "Species name (default: 'Homo sapiens'). Used with query searches.",
                        "default": "Homo sapiens"
                    },
                    "include_participants": {
                        "type": "boolean",
                        "description": "If true and pathway_id is provided, also retrieve participating molecules (genes/proteins). Default: false.",
                        "default": False
                    },
                    "include_hierarchy": {
                        "type": "boolean",
                        "description": "If true and pathway_id is provided, also retrieve ancestor pathways in the hierarchy. Default: false.",
                        "default": False
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of search results to return (default: 10)",
                        "default": 10
                    },
                },
                "required": []
            }
        ),

        # ── Section 19: Gene Ontology Enrichment ─────────────────────────
        Tool(
            name="gene_ontology_enrich",
            description="Perform Gene Ontology (GO) enrichment analysis on a gene list using g:Profiler. Returns statistically enriched GO terms (Biological Process, Molecular Function, Cellular Component) with p-values. Optionally includes KEGG, Reactome, and WikiPathways enrichment. Useful for interpreting gene lists from differential expression, screens, or pathway analyses.",
            inputSchema={
                "type": "object",
                "properties": {
                    "genes": {
                        "type": "string",
                        "description": "Comma-separated list of gene symbols or IDs (e.g., 'TP53,BRCA1,EGFR,KRAS,MYC')"
                    },
                    "organism": {
                        "type": "string",
                        "description": "Organism name in g:Profiler format (default: 'hsapiens'). Common: 'mmusculus', 'rnorvegicus', 'dmelanogaster', 'scerevisiae', 'drerio'",
                        "default": "hsapiens"
                    },
                    "sources": {
                        "type": "string",
                        "description": "Comma-separated data sources to query (default: 'GO:BP,GO:MF,GO:CC'). Options: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP, HP, CORUM, TF, MIRNA",
                        "default": "GO:BP,GO:MF,GO:CC"
                    },
                    "threshold": {
                        "type": "number",
                        "description": "Significance threshold (default: 0.05)",
                        "default": 0.05
                    },
                    "correction_method": {
                        "type": "string",
                        "description": "Multiple testing correction: 'g_SCS' (default, g:Profiler's method), 'fdr', or 'bonferroni'",
                        "default": "g_SCS"
                    },
                    "no_iea": {
                        "type": "boolean",
                        "description": "Exclude electronically inferred GO annotations (IEA). Default: false.",
                        "default": False
                    },
                    "ordered": {
                        "type": "boolean",
                        "description": "If true, treats the gene list as ordered (ranked). Default: false.",
                        "default": False
                    },
                    "background": {
                        "type": "string",
                        "description": "Comma-separated list of background genes (custom statistical background). If not provided, uses all known genes for the organism."
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of enriched terms to return (default: 25)",
                        "default": 25
                    },
                },
                "required": ["genes"]
            }
        ),

        # ── Section 20: COSMIC Somatic Mutations ─────────────────────────
        Tool(
            name="cosmic_search",
            description="Search the COSMIC (Catalogue Of Somatic Mutations In Cancer) database for somatic mutations associated with cancer. Returns mutation IDs, gene names, amino acid changes, genomic positions, primary histology, and tissue site. Data sourced via the NLM Clinical Tables API (COSMIC V89, GRCh37).",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search term — gene name (e.g., 'BRAF'), mutation (e.g., 'V600E'), or keyword"
                    },
                    "gene": {
                        "type": "string",
                        "description": "Filter results to a specific gene (e.g., 'BRAF'). Combined with query if both provided."
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of mutations to return (default: 20, max: 500)",
                        "default": 20
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 21: OMIM Gene-Disease Associations ───────────────────
        Tool(
            name="omim_search",
            description="Search OMIM (Online Mendelian Inheritance in Man) for genetic disease associations. Uses NCBI E-utilities to find OMIM entries linked to a gene or search term, returning MIM numbers, titles, and entry types (gene vs phenotype). For gene input, also retrieves linked disease phenotypes via gene→OMIM cross-references.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search term — gene symbol (e.g., 'BRCA1'), disease name (e.g., 'Marfan syndrome'), or OMIM MIM number (e.g., '191170')"
                    },
                    "gene_id": {
                        "type": "string",
                        "description": "NCBI Gene ID (e.g., '7157' for TP53). If provided, retrieves all OMIM entries linked to this gene."
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of OMIM entries to return (default: 20)",
                        "default": 20
                    },
                },
                "required": []
            }
        ),

        # ── Section 22: Expression Atlas ─────────────────────────────────
        Tool(
            name="expression_atlas_search",
            description="Search EBI Expression Atlas for gene expression experiments — baseline and differential RNA-seq and proteomics datasets across species and conditions. Returns experiment accessions, descriptions, species, assay counts, and experimental factors. Useful for finding public expression datasets relevant to a gene, tissue, disease, or condition.",
            inputSchema={
                "type": "object",
                "properties": {
                    "species": {
                        "type": "string",
                        "description": "Filter by species (e.g., 'homo sapiens', 'mus musculus'). Default: no filter (all species)."
                    },
                    "experiment_type": {
                        "type": "string",
                        "description": "Filter by type: 'baseline' or 'differential'. Default: no filter.",
                        "enum": ["baseline", "differential"]
                    },
                    "keyword": {
                        "type": "string",
                        "description": "Keyword to filter experiments by description (client-side filter). E.g., 'brain', 'cancer', 'immune'."
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of experiments to return (default: 15)",
                        "default": 15
                    },
                },
                "required": []
            }
        ),

        # ── Section 23: NCBI Gene Links ──────────────────────────────────
        Tool(
            name="ncbi_gene_links",
            description="Get gene-gene relationships from NCBI — finds computationally determined gene neighbors (genes with similar sequences, shared protein domains, or co-expression patterns). Uses NCBI E-utilities elink with the gene_gene_neighbors link type.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene_id": {
                        "type": "string",
                        "description": "NCBI Gene ID (e.g., '7157' for TP53)"
                    },
                    "gene_symbol": {
                        "type": "string",
                        "description": "Gene symbol (e.g., 'TP53'). Will be resolved to Gene ID via esearch."
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of related genes to return (default: 15)",
                        "default": 15
                    },
                },
                "required": []
            }
        ),

        # ── Section 24: ENCODE Experiments ───────────────────────────────
        Tool(
            name="encode_get_experiments",
            description="Search the ENCODE Project for functional genomics experiments — ChIP-seq, ATAC-seq, RNA-seq, WGBS, Hi-C, and more. Returns experiment accessions, assay types, biosample information, targets, and file counts. Useful for finding epigenomic and regulatory datasets for a gene or cell type.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search term — gene/target name (e.g., 'TP53'), cell type (e.g., 'K562'), or keyword"
                    },
                    "assay": {
                        "type": "string",
                        "description": "Filter by assay type (e.g., 'ChIP-seq', 'ATAC-seq', 'RNA-seq', 'WGBS', 'Hi-C')"
                    },
                    "organism": {
                        "type": "string",
                        "description": "Filter by organism: 'Homo sapiens' or 'Mus musculus'. Default: no filter."
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of experiments to return (default: 10)",
                        "default": 10
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 25: Primer Design ────────────────────────────────────
        Tool(
            name="primer_design",
            description="Design PCR primers for a target DNA sequence using Primer3 via BioPython. Accepts a nucleotide sequence (or FASTA) and returns primer pairs with Tm, GC%, length, and product size. Useful for cloning, genotyping, and qPCR primer design.",
            inputSchema={
                "type": "object",
                "properties": {
                    "sequence": {
                        "type": "string",
                        "description": "Target nucleotide sequence (raw or FASTA format). Must be at least 100 bp."
                    },
                    "target_start": {
                        "type": "integer",
                        "description": "Start position (0-based) of the region to amplify within the sequence. Default: center of sequence."
                    },
                    "target_length": {
                        "type": "integer",
                        "description": "Length of the target region to amplify (default: 200 bp)",
                        "default": 200
                    },
                    "primer_min_size": {
                        "type": "integer",
                        "description": "Minimum primer length (default: 18)",
                        "default": 18
                    },
                    "primer_opt_size": {
                        "type": "integer",
                        "description": "Optimal primer length (default: 20)",
                        "default": 20
                    },
                    "primer_max_size": {
                        "type": "integer",
                        "description": "Maximum primer length (default: 25)",
                        "default": 25
                    },
                    "primer_min_tm": {
                        "type": "number",
                        "description": "Minimum melting temperature in °C (default: 57.0)",
                        "default": 57.0
                    },
                    "primer_opt_tm": {
                        "type": "number",
                        "description": "Optimal melting temperature in °C (default: 60.0)",
                        "default": 60.0
                    },
                    "primer_max_tm": {
                        "type": "number",
                        "description": "Maximum melting temperature in °C (default: 63.0)",
                        "default": 63.0
                    },
                    "num_primers": {
                        "type": "integer",
                        "description": "Number of primer pairs to return (default: 5)",
                        "default": 5
                    },
                },
                "required": ["sequence"]
            }
        ),

        # ── Section 26: PharmGKB Pharmacogenomics ────────────────────────
        Tool(
            name="pharmgkb_search",
            description="Search PharmGKB for pharmacogenomics annotations — gene-drug-variant associations, clinical annotations, and dosing guidelines. Returns clinical annotations linking genetic variants to drug response (efficacy, toxicity, dosage). Essential for precision medicine and pharmacogenomics research.",
            inputSchema={
                "type": "object",
                "properties": {
                    "gene": {
                        "type": "string",
                        "description": "Gene symbol to search for clinical annotations (e.g., 'CYP2D6', 'VKORC1', 'EGFR')"
                    },
                    "drug": {
                        "type": "string",
                        "description": "Drug name to search for (e.g., 'warfarin', 'tamoxifen', 'imatinib')"
                    },
                    "variant": {
                        "type": "string",
                        "description": "Variant rsID to search for (e.g., 'rs1057910')"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of annotations to return (default: 15)",
                        "default": 15
                    },
                },
                "required": []
            }
        ),

        # ── Section 27: Disease Ontology ─────────────────────────────────
        Tool(
            name="disease_ontology_search",
            description="Search the Disease Ontology (DO) for standardized disease terms, definitions, synonyms, and cross-references to other ontologies (OMIM, MeSH, ICD, SNOMED). Useful for finding canonical disease identifiers and understanding disease classification hierarchy.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search term — disease name (e.g., 'breast cancer'), DOID (e.g., 'DOID:1612'), or keyword"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results to return (default: 10)",
                        "default": 10
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Section 28: Paper Retrieval / Citation Tools ─────────────────
        Tool(
            name="paper_fetch",
            description="Retrieve a biomedical research paper's metadata and full text (when available via PubMed Central or Europe PMC). Accepts a PMID, DOI, PMC ID, or article title. Extracts specific sections (Methods, Results, Discussion) from open-access full-text XML. Falls back to abstract if full text is unavailable.",
            inputSchema={
                "type": "object",
                "properties": {
                    "identifier": {
                        "type": "string",
                        "description": "PMID (e.g., '33057194'), DOI (e.g., '10.1038/s41586-020-2649-2'), PMC ID (e.g., 'PMC7505768'), or article title"
                    },
                    "sections": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Sections to extract from full text: 'methods', 'results', 'discussion', 'introduction', 'conclusions', 'all' (default: ['methods', 'results'])"
                    },
                },
                "required": ["identifier"]
            }
        ),
        Tool(
            name="semantic_scholar_search",
            description="Search Semantic Scholar for paper metadata, citation counts, influential citations, TLDRs, and open-access PDF links. Accepts a paper title, DOI (prefix with 'DOI:'), PMID (prefix with 'PMID:'), or keyword search.",
            inputSchema={
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Paper title, 'DOI:10.1038/...' for DOI lookup, 'PMID:33057194' for PMID lookup, or keyword search"
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Max results for keyword search (default: 5, ignored for DOI/PMID lookup)",
                        "default": 5
                    },
                    "fields_of_study": {
                        "type": "string",
                        "description": "Filter by field: 'Medicine', 'Biology', 'Computer Science', etc."
                    },
                },
                "required": ["query"]
            }
        ),

        # ── Lab Notebook ──────────
        Tool(
            name="lab_notebook_annotate",
            description="Add a manual annotation to the lab notebook session log. Use this to record context, decisions, or observations that tool calls alone don't capture — e.g., 'switching to investigate BRCA1 instead' or 'this approach didn't work because...'.",
            inputSchema={
                "type": "object",
                "properties": {
                    "note": {
                        "type": "string",
                        "description": "The annotation text to add to the lab notebook log"
                    },
                    "session_id": {
                        "type": "string",
                        "description": "Session ID to annotate. If omitted, uses the most recent session."
                    },
                },
                "required": ["note"]
            }
        ),
    ]


def _extract_sequence(text: str) -> str:
    """Extract a raw sequence from text that may be FASTA-formatted."""
    text = text.strip()
    if text.startswith(">"):
        records = list(SeqIO.parse(io.StringIO(text), "fasta"))
        if records:
            return str(records[0].seq)
    return re.sub(r"\s+", "", text).upper()


def _detect_sequence_type(seq: str) -> str:
    """Auto-detect whether a sequence is protein or nucleotide."""
    nuc_chars = set("ATCGURYMKSWHBVDN")
    upper = seq.upper()
    nuc_count = sum(1 for c in upper if c in nuc_chars)
    if len(seq) == 0:
        return "nucleotide"
    if nuc_count / len(seq) > 0.85:
        return "nucleotide"
    return "protein"


def _count_codons(seq: str) -> dict[str, int]:
    """Count codon usage for a nucleotide sequence."""
    seq = seq.upper()
    codons: dict[str, int] = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        if len(codon) == 3:
            codons[codon] = codons.get(codon, 0) + 1
    return dict(sorted(codons.items()))


def _extract_pmc_sections(xml_text: str, requested_sections: list[str]) -> dict[str, Any]:
    """Parse PMC full-text XML and extract requested sections.

    Args:
        xml_text: Full-text XML from PMC/Europe PMC.
        requested_sections: List of section names to extract. Use 'all' for everything.
            Recognized: 'methods', 'results', 'discussion', 'introduction', 'conclusions'.

    Returns:
        Dict with 'sections' (name→text mapping) and 'section_titles_found' (all titles).
    """
    import xml.etree.ElementTree as ET

    MAX_SECTION_CHARS = 15000

    # Map requested names to patterns for matching section titles
    section_patterns: dict[str, list[str]] = {
        "methods": ["method", "materials and method", "experimental", "statistical analysis",
                     "study design", "patients and method", "experimental procedure"],
        "results": ["result", "findings"],
        "discussion": ["discussion"],
        "introduction": ["introduction", "background"],
        "conclusions": ["conclusion", "summary"],
    }

    want_all = "all" in requested_sections
    wanted = set(r.lower() for r in requested_sections)

    def _elem_text(elem) -> str:
        """Recursively extract text from an XML element."""
        parts: list[str] = []
        if elem.text:
            parts.append(elem.text)
        for child in elem:
            if child.tag in ("p", "title", "th", "td", "label", "caption"):
                parts.append(_elem_text(child))
            elif child.tag == "table-wrap":
                parts.append("[Table] " + _elem_text(child))
            elif child.tag == "fig":
                parts.append("[Figure] " + _elem_text(child))
            elif child.tag in ("xref", "ext-link", "italic", "bold", "sup", "sub"):
                parts.append(_elem_text(child))
            if child.tail:
                parts.append(child.tail)
        return " ".join(parts)

    def _match_section(title_text: str) -> str | None:
        """Return canonical section name if title matches any pattern."""
        lower = title_text.lower().strip()
        for canon, patterns in section_patterns.items():
            for pat in patterns:
                if pat in lower:
                    return canon
        return None

    sections: dict[str, str] = {}
    section_titles_found: list[str] = []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return {"sections": {}, "section_titles_found": []}

    # Find body sections — try <body><sec> structure
    body = root.find(".//body")
    if body is None:
        return {"sections": {}, "section_titles_found": []}

    for sec in body.findall(".//sec"):
        title_elem = sec.find("title")
        if title_elem is None:
            continue
        title_text = _elem_text(title_elem).strip()
        if not title_text:
            continue
        section_titles_found.append(title_text)

        canon = _match_section(title_text)
        if canon is None and not want_all:
            continue
        if not want_all and canon not in wanted:
            continue

        key = canon if canon else title_text
        text = _elem_text(sec)
        # Truncate
        if len(text) > MAX_SECTION_CHARS:
            text = text[:MAX_SECTION_CHARS] + "\n... [truncated]"
        if key in sections:
            sections[key] += "\n\n" + text
        else:
            sections[key] = text

    return {"sections": sections, "section_titles_found": section_titles_found}


# Tools that don't require the datasets CLI
BIOPYTHON_TOOLS = {"sequence_align", "sequence_stats", "parse_fasta", "sequence_translate", "primer_design", "lab_notebook_annotate"}
REST_API_TOOLS = {
    "ensembl_lookup_gene", "ensembl_get_sequence", "ensembl_search",
    "ensembl_get_variants", "ensembl_get_homologs",
    "uniprot_search", "uniprot_get_protein", "uniprot_get_features",
    "clinvar_search", "pdb_get_structure", "pdb_search",
    "interpro_get_domains", "string_get_interactions",
    "kegg_get_pathway", "blast_search",
    "pubmed_search", "ensembl_get_regulation", "alphafold_get_prediction",
    "gnomad_get_variant", "gtex_get_expression", "hpo_search",
    "liftover_coordinates",
    "reactome_get_pathway",
    "gene_ontology_enrich",
    "cosmic_search",
    "omim_search",
    "expression_atlas_search",
    "ncbi_gene_links",
    "encode_get_experiments",
    "pharmgkb_search",
    "disease_ontology_search",
    "paper_fetch",
    "semantic_scholar_search",
}
NON_CLI_TOOLS = BIOPYTHON_TOOLS | REST_API_TOOLS


@app.call_tool()
async def call_tool(name: str, arguments: Any) -> list[TextContent]:
    """Handle tool calls."""

    # Only check datasets CLI for tools that need it
    if name not in NON_CLI_TOOLS and not check_datasets_installed():
        return [TextContent(
            type="text",
            text="Error: datasets CLI tool is not installed or not in PATH. Please install it from https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
        )]

    try:
        if name == "datasets_version":
            result = await run_datasets_command(["version"])
            if result["success"]:
                return [TextContent(type="text", text=result["stdout"])]
            else:
                return [TextContent(type="text", text=f"Error: {result['stderr']}")]

        elif name == "datasets_summary_genome":
            args = ["summary", "genome"]

            if arguments.get("accession"):
                args.extend(["accession", arguments["accession"]])
            elif arguments.get("taxon"):
                args.extend(["taxon", arguments["taxon"]])
            else:
                return [TextContent(type="text", text="Error: Either 'accession' or 'taxon' must be provided")]

            if arguments.get("limit"):
                args.extend(["--limit", str(arguments["limit"])])

            result = await run_datasets_command(args)

            if result["success"]:
                try:
                    data = json.loads(result["stdout"])
                    formatted = json.dumps(data, indent=2)
                    return [TextContent(type="text", text=formatted)]
                except json.JSONDecodeError:
                    return [TextContent(type="text", text=result["stdout"])]
            else:
                return [TextContent(type="text", text=f"Error: {result['stderr']}")]

        elif name == "datasets_summary_gene":
            args = ["summary", "gene"]

            if arguments.get("gene_id"):
                args.extend(["gene-id", arguments["gene_id"]])
            elif arguments.get("symbol"):
                args.extend(["symbol", arguments["symbol"]])
            else:
                return [TextContent(type="text", text="Error: Either 'gene_id' or 'symbol' must be provided")]

            if arguments.get("taxon"):
                args.extend(["--taxon", arguments["taxon"]])

            if arguments.get("limit"):
                args.extend(["--limit", str(arguments["limit"])])

            result = await run_datasets_command(args)

            if result["success"]:
                try:
                    data = json.loads(result["stdout"])
                    formatted = json.dumps(data, indent=2)
                    return [TextContent(type="text", text=formatted)]
                except json.JSONDecodeError:
                    return [TextContent(type="text", text=result["stdout"])]
            else:
                return [TextContent(type="text", text=f"Error: {result['stderr']}")]

        elif name == "datasets_download_genome":
            args = ["download", "genome"]

            if arguments.get("accession"):
                args.extend(["accession", arguments["accession"]])
            elif arguments.get("taxon"):
                args.extend(["taxon", arguments["taxon"]])
            else:
                return [TextContent(type="text", text="Error: Either 'accession' or 'taxon' must be provided")]

            if arguments.get("include"):
                for item in arguments["include"]:
                    args.extend(["--include", item])

            if arguments.get("filename"):
                args.extend(["--filename", arguments["filename"]])

            if arguments.get("dry_run"):
                args.append("--dry-run")

            result = await run_datasets_command(args)

            if result["success"]:
                return [TextContent(type="text", text=f"Success:\n{result['stdout']}")]
            else:
                return [TextContent(type="text", text=f"Error: {result['stderr']}")]

        elif name == "datasets_download_gene":
            args = ["download", "gene"]

            if arguments.get("gene_id"):
                args.extend(["gene-id", arguments["gene_id"]])
            elif arguments.get("symbol"):
                args.extend(["symbol", arguments["symbol"]])
            else:
                return [TextContent(type="text", text="Error: Either 'gene_id' or 'symbol' must be provided")]

            if arguments.get("taxon"):
                args.extend(["--taxon", arguments["taxon"]])

            if arguments.get("include"):
                for item in arguments["include"]:
                    args.extend(["--include", item])

            if arguments.get("filename"):
                args.extend(["--filename", arguments["filename"]])

            if arguments.get("dry_run"):
                args.append("--dry-run")

            result = await run_datasets_command(args)

            if result["success"]:
                return [TextContent(type="text", text=f"Success:\n{result['stdout']}")]
            else:
                return [TextContent(type="text", text=f"Error: {result['stderr']}")]

        elif name == "sequence_align":
            seq1 = _extract_sequence(arguments["sequence1"])
            seq2 = _extract_sequence(arguments["sequence2"])
            seq_type = arguments.get("sequence_type") or _detect_sequence_type(seq1)
            mode = arguments.get("mode", "global")

            aligner = PairwiseAligner()
            aligner.mode = mode

            if seq_type == "protein":
                from Bio.Align import substitution_matrices
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                aligner.open_gap_score = -10
                aligner.extend_gap_score = -0.5
            else:
                aligner.match_score = 2
                aligner.mismatch_score = -1
                aligner.open_gap_score = -5
                aligner.extend_gap_score = -0.5

            alignments = aligner.align(seq1, seq2)
            best = alignments[0]

            aligned_str = str(best)
            lines = aligned_str.strip().split("\n")

            # Calculate identity and gaps from the alignment
            aligned_seq1 = lines[0] if len(lines) >= 1 else ""
            aligned_seq2 = lines[2] if len(lines) >= 3 else ""
            identity = 0
            gaps = 0
            align_len = max(len(aligned_seq1), len(aligned_seq2))
            for i in range(min(len(aligned_seq1), len(aligned_seq2))):
                if aligned_seq1[i] == "-" or aligned_seq2[i] == "-":
                    gaps += 1
                elif aligned_seq1[i] == aligned_seq2[i]:
                    identity += 1

            identity_pct = (identity / align_len * 100) if align_len > 0 else 0

            result = {
                "alignment": aligned_str,
                "score": float(best.score),
                "identity_percent": round(identity_pct, 2),
                "identical_positions": identity,
                "alignment_length": align_len,
                "gaps": gaps,
                "sequence_type": seq_type,
                "mode": mode,
                "seq1_length": len(seq1),
                "seq2_length": len(seq2),
            }
            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "sequence_stats":
            raw = arguments["sequence"]
            seq_type = arguments.get("sequence_type")

            # Check if it's a file path
            if os.path.isfile(raw):
                records = list(SeqIO.parse(raw, "fasta"))
                if not records:
                    return [TextContent(type="text", text="Error: No sequences found in FASTA file")]
                seq_str = str(records[0].seq)
                source = f"file:{raw} (first sequence: {records[0].id})"
            else:
                seq_str = _extract_sequence(raw)
                source = "raw input"

            if not seq_type:
                seq_type = _detect_sequence_type(seq_str)

            result: dict[str, Any] = {
                "source": source,
                "length": len(seq_str),
                "sequence_type": seq_type,
            }

            if seq_type == "nucleotide":
                result["gc_content"] = round(gc_fraction(seq_str) * 100, 2)
                result["base_counts"] = {
                    base: seq_str.upper().count(base) for base in "ATCG"
                }
                result["codon_usage"] = _count_codons(seq_str)
            else:
                try:
                    analysis = ProteinAnalysis(seq_str)
                    result["molecular_weight"] = round(analysis.molecular_weight(), 2)
                    result["amino_acid_percent"] = {
                        k: round(v * 100, 2)
                        for k, v in analysis.amino_acids_percent.items()
                    }
                    result["isoelectric_point"] = round(analysis.isoelectric_point(), 2)
                    result["aromaticity"] = round(analysis.aromaticity(), 4)
                except Exception as e:
                    result["analysis_error"] = str(e)

            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "parse_fasta":
            file_path = arguments["file_path"]
            if not os.path.isfile(file_path):
                return [TextContent(type="text", text=f"Error: File not found: {file_path}")]

            filter_pattern = arguments.get("filter_pattern")
            include_sequences = arguments.get("include_sequences", False)
            limit = arguments.get("limit", 100)

            compiled_re = re.compile(filter_pattern, re.IGNORECASE) if filter_pattern else None

            sequences = []
            count = 0
            total_in_file = 0
            for record in SeqIO.parse(file_path, "fasta"):
                total_in_file += 1
                if compiled_re and not (compiled_re.search(record.id) or compiled_re.search(record.description)):
                    continue
                entry: dict[str, Any] = {
                    "id": record.id,
                    "description": record.description,
                    "length": len(record.seq),
                }
                if include_sequences:
                    entry["sequence"] = str(record.seq)
                sequences.append(entry)
                count += 1
                if count >= limit:
                    break

            result = {
                "file": file_path,
                "total_sequences_in_file": total_in_file if count < limit else f"{total_in_file}+",
                "returned": count,
                "filter": filter_pattern,
                "sequences": sequences,
            }
            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "sequence_translate":
            from Bio.Seq import Seq
            from Bio.Data.CodonTable import TranslationError

            raw = arguments["sequence"]
            frames_arg = arguments.get("frames", "all")
            table_id = arguments.get("table", 1)

            # Extract and clean sequence
            if os.path.isfile(raw):
                records = list(SeqIO.parse(raw, "fasta"))
                if not records:
                    return [TextContent(type="text", text="Error: No sequences found in FASTA file")]
                nuc_str = str(records[0].seq).upper()
                source = f"file:{raw} (first sequence: {records[0].id})"
            else:
                nuc_str = _extract_sequence(raw)
                source = "raw input"

            # Validate it's a nucleotide sequence
            if _detect_sequence_type(nuc_str) != "nucleotide":
                return [TextContent(type="text", text="Error: Input does not appear to be a nucleotide sequence")]

            seq_obj = Seq(nuc_str)
            rc_obj = seq_obj.reverse_complement()

            # Determine which frames to translate
            if frames_arg == "all":
                frame_list = [1, 2, 3, -1, -2, -3]
            elif frames_arg == "forward":
                frame_list = [1, 2, 3]
            elif frames_arg == "reverse":
                frame_list = [-1, -2, -3]
            else:
                # Parse specific frame like "+1", "-2", "1", "3"
                try:
                    f = int(frames_arg.replace("+", ""))
                    if f not in (1, 2, 3, -1, -2, -3):
                        return [TextContent(type="text", text="Error: Frame must be one of: +1, +2, +3, -1, -2, -3")]
                    frame_list = [f]
                except ValueError:
                    return [TextContent(type="text", text=f"Error: Invalid frames argument: {frames_arg}")]

            translations = []
            for frame in frame_list:
                if frame > 0:
                    offset = frame - 1
                    sub = seq_obj[offset:]
                    direction = "forward"
                else:
                    offset = abs(frame) - 1
                    sub = rc_obj[offset:]
                    direction = "reverse"

                # Translate (to_stop=False to get full translation with * for stops)
                try:
                    protein = str(sub.translate(table=table_id))
                except Exception:
                    protein = str(sub.translate())

                # Find ORFs (stretches between stop codons)
                orfs = []
                # Find M...* patterns
                for m in re.finditer(r'M[^*]*\*?', protein):
                    orf_seq = m.group()
                    if len(orf_seq) >= 10:  # Only report ORFs >= 10 aa
                        orfs.append({
                            "start_aa": m.start() + 1,
                            "length_aa": len(orf_seq.rstrip("*")),
                            "has_stop": orf_seq.endswith("*"),
                            "sequence": orf_seq[:50] + ("..." if len(orf_seq) > 50 else ""),
                        })

                frame_label = f"+{frame}" if frame > 0 else str(frame)
                translations.append({
                    "frame": frame_label,
                    "direction": direction,
                    "nucleotide_offset": abs(frame) - 1,
                    "protein_length": len(protein),
                    "stop_codons": protein.count("*"),
                    "protein_sequence": protein[:200] + ("..." if len(protein) > 200 else ""),
                    "orfs_10aa_plus": orfs[:10],
                })

            # Find the longest ORF across all frames
            all_orfs = []
            for t in translations:
                for orf in t["orfs_10aa_plus"]:
                    all_orfs.append({"frame": t["frame"], **orf})
            longest_orf = max(all_orfs, key=lambda x: x["length_aa"]) if all_orfs else None

            result = {
                "source": source,
                "nucleotide_length": len(nuc_str),
                "genetic_code_table": table_id,
                "frames_translated": len(frame_list),
                "translations": translations,
                "longest_orf": longest_orf,
            }
            return [TextContent(type="text", text=json.dumps(result, indent=2))]

        elif name == "batch_gene_summary":
            if not arguments.get("gene_ids") and not arguments.get("symbols"):
                return [TextContent(type="text", text="Error: Either 'gene_ids' or 'symbols' must be provided")]

            args = ["summary", "gene"]

            if arguments.get("gene_ids"):
                args.extend(["gene-id", arguments["gene_ids"]])
            elif arguments.get("symbols"):
                args.extend(["symbol", arguments["symbols"]])

            if arguments.get("taxon"):
                args.extend(["--taxon", arguments["taxon"]])

            if arguments.get("limit"):
                args.extend(["--limit", str(arguments["limit"])])

            result = await run_datasets_command(args)

            if result["success"]:
                try:
                    data = json.loads(result["stdout"])
                    reports = data.get("reports", [])
                    summary = {
                        "total_results": len(reports),
                        "query": {
                            "gene_ids": arguments.get("gene_ids"),
                            "symbols": arguments.get("symbols"),
                            "taxon": arguments.get("taxon"),
                        },
                        "results": data,
                    }
                    return [TextContent(type="text", text=json.dumps(summary, indent=2))]
                except json.JSONDecodeError:
                    return [TextContent(type="text", text=result["stdout"])]
            else:
                return [TextContent(type="text", text=f"Error: {result['stderr']}")]

        # ── Section 3: Ensembl REST API Tools ──────────────────────────

        elif name == "ensembl_lookup_gene":
            eid = arguments.get("id")
            symbol = arguments.get("symbol")
            species = arguments.get("species", "homo_sapiens")

            if eid:
                url = f"{ENSEMBL_BASE_URL}/lookup/id/{eid}"
                result = await _rest_get(url, params={"expand": "1"})
            elif symbol:
                url = f"{ENSEMBL_BASE_URL}/lookup/symbol/{species}/{symbol}"
                result = await _rest_get(url, params={"expand": "1"})
            else:
                return [TextContent(type="text", text="Error: Either 'id' or 'symbol' (with 'species') must be provided")]

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            d = result["data"]
            curated = {
                "id": d.get("id"),
                "display_name": d.get("display_name"),
                "description": d.get("description"),
                "species": d.get("species"),
                "biotype": d.get("biotype"),
                "seq_region_name": d.get("seq_region_name"),
                "start": d.get("start"),
                "end": d.get("end"),
                "strand": d.get("strand"),
                "assembly_name": d.get("assembly_name"),
                "source": d.get("source"),
            }
            transcripts = d.get("Transcript", [])
            if transcripts:
                curated["transcript_count"] = len(transcripts)
                curated["canonical_transcript"] = next(
                    (t.get("id") for t in transcripts if t.get("is_canonical")),
                    transcripts[0].get("id") if transcripts else None,
                )
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        elif name == "ensembl_get_sequence":
            eid = arguments["id"]
            seq_type = arguments.get("seq_type", "cdna")
            fmt = arguments.get("format", "json")

            url = f"{ENSEMBL_BASE_URL}/sequence/id/{eid}"
            params: dict[str, Any] = {"type": seq_type}

            if fmt == "fasta":
                result = await _rest_get(url, params=params, headers={"Accept": "text/x-fasta"})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error: {result['error']}")]
                return [TextContent(type="text", text=str(result["data"]))]

            result = await _rest_get(url, params=params)
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            d = result["data"]
            curated = {
                "id": d.get("id"),
                "molecule": d.get("molecule") or seq_type,
                "length": len(d.get("seq", "")),
                "sequence": d.get("seq", ""),
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        elif name == "ensembl_search":
            symbol = arguments["symbol"]
            species = arguments.get("species", "homo_sapiens")
            url = f"{ENSEMBL_BASE_URL}/xrefs/symbol/{species}/{symbol}"

            result = await _rest_get(url)
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            entries = result["data"]
            curated = [
                {
                    "id": e.get("id"),
                    "type": e.get("type"),
                    "db_display_name": e.get("db_display_name"),
                }
                for e in (entries if isinstance(entries, list) else [])
            ]
            return [TextContent(type="text", text=json.dumps({"symbol": symbol, "species": species, "results": curated}, indent=2))]

        elif name == "ensembl_get_variants":
            region = arguments["region"]
            species = arguments.get("species", "homo_sapiens")
            limit = arguments.get("limit", 100)

            url = f"{ENSEMBL_BASE_URL}/overlap/region/{species}/{region}"
            result = await _rest_get(url, params={"feature": "variation", "content-type": "application/json"})

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            variants = result["data"] if isinstance(result["data"], list) else []
            curated = []
            for v in variants[:limit]:
                curated.append({
                    "id": v.get("id"),
                    "start": v.get("start"),
                    "end": v.get("end"),
                    "alleles": v.get("alleles"),
                    "consequence_type": v.get("consequence_type"),
                    "clinical_significance": v.get("clinical_significance", []),
                    "source": v.get("source"),
                })
            return [TextContent(type="text", text=json.dumps({"region": region, "species": species, "count": len(curated), "variants": curated}, indent=2))]

        elif name == "ensembl_get_homologs":
            eid = arguments["id"]
            target = arguments.get("target_species")
            htype = arguments.get("homology_type", "all")

            url = f"{ENSEMBL_BASE_URL}/homology/id/{eid}"
            params = {}
            if target:
                params["target_species"] = target
            if htype != "all":
                params["type"] = htype

            result = await _rest_get(url, params=params)
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            homologies = (
                result["data"]
                .get("data", [{}])[0]
                .get("homologies", [])
            )
            curated = []
            for h in homologies:
                t = h.get("target", {})
                curated.append({
                    "type": h.get("type"),
                    "target_id": t.get("id"),
                    "target_species": t.get("species"),
                    "target_protein_id": t.get("protein_id"),
                    "percent_identity": t.get("perc_id"),
                    "percent_positives": t.get("perc_pos"),
                })
            return [TextContent(type="text", text=json.dumps({"gene_id": eid, "homolog_count": len(curated), "homologs": curated}, indent=2))]

        # ── Section 4: UniProt REST API Tools ──────────────────────────

        elif name == "uniprot_search":
            query = arguments["query"]
            limit = arguments.get("limit", 10)
            reviewed = arguments.get("reviewed")

            if reviewed is True:
                query = f"({query}) AND reviewed:true"

            url = f"{UNIPROT_BASE_URL}/uniprotkb/search"
            params = {"query": query, "size": str(limit), "format": "json"}
            result = await _rest_get(url, params=params)

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            entries = result["data"].get("results", [])
            curated = []
            for e in entries:
                acc = e.get("primaryAccession", "")
                organism = e.get("organism", {})
                genes = e.get("genes", [{}])
                gene_name = genes[0].get("geneName", {}).get("value", "") if genes else ""
                protein_name = (
                    e.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value", "")
                )
                curated.append({
                    "accession": acc,
                    "gene_name": gene_name,
                    "protein_name": protein_name,
                    "organism": organism.get("scientificName", ""),
                    "reviewed": e.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)",
                    "length": e.get("sequence", {}).get("length"),
                })
            return [TextContent(type="text", text=json.dumps({"query": arguments["query"], "count": len(curated), "results": curated}, indent=2))]

        elif name == "uniprot_get_protein":
            accession = arguments["accession"]
            url = f"{UNIPROT_BASE_URL}/uniprotkb/{accession}.json"
            result = await _rest_get(url)

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            d = result["data"]
            genes = d.get("genes", [{}])
            gene_name = genes[0].get("geneName", {}).get("value", "") if genes else ""
            protein_name = (
                d.get("proteinDescription", {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value", "")
            )

            # Extract function comments
            functions = []
            for comment in d.get("comments", []):
                if comment.get("commentType") == "FUNCTION":
                    for txt in comment.get("texts", []):
                        functions.append(txt.get("value", ""))

            # Extract subcellular locations
            locations = []
            for comment in d.get("comments", []):
                if comment.get("commentType") == "SUBCELLULAR LOCATION":
                    for loc in comment.get("subcellularLocations", []):
                        loc_val = loc.get("location", {}).get("value", "")
                        if loc_val:
                            locations.append(loc_val)

            # Extract GO terms from cross-references
            go_terms = []
            for xref in d.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "GO":
                    props = {p["key"]: p["value"] for p in xref.get("properties", [])}
                    go_terms.append({
                        "id": xref.get("id"),
                        "term": props.get("GoTerm", ""),
                    })

            curated = {
                "accession": accession,
                "gene_name": gene_name,
                "protein_name": protein_name,
                "organism": d.get("organism", {}).get("scientificName", ""),
                "length": d.get("sequence", {}).get("length"),
                "functions": functions,
                "subcellular_locations": locations,
                "go_terms": go_terms[:30],
                "entry_type": d.get("entryType", ""),
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        elif name == "uniprot_get_features":
            accession = arguments["accession"]
            feature_types_filter = arguments.get("feature_types")

            url = f"{UNIPROT_BASE_URL}/uniprotkb/{accession}.json"
            result = await _rest_get(url)

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            features_raw = result["data"].get("features", [])
            curated = []
            for f in features_raw:
                ftype = f.get("type", "")
                if feature_types_filter and ftype not in feature_types_filter:
                    continue
                loc = f.get("location", {})
                start = loc.get("start", {}).get("value")
                end = loc.get("end", {}).get("value")
                curated.append({
                    "type": ftype,
                    "description": f.get("description", ""),
                    "start": start,
                    "end": end,
                    "evidences": len(f.get("evidences", [])),
                })

            return [TextContent(type="text", text=json.dumps({"accession": accession, "feature_count": len(curated), "features": curated}, indent=2))]

        # ── Section 5: ClinVar / NCBI E-utilities Tools ───────────────

        elif name == "clinvar_search":
            query = arguments["query"]
            limit = arguments.get("limit", 20)

            # Step 1: esearch to get variant UIDs
            search_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
            search_result = await _rest_get(search_url, params={
                "db": "clinvar",
                "term": query,
                "retmax": str(limit),
                "retmode": "json",
            })

            if not search_result["success"]:
                return [TextContent(type="text", text=f"Error: {search_result['error']}")]

            esearch = search_result["data"].get("esearchresult", {})
            uid_list = esearch.get("idlist", [])

            if not uid_list:
                return [TextContent(type="text", text=json.dumps({
                    "query": query, "count": 0, "results": [],
                    "message": "No ClinVar records found for this query."
                }, indent=2))]

            # Step 2: esummary to get variant details
            summary_url = f"{EUTILS_BASE_URL}/esummary.fcgi"
            summary_result = await _rest_get(summary_url, params={
                "db": "clinvar",
                "id": ",".join(uid_list),
                "retmode": "json",
            }, use_cache=False)

            if not summary_result["success"]:
                return [TextContent(type="text", text=f"Error: {summary_result['error']}")]

            result_data = summary_result["data"].get("result", {})
            curated = []
            for uid in uid_list:
                rec = result_data.get(uid, {})
                if not rec or uid == "uids":
                    continue

                # Extract clinical significance from description or supporting submissions
                clinical_sig = rec.get("clinical_significance", {})
                germline = clinical_sig.get("description", "") if isinstance(clinical_sig, dict) else str(clinical_sig)

                # Extract variation details
                variation_set = rec.get("variation_set", [])
                variant_info = {}
                if variation_set:
                    vs = variation_set[0]
                    variant_info = {
                        "variation_name": vs.get("variation_name", ""),
                        "variation_id": vs.get("variation_xrefs", [{}])[0].get("db_id", "") if vs.get("variation_xrefs") else "",
                        "cdna_change": vs.get("cdna_change", ""),
                    }

                curated.append({
                    "uid": uid,
                    "title": rec.get("title", ""),
                    "accession": rec.get("accession", ""),
                    "clinical_significance": germline,
                    "review_status": rec.get("clinical_significance", {}).get("review_status", "") if isinstance(rec.get("clinical_significance"), dict) else "",
                    "gene_sort": rec.get("gene_sort", ""),
                    "trait_set": [
                        t.get("trait_name", "")
                        for t in rec.get("trait_set", [])
                    ],
                    "variation": variant_info,
                })

            total_count = int(esearch.get("count", len(curated)))
            return [TextContent(type="text", text=json.dumps({
                "query": query,
                "total_in_clinvar": total_count,
                "returned": len(curated),
                "results": curated,
            }, indent=2))]

        # ── Section 6: PDB / RCSB REST API Tools ──────────────────────

        elif name == "pdb_get_structure":
            pdb_id = arguments["pdb_id"].upper()
            url = f"{PDB_DATA_URL}/entry/{pdb_id}"
            result = await _rest_get(url)

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            d = result["data"]
            struct = d.get("struct", {})
            cell = d.get("cell", {})
            exptl = d.get("exptl", [{}])[0] if d.get("exptl") else {}
            refine = d.get("refine", [{}])[0] if d.get("refine") else {}
            rcsb_entry = d.get("rcsb_entry_info", {})

            curated = {
                "pdb_id": pdb_id,
                "title": struct.get("title", ""),
                "experimental_method": exptl.get("method", ""),
                "resolution_angstrom": refine.get("ls_d_res_high"),
                "deposition_date": rcsb_entry.get("deposit_date", ""),
                "release_date": rcsb_entry.get("initial_release_date", ""),
                "polymer_entity_count": rcsb_entry.get("polymer_entity_count"),
                "molecular_weight_kda": rcsb_entry.get("molecular_weight", 0) / 1000 if rcsb_entry.get("molecular_weight") else None,
            }

            # Get polymer entities for organism/gene info
            entity_url = f"{PDB_DATA_URL}/polymer_entity/{pdb_id}/1"
            entity_result = await _rest_get(entity_url)
            if entity_result["success"]:
                ed = entity_result["data"]
                entity_poly = ed.get("entity_poly", {})
                src = ed.get("rcsb_entity_source_organism", [{}])
                src0 = src[0] if src else {}
                gene_names = ed.get("rcsb_gene_name", [])
                curated["entity_1"] = {
                    "description": ed.get("rcsb_polymer_entity", {}).get("pdbx_description", ""),
                    "type": entity_poly.get("type", ""),
                    "organism": src0.get("ncbi_scientific_name", ""),
                    "gene_names": [g.get("value", "") for g in gene_names],
                    "sequence_length": len(entity_poly.get("pdbx_seq_one_letter_code_can", "")),
                }

            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        elif name == "pdb_search":
            query_text = arguments["query"]
            limit = arguments.get("limit", 10)

            search_body = {
                "return_type": "entry",
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": query_text,
                    },
                },
                "request_options": {
                    "paginate": {
                        "start": 0,
                        "rows": limit,
                    },
                    "results_content_type": ["experimental"],
                },
            }

            try:
                async with httpx.AsyncClient(timeout=30) as client:
                    resp = await client.post(
                        PDB_SEARCH_URL,
                        json=search_body,
                        headers={"Content-Type": "application/json"},
                    )

                    if resp.status_code >= 400:
                        return [TextContent(type="text", text=f"Error: PDB search returned HTTP {resp.status_code}: {resp.text[:500]}")]

                    data = resp.json()
            except httpx.TimeoutException:
                return [TextContent(type="text", text="Error: PDB search timed out")]
            except Exception as e:
                return [TextContent(type="text", text=f"Error: {str(e)}")]

            total = data.get("total_count", 0)
            results = data.get("result_set", [])

            # Fetch brief metadata for each hit
            curated = []
            for hit in results:
                pdb_id = hit.get("identifier", "")
                score = hit.get("score", 0)

                entry_url = f"{PDB_DATA_URL}/entry/{pdb_id}"
                entry_result = await _rest_get(entry_url)
                if entry_result["success"]:
                    ed = entry_result["data"]
                    struct = ed.get("struct", {})
                    exptl = ed.get("exptl", [{}])[0] if ed.get("exptl") else {}
                    refine = ed.get("refine", [{}])[0] if ed.get("refine") else {}
                    curated.append({
                        "pdb_id": pdb_id,
                        "title": struct.get("title", ""),
                        "method": exptl.get("method", ""),
                        "resolution": refine.get("ls_d_res_high"),
                        "score": round(score, 2),
                    })
                else:
                    curated.append({"pdb_id": pdb_id, "score": round(score, 2)})

            return [TextContent(type="text", text=json.dumps({
                "query": query_text,
                "total_structures": total,
                "returned": len(curated),
                "results": curated,
            }, indent=2))]

        # ── Section 7: InterPro REST API Tools ────────────────────────

        elif name == "interpro_get_domains":
            accession = arguments["accession"].upper()
            url = f"{INTERPRO_BASE_URL}/protein/uniprot/{accession}"
            result = await _rest_get(url)

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            d = result["data"]
            metadata = d.get("metadata", {})

            # Get InterPro entry matches for this protein
            entries_url = f"{INTERPRO_BASE_URL}/entry/interpro/protein/uniprot/{accession}"
            entries_result = await _rest_get(entries_url)

            domains = []
            if entries_result["success"]:
                for entry in entries_result["data"].get("results", []):
                    emeta = entry.get("metadata", {})
                    # Extract locations from protein match
                    proteins = entry.get("proteins", [])
                    locations = []
                    for prot in proteins:
                        for loc_frag in prot.get("entry_protein_locations", []):
                            for frag in loc_frag.get("fragments", []):
                                locations.append({
                                    "start": frag.get("start"),
                                    "end": frag.get("end"),
                                })
                    domains.append({
                        "accession": emeta.get("accession", ""),
                        "name": emeta.get("name", ""),
                        "type": emeta.get("type", ""),
                        "source_database": emeta.get("source_database", ""),
                        "description": (emeta.get("description", [None]) or [None])[0],
                        "locations": locations,
                    })

            curated = {
                "protein_accession": accession,
                "protein_name": metadata.get("name", ""),
                "protein_length": metadata.get("length"),
                "organism": metadata.get("source_organism", {}).get("scientificName", ""),
                "domain_count": len(domains),
                "domains": domains,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 8: STRING REST API Tools ──────────────────────────

        elif name == "string_get_interactions":
            raw_ids = arguments["identifiers"]
            # Normalize: accept comma-separated or newline-separated
            identifiers = raw_ids.replace(",", "\n").replace("%0d", "\n")
            species = arguments.get("species", 9606)
            limit = arguments.get("limit", 10)
            required_score = arguments.get("required_score", 400)
            network_type = arguments.get("network_type", "functional")

            url = f"{STRING_BASE_URL}/json/interaction_partners"
            result = await _rest_get(url, params={
                "identifiers": identifiers,
                "species": str(species),
                "limit": str(limit),
                "required_score": str(required_score),
                "network_type": network_type,
                "caller_identity": "datasets-mcp-server",
            })

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            interactions = result["data"] if isinstance(result["data"], list) else []
            curated = []
            for ix in interactions:
                curated.append({
                    "protein_a": ix.get("preferredName_A", ix.get("stringId_A", "")),
                    "protein_b": ix.get("preferredName_B", ix.get("stringId_B", "")),
                    "combined_score": ix.get("score", 0),
                    "experimental_score": ix.get("escore", 0),
                    "database_score": ix.get("dscore", 0),
                    "textmining_score": ix.get("tscore", 0),
                    "coexpression_score": ix.get("ascore", 0),
                })

            return [TextContent(type="text", text=json.dumps({
                "query": raw_ids,
                "species": species,
                "interaction_count": len(curated),
                "interactions": curated,
            }, indent=2))]

        # ── Section 9: KEGG REST API Tools ────────────────────────────

        elif name == "kegg_get_pathway":
            query = arguments["query"]
            organism = arguments.get("organism", "hsa")

            # Detect if this is a KEGG pathway ID (e.g., hsa04210, map00010)
            is_pathway_id = bool(re.match(r"^[a-z]{2,4}\d{5}$", query)) or bool(re.match(r"^map\d{5}$", query))

            if is_pathway_id:
                # GET a specific pathway
                url = f"{KEGG_BASE_URL}/get/{query}"
                try:
                    async with httpx.AsyncClient(timeout=30) as client:
                        resp = await client.get(url)
                        if resp.status_code >= 400:
                            return [TextContent(type="text", text=f"Error: KEGG returned HTTP {resp.status_code}")]
                        text = resp.text
                except httpx.TimeoutException:
                    return [TextContent(type="text", text="Error: KEGG request timed out")]
                except Exception as e:
                    return [TextContent(type="text", text=f"Error: {str(e)}")]

                # Parse the KEGG flat file into structured fields
                fields: dict[str, Any] = {"pathway_id": query}
                current_field = ""
                current_value: list[str] = []
                for line in text.split("\n"):
                    if line.startswith("///"):
                        break
                    if line and not line[0].isspace():
                        if current_field:
                            fields[current_field.lower()] = "\n".join(current_value).strip()
                        parts = line.split(None, 1)
                        current_field = parts[0]
                        current_value = [parts[1]] if len(parts) > 1 else []
                    else:
                        current_value.append(line.strip())
                if current_field:
                    fields[current_field.lower()] = "\n".join(current_value).strip()

                # Extract key fields
                curated = {
                    "pathway_id": query,
                    "name": fields.get("name", ""),
                    "description": fields.get("description", ""),
                    "class": fields.get("class", ""),
                    "organism": fields.get("organism", ""),
                }
                # Parse gene list if present
                gene_text = fields.get("gene", "")
                if gene_text:
                    genes = []
                    for gline in gene_text.split("\n"):
                        gline = gline.strip()
                        if gline:
                            parts = gline.split(None, 1)
                            genes.append({
                                "gene_id": parts[0],
                                "description": parts[1] if len(parts) > 1 else "",
                            })
                    curated["gene_count"] = len(genes)
                    curated["genes"] = genes[:50]  # Cap at 50

                return [TextContent(type="text", text=json.dumps(curated, indent=2))]

            else:
                # FIND pathways by keyword
                url = f"{KEGG_BASE_URL}/find/pathway/{query}"
                try:
                    async with httpx.AsyncClient(timeout=30) as client:
                        resp = await client.get(url)
                        if resp.status_code >= 400:
                            return [TextContent(type="text", text=f"Error: KEGG returned HTTP {resp.status_code}")]
                        text = resp.text
                except httpx.TimeoutException:
                    return [TextContent(type="text", text="Error: KEGG request timed out")]
                except Exception as e:
                    return [TextContent(type="text", text=f"Error: {str(e)}")]

                pathways = []
                for line in text.strip().split("\n"):
                    if not line.strip():
                        continue
                    parts = line.split("\t", 1)
                    pathway_id = parts[0].replace("map:", "").replace("path:", "")
                    name = parts[1] if len(parts) > 1 else ""
                    pathways.append({"pathway_id": pathway_id, "name": name})

                    # Also get the organism-specific version
                    if not pathway_id.startswith(organism):
                        org_id = organism + pathway_id.lstrip("map0123456789"[-5:]) if pathway_id.startswith("map") else pathway_id
                        # Just keep the map ID; user can query specific org pathway
                        pass

                return [TextContent(type="text", text=json.dumps({
                    "query": query,
                    "count": len(pathways),
                    "pathways": pathways,
                    "tip": f"Use a specific pathway ID (e.g., '{organism}{pathways[0]['pathway_id'][-5:]}') with this tool to get full details." if pathways else "",
                }, indent=2))]

        # ── Section 10: NCBI BLAST REST API Tools ─────────────────────

        elif name == "blast_search":
            query_seq = arguments["query"]
            program = arguments.get("program", "blastp")
            database = arguments.get("database", "nr")
            evalue = arguments.get("evalue", 0.01)
            max_hits = arguments.get("max_hits", 10)

            # Step 1: Submit BLAST job via PUT
            try:
                async with httpx.AsyncClient(timeout=60) as client:
                    submit_resp = await client.post(
                        BLAST_BASE_URL,
                        data={
                            "CMD": "Put",
                            "QUERY": query_seq,
                            "PROGRAM": program,
                            "DATABASE": database,
                            "EXPECT": str(evalue),
                            "HITLIST_SIZE": str(max_hits),
                            "FORMAT_TYPE": "JSON2",
                        },
                    )
            except Exception as e:
                return [TextContent(type="text", text=f"Error submitting BLAST job: {str(e)}")]

            # Parse RID and RTOE from response
            body = submit_resp.text
            rid_match = re.search(r"RID\s*=\s*(\S+)", body)
            rtoe_match = re.search(r"RTOE\s*=\s*(\d+)", body)

            if not rid_match:
                return [TextContent(type="text", text=f"Error: Could not parse BLAST RID from response: {body[:300]}")]

            rid = rid_match.group(1)
            rtoe = int(rtoe_match.group(1)) if rtoe_match else 15

            # Step 2: Poll for results
            await asyncio.sleep(min(rtoe, 30))  # Wait estimated time, cap at 30s

            max_polls = 20
            poll_interval = 10
            for _ in range(max_polls):
                try:
                    async with httpx.AsyncClient(timeout=60) as client:
                        poll_resp = await client.get(
                            BLAST_BASE_URL,
                            params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"},
                        )
                except Exception as e:
                    return [TextContent(type="text", text=f"Error polling BLAST results: {str(e)}")]

                poll_body = poll_resp.text

                if "Status=WAITING" in poll_body:
                    await asyncio.sleep(poll_interval)
                    continue
                elif "Status=FAILED" in poll_body:
                    return [TextContent(type="text", text="Error: BLAST search failed on the server")]
                elif "Status=UNKNOWN" in poll_body:
                    return [TextContent(type="text", text=f"Error: BLAST RID {rid} expired or is unknown")]
                else:
                    # Results ready — try to parse JSON
                    break
            else:
                return [TextContent(type="text", text=f"Error: BLAST search timed out after polling. RID: {rid} — you can try retrieving results later.")]

            # Step 3: Parse results
            try:
                data = poll_resp.json()
            except Exception:
                # If JSON parse fails, return raw text truncated
                return [TextContent(type="text", text=f"BLAST completed (RID: {rid}) but JSON parsing failed. Raw response (truncated):\n{poll_body[:2000]}")]

            # Extract hits from the BLAST JSON2 format
            search_results = data.get("BlastOutput2", [{}])
            if isinstance(search_results, list) and search_results:
                report = search_results[0].get("report", {})
            else:
                report = search_results.get("report", {})

            search = report.get("results", {}).get("search", {})
            hits = search.get("hits", [])
            query_title = search.get("query_title", "")
            query_len = search.get("query_len", 0)
            stat = search.get("stat", {})

            curated_hits = []
            for h in hits[:max_hits]:
                desc = h.get("description", [{}])[0] if h.get("description") else {}
                hsps = h.get("hsps", [{}])
                best_hsp = hsps[0] if hsps else {}
                curated_hits.append({
                    "accession": desc.get("accession", ""),
                    "title": desc.get("title", ""),
                    "sciname": desc.get("sciname", ""),
                    "evalue": best_hsp.get("evalue"),
                    "bit_score": best_hsp.get("bit_score"),
                    "identity_percent": round(best_hsp.get("identity", 0) / best_hsp.get("align_len", 1) * 100, 1) if best_hsp.get("align_len") else None,
                    "align_length": best_hsp.get("align_len"),
                    "query_coverage": f"{best_hsp.get('query_from', '')}-{best_hsp.get('query_to', '')}",
                    "subject_coverage": f"{best_hsp.get('hit_from', '')}-{best_hsp.get('hit_to', '')}",
                })

            curated = {
                "rid": rid,
                "program": program,
                "database": database,
                "query_title": query_title,
                "query_length": query_len,
                "total_hits": len(hits),
                "returned_hits": len(curated_hits),
                "db_sequences": stat.get("db_num"),
                "hits": curated_hits,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 11: PubMed / NCBI E-utilities Tools ───────────────

        elif name == "pubmed_search":
            query = arguments["query"]
            limit = arguments.get("limit", 10)
            sort = arguments.get("sort", "relevance")

            # Step 1: esearch to get PMIDs
            search_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
            search_result = await _rest_get(search_url, params={
                "db": "pubmed",
                "term": query,
                "retmax": str(limit),
                "retmode": "json",
                "sort": sort,
            })

            if not search_result["success"]:
                return [TextContent(type="text", text=f"Error: {search_result['error']}")]

            esearch = search_result["data"].get("esearchresult", {})
            id_list = esearch.get("idlist", [])

            if not id_list:
                return [TextContent(type="text", text=json.dumps({
                    "query": query, "count": 0, "results": [],
                    "message": "No PubMed articles found for this query."
                }, indent=2))]

            # Step 2: esummary to get article metadata
            summary_url = f"{EUTILS_BASE_URL}/esummary.fcgi"
            summary_result = await _rest_get(summary_url, params={
                "db": "pubmed",
                "id": ",".join(id_list),
                "retmode": "json",
            }, use_cache=False)

            if not summary_result["success"]:
                return [TextContent(type="text", text=f"Error: {summary_result['error']}")]

            result_data = summary_result["data"].get("result", {})

            # Step 3: efetch to get abstracts
            fetch_url = f"{EUTILS_BASE_URL}/efetch.fcgi"
            abstracts: dict[str, str] = {}
            try:
                async with httpx.AsyncClient(timeout=30) as client:
                    fetch_resp = await client.get(fetch_url, params={
                        "db": "pubmed",
                        "id": ",".join(id_list),
                        "rettype": "abstract",
                        "retmode": "xml",
                    })
                    if fetch_resp.status_code == 200:
                        # Parse abstracts from XML
                        xml_text = fetch_resp.text
                        # Simple extraction — find <PMID> and <AbstractText> pairs
                        import xml.etree.ElementTree as ET
                        try:
                            root = ET.fromstring(xml_text)
                            for article in root.findall(".//PubmedArticle"):
                                pmid_elem = article.find(".//PMID")
                                abstract_elems = article.findall(".//AbstractText")
                                if pmid_elem is not None and abstract_elems:
                                    pmid = pmid_elem.text
                                    abstract_parts = []
                                    for ae in abstract_elems:
                                        label = ae.get("Label", "")
                                        text = ae.text or ""
                                        if label:
                                            abstract_parts.append(f"{label}: {text}")
                                        else:
                                            abstract_parts.append(text)
                                    abstracts[pmid] = " ".join(abstract_parts)
                        except ET.ParseError:
                            pass  # Abstracts are optional, continue without them
            except Exception:
                pass  # Abstracts are a bonus; don't fail the whole query

            curated = []
            for pmid in id_list:
                rec = result_data.get(pmid, {})
                if not rec or pmid == "uids":
                    continue

                authors = rec.get("authors", [])
                author_names = [a.get("name", "") for a in authors[:5]]
                if len(authors) > 5:
                    author_names.append(f"... (+{len(authors) - 5} more)")

                elocation = rec.get("elocationid", "")
                doi = ""
                if elocation and elocation.startswith("doi:"):
                    doi = elocation.replace("doi: ", "").replace("doi:", "")

                curated.append({
                    "pmid": pmid,
                    "title": rec.get("title", ""),
                    "authors": author_names,
                    "journal": rec.get("fulljournalname", rec.get("source", "")),
                    "pub_date": rec.get("pubdate", ""),
                    "doi": doi,
                    "abstract": abstracts.get(pmid, ""),
                })

            total_count = int(esearch.get("count", len(curated)))
            return [TextContent(type="text", text=json.dumps({
                "query": query,
                "total_in_pubmed": total_count,
                "returned": len(curated),
                "results": curated,
            }, indent=2))]

        # ── Section 12: Ensembl Regulation Tools ──────────────────────

        elif name == "ensembl_get_regulation":
            region = arguments["region"]
            species = arguments.get("species", "homo_sapiens")

            url = f"{ENSEMBL_BASE_URL}/overlap/region/{species}/{region}"
            result = await _rest_get(url, params={
                "feature": "regulatory",
                "content-type": "application/json",
            })

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            features = result["data"] if isinstance(result["data"], list) else []
            curated = []
            type_counts: dict[str, int] = {}
            for f in features:
                ftype = f.get("feature_type", f.get("description", "unknown"))
                type_counts[ftype] = type_counts.get(ftype, 0) + 1
                curated.append({
                    "id": f.get("id", ""),
                    "feature_type": ftype,
                    "start": f.get("start"),
                    "end": f.get("end"),
                    "strand": f.get("strand"),
                    "description": f.get("description", ""),
                    "source": f.get("source", ""),
                    "biotype": f.get("biotype", ""),
                })

            return [TextContent(type="text", text=json.dumps({
                "region": region,
                "species": species,
                "total_features": len(curated),
                "type_summary": type_counts,
                "features": curated,
            }, indent=2))]

        # ── Section 13: AlphaFold REST API Tools ──────────────────────

        elif name == "alphafold_get_prediction":
            accession = arguments["accession"].upper()
            url = f"{ALPHAFOLD_BASE_URL}/prediction/{accession}"
            result = await _rest_get(url)

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            data = result["data"]
            entries = data if isinstance(data, list) else [data]

            if not entries:
                return [TextContent(type="text", text=json.dumps({
                    "accession": accession,
                    "message": "No AlphaFold prediction found for this accession."
                }, indent=2))]

            # Take the first (canonical) entry
            e = entries[0]
            curated = {
                "accession": accession,
                "entry_id": e.get("entryId", ""),
                "gene": e.get("gene", ""),
                "uniprot_description": e.get("uniprotDescription", ""),
                "organism": e.get("organismScientificName", ""),
                "tax_id": e.get("taxId"),
                "sequence_length": len(e.get("sequence", "")),
                "model_confidence": {
                    "global_plddt": e.get("globalMetricValue"),
                    "fraction_very_high": e.get("fractionPlddtVeryHigh"),
                    "fraction_confident": e.get("fractionPlddtConfident"),
                    "fraction_low": e.get("fractionPlddtLow"),
                    "fraction_very_low": e.get("fractionPlddtVeryLow"),
                },
                "coverage": {
                    "start": e.get("sequenceStart"),
                    "end": e.get("sequenceEnd"),
                },
                "model_version": e.get("latestVersion"),
                "model_created": e.get("modelCreatedDate", ""),
                "urls": {
                    "pdb": e.get("pdbUrl", ""),
                    "cif": e.get("cifUrl", ""),
                    "pae_image": e.get("paeImageUrl", ""),
                },
                "is_reviewed": e.get("isReviewed", False),
                "total_entries": len(entries),
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 14: gnomAD GraphQL API Tools ─────────────────────────

        elif name == "gnomad_get_variant":
            variant_id = arguments["variant_id"].strip()
            dataset = arguments.get("dataset", "gnomad_r4")

            # Determine if this is an rsID or variant ID lookup
            is_rsid = variant_id.lower().startswith("rs")

            if is_rsid:
                query_field = "rsid"
                query_value = variant_id.lower()
            else:
                query_field = "variantId"
                query_value = variant_id

            graphql_query = """
            {
                variant(%s: "%s", dataset: %s) {
                    variantId
                    reference_genome
                    chrom
                    pos
                    ref
                    alt
                    rsids
                    flags
                    exome {
                        ac
                        an
                        ac_hemi
                        ac_hom
                        filters
                        populations {
                            id
                            ac
                            an
                            ac_hemi
                            ac_hom
                        }
                    }
                    genome {
                        ac
                        an
                        ac_hemi
                        ac_hom
                        filters
                        populations {
                            id
                            ac
                            an
                            ac_hemi
                            ac_hom
                        }
                    }
                    sortedTranscriptConsequences {
                        gene_symbol
                        major_consequence
                        hgvsc
                        hgvsp
                        lof
                        canonical
                        transcript_id
                    }
                }
            }
            """ % (query_field, query_value, dataset)

            try:
                async with httpx.AsyncClient(timeout=30) as client:
                    resp = await client.post(
                        GNOMAD_API_URL,
                        json={"query": graphql_query},
                        headers={"Content-Type": "application/json"}
                    )
                    resp.raise_for_status()
                    gql_data = resp.json()
            except Exception as exc:
                return [TextContent(type="text", text=f"Error querying gnomAD: {exc}")]

            errors = gql_data.get("errors")
            if errors:
                return [TextContent(type="text", text=json.dumps({
                    "error": "gnomAD query error",
                    "details": [e.get("message", str(e)) for e in errors]
                }, indent=2))]

            v = gql_data.get("data", {}).get("variant")
            if not v:
                return [TextContent(type="text", text=json.dumps({
                    "variant_id": variant_id,
                    "message": "Variant not found in gnomAD."
                }, indent=2))]

            def _pop_freqs(source_data):
                if not source_data:
                    return None
                ac = source_data.get("ac", 0)
                an = source_data.get("an", 0)
                af = ac / an if an > 0 else 0
                pops = []
                for p in source_data.get("populations", []):
                    p_ac = p.get("ac", 0)
                    p_an = p.get("an", 0)
                    pops.append({
                        "population": p.get("id", ""),
                        "ac": p_ac,
                        "an": p_an,
                        "af": p_ac / p_an if p_an > 0 else 0,
                        "ac_hom": p.get("ac_hom", 0),
                    })
                return {
                    "ac": ac,
                    "an": an,
                    "af": round(af, 8),
                    "ac_hom": source_data.get("ac_hom", 0),
                    "ac_hemi": source_data.get("ac_hemi", 0),
                    "filters": source_data.get("filters", []),
                    "populations": sorted(pops, key=lambda x: x["af"], reverse=True),
                }

            # Get canonical transcript consequence
            consequences = v.get("sortedTranscriptConsequences", []) or []
            canonical = next((c for c in consequences if c.get("canonical")), None)
            top_consequence = canonical or (consequences[0] if consequences else None)

            curated = {
                "variant_id": v.get("variantId", ""),
                "rsids": v.get("rsids", []),
                "reference_genome": v.get("reference_genome", ""),
                "chrom": v.get("chrom", ""),
                "pos": v.get("pos"),
                "ref": v.get("ref", ""),
                "alt": v.get("alt", ""),
                "flags": v.get("flags", []),
                "exome": _pop_freqs(v.get("exome")),
                "genome": _pop_freqs(v.get("genome")),
                "consequence": {
                    "gene_symbol": top_consequence.get("gene_symbol", "") if top_consequence else "",
                    "major_consequence": top_consequence.get("major_consequence", "") if top_consequence else "",
                    "hgvsc": top_consequence.get("hgvsc", "") if top_consequence else "",
                    "hgvsp": top_consequence.get("hgvsp", "") if top_consequence else "",
                    "lof": top_consequence.get("lof", "") if top_consequence else "",
                } if top_consequence else None,
                "dataset": dataset,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 15: GTEx REST API Tools ──────────────────────────────

        elif name == "gtex_get_expression":
            gene = arguments["gene"].strip()
            dataset = arguments.get("dataset", "gtex_v8")

            # Map dataset to gencode version
            gencode_map = {"gtex_v8": "v26", "gtex_v10": "v39", "gtex_v7": "v19"}
            gencode_version = gencode_map.get(dataset, "v26")

            # Step 1: Resolve gene to versioned gencodeId
            gene_url = f"{GTEX_BASE_URL}/reference/gene"
            gene_params = {
                "geneId": gene,
                "gencodeVersion": gencode_version,
                "genomeBuild": "GRCh38/hg38",
            }
            result = await _rest_get(gene_url, params=gene_params)
            if not result["success"]:
                return [TextContent(type="text", text=f"Error looking up gene: {result['error']}")]

            gene_data = result["data"].get("data", [])
            if not gene_data:
                return [TextContent(type="text", text=json.dumps({
                    "gene": gene,
                    "message": f"Gene not found in GTEx ({dataset}). Try a different gene symbol or Ensembl ID."
                }, indent=2))]

            gencode_id = gene_data[0].get("gencodeId", "")
            gene_symbol = gene_data[0].get("geneSymbol", gene)
            gene_description = gene_data[0].get("description", "")

            # Step 2: Get median expression across tissues
            expr_url = f"{GTEX_BASE_URL}/expression/medianGeneExpression"
            expr_params = {
                "gencodeId": gencode_id,
                "datasetId": dataset,
                "itemsPerPage": 250,
            }
            expr_result = await _rest_get(expr_url, params=expr_params)
            if not expr_result["success"]:
                return [TextContent(type="text", text=f"Error fetching expression: {expr_result['error']}")]

            expr_data = expr_result["data"].get("data", [])

            # Sort by expression level (highest first)
            tissues = []
            for t in sorted(expr_data, key=lambda x: x.get("median", 0), reverse=True):
                tissues.append({
                    "tissue": t.get("tissueSiteDetailId", ""),
                    "median_tpm": round(t.get("median", 0), 4),
                    "unit": t.get("unit", "TPM"),
                })

            curated = {
                "gene_symbol": gene_symbol,
                "gencode_id": gencode_id,
                "description": gene_description,
                "dataset": dataset,
                "total_tissues": len(tissues),
                "tissues": tissues,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 16: HPO REST API Tools ───────────────────────────────

        elif name == "hpo_search":
            query = arguments["query"].strip()
            category = arguments.get("category", "search")
            limit = arguments.get("limit", 20)

            if category == "search":
                # General search for terms, genes, diseases
                url = f"{HPO_BASE_URL}/search"
                result = await _rest_get(url, params={"q": query, "max": limit})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error searching HPO: {result['error']}")]

                data = result["data"]
                curated = {
                    "query": query,
                    "terms": [],
                    "genes": [],
                    "diseases": [],
                }
                # Extract terms
                for t in (data.get("terms", []) or [])[:limit]:
                    curated["terms"].append({
                        "id": t.get("ontologyId", ""),
                        "name": t.get("name", ""),
                        "definition": t.get("definition", ""),
                    })
                # Extract genes
                for g in (data.get("genes", []) or [])[:limit]:
                    curated["genes"].append({
                        "gene_symbol": g.get("entrezGeneSymbol", ""),
                        "gene_id": g.get("entrezGeneId"),
                    })
                # Extract diseases
                for d in (data.get("diseases", []) or [])[:limit]:
                    curated["diseases"].append({
                        "disease_id": d.get("diseaseId", ""),
                        "disease_name": d.get("diseaseName", d.get("dbName", "")),
                        "database": d.get("db", ""),
                    })
                curated["total_terms"] = len(curated["terms"])
                curated["total_genes"] = len(curated["genes"])
                curated["total_diseases"] = len(curated["diseases"])
                return [TextContent(type="text", text=json.dumps(curated, indent=2))]

            elif category == "term":
                # Get details for a specific HPO term
                hpo_id = query if query.startswith("HP:") else query
                url = f"{HPO_BASE_URL}/term/{hpo_id}"
                result = await _rest_get(url)
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error fetching HPO term: {result['error']}")]
                data = result["data"]
                details = data.get("details", data)
                relations = data.get("relations", {})
                curated = {
                    "id": details.get("ontologyId", details.get("id", hpo_id)),
                    "name": details.get("name", ""),
                    "definition": details.get("definition", ""),
                    "synonyms": details.get("synonyms", []),
                    "parents": [{"id": p.get("ontologyId", ""), "name": p.get("name", "")}
                                for p in relations.get("parents", [])],
                    "children_count": len(relations.get("children", [])),
                }
                return [TextContent(type="text", text=json.dumps(curated, indent=2))]

            elif category == "genes":
                # Get genes associated with an HPO term
                hpo_id = query if query.startswith("HP:") else query
                url = f"{HPO_BASE_URL}/term/{hpo_id}/genes"
                result = await _rest_get(url, params={"offset": 0, "max": limit})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error fetching genes for HPO term: {result['error']}")]
                data = result["data"]
                genes_raw = data.get("genes", data) if isinstance(data, dict) else data
                if isinstance(genes_raw, dict):
                    genes_raw = genes_raw.get("genes", [])
                genes = []
                for g in (genes_raw or [])[:limit]:
                    genes.append({
                        "gene_symbol": g.get("entrezGeneSymbol", ""),
                        "gene_id": g.get("entrezGeneId"),
                        "disease_count": len(g.get("dbDiseases", [])),
                    })
                curated = {
                    "hpo_term": hpo_id,
                    "total_genes": len(genes),
                    "genes": genes,
                }
                return [TextContent(type="text", text=json.dumps(curated, indent=2))]

            elif category == "diseases":
                # Get diseases associated with an HPO term
                hpo_id = query if query.startswith("HP:") else query
                url = f"{HPO_BASE_URL}/term/{hpo_id}/diseases"
                result = await _rest_get(url, params={"offset": 0, "max": limit})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error fetching diseases for HPO term: {result['error']}")]
                data = result["data"]
                diseases_raw = data.get("diseases", data) if isinstance(data, dict) else data
                if isinstance(diseases_raw, dict):
                    diseases_raw = diseases_raw.get("diseases", [])
                diseases = []
                for d in (diseases_raw or [])[:limit]:
                    diseases.append({
                        "disease_id": d.get("diseaseId", d.get("dbId", "")),
                        "disease_name": d.get("diseaseName", d.get("dbName", "")),
                        "database": d.get("db", ""),
                    })
                curated = {
                    "hpo_term": hpo_id,
                    "total_diseases": len(diseases),
                    "diseases": diseases,
                }
                return [TextContent(type="text", text=json.dumps(curated, indent=2))]

            else:
                return [TextContent(type="text", text=f"Unknown HPO category: {category}. Use 'search', 'term', 'genes', or 'diseases'.")]

        # ── Section 17: Genome Coordinate Tools ─────────────────────────

        elif name == "liftover_coordinates":
            region = arguments["region"].strip()
            source_assembly = arguments.get("source_assembly", "GRCh37")
            target_assembly = arguments.get("target_assembly", "GRCh38")
            species = arguments.get("species", "human")

            # Normalize region format: accept "17:43044295-43170245" or "17:43044295..43170245"
            region = region.replace("-", "..", 1) if ".." not in region else region

            # Ensembl map endpoint: /map/{species}/{asm_one}/{region}/{asm_two}
            url = f"{ENSEMBL_BASE_URL}/map/{species}/{source_assembly}/{region}/{target_assembly}"
            result = await _rest_get(url, headers={"Content-Type": "application/json"})

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            data = result["data"]
            mappings = data.get("mappings", [])

            if not mappings:
                return [TextContent(type="text", text=json.dumps({
                    "region": region,
                    "source_assembly": source_assembly,
                    "target_assembly": target_assembly,
                    "message": "No mapping found. The region may not exist in the target assembly, or the coordinates may be invalid."
                }, indent=2))]

            curated_mappings = []
            for m in mappings:
                original = m.get("original", {})
                mapped = m.get("mapped", {})
                curated_mappings.append({
                    "original": {
                        "assembly": original.get("assembly", source_assembly),
                        "region": f"{original.get('seq_region_name', '')}:{original.get('start', '')}-{original.get('end', '')}",
                        "strand": original.get("strand", 1),
                    },
                    "mapped": {
                        "assembly": mapped.get("assembly", target_assembly),
                        "region": f"{mapped.get('seq_region_name', '')}:{mapped.get('start', '')}-{mapped.get('end', '')}",
                        "strand": mapped.get("strand", 1),
                        "coord_system": mapped.get("coord_system", ""),
                    },
                })

            curated = {
                "source_assembly": source_assembly,
                "target_assembly": target_assembly,
                "species": species,
                "input_region": region,
                "total_mappings": len(curated_mappings),
                "mappings": curated_mappings,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 18: Reactome Pathway Tools ───────────────────────────

        elif name == "reactome_get_pathway":
            pathway_id = arguments.get("pathway_id", "").strip()
            query = arguments.get("query", "").strip()
            species = arguments.get("species", "Homo sapiens")
            include_participants = arguments.get("include_participants", False)
            include_hierarchy = arguments.get("include_hierarchy", False)
            limit = arguments.get("limit", 10)

            if not pathway_id and not query:
                return [TextContent(type="text", text="Error: provide either 'pathway_id' (e.g., 'R-HSA-1640170') or 'query' (e.g., 'TP53' or 'apoptosis')")]

            # --- Direct lookup by Reactome stable ID ---
            if pathway_id:
                url = f"{REACTOME_BASE_URL}/data/query/{pathway_id}"
                result = await _rest_get(url, headers={"Accept": "application/json"})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error: {result['error']}")]

                pw = result["data"]
                curated: dict[str, Any] = {
                    "stableId": pw.get("stId", pathway_id),
                    "name": pw.get("displayName", ""),
                    "species": pw.get("speciesName", ""),
                    "schemaClass": pw.get("schemaClass", ""),
                    "isInDisease": pw.get("isInDisease", False),
                    "compartments": [c.get("displayName", "") for c in pw.get("compartment", [])],
                    "summation": "",
                }
                # Extract summation text if available
                summations = pw.get("summation", [])
                if summations and isinstance(summations, list):
                    texts = [s.get("text", "") for s in summations if s.get("text")]
                    curated["summation"] = " ".join(texts)

                # Sub-events (child pathways/reactions)
                has_event = pw.get("hasEvent", [])
                if has_event:
                    curated["subEvents"] = [
                        {"stableId": e.get("stId", ""), "name": e.get("displayName", ""), "type": e.get("schemaClass", "")}
                        for e in has_event[:50]
                    ]
                    curated["totalSubEvents"] = len(has_event)

                # Optional: participating molecules
                if include_participants:
                    part_url = f"{REACTOME_BASE_URL}/data/participants/{pathway_id}/referenceEntities"
                    part_result = await _rest_get(part_url, headers={"Accept": "application/json"})
                    if part_result["success"] and isinstance(part_result["data"], list):
                        participants = []
                        seen = set()
                        for p in part_result["data"]:
                            display = p.get("displayName", "")
                            if display and display not in seen:
                                seen.add(display)
                                participants.append({
                                    "name": display,
                                    "type": p.get("schemaClass", ""),
                                    "identifier": p.get("identifier", ""),
                                    "database": p.get("databaseName", ""),
                                })
                        curated["participants"] = participants[:100]
                        curated["totalParticipants"] = len(participants)

                # Optional: pathway hierarchy (ancestors)
                if include_hierarchy:
                    anc_url = f"{REACTOME_BASE_URL}/data/event/{pathway_id}/ancestors"
                    anc_result = await _rest_get(anc_url, headers={"Accept": "application/json"})
                    if anc_result["success"] and isinstance(anc_result["data"], list):
                        hierarchy = []
                        for level in anc_result["data"]:
                            if isinstance(level, list):
                                hierarchy.append([{"stableId": a.get("stId", ""), "name": a.get("displayName", "")} for a in level])
                            else:
                                hierarchy.append({"stableId": level.get("stId", ""), "name": level.get("displayName", "")})
                        curated["hierarchy"] = hierarchy

                return [TextContent(type="text", text=json.dumps(curated, indent=2))]

            # --- Search by keyword/gene name ---
            url = f"{REACTOME_BASE_URL}/search/query"
            params = {"query": query, "species": species, "types": "Pathway", "cluster": "true", "rows": str(limit)}
            result = await _rest_get(url, params=params, headers={"Accept": "application/json"})
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            data = result["data"]
            results_list = data.get("results", [])

            pathways = []
            for group in results_list:
                entries = group.get("entries", [])
                for entry in entries:
                    pathways.append({
                        "stableId": entry.get("stId", ""),
                        "name": entry.get("name", ""),
                        "species": entry.get("species", []),
                        "summation": entry.get("summation", ""),
                    })
                    if len(pathways) >= limit:
                        break
                if len(pathways) >= limit:
                    break

            curated_search = {
                "query": query,
                "species": species,
                "totalFound": data.get("found", len(pathways)),
                "returned": len(pathways),
                "pathways": pathways,
            }
            return [TextContent(type="text", text=json.dumps(curated_search, indent=2))]

        # ── Section 19: Gene Ontology Enrichment ─────────────────────────

        elif name == "gene_ontology_enrich":
            genes_str = arguments.get("genes", "").strip()
            if not genes_str:
                return [TextContent(type="text", text="Error: 'genes' is required (comma-separated gene symbols, e.g., 'TP53,BRCA1,EGFR')")]

            gene_list = [g.strip() for g in genes_str.split(",") if g.strip()]
            organism = arguments.get("organism", "hsapiens")
            sources_str = arguments.get("sources", "GO:BP,GO:MF,GO:CC")
            sources = [s.strip() for s in sources_str.split(",") if s.strip()]
            threshold = arguments.get("threshold", 0.05)
            correction = arguments.get("correction_method", "g_SCS")
            no_iea = arguments.get("no_iea", False)
            ordered = arguments.get("ordered", False)
            background_str = arguments.get("background", "")
            limit = arguments.get("limit", 25)

            # Build g:Profiler API payload
            payload: dict[str, Any] = {
                "organism": organism,
                "query": gene_list,
                "sources": sources,
                "user_threshold": threshold,
                "significance_threshold_method": correction,
                "no_iea": no_iea,
                "ordered": ordered,
                "all_results": False,
                "no_evidences": False,
                "combined": False,
                "measure_underrepresentation": False,
                "domain_scope": "annotated",
                "output": "json",
            }
            if background_str:
                payload["background"] = [b.strip() for b in background_str.split(",") if b.strip()]
                payload["domain_scope"] = "custom"

            # POST to g:Profiler gost endpoint
            url = f"{GPROFILER_BASE_URL}/gost/profile/"
            try:
                async with httpx.AsyncClient(timeout=60.0) as client:
                    resp = await client.post(url, json=payload)
                    resp.raise_for_status()
                    data = resp.json()
            except httpx.HTTPStatusError as e:
                return [TextContent(type="text", text=f"Error: g:Profiler returned HTTP {e.response.status_code}")]
            except Exception as e:
                return [TextContent(type="text", text=f"Error calling g:Profiler: {str(e)}")]

            # Parse results
            results = data.get("result", [])
            if not results:
                return [TextContent(type="text", text=json.dumps({
                    "genes": gene_list,
                    "organism": organism,
                    "sources": sources,
                    "message": "No significantly enriched terms found at the given threshold.",
                    "threshold": threshold,
                }, indent=2))]

            # Sort by p-value and limit
            results.sort(key=lambda x: x.get("p_value", 1.0))
            results = results[:limit]

            enriched_terms = []
            for r in results:
                enriched_terms.append({
                    "source": r.get("source", ""),
                    "term_id": r.get("native", ""),
                    "term_name": r.get("name", ""),
                    "p_value": r.get("p_value", None),
                    "term_size": r.get("term_size", 0),
                    "query_size": r.get("query_size", 0),
                    "intersection_size": r.get("intersection_size", 0),
                    "precision": r.get("precision", 0),
                    "recall": r.get("recall", 0),
                    "intersecting_genes": r.get("intersections", []),
                })

            # Metadata from response
            meta = data.get("meta", {})
            query_info = meta.get("genes_metadata", {})

            curated = {
                "query_genes": gene_list,
                "organism": organism,
                "sources": sources,
                "threshold": threshold,
                "correction_method": correction,
                "total_enriched_terms": len(data.get("result", [])),
                "returned": len(enriched_terms),
                "enriched_terms": enriched_terms,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 20: COSMIC Somatic Mutations ─────────────────────────

        elif name == "cosmic_search":
            query = arguments.get("query", "").strip()
            gene_filter = arguments.get("gene", "").strip()
            limit = min(arguments.get("limit", 20), 500)

            if not query:
                return [TextContent(type="text", text="Error: 'query' is required (e.g., 'BRAF', 'V600E', 'lung')")]

            # Build NLM Clinical Tables query
            # ef= requests extra fields in the response
            all_fields = "MutationID,GeneName,MutationCDS,MutationAA,MutationDescription,MutationGenomePosition,MutationStrand,PrimaryHistology,PrimarySite,PubmedPMID,AccessionNumber,HGNC_ID"
            params: dict[str, str] = {
                "terms": query,
                "maxList": str(limit),
                "ef": all_fields,
            }
            # Add gene filter constraint using Elasticsearch query syntax
            if gene_filter:
                params["q"] = f"GeneName:{gene_filter}"

            result = await _rest_get(COSMIC_API_URL, params=params)
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            data = result["data"]
            # Response format: [total_count, code_array, extra_data_hash, display_strings]
            if not isinstance(data, list) or len(data) < 4:
                return [TextContent(type="text", text="Error: unexpected response format from COSMIC API")]

            total_count = data[0]
            mutation_ids = data[1] if data[1] else []
            extra_data = data[2] if data[2] else {}

            mutations = []
            for mid in mutation_ids:
                mid_str = str(mid)
                extras = extra_data.get(mid_str, {})
                mutations.append({
                    "mutation_id": extras.get("MutationID", mid_str),
                    "gene": extras.get("GeneName", ""),
                    "cds_mutation": extras.get("MutationCDS", ""),
                    "aa_mutation": extras.get("MutationAA", ""),
                    "description": extras.get("MutationDescription", ""),
                    "genome_position": extras.get("MutationGenomePosition", ""),
                    "strand": extras.get("MutationStrand", ""),
                    "primary_histology": extras.get("PrimaryHistology", ""),
                    "primary_site": extras.get("PrimarySite", ""),
                    "pubmed_id": extras.get("PubmedPMID", ""),
                })

            curated = {
                "query": query,
                "gene_filter": gene_filter or None,
                "total_found": total_count,
                "returned": len(mutations),
                "note": "Data from COSMIC V89 (GRCh37) via NLM Clinical Tables API",
                "mutations": mutations,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 21: OMIM Gene-Disease Associations ───────────────────

        elif name == "omim_search":
            query = arguments.get("query", "").strip()
            gene_id = arguments.get("gene_id", "").strip()
            limit = arguments.get("limit", 20)

            if not query and not gene_id:
                return [TextContent(type="text", text="Error: provide 'query' (gene symbol, disease name, or MIM number) or 'gene_id' (NCBI Gene ID)")]

            omim_ids: list[str] = []

            # Strategy 1: If gene_id provided, use elink to get gene→omim associations
            if gene_id:
                elink_url = f"{EUTILS_BASE_URL}/elink.fcgi"
                elink_params = {"dbfrom": "gene", "db": "omim", "id": gene_id, "retmode": "json"}
                elink_result = await _rest_get(elink_url, params=elink_params)
                if elink_result["success"]:
                    linksets = elink_result["data"].get("linksets", [])
                    for ls in linksets:
                        for ldb in ls.get("linksetdbs", []):
                            if ldb.get("dbto") == "omim":
                                omim_ids.extend([str(x) for x in ldb.get("links", [])])

            # Strategy 2: If query provided (or elink returned nothing), search OMIM directly
            if query and not omim_ids:
                esearch_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
                esearch_params = {"db": "omim", "term": query, "retmode": "json", "retmax": str(limit)}
                esearch_result = await _rest_get(esearch_url, params=esearch_params)
                if esearch_result["success"]:
                    omim_ids = esearch_result["data"].get("esearchresult", {}).get("idlist", [])
            elif query and omim_ids:
                # Both gene_id and query: also do a search and merge
                esearch_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
                esearch_params = {"db": "omim", "term": query, "retmode": "json", "retmax": str(limit)}
                esearch_result = await _rest_get(esearch_url, params=esearch_params)
                if esearch_result["success"]:
                    search_ids = esearch_result["data"].get("esearchresult", {}).get("idlist", [])
                    # Merge, keeping gene-linked IDs first
                    seen = set(omim_ids)
                    for sid in search_ids:
                        if sid not in seen:
                            omim_ids.append(sid)
                            seen.add(sid)

            if not omim_ids:
                return [TextContent(type="text", text=json.dumps({
                    "query": query or None,
                    "gene_id": gene_id or None,
                    "message": "No OMIM entries found.",
                }, indent=2))]

            # Limit and fetch summaries
            omim_ids = omim_ids[:limit]
            esummary_url = f"{EUTILS_BASE_URL}/esummary.fcgi"
            esummary_params = {"db": "omim", "id": ",".join(omim_ids), "retmode": "json"}
            esummary_result = await _rest_get(esummary_url, params=esummary_params)

            entries = []
            if esummary_result["success"]:
                result_data = esummary_result["data"].get("result", {})
                for oid in omim_ids:
                    entry = result_data.get(oid, {})
                    if not entry or "error" in entry:
                        continue
                    # Determine entry type from OID prefix: * = gene, # = phenotype, + = gene+phenotype, % = phenotype, none = other
                    oid_str = entry.get("oid", "")
                    prefix = oid_str[0] if oid_str and oid_str[0] in "*#+%" else ""
                    type_map = {"*": "gene", "#": "phenotype", "+": "gene_and_phenotype", "%": "phenotype_only"}
                    entry_type = type_map.get(prefix, "other")

                    entries.append({
                        "mim_number": oid,
                        "title": entry.get("title", ""),
                        "type": entry_type,
                        "alternative_titles": entry.get("alttitles", "") or None,
                    })

            curated = {
                "query": query or None,
                "gene_id": gene_id or None,
                "total_entries": len(entries),
                "entries": entries,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 22: Expression Atlas ─────────────────────────────────

        elif name == "expression_atlas_search":
            species = arguments.get("species", "").strip().lower()
            experiment_type = arguments.get("experiment_type", "").strip().lower()
            keyword = arguments.get("keyword", "").strip().lower()
            limit = arguments.get("limit", 15)

            # Fetch experiment listing from Expression Atlas
            result = await _rest_get(EXPRESSION_ATLAS_URL, headers={"Accept": "application/json"})
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            all_experiments = result["data"].get("experiments", [])

            # Apply filters
            filtered = []
            for exp in all_experiments:
                if species and species not in exp.get("species", "").lower():
                    continue
                if experiment_type:
                    et = exp.get("experimentType", "").lower()
                    if experiment_type not in et:
                        continue
                if keyword and keyword not in exp.get("experimentDescription", "").lower():
                    continue
                filtered.append(exp)

            filtered = filtered[:limit]

            experiments = []
            for exp in filtered:
                experiments.append({
                    "accession": exp.get("experimentAccession", ""),
                    "description": exp.get("experimentDescription", ""),
                    "species": exp.get("species", ""),
                    "type": exp.get("experimentType", ""),
                    "technology": exp.get("technologyType", []),
                    "assays": exp.get("numberOfAssays", 0),
                    "factors": exp.get("experimentalFactors", []),
                    "lastUpdate": exp.get("lastUpdate", ""),
                })

            curated = {
                "species_filter": species or None,
                "type_filter": experiment_type or None,
                "keyword_filter": keyword or None,
                "total_matching": len(filtered),
                "returned": len(experiments),
                "experiments": experiments,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 23: NCBI Gene Links ──────────────────────────────────

        elif name == "ncbi_gene_links":
            gene_id = arguments.get("gene_id", "").strip()
            gene_symbol = arguments.get("gene_symbol", "").strip()
            limit = arguments.get("limit", 15)

            if not gene_id and not gene_symbol:
                return [TextContent(type="text", text="Error: provide 'gene_id' (e.g., '7157') or 'gene_symbol' (e.g., 'TP53')")]

            # Resolve gene_symbol to gene_id if needed
            if not gene_id and gene_symbol:
                esearch_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
                esearch_params = {"db": "gene", "term": f"{gene_symbol}[Gene Name] AND Homo sapiens[Organism]", "retmode": "json", "retmax": "1"}
                sr = await _rest_get(esearch_url, params=esearch_params)
                if sr["success"]:
                    ids = sr["data"].get("esearchresult", {}).get("idlist", [])
                    if ids:
                        gene_id = ids[0]
                if not gene_id:
                    return [TextContent(type="text", text=f"Error: could not resolve gene symbol '{gene_symbol}' to a Gene ID")]

            # Get gene neighbors via elink
            elink_url = f"{EUTILS_BASE_URL}/elink.fcgi"
            elink_params = {"dbfrom": "gene", "db": "gene", "id": gene_id, "retmode": "json", "cmd": "neighbor"}
            elink_result = await _rest_get(elink_url, params=elink_params)
            if not elink_result["success"]:
                return [TextContent(type="text", text=f"Error: {elink_result['error']}")]

            neighbor_ids: list[str] = []
            linksets = elink_result["data"].get("linksets", [])
            for ls in linksets:
                for ldb in ls.get("linksetdbs", []):
                    if ldb.get("linkname") == "gene_gene_neighbors":
                        neighbor_ids = [str(x) for x in ldb.get("links", [])]

            if not neighbor_ids:
                return [TextContent(type="text", text=json.dumps({
                    "gene_id": gene_id, "gene_symbol": gene_symbol or None,
                    "message": "No gene neighbors found."
                }, indent=2))]

            neighbor_ids = neighbor_ids[:limit]

            # Get summaries for neighbor genes
            esummary_url = f"{EUTILS_BASE_URL}/esummary.fcgi"
            esummary_params = {"db": "gene", "id": ",".join(neighbor_ids), "retmode": "json"}
            esum_result = await _rest_get(esummary_url, params=esummary_params)

            related_genes = []
            if esum_result["success"]:
                result_data = esum_result["data"].get("result", {})
                for nid in neighbor_ids:
                    entry = result_data.get(nid, {})
                    if not entry or "error" in entry:
                        continue
                    related_genes.append({
                        "gene_id": nid,
                        "symbol": entry.get("nomenclaturesymbol", entry.get("name", "")),
                        "name": entry.get("description", ""),
                        "chromosome": entry.get("chromosome", ""),
                        "map_location": entry.get("maplocation", ""),
                    })

            curated = {
                "query_gene_id": gene_id,
                "query_gene_symbol": gene_symbol or None,
                "total_neighbors": len(related_genes),
                "related_genes": related_genes,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 24: ENCODE Experiments ───────────────────────────────

        elif name == "encode_get_experiments":
            query = arguments.get("query", "").strip()
            assay = arguments.get("assay", "").strip()
            organism = arguments.get("organism", "").strip()
            limit = arguments.get("limit", 10)

            if not query:
                return [TextContent(type="text", text="Error: 'query' is required (gene name, cell type, or keyword)")]

            # Build ENCODE search URL
            params: dict[str, str] = {
                "type": "Experiment",
                "searchTerm": query,
                "format": "json",
                "limit": str(limit),
                "status": "released",
            }
            if assay:
                params["assay_title"] = assay
            if organism:
                params["replicates.library.biosample.donor.organism.scientific_name"] = organism

            url = f"{ENCODE_BASE_URL}/search/"
            result = await _rest_get(url, params=params, headers={"Accept": "application/json"})
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            data = result["data"]
            graph = data.get("@graph", [])

            experiments = []
            for exp in graph[:limit]:
                target_info = exp.get("target", {})
                experiments.append({
                    "accession": exp.get("accession", ""),
                    "assay": exp.get("assay_title", ""),
                    "description": exp.get("description", ""),
                    "target": target_info.get("label", "") if isinstance(target_info, dict) else "",
                    "biosample": exp.get("biosample_summary", ""),
                    "biosample_type": exp.get("biosample_ontology", {}).get("term_name", "") if isinstance(exp.get("biosample_ontology"), dict) else "",
                    "status": exp.get("status", ""),
                    "lab": exp.get("lab", {}).get("title", "") if isinstance(exp.get("lab"), dict) else str(exp.get("lab", "")),
                    "date_released": exp.get("date_released", ""),
                    "files_count": len(exp.get("files", [])),
                })

            curated = {
                "query": query,
                "assay_filter": assay or None,
                "organism_filter": organism or None,
                "total_found": data.get("total", len(experiments)),
                "returned": len(experiments),
                "experiments": experiments,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 25: Primer Design ────────────────────────────────────

        elif name == "primer_design":
            try:
                import primer3
            except ImportError:
                return [TextContent(type="text", text="Error: primer3-py is not installed. Install with: pip install primer3-py")]

            seq_input = arguments.get("sequence", "").strip()
            if not seq_input:
                return [TextContent(type="text", text="Error: 'sequence' is required (nucleotide sequence, at least 100 bp)")]

            # Extract sequence
            nuc_str = _extract_sequence(seq_input)
            if len(nuc_str) < 100:
                return [TextContent(type="text", text=f"Error: sequence is too short ({len(nuc_str)} bp). Minimum 100 bp required for primer design.")]

            target_length = arguments.get("target_length", 200)
            target_start = arguments.get("target_start", None)
            if target_start is None:
                target_start = max(0, (len(nuc_str) // 2) - (target_length // 2))
            target_start = max(0, min(target_start, len(nuc_str) - target_length))

            num_primers = arguments.get("num_primers", 5)
            primer_min = arguments.get("primer_min_size", 18)
            primer_opt = arguments.get("primer_opt_size", 20)
            primer_max = arguments.get("primer_max_size", 25)
            tm_min = arguments.get("primer_min_tm", 57.0)
            tm_opt = arguments.get("primer_opt_tm", 60.0)
            tm_max = arguments.get("primer_max_tm", 63.0)

            design_result = primer3.bindings.design_primers(
                {
                    "SEQUENCE_TEMPLATE": nuc_str,
                    "SEQUENCE_TARGET": [target_start, target_length],
                },
                {
                    "PRIMER_NUM_RETURN": num_primers,
                    "PRIMER_MIN_SIZE": primer_min,
                    "PRIMER_OPT_SIZE": primer_opt,
                    "PRIMER_MAX_SIZE": primer_max,
                    "PRIMER_MIN_TM": tm_min,
                    "PRIMER_OPT_TM": tm_opt,
                    "PRIMER_MAX_TM": tm_max,
                    "PRIMER_PRODUCT_SIZE_RANGE": [[100, 1000]],
                }
            )

            num_returned = design_result.get("PRIMER_PAIR_NUM_RETURNED", 0)

            primer_pairs = []
            for i in range(num_returned):
                left_seq = design_result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
                right_seq = design_result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
                left_pos = design_result.get(f"PRIMER_LEFT_{i}", [0, 0])
                right_pos = design_result.get(f"PRIMER_RIGHT_{i}", [0, 0])
                primer_pairs.append({
                    "pair_index": i,
                    "left_primer": {
                        "sequence": left_seq,
                        "start": left_pos[0],
                        "length": left_pos[1],
                        "tm": round(design_result.get(f"PRIMER_LEFT_{i}_TM", 0), 1),
                        "gc_percent": round(design_result.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0), 1),
                    },
                    "right_primer": {
                        "sequence": right_seq,
                        "start": right_pos[0],
                        "length": right_pos[1],
                        "tm": round(design_result.get(f"PRIMER_RIGHT_{i}_TM", 0), 1),
                        "gc_percent": round(design_result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0), 1),
                    },
                    "product_size": design_result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0),
                    "pair_penalty": round(design_result.get(f"PRIMER_PAIR_{i}_PENALTY", 0), 3),
                })

            curated = {
                "sequence_length": len(nuc_str),
                "target_region": {"start": target_start, "length": target_length},
                "primers_found": num_returned,
                "primer_pairs": primer_pairs,
            }
            if num_returned == 0:
                curated["message"] = "No suitable primer pairs found. Try adjusting parameters (Tm range, primer size, target region)."

            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 26: PharmGKB Pharmacogenomics ────────────────────────

        elif name == "pharmgkb_search":
            gene = arguments.get("gene", "").strip()
            drug = arguments.get("drug", "").strip()
            variant = arguments.get("variant", "").strip()
            limit = arguments.get("limit", 15)

            if not gene and not drug and not variant:
                return [TextContent(type="text", text="Error: provide at least one of 'gene', 'drug', or 'variant'")]

            # Build query for clinical annotations endpoint
            params_dict: dict[str, str] = {"view": "min"}
            if gene:
                params_dict["location.genes.symbol"] = gene
            if variant:
                params_dict["location.variant.name"] = variant

            url = f"{PHARMGKB_BASE_URL}/data/clinicalAnnotation"
            result = await _rest_get(url, params=params_dict, headers={"Accept": "application/json"})

            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            all_annotations = result["data"].get("data", [])

            # Client-side drug filter (API doesn't support direct drug name filter well)
            if drug:
                drug_lower = drug.lower()
                all_annotations = [a for a in all_annotations if any(
                    drug_lower in c.get("name", "").lower() for c in a.get("relatedChemicals", [])
                )]

            annotations = all_annotations[:limit]

            curated_annotations = []
            for ann in annotations:
                location = ann.get("location", {})
                genes_list = [g.get("symbol", "") for g in location.get("genes", [])]
                variant_obj = location.get("variant", {})
                chemicals = [c.get("name", "") for c in ann.get("relatedChemicals", [])]
                evidence = ann.get("levelOfEvidence", {})

                curated_annotations.append({
                    "id": ann.get("accessionId", ""),
                    "name": ann.get("name", ""),
                    "genes": genes_list,
                    "variant": variant_obj.get("name", ""),
                    "chemicals": chemicals,
                    "types": ann.get("types", []),
                    "level_of_evidence": evidence.get("term", "") if isinstance(evidence, dict) else "",
                })

            curated = {
                "gene_filter": gene or None,
                "drug_filter": drug or None,
                "variant_filter": variant or None,
                "total_found": len(all_annotations),
                "returned": len(curated_annotations),
                "clinical_annotations": curated_annotations,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 27: Disease Ontology ─────────────────────────────────

        elif name == "disease_ontology_search":
            query = arguments.get("query", "").strip()
            limit = arguments.get("limit", 10)

            if not query:
                return [TextContent(type="text", text="Error: 'query' is required (disease name, DOID, or keyword)")]

            # Check if query is a DOID — direct lookup
            if query.upper().startswith("DOID:"):
                url = f"{DISEASE_ONTOLOGY_BASE_URL}/term/{query}"
                result = await _rest_get(url, headers={"Accept": "application/json"})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error: {result['error']}")]
                term = result["data"]
                curated = {
                    "query": query,
                    "total_results": 1,
                    "results": [{
                        "doid": term.get("id", query),
                        "name": term.get("name", ""),
                        "definition": term.get("definition", ""),
                        "synonyms": term.get("synonyms", []),
                        "xrefs": term.get("xrefs", []),
                        "parents": [{"id": p.get("id", ""), "name": p.get("name", "")} for p in term.get("parents", [])] if term.get("parents") else [],
                    }],
                }
                return [TextContent(type="text", text=json.dumps(curated, indent=2))]

            # Text search
            url = f"{DISEASE_ONTOLOGY_BASE_URL}/search"
            params_dict2: dict[str, str] = {"query": query}
            result = await _rest_get(url, params=params_dict2, headers={"Accept": "application/json"})
            if not result["success"]:
                return [TextContent(type="text", text=f"Error: {result['error']}")]

            data = result["data"]
            # Response may be a list or an object with results
            results_list = data if isinstance(data, list) else data.get("results", data.get("terms", []))
            if not isinstance(results_list, list):
                results_list = [data] if data else []

            results_list = results_list[:limit]

            terms = []
            for t in results_list:
                terms.append({
                    "doid": t.get("id", t.get("doid", "")),
                    "name": t.get("name", ""),
                    "definition": t.get("definition", ""),
                    "synonyms": t.get("synonyms", []),
                    "xrefs": t.get("xrefs", []),
                })

            curated = {
                "query": query,
                "total_results": len(terms),
                "results": terms,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Section 28: Paper Retrieval / Citation Tools ─────────────────

        elif name == "paper_fetch":
            import xml.etree.ElementTree as ET

            identifier = arguments.get("identifier", "").strip()
            sections_req = arguments.get("sections", ["methods", "results"])

            if not identifier:
                return [TextContent(type="text", text="Error: 'identifier' is required (PMID, DOI, PMC ID, or article title)")]

            # Detect identifier type
            pmid: str | None = None
            pmc_id: str | None = None
            doi: str | None = None

            stripped = identifier.strip()
            if stripped.isdigit():
                pmid = stripped
            elif stripped.upper().startswith("PMC"):
                pmc_id = stripped.upper()
            elif "/" in stripped or stripped.startswith("10."):
                doi = stripped
            else:
                # Title search — use PubMed esearch
                search_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
                search_result = await _rest_get(search_url, params={
                    "db": "pubmed",
                    "term": stripped,
                    "retmax": "1",
                    "retmode": "json",
                })
                if search_result["success"]:
                    id_list = search_result["data"].get("esearchresult", {}).get("idlist", [])
                    if id_list:
                        pmid = id_list[0]
                if not pmid:
                    return [TextContent(type="text", text=json.dumps({
                        "error": "Could not find a paper matching this title",
                        "identifier": identifier,
                    }, indent=2))]

            # Resolve DOI to PMID via PubMed esearch
            if doi and not pmid:
                search_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
                search_result = await _rest_get(search_url, params={
                    "db": "pubmed",
                    "term": f"{doi}[doi]",
                    "retmax": "1",
                    "retmode": "json",
                })
                if search_result["success"]:
                    id_list = search_result["data"].get("esearchresult", {}).get("idlist", [])
                    if id_list:
                        pmid = id_list[0]

            # Get metadata from PubMed esummary (if we have a PMID)
            title = ""
            authors_list: list[str] = []
            journal = ""
            pub_date = ""
            abstract = ""

            if pmid:
                summary_url = f"{EUTILS_BASE_URL}/esummary.fcgi"
                summary_result = await _rest_get(summary_url, params={
                    "db": "pubmed",
                    "id": pmid,
                    "retmode": "json",
                }, use_cache=False)
                if summary_result["success"]:
                    result_data = summary_result["data"].get("result", {})
                    rec = result_data.get(pmid, {})
                    title = rec.get("title", "")
                    authors_raw = rec.get("authors", [])
                    authors_list = [a.get("name", "") for a in authors_raw[:10]]
                    if len(authors_raw) > 10:
                        authors_list.append(f"... (+{len(authors_raw) - 10} more)")
                    journal = rec.get("fulljournalname", rec.get("source", ""))
                    pub_date = rec.get("pubdate", "")
                    elocation = rec.get("elocationid", "")
                    if not doi and elocation and elocation.startswith("doi:"):
                        doi = elocation.replace("doi: ", "").replace("doi:", "")
                    # Check for PMC ID in article IDs
                    for artid in rec.get("articleids", []):
                        if artid.get("idtype") == "pmc" and not pmc_id:
                            pmc_val = artid.get("value", "")
                            if pmc_val:
                                pmc_id = pmc_val if pmc_val.upper().startswith("PMC") else f"PMC{pmc_val}"

            # If we have PMC ID from a pmc_id-only input, try to get PMID
            if pmc_id and not pmid:
                search_url = f"{EUTILS_BASE_URL}/esearch.fcgi"
                search_result = await _rest_get(search_url, params={
                    "db": "pubmed",
                    "term": f"{pmc_id}[pmc]",
                    "retmax": "1",
                    "retmode": "json",
                })
                if search_result["success"]:
                    id_list = search_result["data"].get("esearchresult", {}).get("idlist", [])
                    if id_list:
                        pmid = id_list[0]
                        # Re-fetch metadata
                        summary_url = f"{EUTILS_BASE_URL}/esummary.fcgi"
                        summary_result = await _rest_get(summary_url, params={
                            "db": "pubmed",
                            "id": pmid,
                            "retmode": "json",
                        }, use_cache=False)
                        if summary_result["success"]:
                            result_data = summary_result["data"].get("result", {})
                            rec = result_data.get(pmid, {})
                            title = rec.get("title", "")
                            authors_raw = rec.get("authors", [])
                            authors_list = [a.get("name", "") for a in authors_raw[:10]]
                            journal = rec.get("fulljournalname", rec.get("source", ""))
                            pub_date = rec.get("pubdate", "")

            # Fetch abstract via efetch
            if pmid:
                try:
                    async with httpx.AsyncClient(timeout=30) as client:
                        fetch_resp = await client.get(f"{EUTILS_BASE_URL}/efetch.fcgi", params={
                            "db": "pubmed",
                            "id": pmid,
                            "rettype": "abstract",
                            "retmode": "xml",
                        })
                        if fetch_resp.status_code == 200:
                            try:
                                root = ET.fromstring(fetch_resp.text)
                                for article in root.findall(".//PubmedArticle"):
                                    abstract_elems = article.findall(".//AbstractText")
                                    if abstract_elems:
                                        parts = []
                                        for ae in abstract_elems:
                                            label = ae.get("Label", "")
                                            text_val = ae.text or ""
                                            if label:
                                                parts.append(f"{label}: {text_val}")
                                            else:
                                                parts.append(text_val)
                                        abstract = " ".join(parts)
                            except ET.ParseError:
                                pass
                except Exception:
                    pass

            # Attempt full-text retrieval
            full_text_source = "abstract_only"
            sections_extracted: dict[str, str] = {}
            section_titles_found: list[str] = []

            # Try PMC Open Access first
            if pmc_id:
                pmc_num = pmc_id.replace("PMC", "")
                try:
                    async with httpx.AsyncClient(timeout=60) as client:
                        pmc_resp = await client.get(f"{EUTILS_BASE_URL}/efetch.fcgi", params={
                            "db": "pmc",
                            "id": pmc_num,
                            "rettype": "full",
                            "retmode": "xml",
                        })
                        if pmc_resp.status_code == 200 and "<body>" in pmc_resp.text:
                            parsed = _extract_pmc_sections(pmc_resp.text, sections_req)
                            if parsed["sections"]:
                                full_text_source = "pmc"
                                sections_extracted = parsed["sections"]
                                section_titles_found = parsed["section_titles_found"]
                except Exception:
                    pass

            # Fallback: Europe PMC
            if full_text_source == "abstract_only" and pmid:
                try:
                    async with httpx.AsyncClient(timeout=60) as client:
                        epmc_resp = await client.get(
                            f"{EUROPE_PMC_BASE_URL}/{pmid}/fullTextXML"
                        )
                        if epmc_resp.status_code == 200 and "<body>" in epmc_resp.text:
                            parsed = _extract_pmc_sections(epmc_resp.text, sections_req)
                            if parsed["sections"]:
                                full_text_source = "europe_pmc"
                                sections_extracted = parsed["sections"]
                                section_titles_found = parsed["section_titles_found"]
                except Exception:
                    pass

            curated = {
                "pmid": pmid,
                "pmc_id": pmc_id,
                "doi": doi,
                "title": title,
                "authors": authors_list,
                "journal": journal,
                "pub_date": pub_date,
                "abstract": abstract,
                "full_text_available": full_text_source != "abstract_only",
                "full_text_source": full_text_source,
                "sections": sections_extracted,
                "section_titles_found": section_titles_found,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        elif name == "semantic_scholar_search":
            query = arguments.get("query", "").strip()
            limit = arguments.get("limit", 5)
            fields_of_study = arguments.get("fields_of_study", "")

            if not query:
                return [TextContent(type="text", text="Error: 'query' is required")]

            s2_fields = "paperId,title,abstract,year,authors,venue,publicationDate,citationCount,influentialCitationCount,tldr,openAccessPdf,fieldsOfStudy,externalIds"

            # Detect query type
            is_doi = query.upper().startswith("DOI:") or (query.startswith("10.") and "/" in query)
            is_pmid = query.upper().startswith("PMID:") or (query.isdigit() and len(query) >= 5)

            if is_doi:
                doi_val = query[4:] if query.upper().startswith("DOI:") else query
                url = f"{SEMANTIC_SCHOLAR_BASE_URL}/paper/DOI:{doi_val}"
                result = await _rest_get(url, params={"fields": s2_fields})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error looking up DOI: {result['error']}")]
                papers = [result["data"]]
            elif is_pmid:
                pmid_val = query[5:] if query.upper().startswith("PMID:") else query
                url = f"{SEMANTIC_SCHOLAR_BASE_URL}/paper/PMID:{pmid_val}"
                result = await _rest_get(url, params={"fields": s2_fields})
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error looking up PMID: {result['error']}")]
                papers = [result["data"]]
            else:
                # Keyword search
                url = f"{SEMANTIC_SCHOLAR_BASE_URL}/paper/search"
                params_s2: dict[str, Any] = {
                    "query": query,
                    "limit": str(limit),
                    "fields": s2_fields,
                }
                if fields_of_study:
                    params_s2["fieldsOfStudy"] = fields_of_study
                result = await _rest_get(url, params=params_s2)
                if not result["success"]:
                    return [TextContent(type="text", text=f"Error: {result['error']}")]
                papers = result["data"].get("data", [])
                if not papers:
                    return [TextContent(type="text", text=json.dumps({
                        "query": query,
                        "total_results": 0,
                        "results": [],
                        "message": "No papers found matching this query.",
                    }, indent=2))]

            curated_papers = []
            for paper in papers:
                if not paper:
                    continue
                authors_raw = paper.get("authors", [])
                author_names = [a.get("name", "") for a in (authors_raw or [])[:10]]
                if len(authors_raw or []) > 10:
                    author_names.append(f"... (+{len(authors_raw) - 10} more)")

                tldr_obj = paper.get("tldr")
                tldr_text = tldr_obj.get("text", "") if isinstance(tldr_obj, dict) else ""

                oa_pdf = paper.get("openAccessPdf")
                oa_url = oa_pdf.get("url", "") if isinstance(oa_pdf, dict) else ""

                curated_papers.append({
                    "semantic_scholar_id": paper.get("paperId", ""),
                    "title": paper.get("title", ""),
                    "authors": author_names,
                    "year": paper.get("year"),
                    "venue": paper.get("venue", ""),
                    "pub_date": paper.get("publicationDate", ""),
                    "abstract": (paper.get("abstract") or "")[:2000],
                    "tldr": tldr_text,
                    "citation_count": paper.get("citationCount", 0),
                    "influential_citation_count": paper.get("influentialCitationCount", 0),
                    "open_access_pdf_url": oa_url,
                    "fields_of_study": paper.get("fieldsOfStudy") or [],
                    "external_ids": paper.get("externalIds") or {},
                })

            curated = {
                "query": query,
                "total_results": len(curated_papers),
                "results": curated_papers,
            }
            return [TextContent(type="text", text=json.dumps(curated, indent=2))]

        # ── Lab Notebook ──────────
        elif name == "lab_notebook_annotate":
            note = arguments.get("note", "")
            session_id = arguments.get("session_id")

            notebook_dir = Path(".lab-notebook")
            notebook_dir.mkdir(exist_ok=True)

            if session_id:
                log_file = notebook_dir / f"{session_id}.jsonl"
            else:
                # Find most recent session file
                existing = sorted(notebook_dir.glob("*.jsonl"))
                if not existing:
                    return [TextContent(type="text", text=json.dumps({
                        "error": "No lab notebook session found. Start one with /lab-notebook start <title>"
                    }))]
                log_file = existing[-1]
                session_id = log_file.stem

            from datetime import datetime, timezone
            entry = {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "type": "annotation",
                "note": note,
                "session_id": session_id,
            }

            with open(log_file, "a") as f:
                f.write(json.dumps(entry) + "\n")

            return [TextContent(type="text", text=json.dumps({
                "success": True,
                "session_id": session_id,
                "note": note,
                "message": f"Annotation added to session {session_id}"
            }, indent=2))]

        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]

    except Exception as e:
        return [TextContent(type="text", text=f"Error executing tool: {str(e)}")]


async def main():
    """Main entry point for the server."""
    from mcp.server.stdio import stdio_server

    async with stdio_server() as (read_stream, write_stream):
        await app.run(
            read_stream,
            write_stream,
            app.create_initialization_options()
        )


def run_main():
    """Synchronous entry point for console script."""
    asyncio.run(main())


if __name__ == "__main__":
    run_main()
