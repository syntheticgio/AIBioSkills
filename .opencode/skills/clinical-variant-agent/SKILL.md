---
name: clinical-variant-agent
description: Full clinical variant workup — gnomAD population frequency, ClinVar significance, protein domain impact, AlphaFold structure context, and PubMed literature
---

Perform a comprehensive clinical variant workup for: **$ARGUMENTS**

Use the MCP tools to gather data from all relevant sources, then synthesize a structured clinical variant report. Follow the steps below in order. If a step fails or returns no data, note the gap and continue.

## Input Parsing

Determine the variant format:
1. **rsID** (e.g., `rs11549407`) — use directly for gnomAD lookup
2. **chrom-pos-ref-alt** (e.g., `7-140753336-A-T`) — use directly for gnomAD
3. **HGVS notation** (e.g., `NM_007294.4:c.5266dupC`) — search ClinVar first to get genomic coordinates
4. **Gene + protein change** (e.g., `BRAF V600E`) — search ClinVar first, extract variant details

## Data Gathering Steps

### 1. Population Frequency (gnomAD)

- Call `gnomad_get_variant` with the variant ID or rsID.
- Extract: overall allele frequency (exome + genome), per-population frequencies, homozygote counts, filter flags.
- If the variant is not found in gnomAD, note this — absence from gnomAD is itself informative (suggests very rare).

### 2. Clinical Significance (ClinVar)

- Call `clinvar_search` with the variant identifier (rsID, HGVS, or gene + variant description).
- Extract: clinical significance classification, review status (star rating), associated conditions/diseases, submitter count.
- Note any conflicting interpretations.

### 3. Gene Context

- From the gnomAD or ClinVar result, identify the affected gene symbol.
- Call `datasets_summary_gene` with the gene symbol (taxon: human) for NCBI gene metadata.
- Call `ensembl_lookup_gene` with the gene symbol for genomic coordinates.

### 4. Protein Domain Impact

- Call `uniprot_search` with the gene symbol and `organism_id:9606` (reviewed: true) to find the canonical UniProt entry.
- Call `uniprot_get_features` with the UniProt accession to get domain boundaries, active sites, and functional regions.
- Call `interpro_get_domains` with the UniProt accession for InterPro domain annotations.
- Map the variant position to protein domains — does it fall within a known functional domain?

### 5. Structural Context (AlphaFold)

- Call `alphafold_get_prediction` with the UniProt accession.
- Check the per-residue confidence (pLDDT) at the variant position — high confidence means the structural prediction is reliable for assessing impact.
- Note whether the variant falls in a structurally ordered or disordered region.

### 6. Phenotype Context (HPO)

- If ClinVar returned associated conditions, call `hpo_search` with the condition name to find related HPO phenotype terms.
- Alternatively, call `hpo_search` with the gene symbol (category: search) to find phenotypes associated with the gene.

### 7. Literature (PubMed)

- Call `pubmed_search` with the variant identifier (e.g., "BRAF V600E") and/or the gene symbol + "variant" to find relevant publications.
- Focus on papers discussing the clinical significance, functional impact, or therapeutic relevance of this variant.

## Report Format

```
# Clinical Variant Report: [VARIANT]

## Variant Identity
- Genomic coordinates (GRCh37/GRCh38)
- rsID(s)
- HGVS notation (genomic, coding, protein)
- Gene affected, transcript
- Consequence type (missense, nonsense, frameshift, etc.)

## Population Frequency
- gnomAD overall allele frequency (exome and genome)
- Per-population breakdown table:
  | Population | Allele Count | Allele Number | Frequency | Homozygotes |
  |------------|-------------|---------------|-----------|-------------|
- Interpretation: common (>1%), low-frequency (0.1-1%), rare (<0.1%), ultra-rare (absent)
- Filter flags (if any)

## Clinical Significance
- ClinVar classification and review status
- Associated conditions/diseases
- Number of submitters and any conflicting interpretations
- Date of last review (if available)

## Protein Impact Assessment
- Affected protein domain(s) with positions
- Is the variant in a known functional region (active site, binding site, etc.)?
- Conservation context (is this residue typically conserved?)
- Predicted structural impact based on AlphaFold confidence at this position

## Phenotype Associations
- HPO terms associated with the gene or condition
- Related diseases from HPO

## Relevant Literature
- Key publications discussing this variant (top 3-5)
- Brief summary of findings from the literature

## ACMG Classification Considerations
Based on the evidence gathered above, assess which ACMG/AMP criteria may apply:
- **Population data** (BA1/BS1/PM2): Is the variant too common to be pathogenic, or absent from controls?
- **Computational/predictive** (PP3/BP4): What do consequence predictions suggest?
- **Functional domain** (PM1): Does it fall in a critical domain?
- **Clinical data** (PP5/BP6): What does ClinVar say?
- **Literature** (PS3/BS3): Is there functional evidence?

Note: This is an evidence summary to assist interpretation, NOT a definitive classification. Clinical variant interpretation should always be performed by qualified professionals.

## Data Sources
List which databases were queried and whether each returned data successfully.
```

Keep the report factual — only include data returned by the tools. Do not fabricate allele frequencies, clinical classifications, or literature citations. If a section has no data, write "No data available from [source]."
