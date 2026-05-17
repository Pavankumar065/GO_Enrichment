# Strand-Aware Promoter Extraction Pipeline (v2.0.0)

An automated, strand-sensitive bioinformatics workflow for generating **500 bp upstream promoter intervals** from human genomic features. This pipeline accurately maps single-nucleotide Transcription Start Sites (TSS) and extends them relative to gene orientation using the `bedtools` suite, formatting the output for downstream **Gene Ontology (GO) Enrichment** and **Motif Analysis**.

---

## 📋 Table of Contents
- [Workflow Architecture](#-workflow-architecture)
- [Pipeline Deliverables](#-pipeline-deliverables)
- [Environment Initialization](#-environment-initialization)
- [Input Prerequisites](#-input-prerequisites)
- [Step-by-Step Implementation](#-step-by-step-implementation)
  - [Step 1: Chromosome Indexing](#step-1-chromosome-indexing)
  - [Step 2: Annotation Parsing & TSS Generation](#step-2-annotation-parsing--tss-generation)
  - [Step 3: Strand-Specific Extension](#step-3-strand-specific-extension)
- [Monolithic Automation Script](#-monolithic-automation-script)
- [Technical Notes & Validation](#-technical-notes--validation)

---

## 🧬 Workflow Architecture

The pipeline processes raw tabular transcript annotations and reference genome index shapes to anchor high-confidence promoter coordinates. By evaluating transcript orientation, it circumvents coordinate-shift errors typical of strand-agnostic flanking operations.
[human_gene_annotation.tsv.gz] ---> (AWK Parser) ---> [genes_tss_clean.bed]
|
[hg38.fa] ---> (samtools faidx) ---> [hg38.genome] --------> (grep Filter)
|
[genes_tss_final.bed]
|
(bedtools slop -s)
|
[promoters_500bp.bed]


---

## 📦 Pipeline Deliverables

The execution produces two foundational data structures:

| Target File | Format | Description |
| :--- | :---: | :--- |
| `genes_tss_final.bed` | BED6 | Quality-filtered, single-nucleotide resolution coordinates mapping valid TSS positions to the hg38 assembly. |
| `promoters_500bp.bed` | BED6 | Final strand-adjusted 500 bp upstream windows optimized for motif scanning matrices and enrichment tools. |

---

## 🛠️ Environment Initialization

Deploy an isolated package environment running **Python 3.12** through the `Mamba` package manager, and pull core binaries from the `Bioconda` channel:

```bash
# Create and initialize the environment
mamba create -n go_enrichment python=3.12 -y
mamba activate go_enrichment

# Install pipeline dependencies
mamba install -c bioconda bedtools samtools emboss -y


🗂️ Input Prerequisites
1. Reference Genome Assembly (hg38.fa)
Retrieve the standard GRCh38/hg38 human primary genome assembly directly from the UCSC GoldenPath archive:

URL: UCSC hg38 BigZips Repository


2. Gene Annotation Matrix (human_gene_annotation.tsv.gz)
A gzipped, tab-separated metadata file containing feature boundaries, identifiers, and orientation markers.


💻 Step-by-Step Implementation
Step 1: Chromosome Indexing
Generate a standard FASTA index (.fai) to fetch sequence structures and map lengths into a dedicated coordinate bounds file.


Bash
# Index the FASTA sequence
samtools faidx hg38.fa


# Extract identifiers and sequence lengths for bedtools compatibility
cut -f1,2 hg38.fa.fai > hg38.genome
Step 2: Annotation Parsing & TSS Generation
Filter and map raw annotations into standard 6-column BED records representing structural transcription start nodes.


BED6 Output Schema Mapping
Column 1 (Chrom): Chromosome name (chrN format; MT mapped natively to M)

Column 2 (Start): 0-based start coordinate

Column 3 (End): 1-based exclusive end coordinate, equal to start + 1

Column 4 (Name): chr@start-end|gene_symbol structural identifier

Column 5 (Score): Null placeholder string .

Column 6 (Strand): + or - structural orientation


Bash
# Decompress, validate fields, and map numeric metadata states to BED fields
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR>1{
    # Evict rows missing defined TSS markers or locus names
    if($8 == -1 || $8 == "" || $7 == "")
        next

    # Construct uniform UCSC chromosome names
    chrom = "chr"$5
    if(chrom == "chrMT")
        chrom = "chrM"

    tss = $8 + 0

    # Translate orientation integers to formal string states
    strand = "+"
    if($6 == -1)
        strand = "-"

    # Write mapped record to disk
    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand
}' > genes_tss_clean.bed


# Purge coordinates targeting orphaned or unindexed structural contigs
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed

Step 3: Strand-Specific ExtensionWiden the filtered 1 bp TSS markers into 500 bp promoter windows upstream of the transcription path. Using the -s parameter ensures calculations shift correctly around features based on the transcript direction.Bashbedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed

Strand Translation MatrixWithout passing the orientation flag (-s), coordinate shifts default strictly to lower linear numerical spaces. The pipeline tracks directionality as follows:StrandContextTarget DirectionCoordinate Modification+Sense Feature5′ Regulatory UpstreamShifts toward smaller coordinates-Antisense Feature5′ Regulatory UpstreamShifts toward larger coordinatesTrace Validation Example (Reverse/Minus Strand Locus):Input Feature: chrM  4400  4401  chrM@4400-4401|MT-TQ  .  -Calculated Upstream Interval: Because the transcript is on the minus strand, the upstream sequence occupies the larger coordinate direction. The final interval is generated seamlessly from 4400 to 4901, safely containing the 500 bp promoter window.🚀 Monolithic Automation ScriptSave this workflow as a bash script (e.g., run_pipeline.sh) to execute the sequence automatically:Bash#!/usr/bin/env bash



# Configure rigorous error tracing
set -euo pipefail


echo "[INFO] Activating target runtime environment..."
# Ensure shell profiles can access mamba hooks
if info="$(type -p mamba)" && [ -n "$info" ]; then
    source "$(dirname "$info")/../etc/profile.d/conda.sh"
fi
conda activate go_enrichment


echo "[INFO] Indexing assembly and preparing genome space structures..."
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.genome


echo "[INFO] Generating single-nucleotide TSS boundaries..."
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR>1{
    if($8 == -1 || $8 == "" || $7 == "")
        next

    chrom = "chr"$5
    if(chrom == "chrMT")
        chrom = "chrM"

    tss    = $8 + 0
    strand = "+"
    if($6 == -1)
        strand = "-"

    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand
}' > genes_tss_clean.bed

echo "[INFO] Sanitizing genomic loci names against chromosome references..."
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed

echo "[INFO] Executing strand-sensitive upstream expansion..."
bedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed

# Cleanup intermediate step files
rm genes_tss_clean.bed

echo "[SUCCESS] Execution finalized. Destination output saved to: promoters_500bp.bed"


📊 Technical Notes & ValidationAssembly Identity: Ensure the input annotation coordinates perfectly match the assembly coordinates of the hg38.fa file to avoid offset inaccuracies.Compatibility: The generated promoters_500bp.bed can be used directly with tools like MEME Suite, Homer, or converted to FASTA format using bedtools getfasta for downstream motif discovery algorithms.
