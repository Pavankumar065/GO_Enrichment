# Pipeline: Strand-Aware Promoter Core Extraction Engine

> Architectural bioinformatic workflow for generating high-fidelity, 500 bp promoter intervals anchored to the human (hg38) Transcription Start Site (TSS). Designed specifically to format raw coordinates for downstream Motif Discovery, Scanning, and Gene Ontology (GO) enrichment pipelines.

---

## 🛠️ Compute Environment Specification

To guarantee cross-platform reproducibility, execution environments are strictly isolated via virtualized micro-environments using `mamba`/`conda` packages distributed through the **Bioconda** channel.

```bash
# 1. Initialize structural virtual environment
mamba create --name genomic_extraction_env python=3.12 --yes
mamba activate genomic_extraction_env

# 2. Deploy binary utilities via Bioconda channel
mamba install --channel bioconda bedtools samtools emboss --yes
💾 Core Pipeline Assets & ManifestRequired Data Ingestion Arrayshg38.fa: Complete reference sequence assemblies for Homo sapiens (UCSC Genome Browser Distribution).human_gene_annotation.tsv.gz: Raw genomic annotation index matrix containing explicit gene coordinates, strand orientation factors, and assigned TSS offsets.Pipeline Deliverables Summarypipeline_workspace/
├── hg38.genome           <- Calculated chromosome-level coordinate boundary map
├── genes_tss_final.bed   <- Filtered single-nucleotide TSS seed coordinates (BED6 format)
└── promoters_500bp.bed   <- Final 500 bp strand-corrected promoter intervals
📐 Algorithmic Pipeline Architecture[human_gene_annotation.tsv.gz] ---> (awk coordinate parsing engine)
                                                 │
                                                 ▼
[hg38.fa] ---> [samtools faidx] ---> [genes_tss_final.bed] (Valid Chroms Only)
                                                 │
                                                 ▼
                                     [bedtools slop -s -l 500 -r 0]
                                                 │
                                                 ▼
                                      [promoters_500bp.bed]
Stage 1: Structural Genome Space IndexingBefore manipulating interval spans, we create an absolute coordinate boundary engine using the genome index .fai. This prevents structural out-of-bounds overflow exceptions during later interval expansions.Bash# Compile structural binary index mapping file
samtools faidx hg38.fa

# Isolate sequence identity fields and max coordinate spaces
cut -f1,2 hg38.fa.fai > hg38.genome
Stage 2: Single-Nucleotide TSS VectorizationThis process parses raw genomic matrix elements, scales numeric representations, fixes mitochondrial syntax structural naming conflicts, and creates a strict BED6-compliant data structure.Bashzcat human_gene_annotation.tsv.gz | awk '
BEGIN { 
    OFS = "\t" 
}
NR > 1 {
    # Data Sanity Validation: Drop rows containing missing target matrices
    if ($8 == -1 || $8 == "" || $7 == "") next;

    # Normalize Chromosome names to UCSC-Standard Syntax (e.g., 1 -> chr1, MT -> chrM)
    chrom = "chr"$5;
    if (chrom == "chrMT") chrom = "chrM";

    # Enforce strict float/integer conversion for coordinates
    tss_start = $8 + 0;
    tss_end   = tss_start + 1;

    # Map binary evaluation coordinates to standard string characters (+ / -)
    strand = ($6 == -1) ? "-" : "+";

    # Compile explicit unique identifiers: chr@start-end|GeneSymbol
    record_id = chrom "@" tss_start "-" tss_end "|" $7;

    # Output standard 6-column BED schema
    print chrom, tss_start, tss_end, record_id, ".", strand;
}' > genes_tss_clean.bed

# Boundary Filter Checkpoint: Discard scaffold fragments not represented inside hg38.genome
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed
Output BED6 Layout Specifications:Output BED6 ColumnData TypeInternal Attribute Value MappingCol 1: ChromosomeStringUCSC nomenclature format prefix (chr1 ... chrX, chrM)Col 2: StartIntegerExact genomic position of the Transcription Start SiteCol 3: EndIntegerClosed-interval tracking coordinate limit ($Start + 1$)Col 4: NameStringPrimary lookup index key compound identifier stringCol 5: ScoreCharacterUnutilized structural filler attribute (.)Col 6: StrandCharacterDirectional feature orientation tracker (+ or -)Stage 3: Strand-Sensitive Directional Window ScalingTo isolate the upstream promoter environment, calculations must change behavior based on the strand orientation flag. Standard coordinate expansion expands coordinates symmetrically or blind to direction. Here, we use directional parameters via bedtools slop.Bashbedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed
Directional Math Transformations Applied:The -s runtime parameter shifts coordinate transformations using the following mathematical logic:Sense Orientation ($+$ Strand): The upstream direction moves left toward lower index coordinates.$$\text{New Start} = \max(0, \text{Original Start} - 500)$$$$\text{New End} = \text{Original End}$$Antisense Orientation ($-$ Strand): The upstream direction moves right toward higher index coordinates.$$\text{New Start} = \text{Original Start}$$$$\text{New End} = \min(\text{Chromosome Limit}, \text{Original End} + 500)$$🚀 Unified Core Bash Script Pipeline (run_pipeline.sh)This self-contained executable contains safety locks (set -euo pipefail) to ensure immediate process termination if any intermediate pipeline subprocess returns an error code.Bash#!/usr/bin/env bash
set -euo pipefail

# Confirm presence of runtime data inputs
for asset in hg38.fa human_gene_annotation.tsv.gz; do
    if [[ ! -f "$asset" ]]; then
        echo "[FATAL ERROR]: Missing critical genomic asset dependency -> $asset" >&2
        exit 1
    fi
done

echo "[STAGE 1/3]: Compiling Reference Genome Spatial Coordinates..."
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.genome

echo "[STAGE 2/3]: Converting Transcriptome Annotations to Strict BED6 Formats..."
zcat human_gene_annotation.tsv.gz | awk '
BEGIN { OFS="\t" } 
NR>1 {
    if($8==-1 || $8=="" || $7=="") next;
    chrom = "chr"$5;
    if(chrom == "chrMT") chrom = "chrM";
    tss = $8 + 0;
    strand = ($6 == -1) ? "-" : "+";
    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand;
}' > genes_tss_clean.bed

# Drop out-of-bounds references
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed
rm genes_tss_clean.bed

echo "[STAGE 3/3]: Running Strand-Aware Promoter Window Calculations..."
bedtools slop -i genes_tss_final.bed -g hg38.genome -l 500 -r 0 -s > promoters_500bp.bed

echo "[SUCCESS]: Promoter Core Matrix exported safely to -> promoters_500bp.bed"
