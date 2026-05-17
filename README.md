# 🧬 Promoter Region Extraction Pipeline
### TSS-anchored 500 bp windows for motif & GO enrichment analysis

---

## What This Does

Builds strand-aware **500 bp promoter intervals** for every human gene in the hg38 annotation, anchored at each gene's Transcription Start Site (TSS). Output is a BED file ready for motif scanning and Gene Ontology enrichment tools.

```
Gene annotation (TSV.gz)  ──▶  TSS coordinates (BED)  ──▶  500 bp promoter windows (BED)
```

---

## Output Files

| File | Description |
|---|---|
| `genes_tss_final.bed` | Quality-filtered, single-nucleotide TSS positions |
| `promoters_500bp.bed` | Final 500 bp upstream promoter windows (strand-corrected) |

---

## Prerequisites

### 1 · Create the Conda environment

```bash
mamba create -n go_enrichment python=3.12
mamba activate go_enrichment
mamba install -c bioconda bedtools samtools emboss
```

### 2 · Obtain input files

| File | Source |
|---|---|
| `hg38.fa` | [UCSC hg38 bigZips](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) |
| `human_gene_annotation.tsv.gz` | Provided separately |

---

## Pipeline Steps

### Step 1 — Chromosome size reference

Index the genome FASTA and extract chromosome lengths:

```bash
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.genome
```

---

### Step 2 — Convert annotation to TSS BED

Parse the gzipped annotation into a 6-column BED file. Each row is a single-nucleotide TSS:

| # | Column | Content |
|---|---|---|
| 1 | chrom | `chr`-prefixed chromosome name (`chrMT` → `chrM`) |
| 2 | start | TSS position |
| 3 | end | TSS + 1 |
| 4 | name | `chr@start-end\|gene` |
| 5 | score | `.` (placeholder) |
| 6 | strand | `+` or `-` |

```bash
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR>1{

    # Skip entries with undefined TSS or gene name
    if($8==-1 || $8=="" || $7=="") next

    chrom = "chr"$5
    if(chrom == "chrMT") chrom = "chrM"

    tss    = $8 + 0
    strand = ($6 == -1) ? "-" : "+"

    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand

}' > genes_tss_clean.bed
```

Then filter out any chromosomes absent from the hg38 genome file:

```bash
grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed
```

---

### Step 3 — Expand TSS → 500 bp promoter windows

Use `bedtools slop` with the `-s` flag to expand each TSS point **500 bp upstream**, respecting strand orientation:

```bash
bedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed
```

> **Why `-s` matters**
>
> | Strand | Upstream direction | Window extends toward |
> |---|---|---|
> | `+` (sense) | Lower coordinates | ← smaller coords |
> | `−` (antisense) | Higher coordinates | → larger coords |
>
> Without `-s`, upstream would always mean lower coordinates — incorrect for minus-strand genes.
>
> **Example** — reverse-strand gene on chrM:
> ```
> Input:   chrM  4400  4401  chrM@4400-4401|MT-TQ  .  -
> Output:  chrM  4400  4901  ...                    .  -
> ```
> The window expands rightward (toward higher coords) because the gene is on the minus strand.

---

## Run the Full Pipeline

Save the following as `run_pipeline.sh` and execute it:

```bash
#!/usr/bin/env bash
set -euo pipefail

mamba activate go_enrichment

# Step 1: Chromosome sizes
samtools faidx hg38.fa
cut -f1,2 hg38.fa.fai > hg38.genome

# Step 2: TSS BED
zcat human_gene_annotation.tsv.gz | \
awk 'BEGIN{OFS="\t"} NR>1{
    if($8==-1 || $8=="" || $7=="") next
    chrom = "chr"$5
    if(chrom == "chrMT") chrom = "chrM"
    tss    = $8 + 0
    strand = ($6 == -1) ? "-" : "+"
    print chrom, tss, tss+1, chrom"@"tss"-"(tss+1)"|"$7, ".", strand
}' > genes_tss_clean.bed

grep -Fwf <(cut -f1 hg38.genome) genes_tss_clean.bed > genes_tss_final.bed

# Step 3: Promoter windows
bedtools slop \
    -i genes_tss_final.bed \
    -g hg38.genome \
    -l 500 \
    -r 0 \
    -s \
    > promoters_500bp.bed

echo "✅  Done — output: promoters_500bp.bed"
```

```bash
bash run_pipeline.sh
```

---

## Result

`promoters_500bp.bed` contains a strand-corrected 500 bp upstream promoter window for every annotated human gene. It is directly compatible with motif scanning tools (e.g. MEME, FIMO) and GO enrichment pipelines.
