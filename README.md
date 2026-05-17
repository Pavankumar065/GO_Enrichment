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
