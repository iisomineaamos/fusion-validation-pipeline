# fusion-validation-pipeline
=======
Long-Read Fusion Validation Pipeline
A modular R-based pipeline for validating gene fusions in long-read RNA-seq datasets using supplementary alignments, split-reads, soft-clipped re-alignment, and full-length fusion spanning detection.
[Flowchart diagram placeholder]

📂 Overview
This pipeline performs multi-stage validation of gene fusions detected in long-read RNA sequencing data using BAM files and a master fusion list. It uses samtools, minimap2, and R packages such as tidyverse, optparse, and IRanges.

Validation Stages
1. Supplementary Alignment Detection: Find supporting alignments flagged as supplementary near fusion breakpoints.
2. Split-Read Detection: Detect reads with parts mapped to both 5′ and 3′ fusion regions.
3. Soft-Clipped Read Extraction (New):
Extracts reads with soft clips (S in CIGAR) near fusion breakpoints
Generates FASTA file from clipped portion
3. Soft-Clip Re-alignment: Re-align soft-clipped reads across fusion breakpoints using minimap2.
4. Full-Length Fusion Span Detection: Identify long reads that span both fusion partners completely.

📁 Folder Structure

fusion-validation/
├── scripts/
│   ├── stage1_supplementary_alignments.R
│   ├── stage2_split_reads.R
│   ├── stage3_softclip_realignment.R
│   ├── stage4_fullspan_reads.R
├── test_data/
│   ├── fusion_master.tsv
│   ├── A549.bam
│   ├── Homo_sapiens.GRCh38.91.gtf
├── output/
│   ├── supplementary_alignments.tsv
│   ├── split_read_support.tsv
│   ├── softclip_realignments_filtered.tsv
│   ├── full_length_fusion_reads.tsv
├── environment.yml
└── README.md

Setup
Install conda environment:
conda env create -f environment.yml
conda activate fusion-validation

Usage
# Run all stages
Rscript scripts/fusion_validation_longread_advanced.R \
  -i test_data/fusion_input.tsv \
  -r reference/hg38.fa

# Run only Stage 3a (Soft-clip extraction)
Rscript scripts/fusion_validation_longread_advanced.R \
  -i test_data/fusion_input.tsv \
  --run_stage stage3a

# Run only Stage 3b (Soft-clip re-alignment)
Rscript scripts/fusion_validation_longread_advanced.R \
  -i test_data/fusion_input.tsv \
  -r reference/hg38.fa \
  --run_stage stage3b

Features
- Long-read specific enhancements
- Handles large BAMs via samtools region fetches
- Modular design (run stages independently or in sequence)
- Human-readable summaries for downstream scoring/ML

Future Improvements
- Nextflow integration
- Docker/Singularity containerization
- Integration with ML scoring modules
- Interactive HTML reports (via RMarkdown)

Citation
If you use this pipeline, please cite:
Your Name (2025). Long-Read Fusion Validation Pipeline. [GitHub Repo URL]

Requirements
- samtools ≥ 1.10
- minimap2
- R ≥ 4.1
- Packages: tidyverse, optparse, IRanges
=======
# fusion-validation-pipeline
