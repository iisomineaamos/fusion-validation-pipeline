# Transcript-Centric Gene Fusion Validation Pipeline.

This repository contains the official R-based pipeline for the paper: "A WGS-Independent Framework for Gene Fusion Validation in Long-Read RNA-Seq Data."

## 📂 Overview

This is a modular, automated pipeline designed to validate gene fusion candidates from long-read RNA-Seq data. Its key contribution is that it operates **without requiring matched Whole-Genome Sequencing (WGS) data**.

The pipeline works in two stages:
1.  **Evidence Gathering:** It systematically quantifies three key types of transcript-centric evidence from a coordinate-sorted BAM file:
    * Full-Length Spanning Reads
    * Supplementary Alignments
    * Realigned Soft-Clipped Reads
2.  **ML Classification:** It uses the quantified evidence as features for a pre-trained Random Forest model to generate a final, high-confidence classification (TRUE vs. FALSE) for each candidate.

---

## 🚀 Reproducibility with Singularity

This pipeline is designed to be fully reproducible. The recommended method for running the pipeline is by using the provided Singularity container, which includes all required software dependencies (R, Samtools, Minimap2, and all R packages).

The recipe file (`singularity/fusion_pipeline.def`) is included in this repository.

### 1. Build the Container

You must have Singularity installed on your Linux system.

```bash
# Clone this repository
git clone [https://github.com/iisomineaamos/fusion-validation-pipeline.git](https://github.com/iisomineaamos/fusion-validation-pipeline.git)
cd fusion-validation-pipeline

# Build the container (requires sudo/root privileges)
sudo singularity build singularity/fusion_pipeline.sif singularity/fusion_pipeline.def

This will create a single container file named fusion_pipeline.sif inside the singularity/ directory.

### 2. Run the Pipeline

To run the pipeline, use the singularity exec command. You must use the --bind (-B) flag to mount your data directory (containing BAMs, FASTA, etc.) and your code directory (this repository) into the container's environment.
# --- 1. Define your paths ---

# (e.g., /mnt/projects/my_project/data)
DATA_DIR=/path/to/your/data

# (e.g., /home/user/code/fusion-validation-pipeline)
PIPELINE_REPO_DIR=/path/to/where/you/cloned/this/repo

singularity exec \
    -B ${DATA_DIR}:/data \
    -B ${PIPELINE_REPO_DIR}:/pipeline \
    singularity/fusion_pipeline.sif \
    Rscript /pipeline/scripts/fusion_validation_pipeline.R \
        --input_file /data/your_fusion_candidates.tsv \
        --bam_file /data/your_alignment.bam \
        --ref_genome /data/your_genome.fa

Explanation of the Command:
--bind /path/to/your/data:/data: This makes your local data folder available inside the container at the path /data.

--bind /path/to/your/code:/pipeline: This makes the cloned GitHub folder available inside the container at the path /pipeline.

singularity/fusion_pipeline.sif: This is the container file you built.

Rscript /pipeline/scripts/fusion_validation_pipeline.R: This tells R to run your pipeline script from its path inside the container.

--input_file /data/your_fusion_candidates.tsv: This is critical. All paths to your data files must now use the container's internal path (e.g., /data/your_file.bam).

📁 Repository Structure

fusion-validation-pipeline/
├── scripts/
│   └── fusion_validation_pipeline.R  # The main pipeline script
├── singularity/
│   └── fusion_pipeline.def           # The Singularity recipe file
├── test_data/
│   ├── fusion_master.tsv             # Example input
│   └── ...
├── output/                           # Directory created by the pipeline
│   └── ...
└── README.md

Features
WGS-Independent: Validates fusions using only long-read RNA-Seq data.

Multi-Modal Evidence: Quantifies full-length, supplementary, and soft-clip support.

ML-Powered: Integrates a Random Forest model for robust classification.
roducible: Bundled in a Singularity container for ease of use.

