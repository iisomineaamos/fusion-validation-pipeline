#!/usr/bin/env Rscript
# ================================
# Fusion Detection Pipeline (Stages 1-4 Combined)
# ================================
library(tidyverse)
library(optparse)
library(stringr)
library(glue)
library(tibble)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# ---- CLI arguments ----
# Define options
option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", help = "Path to fusion TSV"),
  make_option(c("-b", "--bam_file"), type = "character", default = "", help = "Path to BAM file (optional override)"),
  make_option(c("-f", "--fusions"), type = "character", help = "Fusion TSV file path"),
  make_option(c("-r", "--ref_genome"), type = "character", help = "Path to reference genome FASTA"),
  make_option(c("-m", "--minimap2"), type = "character", default = "minimap2", help = "Path to minimap2"),
  make_option(c("-w", "--window"), type = "integer", default = 1000, help = "Window around fusion breakpoint"),
  make_option(c("-q", "--mapq_cutoff"), type = "integer", default = 20, help = "Minimum MAPQ"),
  make_option(c("-s", "--softclip_threshold"), type = "integer", default = 20, help = "Minimum soft-clip length to consider (bp) [default %default]", metavar = "number")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# ✅ ADD THIS CHECK after parsing
if (is.null(opt$fusions)) {
  stop("❌ Fusion TSV file not provided. Use -f to specify the path.")
}

# Assign options to variables
fusion_file <- opt$input_file
bam_override <- opt$bam_file
ref_genome <- opt$ref_genome
minimap2 <- opt$minimap2
window <- opt$window
mapq_cutoff <- opt$mapq_cutoff
softclip_threshold <- opt$softclip_threshold  # New parameter

# ---- Generate output prefix from input filename ----
output_prefix <- tools::file_path_sans_ext(basename(fusion_file))

# ---- Parameters ----
softclip_window <- 1000
sam_fields <- c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","ISIZE","SEQ","QUAL","OPT")

# ================================
# Stage 1: Supplementary Alignments
# ================================
print("🔍 Stage 1: Extracting supplementary alignments...")

# ---- Load input TSV ----
fusion_input <- read_tsv(fusion_file, show_col_types = FALSE)

# ---- Detect if BAM uses chr-prefixed references ----
bam_to_check <- if (bam_override != "") bam_override else fusion_input$path[1]
bam_header <- system(paste("samtools view -H", bam_to_check, "| grep '^@SQ'"), intern = TRUE)
bam_uses_chr <- any(grepl("SN:chr", bam_header))

# ---- Apply chr-fix accordingly ----
fusion_input <- fusion_input %>%
  mutate(
    fiveprime_chr = if (bam_uses_chr) {
      ifelse(!grepl("^chr", fiveprime_chr), paste0("chr", fiveprime_chr), fiveprime_chr)
    } else {
      gsub("^chr", "", fiveprime_chr)
    },
    threeprime_chr = if (bam_uses_chr) {
      ifelse(!grepl("^chr", threeprime_chr), paste0("chr", threeprime_chr), threeprime_chr)
    } else {
      gsub("^chr", "", threeprime_chr)
    }
  )

# ---- Build regions and extract supplementary alignments ----
dir.create("output", showWarnings = FALSE)
supp_aligns <- fusion_input %>%
  mutate(region_5p = paste0(fiveprime_chr, ":", fiveprime_search_start - softclip_window, "-", fiveprime_search_end + softclip_window),
         region_3p = paste0(threeprime_chr, ":", threeprime_search_start - softclip_window, "-", threeprime_search_end + softclip_window)) %>%
  pivot_longer(cols = c(region_5p, region_3p), names_to = "region_type", values_to = "region") %>%
  mutate(bam = ifelse(bam_override != "", bam_override, path)) %>%
  mutate(cmd = paste0("samtools view -f 2048 ", bam, " ", region)) %>%
  mutate(results = map(cmd, ~system2("bash", c("-c", .x), stdout = TRUE))) %>%
  select(fusion_id, sample_id, path, region_type, region, results) %>%
  unnest(results) %>%
  separate(results, into = sam_fields, sep = "\t", extra = "merge") %>%
  filter(!is.na(POS))

supp_output <- paste0("output/", output_prefix, "_supplementary_alignments.tsv")
write_tsv(supp_aligns, supp_output)
print(paste("✅ Stage 1 complete. Output:", supp_output))


# Stage 2: Extract all supplementary alignments once
message("🔎 Extracting supplementary alignments...")

# Get all supplementary alignments into memory
supp_reads <- system2("samtools", c("view", "-f", "2048", opt$bam), stdout = TRUE)

# Parse the results
supp_df <- tibble::tibble(raw = supp_reads) %>%
  separate(raw, into = c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", 
                         "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "TAGS"), 
           sep = "\t", fill = "right", extra = "merge")

# Ensure numeric
supp_df <- supp_df %>% mutate(POS = as.numeric(POS))

# Match to fusion regions (allow ± window)
supp_hits <- fusion_input %>%
  rowwise() %>%
  mutate(
    matched_reads_5p = list(
      supp_df %>%
        filter(RNAME == fiveprime_chr & abs(POS - fiveprime_search_start) < window)
    ),
    matched_reads_3p = list(
      supp_df %>%
        filter(RNAME == threeprime_chr & abs(POS - threeprime_search_start) < window)
    ),
    matched_reads = list(bind_rows(matched_reads_5p, matched_reads_3p))
  ) %>%
  unnest(matched_reads) %>%
  ungroup()

# Write output
write_tsv(supp_hits, file.path("output", paste0(tools::file_path_sans_ext(basename(opt$fusions)), "_supplementary_alignments.tsv")))
message("✅ Stage 2 complete.")

# ================================
# Stage 3a: Soft-Clipped Read Extraction with fusion_id tagging
# ================================

message(glue::glue("🔍 Stage 3a: Extracting soft-clipped reads (threshold: {softclip_threshold} bp)..."))

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(glue)
})

bam_file <- if (bam_override != "") bam_override else fusion_input$path[1]
if (!file.exists(paste0(bam_file, ".bai"))) {
  message("⏳ Indexing BAM file...")
  system(paste("samtools index", shQuote(bam_file)))
}

output_prefix <- tools::file_path_sans_ext(basename(fusion_file))
softclip_output <- glue("output/{output_prefix}_softclipped_reads_{softclip_threshold}bp.tsv")
softclip_fastq <- glue("output/{output_prefix}_softclipped_reads_{softclip_threshold}bp.fastq")

fusion_chrs <- unique(c(fusion_input$fiveprime_chr, fusion_input$threeprime_chr))
regions <- paste(fusion_chrs, collapse = " ")
cmd <- glue("samtools view -F 0 {shQuote(bam_file)} {regions}")
message("🔍 Running: {cmd}")

sam_lines <- tryCatch(system(cmd, intern = TRUE), error = function(e) character())

if (length(sam_lines) == 0) {
  message("⚠️ No reads found.")
  write_tsv(tibble(), softclip_output)
} else {
  softclip_df <- tibble(raw = sam_lines) %>%
    mutate(fields = str_split(raw, "\t")) %>%
    filter(map_int(fields, length) >= 11) %>%
    transmute(
      QNAME = map_chr(fields, 1),
      FLAG = map_chr(fields, 2),
      RNAME = map_chr(fields, 3),
      POS = as.integer(map_chr(fields, 4)),
      MAPQ = as.integer(map_chr(fields, 5)),
      CIGAR = map_chr(fields, 6),
      SEQ = map_chr(fields, 10),
      QUAL = map_chr(fields, 11)
    ) %>%
    mutate(
      left_softclip = as.integer(str_match(CIGAR, "^(\\d+)S")[,2]),
      right_softclip = as.integer(str_match(CIGAR, "(\\d+)S$")[,2]),
      total_softclip = coalesce(left_softclip, 0) + coalesce(right_softclip, 0)
    ) %>%
    filter(total_softclip >= softclip_threshold)

# Assign fusion_id by checking overlap with fusion_input regions
message("🔗 Assigning fusion_id based on region overlap...")

softclip_df <- softclip_df %>%
  mutate(
    fusion_id = map_chr(seq_along(POS), function(i) {
      hit <- fusion_input %>%
        filter(
          (fiveprime_chr == RNAME[i] & POS[i] >= (fiveprime_search_start - window) & POS[i] <= (fiveprime_search_end + window)) |
          (threeprime_chr == RNAME[i] & POS[i] >= (threeprime_search_start - window) & POS[i] <= (threeprime_search_end + window))
        ) %>%
        pull(fusion_id)

      if (length(hit) == 0) NA_character_ else hit[1]
    })
  ) %>%
  filter(!is.na(fusion_id))

message(glue("✅ Found {nrow(softclip_df)} soft-clipped reads ≥{softclip_threshold}bp."))

write_tsv(softclip_df, softclip_output)
message("📁 Output saved to: {softclip_output}")

  # Write FASTQ
  if (nrow(softclip_df) > 0) {
    fastq_lines <- purrr::map_chr(1:nrow(softclip_df), function(i) {
      paste0("@", softclip_df$QNAME[i], "\n", softclip_df$SEQ[i], "\n+\n", softclip_df$QUAL[i])
    })
    writeLines(fastq_lines, softclip_fastq)
    message("📁 Soft-clipped reads saved to FASTQ format: ", softclip_fastq)
  } else {
    message("⚠️ No soft-clipped reads to write to FASTQ.")
  }
}

# ================================
# Stage 3b: Realign Soft-Clipped Reads
# ================================
message("🔍 Stage 3b: Realigning soft-clipped reads...")

# Ensure that the softclip_fastq is available and not empty
if (file.exists(softclip_fastq) && file.info(softclip_fastq)$size > 0) {
  
  # Define the output for the realigned reads
  realign_output <- glue("output/{output_prefix}_softclipped_realigned.sam")
  
  # Command for minimap2 to align the soft-clipped reads to the reference genome
  align_cmd <- glue("minimap2 -a {ref_genome} {softclip_fastq} > {realign_output}")
  
  message("⏳ Aligning soft-clipped reads with minimap2...")
  
  # Run the alignment command
  system(align_cmd)
  
  message(glue("✅ Stage 3b complete. Output: {realign_output}"))
} else {
  message("⚠️ Skipping Stage 3b: No soft-clipped reads available or file is empty")
}

# ================================
# Stage 4: Full-Length Spanning Fusion Reads
# ================================
print("🔍 Stage 4: Detecting full-length spanning fusion reads...")

fullspan_reads <- list()

# ---- Loop through each fusion ----
for (i in 1:nrow(fusion_input)) {
  fusion <- fusion_input[i, ]
  fusion_id <- fusion$fusion_id
  sample_id <- fusion$sample_id
  
  region_5p <- paste0(fusion$fiveprime_chr, ":", fusion$fiveprime_search_start, "-", fusion$fiveprime_search_end)
  region_3p <- paste0(fusion$threeprime_chr, ":", fusion$threeprime_search_start, "-", fusion$threeprime_search_end)

  # Extract reads in both regions (primary only, MAPQ filtered)
  reads_5p <- system(paste0("samtools view -F 2048 -q ", mapq_cutoff, " ", bam_file, " ", region_5p), intern = TRUE)
  reads_3p <- system(paste0("samtools view -F 2048 -q ", mapq_cutoff, " ", bam_file, " ", region_3p), intern = TRUE)

  # Extract read names (QNAMEs)
  qnames_5p <- unique(str_extract(reads_5p, "^[^\\t]+"))
  qnames_3p <- unique(str_extract(reads_3p, "^[^\\t]+"))

  # Detect intersection
  spanning_reads <- intersect(qnames_5p, qnames_3p)

  # Save if found
  if (length(spanning_reads) > 0) {
    fullspan_reads[[length(fullspan_reads) + 1]] <- tibble(
      fusion_id = fusion_id,
      sample_id = sample_id,
      read_name = spanning_reads
    )
  }
}

# ---- Combine and Save ----
if (length(fullspan_reads) > 0) {
  fullspan_output <- paste0("output/", output_prefix, "_full_length_fusion_reads.tsv")
  fullspan_df <- bind_rows(fullspan_reads)
  write_tsv(fullspan_df, fullspan_output)
  print(paste("✅ Stage 4 complete. Output:", fullspan_output))
} else {
  print("⚠️ No full-length spanning reads detected above MAPQ cutoff.")
}

print("🎉 Pipeline completed successfully!")


# ================================
# Stage 5: Final Fusion Validation
# ================================
message("🔍 Stage 5: Validating fusion events using all evidence...")

# Construct paths using output_prefix
supp_file <- file.path("output", paste0(output_prefix, "_supplementary_alignments.tsv"))
softclip_file <- file.path("output", paste0(output_prefix, "_softclipped_reads_", softclip_threshold, "bp.tsv"))
fullspan_file <- file.path("output", paste0(output_prefix, "_full_length_fusion_reads.tsv"))

# Check existence
required_files <- c(supp_file, softclip_file, fullspan_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("❌ One or more required output files are missing:\n", paste(missing_files, collapse = "\n"))
}

# Load all outputs
supp_df <- read_tsv(supp_file, show_col_types = FALSE)
soft_df <- read_tsv(softclip_file, show_col_types = FALSE)
full_df <- read_tsv(fullspan_file, show_col_types = FALSE)

# If any of them are completely empty, warn and create empty frames
if (nrow(supp_df) == 0) message("⚠️ Supplementary file is empty.")
if (nrow(soft_df) == 0) message("⚠️ Soft-clipped file is empty.")
if (nrow(full_df) == 0) message("⚠️ Full-span file is empty.")

# Combine by fusion_id (adjust if different identifier used)
combined_validated <- fusion_input %>%
  mutate(
    supplementary = fusion_id %in% supp_df$fusion_id,
    softclipped   = fusion_id %in% soft_df$fusion_id,
    fullspan      = fusion_id %in% full_df$fusion_id
  ) %>%
  mutate(
    total_support = supplementary + softclipped + fullspan,
    validated = total_support > 0
  )

# Save final validated fusions
validated_output <- file.path("output", paste0(output_prefix, "_validated_fusions.tsv"))
write_tsv(combined_validated, validated_output)

message("✅ Stage 5 complete. Validated fusion output: ", validated_output)
