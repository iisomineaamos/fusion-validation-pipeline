#!/usr/bin/env Rscript
# ================================
# Fusion Detection Pipeline (Stages 1-4 Combined)
# ================================
library(tidyverse)
library(optparse)
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
# Stage 3a: Soft-Clipped Read Extraction
# ================================
print(paste("🔍 Stage 3a: Extracting soft-clipped reads (threshold:", softclip_threshold, "bp)..."))

# Check if BAM file exists and is indexed
bam_file <- if (bam_override != "") bam_override else fusion_input$path[1]
if (!file.exists(bam_file)) {
  stop("❌ BAM file not found: ", bam_file)
}

# Check/create BAM index
if (!file.exists(paste0(bam_file, ".bai"))) {
  message("⏳ Indexing BAM file...")
  system(paste("samtools index", bam_file))
}

# --- Extract chromosome lengths from BAM header ---
bam_header <- system(paste("samtools view -H", bam_file, "| grep '^@SQ'"), intern = TRUE)

chr_info <- do.call(rbind, lapply(bam_header, function(line) {
  parts <- unlist(strsplit(line, "\t"))
  chr <- sub("SN:", "", parts[2])
  len <- as.numeric(sub("LN:", "", parts[3]))
  return(data.frame(chr = chr, len = len))
}))

chr_lengths <- setNames(chr_info$len, chr_info$chr)

# --- Define window regions with safe bounds ---
fusion_regions <- fusion_input %>%
  rowwise() %>%
  mutate(
    region_5p_chr = fiveprime_chr,
    region_5p_start = max(1, fiveprime_search_start - window),
    region_5p_end = min(chr_lengths[[fiveprime_chr]], fiveprime_search_end + window),
    region_3p_chr = threeprime_chr,
    region_3p_start = max(1, threeprime_search_start - window),
    region_3p_end = min(chr_lengths[[threeprime_chr]], threeprime_search_end + window)
  ) %>%
  ungroup() %>%
  mutate(
    region_5p = ifelse(!is.na(chr_lengths[region_5p_chr]),
                       paste0(region_5p_chr, ":", region_5p_start, "-", region_5p_end), NA),
    region_3p = ifelse(!is.na(chr_lengths[region_3p_chr]),
                       paste0(region_3p_chr, ":", region_3p_start, "-", region_3p_end), NA)
  ) %>%
  pivot_longer(cols = c(region_5p, region_3p),
               names_to = "region_type", values_to = "region") %>%
  filter(!is.na(region))

message("ℹ Searching in ", nrow(fusion_regions), " regions for soft-clipped reads (≥", softclip_threshold, "bp)")

# --- Extract SAM lines using samtools ---
softclip_results <- fusion_regions %>%
  mutate(cmd = paste("samtools view", bam_file, region, "2>&1")) %>%
  mutate(output = map(cmd, function(x) {
    out <- tryCatch({
      result <- system(x, intern = TRUE, ignore.stderr = FALSE)
      if (any(grepl("error|failed|invalid", result, ignore.case = TRUE))) return(NA)
      if (length(result) == 0) return(NA)
      result
    }, error = function(e) NA)
    out
  })) %>%
  filter(!is.na(output)) %>%
  select(fusion_id, region_type, region, output) %>%
  unnest(output)

# --- Process reads ---
if (nrow(softclip_results) == 0) {
  message("⚠️ No reads found in any regions.")
} else {
  softclip_results <- softclip_results %>%
    separate(output, into = sam_fields, sep = "\t", extra = "merge", fill = "right") %>%
    mutate(
      left_softclip = as.numeric(str_extract(CIGAR, "^(\\d+)S")),
      right_softclip = as.numeric(str_extract(CIGAR, "(\\d+)S$")),
      total_softclip = coalesce(left_softclip, 0) + coalesce(right_softclip, 0),
      has_softclip = total_softclip >= softclip_threshold
    )
  
  message("ℹ Found ", sum(softclip_results$has_softclip), 
          " reads with soft-clips ≥", softclip_threshold, "bp")

  if (sum(softclip_results$has_softclip) == 0) {
    message("⚠️ No reads met soft-clip threshold.")
  } else {
    softclip_output <- paste0("output/", output_prefix, 
                             "_softclipped_reads_", softclip_threshold, "bp.tsv")
    write_tsv(filter(softclip_results, has_softclip), softclip_output)
    message("✅ Saved ", sum(softclip_results$has_softclip), 
            " soft-clipped reads to ", softclip_output)
  }
}
# Modified Final Check in Stage 3a
if (exists("softclip_results") && 
    nrow(softclip_results) > 0 && 
    sum(softclip_results$has_softclip, na.rm = TRUE) > 0) {
  
  softclip_output <- paste0("output/", output_prefix, "_softclipped_reads_", softclip_threshold, "bp.tsv")
  write_tsv(softclip_results %>% filter(has_softclip), softclip_output)
  message("✅ Saved ", sum(softclip_results$has_softclip), " soft-clipped reads")
  
} else {
  message("⚠️ No soft-clipped reads ≥", softclip_threshold, "bp found")
  # Create empty file to prevent downstream errors
  softclip_output <- paste0("output/", output_prefix, "_softclipped_reads_", softclip_threshold, "bp.tsv")
  write_tsv(tibble(), softclip_output)
}

# ================================
# Stage 3b: Realignment of Soft-Clipped Reads
# ================================
# Add this at the start

print("🔍 Stage 3b: Realigning soft-clipped reads...")

# Check if we have valid soft-clipped reads
if (file.size(paste0("output/", output_prefix, "_softclipped_reads_", softclip_threshold, "bp.tsv")) == 0) {
  message("⚠️ Skipping Stage 3b: No soft-clipped reads available")
} else {
  # Read the filtered soft-clipped reads
  softclipped <- read_tsv(paste0("output/", output_prefix, "_softclipped_reads_", softclip_threshold, "bp.tsv"),
                         show_col_types = FALSE)
  
  # ---- Extract soft-clipped sequences ----
  softclip_fasta <- paste0("output/", output_prefix, "_softclip_sequences.fa")
  
  # Create FASTA only from reads with valid soft-clips
  fasta_lines <- softclipped %>%
    filter(has_softclip) %>%
    mutate(
      fasta_entry = ifelse(
        !is.na(left_softclip),
        paste0(">", QNAME, "_left\n", str_sub(SEQ, 1, left_softclip), "\n",
               ">", QNAME, "_right\n", str_sub(SEQ, nchar(SEQ) - right_softclip + 1, nchar(SEQ))),
        paste0(">", QNAME, "\n", SEQ)
      )
    ) %>%
    pull(fasta_entry)
  
  writeLines(fasta_lines, softclip_fasta)
  
  # ---- Run minimap2 ----
  realign_output <- paste0("output/", output_prefix, "_softclip_realignments.sam")
  realign_cmd <- paste(minimap2, "-a -x map-ont", ref_genome, softclip_fasta, ">", realign_output)
  system(realign_cmd)
  
  # ---- Parse and filter alignments ----
  if (file.exists(realign_output) && file.size(realign_output) > 0) {
    realigned <- read_tsv(
      realign_output,
      comment = "@",
      col_names = c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", 
                    "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "TAGS"),
      show_col_types = FALSE,
      na = c("", "NA", "*")
    ) %>%
      mutate(MAPQ = as.integer(MAPQ)) %>%
      filter(MAPQ >= mapq_cutoff)
    
    if (nrow(realigned) > 0) {
      realign_output_tsv <- paste0("output/", output_prefix, "_softclip_realignments_filtered.tsv")
      write_tsv(realigned, realign_output_tsv)
      message("✅ Saved ", nrow(realigned), " high-quality realignments")
    } else {
      message("⚠️ No realignments passed MAPQ cutoff (", mapq_cutoff, ")")
    }
  } else {
    message("⚠️ Realignment failed - no output generated")
  }
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
