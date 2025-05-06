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
  make_option(c("-q", "--mapq_cutoff"), type = "integer", default = 20, help = "Minimum MAPQ")
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
print("🔍 Stage 3a: Extracting soft-clipped reads...")

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

chr_lengths <- setNames(chr_info$len, chr_info$chr)  # Named vector: chr -> length

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
  filter(!is.na(region))  # Drop regions with invalid chromosomes

message("ℹ Searching in ", nrow(fusion_regions), " regions for soft-clipped reads")

# --- Extract SAM lines using samtools with robust error handling ---
softclip_results <- fusion_regions %>%
  mutate(cmd = paste("samtools view", bam_file, region, "2>&1")) %>%
  mutate(output = map(cmd, function(x) {
    message("  Processing: ", x)
    out <- tryCatch({
      result <- system(x, intern = TRUE, ignore.stderr = FALSE)
      
      # Check for samtools errors
      if (any(grepl("error|failed|invalid", result, ignore.case = TRUE))) {
        message("⚠️ Samtools error in region: ", x)
        message(paste(head(result, 3), collapse = "\n"))
        return(NA)
      }
      
      # Return results if we have any
      if (length(result) > 0) result else NA
    }, error = function(e) {
      message("⚠️ Command failed: ", x)
      message("Error: ", e$message)
      return(NA)
    })
    out
  })) %>%
  filter(!is.na(output)) %>%  # Remove failed regions
  select(fusion_id, region_type, region, output) %>%
  unnest(output)

# --- Check if any reads returned ---
if (nrow(softclip_results) == 0) {
  message("⚠️ No reads found in any regions.")
} else {
  message("ℹ Found ", nrow(softclip_results), " reads before soft-clip filtering")
  
  # --- Parse SAM fields and filter soft-clips ---
  softclip_results <- softclip_results %>%
    separate(output, into = c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "MRNM", "MPOS", "ISIZE", "SEQ", "QUAL"), 
             sep = "\t", extra = "merge", fill = "right") %>%
    mutate(
      # Parse CIGAR for detailed soft-clip info
      left_softclip = as.numeric(str_extract(CIGAR, "^(\\d+)S")),
      right_softclip = as.numeric(str_extract(CIGAR, "(\\d+)S$")),
      total_softclip = coalesce(left_softclip, 0) + coalesce(right_softclip, 0),
      has_softclip = str_detect(CIGAR, "S")
    )
  
  # Debug output
  message("ℹ CIGAR string examples:")
  print(head(softclip_results$CIGAR, 5))
  message("ℹ Soft-clip summary:")
  print(table(softclip_results$has_softclip))
  
  # Filter for significant soft-clipping (≥20bp)
  softclip_results <- softclip_results %>%
    filter(has_softclip, total_softclip >= 10)
  
  if (nrow(softclip_results) == 0) {
    message("⚠️ Reads found, but none with significant soft-clipping (≥20bp).")
  } else {
    softclip_output <- paste0("output/", output_prefix, "_softclipped_reads_longread.tsv")
    write_tsv(softclip_results, softclip_output)
    print(paste("✅ Stage 3a complete. Found", nrow(softclip_results), 
                "soft-clipped reads. Output:", softclip_output))
  }
}
# ================================
# Stage 3b: Realignment of Soft-Clipped Reads
# ================================
if (exists("softclip_results") && nrow(softclip_results) > 0) {
  print("🔍 Stage 3b: Realigning soft-clipped reads...")
  
  # ---- Extract soft-clipped sequences to FASTA ----
  softclip_fasta <- paste0("output/", output_prefix, "_softclip_sequences.fa")
  
  writeLines(
    unlist(
      lapply(1:nrow(softclip_results), function(i) {
        paste0(">", softclip_results$QNAME[i], "\n", softclip_results$SEQ[i])
      })
    ),
    con = softclip_fasta
  )
  
  # ---- Run minimap2 to re-align softclips ----
  realign_output <- paste0("output/", output_prefix, "_softclip_realignments.sam")
  realign_cmd <- paste(
    minimap2, "-a", ref_genome, softclip_fasta,
    "| samtools view -h -q 20 -F 4 - >", realign_output
  )
  
  status <- system(realign_cmd)
  
  if (status != 0) {
    message("❌ minimap2 alignment failed. Check if minimap2 is installed.")
  } else {
    # ---- Parse aligned reads ----
    realigned <- read_tsv(realign_output, comment = "@", col_names = FALSE, show_col_types = FALSE) %>%
      select(QNAME = X1, FLAG = X2, RNAME = X3, POS = X4, MAPQ = X5, CIGAR = X6, SEQ = X10)
    
    # ---- Filter by MAPQ ----
    realigned_summary <- softclip_results %>%
      inner_join(realigned, by = "QNAME") %>%
      filter(MAPQ >= 20)
    
    realign_output_tsv <- paste0("output/", output_prefix, "_softclip_realignments_filtered.tsv")
    write_tsv(realigned_summary, realign_output_tsv)
    print(paste("✅ Stage 3b complete. Output:", realign_output_tsv))
  }
} else {
  message("⚠️ Skipping Stage 3b: No soft-clipped reads available for realignment.")
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
