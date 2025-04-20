
# ============================================================
# Advanced Long-Read Fusion Validation Pipeline (All Stages)
# ============================================================

library(tidyverse)
library(IRanges)
library(optparse)

# ---------------------------
# Command-line arguments
# ---------------------------
option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", help = "Path to fusion TSV"),
  make_option(c("-r", "--ref_genome"), type = "character", help = "Path to reference genome FASTA"),
  make_option(c("-b", "--bam_override"), type = "character", default = "", help = "Path to BAM file (optional override)"),
  make_option(c("-q", "--mapq_cutoff"), type = "integer", default = 20, help = "Minimum MAPQ"),
  make_option("--run_stage", type = "character", default = "all", help = "Stage to run: all, stage1, stage2, stage3, stage4")
)

opt <- parse_args(OptionParser(option_list = option_list))
fusion_file <- opt$input_file
ref_genome <- opt$ref_genome
bam_override <- opt$bam_override
mapq_cutoff <- opt$mapq_cutoff
run_stage <- opt$run_stage

dir.create("output", showWarnings = FALSE)

# ---------------------------
# Load input
# ---------------------------
fusion_input <- read_tsv(fusion_file)
sam_fields <- c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","ISIZE","SEQ","QUAL","OPT")
softclip_window <- 1000
min_clip_length <- 10
min_clip_quality <- 15

# ================================
# Stage 1: Supplementary Alignments
# ================================
if (run_stage %in% c("all", "stage1")) {
  print("🔍 Stage 1: Extracting supplementary alignments...")

  supp_aligns <- fusion_input %>%
    mutate(region_5p = paste0(fiveprime_chr, ":", fiveprime_search_start - softclip_window, "-", fiveprime_search_end + softclip_window),
           region_3p = paste0(threeprime_chr, ":", threeprime_search_start - softclip_window, "-", threeprime_search_end + softclip_window)) %>%
    pivot_longer(cols = c(region_5p, region_3p), names_to = "region_type", values_to = "region") %>%
    mutate(bam = ifelse(bam_override != "", bam_override, path)) %>%
    mutate(cmd = paste0("samtools view -f 2048 ", bam, " ", region)) %>%
    mutate(results = map(cmd, ~system2("bash", c("-c", .x), stdout = TRUE))) %>%
    select(fusion_id, sample_id, region_type, region, results) %>%
    unnest(results) %>%
    separate(results, into = sam_fields, sep = "\t", extra = "merge") %>%
    filter(!is.na(POS))

  write_tsv(supp_aligns, "output/supplementary_alignments.tsv")
  print("✅ Stage 1 complete.")
}

# ================================
# Stage 2: Split-read Detection
# ================================
if (run_stage %in% c("all", "stage2")) {
  print("🔍 Stage 2: Detecting split-reads...")

  supp <- read_tsv("output/supplementary_alignments.tsv")

  split_candidates <- supp %>%
    group_by(fusion_id, QNAME) %>%
    summarise(n_regions = n_distinct(region_type), .groups = "drop") %>%
    filter(n_regions >= 2)

  split_reads <- supp %>%
    semi_join(split_candidates, by = c("fusion_id", "QNAME")) %>%
    arrange(fusion_id, QNAME)

  write_tsv(split_reads, "output/split_read_support.tsv")

  summary_split <- split_reads %>%
    group_by(fusion_id) %>%
    summarise(n_split_reads = n_distinct(QNAME), .groups = "drop")

  write_tsv(summary_split, "output/split_read_summary.tsv")
  print("✅ Stage 2 complete.")
}


# ================================
# Stage 3a: Soft-Clipped Read Extraction
# ================================
if (run_stage %in% c("all", "stage3a")) {
  print("🔍 Stage 3a: Extracting soft-clipped reads...")

  softclip_window <- 1000
  softclip_reads <- list()

  for (i in 1:nrow(fusion_input)) {
    fusion <- fusion_input[i, ]
    fusion_id <- fusion$fusion_id
    sample_bam <- ifelse(bam_override != "", bam_override, fusion$path)

    if (!file.exists(sample_bam)) next

    region5 <- paste0(fusion$fiveprime_chr, ":", fusion$fiveprime_search_start - softclip_window, "-", fusion$fiveprime_search_end + softclip_window)
    region3 <- paste0(fusion$threeprime_chr, ":", fusion$threeprime_search_start - softclip_window, "-", fusion$threeprime_search_end + softclip_window)

    cmd5 <- paste0("samtools view -F 2048 -q ", mapq_cutoff, " ", sample_bam, " ", region5)
    cmd3 <- paste0("samtools view -F 2048 -q ", mapq_cutoff, " ", sample_bam, " ", region3)

    reads5 <- system(cmd5, intern = TRUE)
    reads3 <- system(cmd3, intern = TRUE)

    extract_softclip <- function(reads) {
      parsed <- strsplit(reads, "	")
      tibble(
        QNAME = sapply(parsed, `[[`, 1),
        CIGAR = sapply(parsed, `[[`, 6),
        SEQ = sapply(parsed, `[[`, 10)
      ) %>%
      filter(grepl("^[0-9]+S", CIGAR) | grepl("[0-9]+S$", CIGAR)) %>%
      transmute(QNAME, softclip_seq = SEQ)
    }

    sc5 <- extract_softclip(reads5)
    sc3 <- extract_softclip(reads3)

    combined <- bind_rows(sc5, sc3)
    combined <- distinct(combined)

    if (nrow(combined) > 0) {
      softclip_reads[[length(softclip_reads) + 1]] <- combined
    }
  }

  if (length(softclip_reads) > 0) {
    final_softclips <- bind_rows(softclip_reads)
    write_tsv(final_softclips, "output/softclipped_reads_longread.tsv")
    print("✅ Stage 3a complete: Soft-clipped reads extracted.")
  } else {
    print("⚠️ No soft-clipped reads found in regions.")
  }
}

# ================================
# Stage 3b: Softclip Re-alignment
# ================================
if (run_stage %in% c("all", "stage3")) {
  print("🔍 Stage 3: Re-aligning soft-clipped reads...")

  softclips <- read_tsv("output/softclipped_reads_longread.tsv")

  softclip_fasta <- "output/softclip_sequences.fa"
  writeLines(
    unlist(lapply(1:nrow(softclips), function(i) {
      paste0(">", softclips$QNAME[i], "\n", softclips$softclip_seq[i])
    })),
    con = softclip_fasta
  )

  system(paste(
    "minimap2 -a", ref_genome, softclip_fasta,
    "| samtools view -h -q", mapq_cutoff, "-F 4 - > output/softclip_realignments.sam"
  ))

  realigned <- read_tsv("output/softclip_realignments.sam", comment = "@", col_names = FALSE) %>%
    select(QNAME = X1, FLAG = X2, RNAME = X3, POS = X4, MAPQ = X5, CIGAR = X6, SEQ = X10)

  realigned_summary <- softclips %>%
    inner_join(realigned, by = "QNAME") %>%
    filter(MAPQ >= mapq_cutoff)

  write_tsv(realigned_summary, "output/softclip_realignments_filtered.tsv")
  print("✅ Stage 3 complete.")
}

# ================================
# Stage 4: Full-Length Fusion Spans
# ================================
if (run_stage %in% c("all", "stage4")) {
  print("🔍 Stage 4: Detecting full-length spanning reads...")

  fullspan_reads <- list()

  for (i in 1:nrow(fusion_input)) {
    fusion <- fusion_input[i, ]
    fusion_id <- fusion$fusion_id
    bam_file <- ifelse(bam_override != "", bam_override, fusion$path)

    if (!file.exists(bam_file)) next

    region5 <- paste0(fusion$fiveprime_chr, ":", fusion$fiveprime_search_start, "-", fusion$fiveprime_search_end)
    region3 <- paste0(fusion$threeprime_chr, ":", fusion$threeprime_search_start, "-", fusion$threeprime_search_end)

    reads_5p <- system(paste0("samtools view -F 2048 -q ", mapq_cutoff, " ", bam_file, " ", region5), intern = TRUE)
    reads_3p <- system(paste0("samtools view -F 2048 -q ", mapq_cutoff, " ", bam_file, " ", region3), intern = TRUE)

    qnames_5p <- unique(str_extract(reads_5p, "^[^\t]+"))
    qnames_3p <- unique(str_extract(reads_3p, "^[^\t]+"))

    fullspan <- intersect(qnames_5p, qnames_3p)

    if (length(fullspan) > 0) {
      fullspan_reads[[length(fullspan_reads) + 1]] <- tibble(
        fusion_id = fusion_id,
        sample_id = fusion$sample_id,
        fullspan_read = fullspan
      )
    }
  }

  if (length(fullspan_reads) > 0) {
    final_fullspan <- bind_rows(fullspan_reads)
    write_tsv(final_fullspan, "output/full_length_fusion_reads.tsv")
    print("✅ Stage 4 complete.")
  } else {
    print("No full-length spanning reads found.")
  }
}
