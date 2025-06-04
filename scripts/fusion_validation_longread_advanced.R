#!/usr/bin/env Rscript
# ================================
# Fusion Detection Pipeline (Consolidated & Revised)
# ================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(stringr)
  library(glue)
  library(tibble)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

# ================================
# Stage 0: Setup - CLI Arguments and Parameters
# ================================
print("🚀 Stage 0: Initializing script and parameters...")

option_list <- list(
  make_option(c("-i", "--input_file"), type = "character", 
              help = "Path to fusion TSV (used for main fusion list and output prefix). This is the primary input for fusion candidates."),
  make_option(c("-b", "--bam_file"), type = "character", 
              help = "Path to the coordinate-sorted and indexed BAM file (required)."),
  make_option(c("-f", "--fusions_legacy_stage2_name"), type = "character", default = NULL, # Changed option name
              help = "Legacy: Fusion TSV file path used for Stage 2's output naming convention. If not provided, will use --input_file for naming."),
  make_option(c("-r", "--ref_genome"), type = "character", 
              help = "Path to reference genome FASTA file (required for Stage 3b realignment)."),
  make_option(c("-m", "--minimap2_path"), type = "character", default = "minimap2", 
              help = "Path to the minimap2 executable."),
  make_option(c("-w", "--window"), type = "integer", default = 1000, 
              help = "Window size (bp) around fusion breakpoints for searching supporting reads [default %default]."),
  make_option(c("-q", "--mapq_cutoff"), type = "integer", default = 20, 
              help = "Minimum mapping quality (MAPQ) for filtering reads [default %default]."),
  make_option(c("-s", "--softclip_threshold"), type = "integer", default = 20, 
              help = "Minimum soft-clip length (bp) to consider for extraction [default %default].")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)) {
  print_help(opt_parser); stop("❌ Fusion TSV file (--input_file) must be provided.", call. = FALSE)
}
if (is.null(opt$bam_file)) {
  print_help(opt_parser); stop("❌ BAM file (--bam_file) must be provided.", call. = FALSE)
}
if (is.null(opt$ref_genome)) {
  print_help(opt_parser); stop("❌ Reference genome FASTA (--ref_genome) must be provided for Stage 3b.", call. = FALSE)
}

# Use consistent variable names
fusions_input_file_path <- opt$input_file
bam_file_path <- opt$bam_file
ref_genome_path <- opt$ref_genome
minimap2_executable <- opt$minimap2_path
search_window <- opt$window
min_mapq <- opt$mapq_cutoff
min_softclip_len <- opt$softclip_threshold

# For Stage 2 output naming (from original script)
fusions_file_for_stage2_naming <- if (!is.null(opt$fusions_legacy_stage2_name)) opt$fusions_legacy_stage2_name else fusions_input_file_path
if (!is.null(opt$fusions_legacy_stage2_name) && opt$fusions_legacy_stage2_name != fusions_input_file_path) {
    message(paste0("⚠️ Warning: Different files specified for --input_file (", fusions_input_file_path, 
                   ") and --fusions_legacy_stage2_name (", opt$fusions_legacy_stage2_name, 
                   "). Using --input_file for main fusion list. Stage 2 output name uses --fusions_legacy_stage2_name value."))
}


output_prefix <- tools::file_path_sans_ext(basename(fusions_input_file_path))
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)
print(paste("📂 Output prefix:", output_prefix))
print(paste("📂 Output directory:", output_dir))

# Standard SAM fields. Original script used "OPT" for the 12th field.
# This definition will be used when parsing SAM strings from system calls.
sam_fields_definition <- c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","ISIZE","SEQ","QUAL","OPT")

# ---- Load main fusion candidates list ----
# Ensure you address any issues with extra header lines like "" here,
# for example, by using skip=1 if that line is present and not part of the actual header names.
# fusion_input <- read_tsv(fusions_input_file_path, show_col_types = FALSE, skip = X) # if skipping lines
fusion_input <- read_tsv(fusions_input_file_path, show_col_types = FALSE)
if (nrow(fusion_input) == 0) stop("No fusions in input file.", call. = FALSE)

# MODIFICATION START: Rename partner columns to gene columns if they exist
if ("fiveprimepartner" %in% names(fusion_input)) {
  message("ℹ️ Renaming 'fiveprimepartner' to 'fiveprime_gene'.")
  fusion_input <- fusion_input %>% rename(fiveprime_gene = fiveprimepartner)
}
if ("threeprimepartner" %in% names(fusion_input)) {
  message("ℹ️ Renaming 'threeprimepartner' to 'threeprime_gene'.")
  fusion_input <- fusion_input %>% rename(threeprime_gene = threeprimepartner)
}
# MODIFICATION END

# Add optional gene name columns if they don't exist in fusion_input AFTER potential rename
# These specific names ("fiveprime_gene", "threeprime_gene") are now expected by Stage 5
optional_gene_cols <- c("fiveprime_gene", "threeprime_gene")
for (col_name in optional_gene_cols) {
  if (! (col_name %in% names(fusion_input))) {
    message(paste0("ℹ️ Adding missing optional column '", col_name, "' to fusion_input with NA values, as it wasn't found or created by renaming."))
    fusion_input[[col_name]] <- NA_character_
  }
}

# ---- Adjust chromosome names based on BAM ----
print(paste("🔬 Checking BAM for 'chr' prefix in chromosome names using:", bam_file_path))
bam_header_sq_lines <- system(paste("samtools view -H", shQuote(bam_file_path), "| grep '^@SQ'"), intern = TRUE)
if (length(bam_header_sq_lines) == 0) {
    stop(paste("❌ Could not read header or find @SQ lines in BAM file:", bam_file_path), call. = FALSE)
}
bam_uses_chr_prefix <- any(grepl("SN:chr", bam_header_sq_lines))
fusion_input <- fusion_input %>%
  mutate(
    fiveprime_chr = if (bam_uses_chr_prefix) ifelse(!grepl("^chr", fiveprime_chr), paste0("chr", fiveprime_chr), fiveprime_chr) else gsub("^chr", "", fiveprime_chr),
    threeprime_chr = if (bam_uses_chr_prefix) ifelse(!grepl("^chr", threeprime_chr), paste0("chr", threeprime_chr), threeprime_chr) else gsub("^chr", "", threeprime_chr)
  )
print("✅ Chromosome name adjustment complete based on BAM header.")

# ================================
# Stage 1: Supplementary Alignments (Original script's per-region approach)
# ================================
print("🔍 Stage 1: Extracting supplementary alignments (per-region)...")
# This stage is kept for fidelity to original script, though Stage 2 is generally more efficient.
# Its output is not directly used by the refined Stage 5, which prefers Stage 2's output.
supp_aligns_stage1 <- fusion_input %>%
  mutate(region_5p = paste0(fiveprime_chr, ":", fiveprime_search_start - search_window, "-", fiveprime_search_end + search_window),
         region_3p = paste0(threeprime_chr, ":", threeprime_search_start - search_window, "-", threeprime_search_end + search_window)) %>%
  pivot_longer(cols = c(region_5p, region_3p), names_to = "region_type", values_to = "region") %>%
  mutate(cmd = paste0("samtools view -f 2048 ", shQuote(bam_file_path), " ", shQuote(region))) %>%
  mutate(results = map(cmd, ~tryCatch(system(.x, intern = TRUE, ignore.stderr = TRUE), 
                                      warning = function(w) {message(paste("Warning for Stage 1 cmd:",.x, w$message)); character(0)}, 
                                      error = function(e) {message(paste("Error for Stage 1 cmd:",.x, e$message)); character(0)}
                                     ))) %>%
  # Original script had 'path' column, assuming it might be sample-specific BAM path from TSV.
  # Since we now mandate a single -b BAM, 'path' from fusion_input might be less relevant here.
  # For safety, if 'path' column exists in fusion_input, keep it, otherwise it won't be selected.
  select(any_of(c("fusion_id", "sample_id", "path")), region_type, region, results) %>% 
  unnest(results) %>%
  filter(results != "" & !is.na(results)) %>% 
  separate(results, into = sam_fields_definition, sep = "\t", extra = "merge", fill = "right") %>% 
  filter(!is.na(POS) & POS != "0") # Basic check on POS

supp_output_stage1_path <- file.path(output_dir, paste0(output_prefix, "_stage1_supplementary_alignments.tsv"))
write_tsv(supp_aligns_stage1, supp_output_stage1_path)
print(paste("✅ Stage 1 complete. Output:", supp_output_stage1_path, "(", nrow(supp_aligns_stage1), "records )"))

# ================================
# Stage 2: Extract all supplementary alignments once (Original script's Stage 2)
# ================================
message("🔎 Stage 2: Extracting all supplementary alignments once from BAM...")
supp_reads_stage2_raw <- system2("samtools", c("view", "-f", "2048", shQuote(bam_file_path)), stdout = TRUE, stderr = FALSE)
supp_df_stage2 <- tibble::tibble(raw = supp_reads_stage2_raw) %>%
  filter(raw != "") %>%
  # Your original script used specific names here, then "TAGS". Let's use sam_fields_definition
  separate(raw, into = sam_fields_definition, sep = "\t", fill = "right", extra = "merge") %>%
  rename(TAGS = OPT) %>% # Rename OPT to TAGS for consistency with SA parsing function
  mutate(POS = as.integer(POS), MAPQ = as.integer(MAPQ)) %>%
  filter(!is.na(POS) & !is.na(QNAME))

supp_hits_stage2 <- fusion_input %>%
  rowwise() %>%
  mutate(
    matched_reads_5p = list(
      supp_df_stage2 %>%
        filter(RNAME == fiveprime_chr &
               POS >= (fiveprime_search_start - search_window) &
               POS <= (fiveprime_search_end + search_window) &
               !is.na(MAPQ) & MAPQ >= min_mapq
              )
    ),
    matched_reads_3p = list(
      supp_df_stage2 %>%
        filter(RNAME == threeprime_chr &
               POS >= (threeprime_search_start - search_window) &
               POS <= (threeprime_search_end + search_window) &
               !is.na(MAPQ) & MAPQ >= min_mapq
              )
    )
  ) %>%
  mutate(
      all_hits_for_fusion = list(bind_rows(
          matched_reads_5p %>% mutate(matched_partner_region = "5prime"),
          matched_reads_3p %>% mutate(matched_partner_region = "3prime")
      ) %>% distinct(QNAME, RNAME, POS, .keep_all = TRUE)
    )
  ) %>%
  select(any_of(c("fusion_id", "sample_id")), all_hits_for_fusion) %>% # Select only existing columns
  unnest(all_hits_for_fusion) %>%
  ungroup() %>%
  filter(!is.na(QNAME))

# MODIFICATION: Clean the TAGS column before writing to TSV
if ("TAGS" %in% names(supp_hits_stage2)) {
  message("ℹ️ Stage 2: Cleaning internal tabs from TAGS column by replacing with semicolons.")
  supp_hits_stage2 <- supp_hits_stage2 %>%
    mutate(TAGS = str_replace_all(TAGS, "\t", ";"))
}
# END OF MODIFICATION

stage2_supp_output_basename <- tools::file_path_sans_ext(basename(fusions_file_for_stage2_naming))
stage2_supp_output_path <- file.path(output_dir, paste0(stage2_supp_output_basename, "_stage2_supplementary_alignments.tsv"))
write_tsv(supp_hits_stage2, stage2_supp_output_path)
message(paste("✅ Stage 2 complete. Candidate supplementary alignments written to:", stage2_supp_output_path, "(", nrow(supp_hits_stage2), "records )"))

#=================================
# Stage 3a: Soft-Clipped Read Extraction (Modified for Chunking)
# ================================
message(glue("🔍 Stage 3a: Extracting soft-clipped reads (threshold: {min_softclip_len} bp)..."))
# ... (BAM indexing code: if (!file.exists(paste0(bam_file_path, ".bai"))) { ... } ) ...
# Ensure your BAM indexing code is here if it's not in the ellipsis above.

softclip_tsv_path <- file.path(output_dir, glue("{output_prefix}_softclipped_reads_min{min_softclip_len}bp.tsv"))
# We will not write a single large FASTQ path here anymore, but will generate paths for chunks.
# softclip_fastq_path <- file.path(output_dir, glue("{output_prefix}_softclipped_reads_min{min_softclip_len}bp.fastq"))

fusion_chrs_stage3 <- unique(c(fusion_input$fiveprime_chr, fusion_input$threeprime_chr))
if (length(fusion_chrs_stage3) > 0) {
    regions_stage3 <- paste(shQuote(fusion_chrs_stage3), collapse = " ") # Define regions_stage3

    # ================================================================================
    # THIS IS THE MISSING BLOCK TO RE-INSERT OR UNCOMMENT
    # ================================================================================
    # Using -F 2304 for primary, non-supp, non-secondary, and min_mapq for the original read
    cmd_stage3a <- glue("samtools view -F 2304 -q {min_mapq} {shQuote(bam_file_path)} {regions_stage3}")
    message(glue("🔍 Running for Stage 3a: {cmd_stage3a}"))
    sam_lines_stage3a <- tryCatch(system(cmd_stage3a, intern = TRUE, ignore.stderr = TRUE),
                                  error = function(e) {message(paste("Error in samtools view (Stage 3a):",e$message)); character()})
    # ================================================================================

    # --- ADD THIS FOR DEBUGGING (Optional, but recommended for now) ---
    if (exists("sam_lines_stage3a")) {
      message(paste("DEBUG: sam_lines_stage3a object exists. Length:", length(sam_lines_stage3a)))
      if (length(sam_lines_stage3a) > 0 && length(sam_lines_stage3a) < 10) { # Print first few lines only if not too many
        message(paste("DEBUG: First few lines of sam_lines_stage3a:", paste(head(sam_lines_stage3a), collapse=" | ")))
      } else if (length(sam_lines_stage3a) >= 10) {
        message(paste("DEBUG: First line of sam_lines_stage3a:", head(sam_lines_stage3a, 1)))
      } else {
        message("DEBUG: sam_lines_stage3a is empty.")
      }
    } else {
      message("DEBUG: sam_lines_stage3a object DOES NOT EXIST here.")
    }
    # --- END OF DEBUGGING BLOCK ---

    if (length(sam_lines_stage3a) == 0) {
      message("⚠️ No reads found by samtools in Stage 3a for soft-clip analysis, or command failed.")
      write_tsv(tibble(), softclip_tsv_path) # Still write empty TSV for consistency
      # No FASTQ files to write. Ensure Stage 3b knows this.
    } else {
      message(paste("📄 Stage 3a: Retrieved", length(sam_lines_stage3a), "records from samtools view."))

      softclip_df_stage3a <- tibble(raw = sam_lines_stage3a) %>%
        # ... (the rest of your existing separate, filter, mutate chain) ...
        separate(raw, into = sam_fields_definition, sep = "\t", fill = "right", extra = "merge") %>%
        filter(!is.na(SEQ) & SEQ != "*" & !is.na(QUAL) & QUAL != "*" & QUAL != "0" & nchar(SEQ) == nchar(QUAL)) %>% 
        mutate(
          POS = as.integer(POS), MAPQ = as.integer(MAPQ),
          left_softclip = as.integer(str_match(CIGAR, "^(\\d+)S")[,2]),
          right_softclip = as.integer(str_match(CIGAR, "(\\d+)S$")[,2]),
          total_softclip = coalesce(left_softclip, 0L) + coalesce(right_softclip, 0L)
        ) %>%
        filter(total_softclip >= min_softclip_len)
      message(glue("✅ Stage 3a: Found {nrow(softclip_df_stage3a)} reads with valid SEQ/QUAL and soft-clips ≥{min_softclip_len}bp."))
      
      if ("OPT" %in% names(softclip_df_stage3a)) {
        message("ℹ️ Stage 3a: Cleaning internal tabs from OPT (tags) column by replacing with semicolons for TSV.")
        softclip_df_stage3a_for_tsv <- softclip_df_stage3a %>%
          mutate(OPT = str_replace_all(OPT, "\t", ";"))
        write_tsv(softclip_df_stage3a_for_tsv, softclip_tsv_path)
      } else {
        write_tsv(softclip_df_stage3a, softclip_tsv_path)
      }
      message(paste("📁 Stage 3a TSV output saved to:", softclip_tsv_path))

      # --- CHUNKING LOGIC FOR FASTQ ---
      # ... (rest of your chunking logic for FASTQ remains the same) ...
      if (nrow(softclip_df_stage3a) > 0) {
        reads_per_chunk <- 50000 
        num_chunks <- ceiling(nrow(softclip_df_stage3a) / reads_per_chunk)
        message(glue("ℹ️ Stage 3a: Splitting {nrow(softclip_df_stage3a)} reads into {num_chunks} FASTQ chunks of up to {reads_per_chunk} reads each."))
        for (chunk_num in 1:num_chunks) {
          start_row <- (chunk_num - 1) * reads_per_chunk + 1
          end_row <- min(chunk_num * reads_per_chunk, nrow(softclip_df_stage3a))
          current_chunk_df <- softclip_df_stage3a[start_row:end_row, ]
          chunk_fastq_path <- file.path(output_dir, glue("{output_prefix}_softclipped_reads_min{min_softclip_len}bp_chunk{chunk_num}.fastq"))
          fastq_lines_chunk <- purrr::map_chr(1:nrow(current_chunk_df), function(i) {
            paste0("@", current_chunk_df$QNAME[i], "\n",
                   current_chunk_df$SEQ[i], "\n", "+\n", current_chunk_df$QUAL[i])
          })
          writeLines(fastq_lines_chunk, chunk_fastq_path)
          message(paste("📁 Stage 3a FASTQ chunk", chunk_num, "saved to:", chunk_fastq_path))
        }
      } else {
        message("⚠️ No soft-clipped reads with valid SEQ/QUAL to write to FASTQ for Stage 3a.")
      }
    }
} else {
    message("⚠️ No chromosomes identified for Stage 3a based on fusion input.")
    write_tsv(tibble(), softclip_tsv_path) 
}
message("✅ Stage 3a complete.")

# ================================
# Stage 3b: Realign Reads that Originally had Soft-Clips (Modified for Chunking)
# ================================
message("🔍 Stage 3b: Realigning full reads that originally had soft-clips (in chunks)...")

# Define the final combined SAM output path
realign_output_sam_path <- file.path(output_dir, glue("{output_prefix}_softclipped_reads_realigned.sam")) # This is the final combined file

# Find the FASTQ chunk files
chunk_fastq_files <- list.files(path = output_dir, 
                                pattern = glue("^{output_prefix}_softclipped_reads_min{min_softclip_len}bp_chunk\\d+\\.fastq$"), 
                                full.names = TRUE)

if (length(chunk_fastq_files) > 0) {
  minimap2_preset_3b <- "-ax map-ont" # Or your chosen long-read preset
  message(paste("ℹ️ Using minimap2 preset:", minimap2_preset_3b, "for Stage 3b realignment (long reads)."))
  
  all_chunk_sam_outputs <- c() # To store paths of individual SAM outputs

  for (i in 1:length(chunk_fastq_files)) {
    current_fastq_chunk <- chunk_fastq_files[i]
    chunk_sam_output_path <- sub("\\.fastq$", ".sam", current_fastq_chunk) # Name SAM chunk based on FASTQ chunk
    all_chunk_sam_outputs <- c(all_chunk_sam_outputs, chunk_sam_output_path)

    message(paste("⏳ Aligning chunk", i, "/", length(chunk_fastq_files), ":", basename(current_fastq_chunk)))
    align_cmd_3b_chunk <- glue("{shQuote(minimap2_executable)} {minimap2_preset_3b} {shQuote(ref_genome_path)} {shQuote(current_fastq_chunk)} > {shQuote(chunk_sam_output_path)}")
    message(paste("   Running:", align_cmd_3b_chunk))
    
    system_status_3b_chunk <- system(align_cmd_3b_chunk)
    
    if(system_status_3b_chunk != 0) {
      message(paste("⚠️ Minimap2 failed for chunk", i, "(status:", system_status_3b_chunk, "). Output file:", chunk_sam_output_path))
      # Decide how to handle partial failure: stop, or continue and try to merge what's available?
      # For now, let's assume we might want to merge successful ones, or just note the failure.
      # If a chunk fails, its SAM file might be empty or incomplete.
    } else {
      message(paste("✅ Chunk", i, "alignment complete. Output:", chunk_sam_output_path))
    }
  }

  # --- Combine SAM files ---
  message("⏳ Combining individual SAM chunk outputs...")
  first_sam_written <- FALSE
  final_sam_con <- file(realign_output_sam_path, "w")

  for (sam_file_path in all_chunk_sam_outputs) {
    if (file.exists(sam_file_path) && file.info(sam_file_path)$size > 0) {
      sam_lines <- readLines(sam_file_path)
      header_lines <- sam_lines[startsWith(sam_lines, "@")]
      alignment_lines <- sam_lines[!startsWith(sam_lines, "@")]
      alignment_lines <- alignment_lines[alignment_lines != ""] # Remove empty lines

      if (!first_sam_written && length(header_lines) > 0) {
        writeLines(header_lines, final_sam_con)
        first_sam_written <- TRUE
      }
      if (length(alignment_lines) > 0) {
        writeLines(alignment_lines, final_sam_con)
      }
    } else {
      message(paste("⚠️ Skipping empty or non-existent SAM chunk:", sam_file_path))
    }
  }
  close(final_sam_con)

  if (file.exists(realign_output_sam_path) && file.info(realign_output_sam_path)$size > 0) {
    message(paste("✅ Stage 3b complete. Combined realigned SAM output:", realign_output_sam_path))
    # Optionally, delete individual SAM chunks after successful merge
    # for (sam_chunk in all_chunk_sam_outputs) { if(file.exists(sam_chunk)) file.remove(sam_chunk) }
    # message("ℹ️ Individual SAM chunks deleted.")
  } else {
    message(paste("⚠️ Stage 3b resulted in an empty combined SAM file:", realign_output_sam_path))
    # Ensure an empty file exists if it's expected by later stages, even if processing failed.
    if(!file.exists(realign_output_sam_path)) writeLines(character(), realign_output_sam_path)
  }

} else {
  message("⚠️ Skipping Stage 3b: No FASTQ chunks from Stage 3a available or files are empty.")
  writeLines(character(), realign_output_sam_path) # Create empty SAM if no FASTQ
}
# ================================
# Stage 4: Full-Length Spanning Fusion Reads (Original script's definition)
# ================================
print("🔍 Stage 4: Detecting reads with primary alignments in both partner regions...")
fullspan_reads_list_stage4 <- list()
# Original script named this output "_full_length_fusion_reads.tsv"
stage4_output_path <- file.path(output_dir, paste0(output_prefix, "_full_length_fusion_reads.tsv")) 
for (i in 1:nrow(fusion_input)) {
  current_fusion_s4 <- fusion_input[i, ]
  fusion_id_s4 <- current_fusion_s4$fusion_id
  sample_id_s4 <- current_fusion_s4$sample_id 
  region_5p_s4 <- paste0(current_fusion_s4$fiveprime_chr, ":", current_fusion_s4$fiveprime_search_start, "-", current_fusion_s4$fiveprime_search_end)
  region_3p_s4 <- paste0(current_fusion_s4$threeprime_chr, ":", current_fusion_s4$threeprime_search_start, "-", current_fusion_s4$threeprime_search_end)
  # Original script: -F 2048. Using -F 2304 (not secondary, not supplementary) for more specificity.
  cmd_5p_s4 <- paste0("samtools view -F 2304 -q ", min_mapq, " ", shQuote(bam_file_path), " ", shQuote(region_5p_s4))
  cmd_3p_s4 <- paste0("samtools view -F 2304 -q ", min_mapq, " ", shQuote(bam_file_path), " ", shQuote(region_3p_s4))
  reads_5p_s4_raw <- tryCatch(system(cmd_5p_s4, intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
  reads_3p_s4_raw <- tryCatch(system(cmd_3p_s4, intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
  if (length(reads_5p_s4_raw) > 0 && length(reads_3p_s4_raw) > 0) {
    qnames_5p_s4 <- tibble(raw=reads_5p_s4_raw) %>% separate(raw, into="QNAME", sep="\t", extra="drop", fill="right") %>% pull(QNAME) %>% unique()
    qnames_3p_s4 <- tibble(raw=reads_3p_s4_raw) %>% separate(raw, into="QNAME", sep="\t", extra="drop", fill="right") %>% pull(QNAME) %>% unique()
    spanning_qnames_s4 <- intersect(qnames_5p_s4, qnames_3p_s4)
    if (length(spanning_qnames_s4) > 0) {
      fullspan_reads_list_stage4[[length(fullspan_reads_list_stage4) + 1]] <- tibble(
        fusion_id = fusion_id_s4, sample_id = sample_id_s4, read_name = spanning_qnames_s4) # Keeping 'read_name' as per original script's output
    }
  }
}
if (length(fullspan_reads_list_stage4) > 0) {
  fullspan_df_stage4 <- bind_rows(fullspan_reads_list_stage4)
  write_tsv(fullspan_df_stage4, stage4_output_path) # Use stage4_output_path
  print(paste("✅ Stage 4 complete. Output:", stage4_output_path, "(", nrow(fullspan_df_stage4), "events )"))
} else {
  print("⚠️ No Stage 4 type 'full-length spanning' reads detected.")
  write_tsv(tibble(), stage4_output_path) 
}
print("🎉 Original pipeline Stages 1-4 completed!")

# ===============================================================
# Stage 5: Consolidate Evidence and Create Validation Summary Table (Optimized v2)
# ===============================================================
print("🔍 Stage 5 (Optimized v2): Consolidating evidence and creating validation summary...")

search_window_s5 <- opt$window # Already defined, but for clarity if block is moved
min_mapq_s5 <- opt$mapq_cutoff   # Already defined

# ---- Define paths to input evidence files from THIS SCRIPT's previous stages ----
supp_align_evidence_file_s5 <- stage2_supp_output_path # From Stage 2
realigned_softclip_origin_sam_file_s5 <- realign_output_sam_path # From Stage 3b
original_softclip_info_tsv_s5 <- softclip_tsv_path # From Stage 3a
full_length_spanning_evidence_file_s5 <- stage4_output_path # From Stage 4

final_summary_table_path_s5 <- file.path(output_dir, paste0(output_prefix, "_fusions_validation_summary_consolidated.tsv"))

parse_sa_tag_from_sam_tags_s5 <- function(tags_string) {
  if (is.na(tags_string) || tags_string == "") return(tibble())
  sa_tag_value <- str_extract(tags_string, "SA:Z:[^\\t;]+")
  if (is.na(sa_tag_value)) return(tibble())
  sa_string_clean <- sub("^SA:Z:", "", sa_tag_value)
  alignments <- str_split(sa_string_clean, ";")[[1]]
  alignments <- alignments[alignments != "" & !is.na(alignments)]
  if(length(alignments) == 0) return(tibble())
  map_dfr(alignments, function(aln) {
    parts <- str_split(aln, ",")[[1]]
    if (length(parts) >= 6) {
        tryCatch({
            tibble(sa_rname = parts[1], sa_pos = as.integer(parts[2]), sa_strand = parts[3],
                   sa_cigar = parts[4], sa_mapq = as.integer(parts[5]), sa_nm = as.integer(parts[6]))
        }, error = function(e) tibble())
    } else { tibble() }
  })
}

# ---- Load and Pre-process Evidence Data ----
print("⏳ Stage 5: Loading and pre-processing evidence data...")

loaded_supp_align_evidence_s5 <- tibble()
if (file.exists(supp_align_evidence_file_s5) && file.info(supp_align_evidence_file_s5)$size > 0) {
  loaded_supp_align_evidence_s5 <- read_tsv(supp_align_evidence_file_s5, show_col_types = FALSE,
                                         col_types = cols(MAPQ=col_integer(), POS=col_integer(), .default=col_character()))
  print(problems(loaded_supp_align_evidence_s5)) # Keep this to check for parsing issues
  if ("TAGS" %in% names(loaded_supp_align_evidence_s5)) {
    loaded_supp_align_evidence_s5 <- loaded_supp_align_evidence_s5 %>%
      mutate(SA_DETAILS_PARSED = map(TAGS, ~parse_sa_tag_from_sam_tags_s5(.x)))
  } else {
    message("⚠️ Warning: TAGS column not found in Stage 2 supplementary alignment output. SA tag parsing from this source will be limited.")
    loaded_supp_align_evidence_s5$SA_DETAILS_PARSED <- list(tibble())
  }
  print(paste("Loaded and pre-parsed", nrow(loaded_supp_align_evidence_s5), "records from Stage 2 supp. align. output:", supp_align_evidence_file_s5))
} else { print(paste("⚠️ Warning: Stage 2 supp. align. file not found/empty:", supp_align_evidence_file_s5)) }

loaded_realigned_sc_sam_s5 <- tibble()
if (file.exists(realigned_softclip_origin_sam_file_s5) && file.info(realigned_softclip_origin_sam_file_s5)$size > 0) {
  realigned_sam_lines_s5 <- readLines(realigned_softclip_origin_sam_file_s5)
  realigned_sam_data_lines_s5 <- realigned_sam_lines_s5[!startsWith(realigned_sam_lines_s5, "@")]
  if(length(realigned_sam_data_lines_s5) > 0){
    loaded_realigned_sc_sam_s5 <- tibble(raw_sam = realigned_sam_data_lines_s5) %>%
      separate(raw_sam, into = sam_fields_definition, sep = "\t", fill = "right", extra = "merge") %>%
      rename(TAGS = OPT) %>%
      mutate(POS = as.integer(POS), MAPQ = as.integer(MAPQ), FLAG = as.integer(FLAG)) %>%
      rename(ORIG_QNAME = QNAME) %>%
      filter(!is.na(RNAME) & RNAME != "*" & !is.na(POS) & !is.na(MAPQ) & MAPQ >= min_mapq_s5)

    if (nrow(loaded_realigned_sc_sam_s5) > 0 && "TAGS" %in% names(loaded_realigned_sc_sam_s5)) {
        loaded_realigned_sc_sam_s5 <- loaded_realigned_sc_sam_s5 %>%
            mutate(SA_DETAILS_PARSED = map(TAGS, ~parse_sa_tag_from_sam_tags_s5(.x)))
    } else {
        message("⚠️ Warning: TAGS column not found in realigned SAM or no records after filter. SA tag parsing from this source limited.")
        if(nrow(loaded_realigned_sc_sam_s5) > 0) loaded_realigned_sc_sam_s5$SA_DETAILS_PARSED <- list(tibble())
    }
    print(paste("Loaded and pre-parsed", nrow(loaded_realigned_sc_sam_s5), "mapped realigned records from Stage 3b SAM (MAPQ >=", min_mapq_s5, ")."))
  } else { print(paste("⚠️ Warning: Stage 3b realigned SAM has no data lines:", realigned_softclip_origin_sam_file_s5)) }
} else { print(paste("⚠️ Warning: Stage 3b realigned SAM file not found/empty:", realigned_softclip_origin_sam_file_s5)) }

original_sc_qnames_s5_df <- tibble()
if(file.exists(original_softclip_info_tsv_s5) && file.info(original_softclip_info_tsv_s5)$size > 0){
    original_sc_info_s5 <- read_tsv(original_softclip_info_tsv_s5, show_col_types = FALSE, col_types = cols(QNAME=col_character()))
    print(problems(original_sc_info_s5)) # Keep this to check for parsing issues
    if(nrow(original_sc_info_s5) > 0 && "QNAME" %in% names(original_sc_info_s5)){
        original_sc_qnames_s5_df <- tibble(ORIG_QNAME = unique(original_sc_info_s5$QNAME))
        print(paste("Loaded", nrow(original_sc_qnames_s5_df), "unique QNAMEs from Stage 3a TSV."))
    }
} else { print(paste("⚠️ Warning: Stage 3a original soft-clip TSV not found/empty:", original_softclip_info_tsv_s5)) }

loaded_full_length_spanning_evidence_s5 <- tibble()
if (file.exists(full_length_spanning_evidence_file_s5) && file.info(full_length_spanning_evidence_file_s5)$size > 0) {
  loaded_full_length_spanning_evidence_s5 <- read_tsv(full_length_spanning_evidence_file_s5, show_col_types = FALSE)
  print(problems(loaded_full_length_spanning_evidence_s5)) # Keep this to check for parsing issues
  if(!("QNAME" %in% names(loaded_full_length_spanning_evidence_s5)) && ("read_name" %in% names(loaded_full_length_spanning_evidence_s5))){
      loaded_full_length_spanning_evidence_s5 <- loaded_full_length_spanning_evidence_s5 %>% rename(QNAME = read_name)
  }
  print(paste("Loaded", nrow(loaded_full_length_spanning_evidence_s5), "records from Stage 4 full-length spanning output."))
} else { print(paste("⚠️ Warning: Stage 4 full-length spanning file not found/empty:", full_length_spanning_evidence_file_s5)) }

validation_summary_list_s5 <- list()
if (!exists("fusion_input") || nrow(fusion_input) == 0) {
    stop("❌ Stage 5: fusion_input data frame not found or empty.")
}
print(paste("🚀 Stage 5: Processing", nrow(fusion_input), "fusion candidates for final validation summary..."))

realigned_sc_sam_filtered_once_s5 <- tibble()
if(nrow(loaded_realigned_sc_sam_s5) > 0 && nrow(original_sc_qnames_s5_df) > 0) {
    realigned_sc_sam_filtered_once_s5 <- loaded_realigned_sc_sam_s5 %>%
        semi_join(original_sc_qnames_s5_df, by = "ORIG_QNAME")
    print(paste("Pre-filtered realigned SAM to", nrow(realigned_sc_sam_filtered_once_s5), "records matching original soft-clipped QNAMEs."))
}

for (idx in 1:nrow(fusion_input)) {
  current_fusion <- fusion_input[idx, ]; fusion_id_val <- current_fusion$fusion_id; sample_id_val <- current_fusion$sample_id
  fp_gene <- current_fusion[["fiveprime_gene"]]; tp_gene <- current_fusion[["threeprime_gene"]]
  fp_chr <- current_fusion$fiveprime_chr; fp_start_win <- current_fusion$fiveprime_search_start - search_window_s5; fp_end_win <- current_fusion$fiveprime_search_end + search_window_s5
  tp_chr <- current_fusion$threeprime_chr; tp_start_win <- current_fusion$threeprime_search_start - search_window_s5; tp_end_win <- current_fusion$threeprime_search_end + search_window_s5
  num_supp_reads <- 0; qnames_supp <- character(0); avg_mapq_supp <- NA_real_
  num_realigned_sc_reads <- 0; qnames_realigned_sc <- character(0); avg_mapq_realigned_sc <- NA_real_
  num_full_length_spanning_s4 <- 0; qnames_full_length_spanning_s4 <- character(0)

  # --- 1. Process Supplementary Alignment Evidence ---
  if (nrow(loaded_supp_align_evidence_s5) > 0 && "fusion_id" %in% names(loaded_supp_align_evidence_s5)) {
    # Pre-filter supplementary alignments for the current fusion (OPTIONAL but good practice if loaded_supp_align_evidence_s5 is huge)
    # This step is slightly different because loaded_supp_align_evidence_s5 already has fusion_id.
    sa_for_this_fusion <- loaded_supp_align_evidence_s5 %>% filter(fusion_id == fusion_id_val)
    
    # The rest of the supplementary processing remains the same...
    if(nrow(sa_for_this_fusion) > 0 && "SA_DETAILS_PARSED" %in% names(sa_for_this_fusion)){
        qname_sa_summary <- sa_for_this_fusion %>%
            filter(!is.na(MAPQ) & MAPQ >= min_mapq_s5) %>%
            group_by(QNAME) %>%
            summarise(
                has_sa_near_5p = any(!is.na(RNAME)&RNAME==fp_chr & !is.na(POS)&POS>=fp_start_win&POS<=fp_end_win, na.rm=TRUE),
                has_sa_near_3p = any(!is.na(RNAME)&RNAME==tp_chr & !is.na(POS)&POS>=tp_start_win&POS<=tp_end_win, na.rm=TRUE),
                sa_tag_points_to_5p_calc = any(sapply(SA_DETAILS_PARSED, function(sa_df) {
                    if(nrow(sa_df) > 0 && "sa_rname" %in% names(sa_df)) {
                        any(sa_df$sa_rname == fp_chr & !is.na(sa_df$sa_pos) & sa_df$sa_pos >= fp_start_win & sa_df$sa_pos <= fp_end_win & !is.na(sa_df$sa_mapq) & sa_df$sa_mapq >= min_mapq_s5, na.rm = TRUE)
                    } else { FALSE }
                })),
                sa_tag_points_to_3p_calc = any(sapply(SA_DETAILS_PARSED, function(sa_df) {
                    if(nrow(sa_df) > 0 && "sa_rname" %in% names(sa_df)) {
                        any(sa_df$sa_rname == tp_chr & !is.na(sa_df$sa_pos) & sa_df$sa_pos >= tp_start_win & sa_df$sa_pos <= tp_end_win & !is.na(sa_df$sa_mapq) & sa_df$sa_mapq >= min_mapq_s5, na.rm = TRUE)
                    } else { FALSE }
                })),
                mapq_scores_of_relevant_sa_segments = list(MAPQ[
                    (!is.na(RNAME)&RNAME==fp_chr & !is.na(POS)&POS>=fp_start_win&POS<=fp_end_win) |
                    (!is.na(RNAME)&RNAME==tp_chr & !is.na(POS)&POS>=tp_start_win&POS<=tp_end_win)
                ]),
                .groups = 'drop'
            ) %>%
            filter( (has_sa_near_5p & (has_sa_near_3p | sa_tag_points_to_3p_calc)) |
                    (has_sa_near_3p & (has_sa_near_5p | sa_tag_points_to_5p_calc)) )

        num_supp_reads <- nrow(qname_sa_summary)
        if(num_supp_reads > 0) {
           qnames_supp <- qname_sa_summary %>% pull(QNAME)
            all_mapqs <- unlist(qname_sa_summary$mapq_scores_of_relevant_sa_segments)
            valid_mapqs <- as.numeric(all_mapqs[!is.na(all_mapqs)])
            if(length(valid_mapqs) > 0) avg_mapq_supp <- mean(valid_mapqs, na.rm = TRUE)
        }
    }
  }

  # --- 2. Process Realigned Soft-Clip Origin Read Evidence ---
  if (nrow(realigned_sc_sam_filtered_once_s5) > 0 && "SA_DETAILS_PARSED" %in% names(realigned_sc_sam_filtered_once_s5)) {

    # ====================================================================================
    # START OF PRE-FILTERING LOGIC FOR REALIGNED SOFT-CLIP EVIDENCE
    # ====================================================================================
    message(glue("   Pre-filtering realigned reads for fusion_id: {fusion_id_val} ({fp_gene} --- {tp_gene}). Initial set: {nrow(realigned_sc_sam_filtered_once_s5)} reads."))
    
    # Create a temporary logical vector for filtering to handle potential NA from sapply if SA_DETAILS_PARSED is missing for a row
    # This ensures that filter() gets a logical vector of the correct length.
    
    # First, create logical vectors for primary alignment hits
    primary_hits_5p_vec <- (!is.na(realigned_sc_sam_filtered_once_s5$RNAME) & realigned_sc_sam_filtered_once_s5$RNAME == fp_chr & 
                           !is.na(realigned_sc_sam_filtered_once_s5$POS) & realigned_sc_sam_filtered_once_s5$POS >= fp_start_win & 
                           realigned_sc_sam_filtered_once_s5$POS <= fp_end_win)
    
    primary_hits_3p_vec <- (!is.na(realigned_sc_sam_filtered_once_s5$RNAME) & realigned_sc_sam_filtered_once_s5$RNAME == tp_chr & 
                           !is.na(realigned_sc_sam_filtered_once_s5$POS) & realigned_sc_sam_filtered_once_s5$POS >= tp_start_win & 
                           realigned_sc_sam_filtered_once_s5$POS <= tp_end_win)

    # Next, create a logical vector for SA tag hits
    # Ensure SA_DETAILS_PARSED exists and is a list before applying sapply
    sa_tag_hits_partner_vec <- rep(FALSE, nrow(realigned_sc_sam_filtered_once_s5)) # Default to FALSE
    if(is.list(realigned_sc_sam_filtered_once_s5$SA_DETAILS_PARSED)) {
        sa_tag_hits_partner_vec <- sapply(realigned_sc_sam_filtered_once_s5$SA_DETAILS_PARSED, function(sa_df_item) {
            if(!is.null(sa_df_item) && nrow(sa_df_item) > 0 && "sa_rname" %in% names(sa_df_item)) {
                any(
                    (sa_df_item$sa_rname == fp_chr & !is.na(sa_df_item$sa_pos) & sa_df_item$sa_pos >= fp_start_win & sa_df_item$sa_pos <= fp_end_win & !is.na(sa_df_item$sa_mapq) & sa_df_item$sa_mapq >= min_mapq_s5) |
                    (sa_df_item$sa_rname == tp_chr & !is.na(sa_df_item$sa_pos) & sa_df_item$sa_pos >= tp_start_win & sa_df_item$sa_pos <= tp_end_win & !is.na(sa_df_item$sa_mapq) & sa_df_item$sa_mapq >= min_mapq_s5),
                    na.rm = TRUE
                )
            } else {
                FALSE
            }
        })
    }
    # Ensure sa_tag_hits_partner_vec has the correct length and handles cases where SA_DETAILS_PARSED might be NULL for some rows
    if(length(sa_tag_hits_partner_vec) != nrow(realigned_sc_sam_filtered_once_s5)) {
        message(glue("Warning: Length mismatch in sa_tag_hits_partner_vec for fusion {fusion_id_val}. Defaulting to no SA hits for pre-filter."))
        sa_tag_hits_partner_vec <- rep(FALSE, nrow(realigned_sc_sam_filtered_once_s5))
    }


    realigned_reads_to_check <- realigned_sc_sam_filtered_once_s5[
      primary_hits_5p_vec | primary_hits_3p_vec | sa_tag_hits_partner_vec, 
    ]
    
    message(glue("   After pre-filtering: {nrow(realigned_reads_to_check)} reads to check for this fusion."))
    # ====================================================================================
    # END OF PRE-FILTERING LOGIC
    # ====================================================================================
    
    if(nrow(realigned_reads_to_check) > 0){ # Check if any reads remain after pre-filtering
        qname_realigned_summary <- realigned_reads_to_check %>%
            group_by(ORIG_QNAME) %>%
            summarise(
                realigned_aln_hits_5p = any(!is.na(RNAME)&RNAME==fp_chr & !is.na(POS)&POS>=fp_start_win&POS<=fp_end_win, na.rm=TRUE),
                realigned_aln_hits_3p = any(!is.na(RNAME)&RNAME==tp_chr & !is.na(POS)&POS>=tp_start_win&POS<=tp_end_win, na.rm=TRUE),
                sa_tag_points_to_5p_realigned_calc = any(sapply(SA_DETAILS_PARSED, function(sa_df) {
                    if(nrow(sa_df) > 0 && "sa_rname" %in% names(sa_df)) {
                        any(sa_df$sa_rname == fp_chr & !is.na(sa_df$sa_pos) & sa_df$sa_pos >= fp_start_win & sa_df$sa_pos <= fp_end_win & !is.na(sa_df$sa_mapq) & sa_df$sa_mapq >= min_mapq_s5, na.rm = TRUE)
                    } else { FALSE }
                })),
                sa_tag_points_to_3p_realigned_calc = any(sapply(SA_DETAILS_PARSED, function(sa_df) {
                    if(nrow(sa_df) > 0 && "sa_rname" %in% names(sa_df)) {
                        any(sa_df$sa_rname == tp_chr & !is.na(sa_df$sa_pos) & sa_df$sa_pos >= tp_start_win & sa_df$sa_pos <= tp_end_win & !is.na(sa_df$sa_mapq) & sa_df$sa_mapq >= min_mapq_s5, na.rm = TRUE)
                    } else { FALSE }
                })),
                mapq_scores_of_realigned_segments = list(MAPQ[!is.na(MAPQ)]), # MAPQ from the primary alignment of the realigned read
                .groups = 'drop'
            ) %>%
            filter( (realigned_aln_hits_5p & (realigned_aln_hits_3p | sa_tag_points_to_3p_realigned_calc)) |
                    (realigned_aln_hits_3p & (realigned_aln_hits_5p | sa_tag_points_to_5p_realigned_calc)) )
        num_realigned_sc_reads <- nrow(qname_realigned_summary)
        if(num_realigned_sc_reads > 0) {
            qnames_realigned_sc <- qname_realigned_summary %>% pull(ORIG_QNAME)
            all_mapqs_realigned <- unlist(qname_realigned_summary$mapq_scores_of_realigned_segments)
            valid_mapqs_realigned <- as.numeric(all_mapqs_realigned[!is.na(all_mapqs_realigned)])
            if(length(valid_mapqs_realigned) > 0) avg_mapq_realigned_sc <- mean(valid_mapqs_realigned, na.rm = TRUE)
        }
    }
  }

  # --- 3. Process Full-Length Spanning Reads ---
  # ... (rest of Stage 5 remains the same) ...
  if (nrow(loaded_full_length_spanning_evidence_s5) > 0 && "fusion_id" %in% names(loaded_full_length_spanning_evidence_s5)) {
    if ("QNAME" %in% names(loaded_full_length_spanning_evidence_s5)) {
        full_length_spanning_for_fusion <- loaded_full_length_spanning_evidence_s5 %>% filter(fusion_id == fusion_id_val)
        num_full_length_spanning_s4 <- nrow(full_length_spanning_for_fusion)
        if (num_full_length_spanning_s4 > 0) qnames_full_length_spanning_s4 <- full_length_spanning_for_fusion %>% distinct(QNAME) %>% pull(QNAME)
    } else {
        message("⚠️ Warning: QNAME column missing in loaded Stage 4 data for fusion_id: ", fusion_id_val)
    }
  }
  
  all_supporting_qnames <- unique(c(qnames_supp, qnames_realigned_sc, qnames_full_length_spanning_s4))
  total_unique_reads <- length(all_supporting_qnames)
  evidence_cats <- c()
  if (num_supp_reads > 0) evidence_cats <- c(evidence_cats, "Supplementary")
  if (num_realigned_sc_reads > 0) evidence_cats <- c(evidence_cats, "RealignedSoftClipOrigin")
  if (num_full_length_spanning_s4 > 0) evidence_cats <- c(evidence_cats, "FullLengthSpanning_S4")
  is_validated <- total_unique_reads > 0
  

  validation_summary_list_s5[[length(validation_summary_list_s5) + 1]] <- tibble(
    fusion_id = fusion_id_val, 
    sample_id = sample_id_val, 
    fiveprime_gene = fp_gene, 
    threeprime_gene = tp_gene,
    # ADDING NEW COLUMNS HERE:
    fiveprime_chr = fp_chr, # fp_chr is current_fusion$fiveprime_chr
    fiveprime_junction = current_fusion$fiveprime_junction, # Assuming this column exists in current_fusion
    fiveprime_strand = current_fusion$fiveprime_strand,     # Assuming this column exists in current_fusion
    threeprime_chr = tp_chr, # tp_chr is current_fusion$threeprime_chr
    threeprime_junction = current_fusion$threeprime_junction, # Assuming this column exists in current_fusion
    threeprime_strand = current_fusion$threeprime_strand,   # Assuming this column exists in current_fusion
    # End of new columns
    fiveprime_locus = paste0(fp_chr, ":", current_fusion$fiveprime_search_start, "-", current_fusion$fiveprime_search_end), # Keeping this for context
    threeprime_locus = paste0(tp_chr, ":", current_fusion$threeprime_search_start, "-", current_fusion$threeprime_search_end), # Keeping this for context
    num_supplementary_reads = num_supp_reads, 
    avg_mapq_supplementary = round(avg_mapq_supp, 2),
    num_realigned_softclip_origin_reads = num_realigned_sc_reads, 
    avg_mapq_realigned_sc = round(avg_mapq_realigned_sc, 2),
    num_full_length_spanning_reads_S4 = num_full_length_spanning_s4, 
    total_unique_supporting_reads = total_unique_reads, 
    evidence_types = paste(evidence_cats, collapse = "; "),
    is_validated_by_pipeline = is_validated
    # Potentially other columns here
  )

}

if (length(validation_summary_list_s5) > 0) {
  final_df_s5 <- bind_rows(validation_summary_list_s5)
  write_tsv(final_df_s5, final_summary_table_path_s5)
  print(paste("✅ Stage 5 (Validation Summary) complete. Final summary table written to:", final_summary_table_path_s5))
  print(paste("Total fusions validated by at least one evidence type in Stage 5:", sum(final_df_s5$is_validated_by_pipeline, na.rm=TRUE)))
} else {
  print("⚠️ No fusion candidates processed in Stage 5 or no evidence found. Final validation table not generated.")
  write_tsv(tibble(), final_summary_table_path_s5)
}

print("🎉🎉🎉 Full pipeline (Consolidated & Revised) completed successfully! 🎉🎉🎉")
