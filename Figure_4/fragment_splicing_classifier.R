#!/usr/bin/env Rscript

# Fragment-level classifier for paired-end amplicon sequencing.
#
# Each paired-end fragment is one observation. The two mates are paired by QNAME
# and classified together.
#
# Standard classes:
#   FULL          = E1-E2 AND E2-E3 junctions on the same fragment
#   SKIP          = E1-E3 junction on the fragment
#   UNSPLICED     = no splice junction + fragment covers Exon 1 and Exon 3
#   PARTIAL       = exactly one canonical junction (E1-E2 or E2-E3)
#   OTHER_SPLICED = one or more non-canonical/other splice junctions
#   AMBIGUOUS     = not enough evidence to classify confidently
#
# Inputs for each sample:
#   <sample>.gtf or <sample>_full.gtf[.gz]
#   <sample>_namesorted.bam (recommended)
#
# The GTF must describe the three-exon construct. The first three exons in
# transcript order are used as Exon 1, Exon 2, Exon 3.
#
# Usage:
#   Rscript fragment_splicing_classifier.R <input_dir> <output_dir>
#
# Example:
#   Rscript fragment_splicing_classifier.R D341_rep1 D341_rep1_fragment_classification
#
# Outputs:
#   <sample>_fragment_classification.tsv
#   <sample>_summary.tsv
#   all_samples_fragment_summary.tsv
#   processing_log.tsv
#
# IMPORTANT:
#   BAMs should be NAME-SORTED so the two mates of each fragment are adjacent.
#   For example:
#     samtools sort -n -o sample_namesorted.bam sample_Aligned.out_sorted.bam

suppressPackageStartupMessages({
  req <- c("Rsamtools", "rtracklayer", "GenomicRanges")
  miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(
      "Missing R/Bioconductor packages: ", paste(miss, collapse = ", "),
      "\nInstall with:\n",
      "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "BiocManager::install(c('Rsamtools','rtracklayer','GenomicRanges'))"
    )
  }
  library(Rsamtools)
  library(rtracklayer)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat(
    "Usage:\n",
    "  Rscript fragment_splicing_classifier.R <input_dir> <output_dir>\n\n",
    "Example:\n",
    "  Rscript fragment_splicing_classifier.R D341_rep1 D341_fragment_results\n"
  )
  quit(status = 1)
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- normalizePath(args[[2]], mustWork = FALSE)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# File matching
# -----------------------------
strip_suffix <- function(path) {
  x <- basename(path)
  x <- sub("_namesorted\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("_Aligned\\.out_sorted\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("_full\\.gtf(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- sub("\\.gtf(\\.gz)?$", "", x, ignore.case = TRUE)
  x
}

find_matching_file <- function(files, sample) {
  if (!length(files)) return(NA_character_)
  stems <- vapply(files, strip_suffix, character(1))
  exact <- which(stems == sample)
  if (length(exact)) return(files[[exact[1]]])
  hit <- which(grepl(sample, stems, fixed = TRUE))
  if (length(hit)) {
    scores <- nchar(stems[hit])
    return(files[[hit[which.min(scores)]]])
  }
  NA_character_
}

gtf_files <- list.files(
  input_dir,
  pattern = "\\.gtf(\\.gz)?$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

bam_files <- list.files(
  input_dir,
  pattern = "\\.bam$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

if (!length(gtf_files)) stop("No GTF files found in: ", input_dir)
if (!length(bam_files)) stop("No BAM files found in: ", input_dir)

samples <- unique(c(vapply(gtf_files, strip_suffix, character(1)),
                    vapply(bam_files, strip_suffix, character(1))))
samples <- samples[nzchar(samples)]

message("Found ", length(samples), " sample name(s).")

# -----------------------------
# Read the first three exons in transcript order
# -----------------------------
read_three_exons <- function(gtf_file) {
  gr <- rtracklayer::import(gtf_file)
  ex <- gr[as.character(mcols(gr)$type) == "exon"]
  if (!length(ex)) stop("No exon features found in ", basename(gtf_file))
  
  df <- data.frame(
    seqname = as.character(seqnames(ex)),
    start = as.integer(start(ex)),
    end = as.integer(end(ex)),
    strand = as.character(strand(ex)),
    transcript_id = if ("transcript_id" %in% names(mcols(ex))) as.character(mcols(ex)$transcript_id) else NA_character_,
    exon_number = if ("exon_number" %in% names(mcols(ex))) suppressWarnings(as.integer(as.character(mcols(ex)$exon_number))) else NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # Select one transcript if multiple are present: prefer the transcript
  # with the largest number of exon rows.
  tx <- unique(df$transcript_id[!is.na(df$transcript_id) & nzchar(df$transcript_id)])
  if (length(tx) > 1) {
    tx_counts <- sort(table(df$transcript_id), decreasing = TRUE)
    df <- df[df$transcript_id == names(tx_counts)[1], , drop = FALSE]
  }
  
  if (length(unique(df$seqname)) > 1) {
    seq_counts <- sort(table(df$seqname), decreasing = TRUE)
    df <- df[df$seqname == names(seq_counts)[1], , drop = FALSE]
  }
  
  # Biological/transcript order.
  strand_val <- unique(df$strand[df$strand %in% c("+", "-")])[1]
  if (is.na(strand_val)) strand_val <- "+"
  if (strand_val == "+") {
    df <- df[order(df$start, df$end), , drop = FALSE]
  } else {
    df <- df[order(-df$start, -df$end), , drop = FALSE]
  }
  
  # Remove duplicate exon intervals and keep first 3.
  df <- unique(df[, c("seqname", "start", "end", "strand")])
  if (nrow(df) < 3) stop("Fewer than 3 exons found in ", basename(gtf_file))
  df <- df[seq_len(3), , drop = FALSE]
  df$exon <- 1:3
  rownames(df) <- NULL
  df
}

# -----------------------------
# CIGAR helpers
# -----------------------------
get_junctions <- function(cigar, pos1) {
  if (is.na(cigar) || !nzchar(cigar)) return(character())
  
  pieces <- regmatches(
    cigar,
    gregexpr("([0-9]+)([MIDNSHP=X])", cigar, perl = TRUE)
  )[[1]]
  if (!length(pieces)) return(character())
  
  ref_pos <- as.integer(pos1)
  out <- character()
  
  for (piece in pieces) {
    n <- as.integer(sub("([0-9]+).*", "\\1", piece))
    op <- sub("[0-9]+", "", piece)
    
    if (op == "N") {
      js <- ref_pos
      je <- ref_pos + n - 1L
      out <- c(out, paste0(js, "-", je))
      ref_pos <- ref_pos + n
    } else if (op %in% c("M", "D", "=", "X")) {
      ref_pos <- ref_pos + n
    }
  }
  
  unique(out)
}

aligned_reference_span <- function(pos1, cigar) {
  if (is.na(pos1) || is.na(cigar) || !nzchar(cigar)) return(c(NA_integer_, NA_integer_))
  
  pieces <- regmatches(
    cigar,
    gregexpr("([0-9]+)([MIDNSHP=X])", cigar, perl = TRUE)
  )[[1]]
  if (!length(pieces)) return(c(as.integer(pos1), as.integer(pos1)))
  
  ref_len <- 0L
  for (piece in pieces) {
    n <- as.integer(sub("([0-9]+).*", "\\1", piece))
    op <- sub("[0-9]+", "", piece)
    if (op %in% c("M", "D", "N", "=", "X")) ref_len <- ref_len + n
  }
  
  c(as.integer(pos1), as.integer(pos1) + ref_len - 1L)
}

interval_overlap <- function(a, b, x, y) {
  !anyNA(c(a, b, x, y)) && a <= y && b >= x
}

# -----------------------------
# Classify one fragment
# -----------------------------
classify_fragment <- function(r1, r2, exons) {
  reads <- list(r1, r2)
  reads <- reads[!vapply(reads, is.null, logical(1))]
  
  if (!length(reads)) {
    return(list(classification = "AMBIGUOUS", junctions = character(), evidence = "no_reads"))
  }
  
  all_junctions <- character()
  spans <- list()
  
  for (r in reads) {
    all_junctions <- c(all_junctions, get_junctions(r$cigar, r$pos))
    spans[[length(spans) + 1L]] <- aligned_reference_span(r$pos, r$cigar)
  }
  all_junctions <- unique(all_junctions)
  
  # Expected junctions from exon boundaries.
  e1e2 <- paste0(exons$end[1] + 1L, "-", exons$start[2] - 1L)
  e2e3 <- paste0(exons$end[2] + 1L, "-", exons$start[3] - 1L)
  e1e3 <- paste0(exons$end[1] + 1L, "-", exons$start[3] - 1L)
  
  has_e1e2 <- e1e2 %in% all_junctions
  has_e2e3 <- e2e3 %in% all_junctions
  has_e1e3 <- e1e3 %in% all_junctions
  
  if (has_e1e2 && has_e2e3) {
    return(list("classification" = "FULL",
                "junctions" = all_junctions,
                "evidence" = "E1-E2_and_E2-E3"))
  }
  
  if (has_e1e3) {
    return(list("classification" = "SKIP",
                "junctions" = all_junctions,
                "evidence" = "E1-E3"))
  }
  
  if (has_e1e2 || has_e2e3) {
    return(list("classification" = "PARTIAL",
                "junctions" = all_junctions,
                "evidence" = if (has_e1e2) "E1-E2_only" else "E2-E3_only"))
  }
  
  # No N operations. A fragment is called unspliced only if its two mates
  # collectively cover Exon 1 and Exon 3.
  if (!length(all_junctions)) {
    has_e1 <- any(vapply(spans, function(s) {
      interval_overlap(s[1], s[2], exons$start[1], exons$end[1])
    }, logical(1)))
    
    has_e3 <- any(vapply(spans, function(s) {
      interval_overlap(s[1], s[2], exons$start[3], exons$end[3])
    }, logical(1)))
    
    if (has_e1 && has_e3) {
      return(list("classification" = "UNSPLICED",
                  "junctions" = character(),
                  "evidence" = "no_N_and_E1_E3_coverage"))
    }
  }
  
  if (length(all_junctions)) {
    return(list("classification" = "OTHER_SPLICED",
                "junctions" = all_junctions,
                "evidence" = paste(all_junctions, collapse = ";")))
  }
  
  list("classification" = "AMBIGUOUS",
       "junctions" = character(),
       "evidence" = "insufficient_E1_E3_coverage")
}

# -----------------------------
# Process one BAM
# -----------------------------
process_bam <- function(bam_file, gtf_file, sample, out_dir) {
  exons <- read_three_exons(gtf_file)
  
  message("Processing ", sample, "...")
  message("  reference = ", exons$seqname[1])
  message("  E1 = ", exons$start[1], "-", exons$end[1])
  message("  E2 = ", exons$start[2], "-", exons$end[2])
  message("  E3 = ", exons$start[3], "-", exons$end[3])
  
  bf <- BamFile(bam_file, yieldSize = 100000L)
  open(bf)
  on.exit(close(bf), add = TRUE)
  
  param <- ScanBamParam(
    what = c("qname", "flag", "rname", "pos", "cigar"),
    flag = scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      isSupplementaryAlignment = FALSE
    )
  )
  
  rows <- list()
  row_i <- 0L
  pending <- NULL
  pending_qname <- NULL
  non_adjacent_mates <- 0L
  
  make_record <- function(x, i) {
    lapply(x, function(z) if (length(z)) z[[i]] else NULL)
  }
  
  process_pair <- function(a, b) {
    reads <- list(a, b)
    reads <- reads[!vapply(reads, is.null, logical(1))]
    
    r1 <- NULL
    r2 <- NULL
    for (r in reads) {
      flag <- as.integer(r$flag)
      if (bitwAnd(flag, 64L) != 0L) r1 <- r
      if (bitwAnd(flag, 128L) != 0L) r2 <- r
    }
    
    if (is.null(r1) || is.null(r2)) {
      qn <- if (!is.null(r1)) r1$qname else r2$qname
      return(data.frame(
        sample = sample,
        qname = qn,
        classification = "AMBIGUOUS",
        junctions = NA_character_,
        evidence = "missing_mate",
        R1_pos = if (!is.null(r1)) r1$pos else NA_integer_,
        R1_cigar = if (!is.null(r1)) r1$cigar else NA_character_,
        R2_pos = if (!is.null(r2)) r2$pos else NA_integer_,
        R2_cigar = if (!is.null(r2)) r2$cigar else NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    
    cls <- classify_fragment(r1, r2, exons)
    
    data.frame(
      sample = sample,
      qname = r1$qname,
      classification = cls$classification,
      junctions = if (length(cls$junctions)) paste(cls$junctions, collapse = ";") else NA_character_,
      evidence = cls$evidence,
      R1_pos = r1$pos,
      R1_cigar = r1$cigar,
      R2_pos = r2$pos,
      R2_cigar = r2$cigar,
      stringsAsFactors = FALSE
    )
  }
  
  repeat {
    chunk <- scanBam(bf, param = param)[[1]]
    if (!length(chunk$qname)) break
    
    n <- length(chunk$qname)
    for (i in seq_len(n)) {
      rec <- make_record(chunk, i)
      qn <- rec$qname
      
      if (is.null(pending)) {
        pending <- rec
        pending_qname <- qn
      } else if (identical(qn, pending_qname)) {
        row_i <- row_i + 1L
        rows[[row_i]] <- process_pair(pending, rec)
        pending <- NULL
        pending_qname <- NULL
      } else {
        # If mates are not adjacent, do not guess. The BAM was not name-sorted
        # or contains an unusual ordering. The orphan is counted as ambiguous.
        non_adjacent_mates <- non_adjacent_mates + 1L
        row_i <- row_i + 1L
        rows[[row_i]] <- process_pair(pending, NULL)
        pending <- rec
        pending_qname <- qn
      }
    }
  }
  
  if (!is.null(pending)) {
    row_i <- row_i + 1L
    rows[[row_i]] <- process_pair(pending, NULL)
  }
  
  if (!length(rows)) stop("No mapped primary alignments found in ", basename(bam_file))
  
  fragments <- do.call(rbind, rows)
  
  # Standard categories, including zeroes.
  standard <- c("FULL", "SKIP", "UNSPLICED", "PARTIAL", "OTHER_SPLICED", "AMBIGUOUS")
  counts <- as.data.frame(table(factor(fragments$classification, levels = standard)), stringsAsFactors = FALSE)
  names(counts) <- c("classification", "fragments")
  
  informative <- sum(counts$fragments[counts$classification %in% c("FULL", "SKIP", "UNSPLICED")])
  counts$informative_fragments <- informative
  counts$fraction_informative <- if (informative > 0) counts$fragments / informative else NA_real_
  counts$sample <- sample
  counts <- counts[, c("sample", "classification", "fragments", "informative_fragments", "fraction_informative")]
  
  # Useful warning for this workflow.
  if (non_adjacent_mates > 0L) {
    warning(sample, ": ", non_adjacent_mates,
            " fragment(s) were not adjacent in the BAM and were treated as ambiguous. Name-sort the BAM with samtools sort -n.")
  }
  
  frag_out <- file.path(out_dir, paste0(sample, "_fragment_classification.tsv"))
  sum_out <- file.path(out_dir, paste0(sample, "_summary.tsv"))
  
  write.table(fragments, frag_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write.table(counts, sum_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  
  list(
    summary = counts,
    fragment_file = frag_out,
    summary_file = sum_out,
    non_adjacent = non_adjacent_mates
  )
}

# -----------------------------
# Batch processing
# -----------------------------
log_rows <- list()
summary_rows <- list()
k <- 0L

for (sample in samples) {
  gtf <- find_matching_file(gtf_files, sample)
  bam <- find_matching_file(bam_files, sample)
  
  if (is.na(gtf) || is.na(bam)) {
    msg <- paste0(
      "Missing ",
      if (is.na(gtf)) "GTF" else "",
      if (is.na(gtf) && is.na(bam)) " and " else "",
      if (is.na(bam)) "BAM" else ""
    )
    warning(sample, ": ", msg)
    k <- k + 1L
    log_rows[[k]] <- data.frame(sample = sample, status = "SKIPPED", message = msg, stringsAsFactors = FALSE)
    next
  }
  
  res <- tryCatch(
    process_bam(bam, gtf, sample, out_dir),
    error = function(e) {
      warning(sample, ": FAILED: ", conditionMessage(e))
      NULL
    }
  )
  
  k <- k + 1L
  if (is.null(res)) {
    log_rows[[k]] <- data.frame(sample = sample, status = "FAILED", message = "See warning above", stringsAsFactors = FALSE)
  } else {
    log_rows[[k]] <- data.frame(sample = sample, status = "OK", message = "", stringsAsFactors = FALSE)
    summary_rows[[length(summary_rows) + 1L]] <- res$summary
  }
}

log_df <- do.call(rbind, log_rows)
write.table(log_df, file.path(out_dir, "processing_log.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (length(summary_rows)) {
  all_summary <- do.call(rbind, summary_rows)
  rownames(all_summary) <- NULL
  write.table(
    all_summary,
    file.path(out_dir, "all_samples_fragment_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

cat("\nFinished. Results written to: ", out_dir, "\n", sep = "")