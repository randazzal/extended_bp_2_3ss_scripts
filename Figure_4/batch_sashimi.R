#!/usr/bin/env Rscript

# Publication-style batch Sashimi plots for targeted minigene/splicing assays.
#
# Usage:
#   Rscript batch_sashimi_v3.R <input_dir> <output_dir> [min_junction_reads]
#
# Expected files for each mutant/sample:
#   <sample>.gtf
#   <sample>_SJ.out.tab
#   <sample>_Aligned.sortedByCoord.out.bam   (recommended, optional)
#
# The script will also accept .gtf.gz and BAM/SJ files in subdirectories.
#
# Outputs:
#   <sample>_sashimi.pdf
#   <sample>_sashimi.png
#   <sample>_junctions.tsv
#   all_mutants_junction_summary.tsv
#   all_mutants_junction_heatmap.pdf/png
#   processing_log.tsv

suppressPackageStartupMessages({
  required <- c(
    "ggplot2", "rtracklayer", "Rsamtools", "GenomicAlignments",
    "GenomicRanges", "IRanges", "S4Vectors"
  )
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      "\n\nInstall CRAN packages with:\n",
      "install.packages('ggplot2')\n\n",
      "Install Bioconductor packages with:\n",
      "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "BiocManager::install(c('rtracklayer','Rsamtools','GenomicAlignments','GenomicRanges','IRanges','S4Vectors'))"
    )
  }
})

library(S4Vectors)
# -----------------------------
# Defaults
# -----------------------------
min_junction_reads <- 3
coverage_quantile <- 0.995
pdf_width <- 8.5
pdf_height <- 6.4

# -----------------------------
# CLI
# -----------------------------
usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript batch_sashimi_v3.R <input_dir> <output_dir> [min_junction_reads]\n\n",
    "Example:\n",
    "  Rscript batch_sashimi_v3.R D341_rep1 D341_1_plots 10\n\n",
    "Notes:\n",
    "  - GTF and SJ.out.tab are required for every plotted mutant.\n",
    "  - BAM files are optional; if present they add the coverage track.\n",
    "  - Junctions below min_junction_reads are retained in TSV output but not drawn.\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  usage()
  quit(status = 1)
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (length(args) == 3) {
  min_junction_reads <- suppressWarnings(as.numeric(args[[3]]))
}
if (!is.finite(min_junction_reads) || min_junction_reads < 1) {
  stop("min_junction_reads must be a single number >= 1")
}
min_junction_reads <- as.integer(ceiling(min_junction_reads))

# -----------------------------
# File naming helpers
# -----------------------------
strip_suffix <- function(path) {
  x <- basename(path)
  x <- sub("_full\\.gtf(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- sub("_?SJ\\.out\\.tab$", "", x, ignore.case = TRUE)
  x <- sub("_?Aligned\\.out_sorted\\.bam$", "", x, ignore.case = TRUE)
  
  x
}

find_matching_file <- function(files, sample) {
  if (!length(files)) return(NA_character_)
  stems <- vapply(files, strip_suffix, character(1))
  exact <- which(stems == sample)
  if (length(exact)) return(files[[exact[1]]])
  # Fallback for small naming differences: longest filename containing sample.
  contains <- which(grepl(sample, stems, fixed = TRUE))
  if (length(contains)) {
    scores <- nchar(stems[contains])
    return(files[[contains[which.min(scores)]]])
  }
  NA_character_
}

# -----------------------------
# Discover files
# -----------------------------
gtf_files <- list.files(
  input_dir, pattern = "_full\\.gtf(\\.gz)?$", recursive = TRUE,
  full.names = TRUE, ignore.case = TRUE
)
sj_files <- list.files(
  input_dir, pattern = "_?SJ\\.out\\.tab$", recursive = FALSE,
  full.names = TRUE, ignore.case = TRUE
)
bam_files <- list.files(
  input_dir, pattern = "\\.bam$", recursive = TRUE,
  full.names = TRUE, ignore.case = TRUE
)

if (!length(gtf_files)) stop("No .gtf or .gtf.gz files found in: ", input_dir)
if (!length(sj_files)) stop("No SJ.out.tab files found in: ", input_dir)

# Use GTF sample names as the primary sample list.
samples <- unique(vapply(gtf_files, strip_suffix, character(1)))
samples <- samples[nzchar(samples)]

message("Found ", length(samples), " GTF/sample(s).")
message("Found ", length(sj_files), " SJ.out.tab file(s) and ", length(bam_files), " BAM file(s).")

# -----------------------------
# GTF parser
# -----------------------------
read_exons <- function(gtf_file) {
  gr <- rtracklayer::import(gtf_file)
  ex <- gr[as.character(mcols(gr)$type) == "exon"]
  if (!length(ex)) stop("No exon features found in ", gtf_file)
  
  df <- data.frame(
    seqname = as.character(GenomicRanges::seqnames(ex)),
    start = as.integer(GenomicRanges::start(ex)),
    end = as.integer(GenomicRanges::end(ex)),
    strand = as.character(GenomicRanges::strand(ex)),
    transcript_id = if ("transcript_id" %in% names(mcols(ex))) as.character(mcols(ex)$transcript_id) else NA_character_,
    exon_number = if ("exon_number" %in% names(mcols(ex))) suppressWarnings(as.integer(as.character(mcols(ex)$exon_number))) else NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # Most custom minigene GTFs contain one transcript. If more than one is present,
  # select the transcript with the largest number of exons and report it.
  valid_tx <- unique(df$transcript_id[!is.na(df$transcript_id) & nzchar(df$transcript_id)])
  if (length(valid_tx) > 1) {
    tx_counts <- table(df$transcript_id)
    selected <- names(tx_counts)[which.max(tx_counts)]
    message("  Multiple transcripts detected; using transcript ", selected)
    df <- df[df$transcript_id == selected, , drop = FALSE]
  }
  
  if (length(unique(df$seqname)) > 1) {
    # A custom construct should normally occupy one reference sequence.
    seq_counts <- sort(table(df$seqname), decreasing = TRUE)
    selected_seq <- names(seq_counts)[1]
    message("  Multiple reference sequences detected; using ", selected_seq)
    df <- df[df$seqname == selected_seq, , drop = FALSE]
  }
  
  # Remove accidental duplicate exon rows.
  df <- unique(df[, c("seqname", "start", "end", "strand", "transcript_id", "exon_number")])
  
  # Order biologically by transcript direction.
  strand <- unique(df$strand)
  strand <- strand[strand %in% c("+", "-")][1]
  if (is.na(strand) || !length(strand)) strand <- "+"
  
  if (strand == "+") {
    df <- df[order(df$start, df$end), , drop = FALSE]
  } else {
    df <- df[order(-df$start, -df$end), , drop = FALSE]
  }
  
  if (all(is.na(df$exon_number))) {
    df$exon_number <- seq_len(nrow(df))
  } else {
    # For plotting/classification, transcript order is more important than GTF numbering.
    df$exon_number <- seq_len(nrow(df))
  }
  
  rownames(df) <- NULL
  df
}

# -----------------------------
# STAR junction parser
# -----------------------------
read_sj <- function(sj_file) {
  x <- utils::read.table(
    sj_file, header = FALSE, sep = "\t", quote = "", comment.char = "",
    stringsAsFactors = FALSE, fill = TRUE, check.names = FALSE
  )
  if (ncol(x) < 9) stop("SJ.out.tab has fewer than 9 columns: ", sj_file)
  x <- x[, 1:9, drop = FALSE]
  names(x) <- c(
    "chrom", "intron_start", "intron_end", "strand", "motif",
    "annotated", "unique_reads", "multi_reads", "max_overhang"
  )
  
  x$chrom <- as.character(x$chrom)
  x$intron_start <- as.integer(x$intron_start)
  x$intron_end <- as.integer(x$intron_end)
  x$strand <- as.integer(x$strand)
  x$motif <- as.integer(x$motif)
  x$annotated <- as.integer(x$annotated)
  x$unique_reads <- as.integer(x$unique_reads)
  x$multi_reads <- as.integer(x$multi_reads)
  x$max_overhang <- as.integer(x$max_overhang)
  
  x <- x[(x$unique_reads + x$multi_reads) > 0, , drop = FALSE]
  rownames(x) <- NULL
  x
}

# -----------------------------
# Junction classification
# -----------------------------
make_reference_junctions <- function(exons) {
  if (nrow(exons) < 2) return(data.frame())
  
  # Coordinates in SJ.out.tab are intron start/end, so for a plus-strand transcript:
  # intron start = upstream exon end + 1; intron end = downstream exon start - 1.
  # For the minus strand, the same genomic boundaries apply; only transcript direction differs.
  out <- list()
  k <- 1L
  
  for (i in seq_len(nrow(exons) - 1L)) {
    for (j in (i + 1L):nrow(exons)) {
      event <- if (j == i + 1L) "canonical" else "exon_skipping"
      out[[k]] <- data.frame(
        intron_start = min(exons$end[i] + 1L, exons$end[j] + 1L),
        intron_end = max(exons$start[i] - 1L, exons$start[j] - 1L),
        upstream_exon = i,
        downstream_exon = j,
        event = event,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  
  # The min/max expression above is intentionally simplified by genomic order, but
  # the exon ordering can be reverse-strand. Recalculate each boundary from genomic span.
  ref <- do.call(rbind, out)
  for (r in seq_len(nrow(ref))) {
    i <- ref$upstream_exon[r]
    j <- ref$downstream_exon[r]
    ref$intron_start[r] <- min(exons$end[c(i, j)]) + 1L
    ref$intron_end[r] <- max(exons$start[c(i, j)]) - 1L
  }
  ref
}

classify_junctions <- function(sj, exons) {
  ref <- make_reference_junctions(exons)
  if (!nrow(sj)) return(sj)
  
  sj$junction <- paste0(sj$intron_start, "-", sj$intron_end)
  sj$span <- sj$intron_end - sj$intron_start + 1L
  
  if (!nrow(ref)) {
    sj$event <- "cryptic"
    sj$reference_match <- FALSE
    return(sj)
  }
  
  keys_ref <- paste(ref$intron_start, ref$intron_end, sep = "-")
  keys_sj <- sj$junction
  match_idx <- match(keys_sj, keys_ref)
  
  sj$reference_match <- !is.na(match_idx)
  sj$event <- "cryptic"
  sj$event[sj$reference_match] <- ref$event[match_idx[sj$reference_match]]
  sj$reference_event_exons <- NA_character_
  sj$reference_event_exons[sj$reference_match] <- paste0(
    "Exon ", ref$upstream_exon[match_idx[sj$reference_match]],
    "-Exon ", ref$downstream_exon[match_idx[sj$reference_match]]
  )
  sj
}

# -----------------------------
# Arc geometry
# -----------------------------
arc_data <- function(start, end, reads, max_reads, max_span, baseline = 0.15) {
  span <- max(1, end - start)
  n <- max(50L, min(300L, as.integer(span / 2L) + 50L))
  t <- seq(0, 1, length.out = n)
  
  # Height is primarily determined by genomic span, with a modest read-count effect.
  height <- 0.55 +
    1.15 * sqrt(span / max(1, max_span)) +
    0.45 * sqrt(reads / max(1, max_reads))
  
  data.frame(
    x = start + (end - start) * t,
    y = baseline + height * sin(pi * t),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# BAM coverage
# -----------------------------
get_coverage <- function(bam_file, seqname, start_pos, end_pos) {
  if (is.na(bam_file) || !nzchar(bam_file) || !file.exists(bam_file)) return(NULL)
  
  bf <- Rsamtools::BamFile(bam_file, yieldSize = 100000L)
  gr <- GenomicRanges::GRanges(
    seqnames = seqname,
    ranges = IRanges::IRanges(start = start_pos, end = end_pos)
  )
  param <- Rsamtools::ScanBamParam(which = gr)
  
  ga <- tryCatch(
    GenomicAlignments::readGAlignments(bf, param = param, use.names = FALSE),
    error = function(e) {
      warning("Could not read BAM coverage from ", basename(bam_file), ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(ga) || !length(ga)) {
    return(data.frame(position = seq.int(start_pos, end_pos), coverage = 0))
  }
  
  cv <- GenomicAlignments::coverage(ga)
  if (!(seqname %in% names(cv))) {
    return(data.frame(position = seq.int(start_pos, end_pos), coverage = 0))
  }
  
  v <- as.numeric(cv[[seqname]])
  needed <- end_pos - start_pos + 1L
  if (length(v) < end_pos) v <- c(v, rep(0, end_pos - length(v)))
  
  data.frame(
    position = seq.int(start_pos, end_pos),
    coverage = v[start_pos:end_pos],
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Plotting
# -----------------------------
event_colors <- c(
  canonical = "#4C78A8",
  exon_skipping = "#7B3294",
  cryptic = "#D55E00"
)

event_labels <- c(
  canonical = "Canonical",
  exon_skipping = "Exon skipping",
  cryptic = "Cryptic / alternative"
)

plot_sample <- function(sample, gtf_file, sj_file, bam_file, out_dir) {
  message("Processing: ", sample)
  
  exons <- read_exons(gtf_file)
  sj <- read_sj(sj_file)
  
  seqname <- exons$seqname[1]
  sj <- sj[sj$chrom == seqname, , drop = FALSE]
  if (!nrow(sj)) stop("No SJ.out.tab junctions match GTF sequence name '", seqname, "'")
  
  classified <- classify_junctions(sj, exons)
  classified$sample <- sample
  classified$plot_included <- classified$unique_reads >= min_junction_reads
  
  # Save the complete junction table before filtering.
  out_tsv <- file.path(out_dir, paste0(sample, "_junctions.tsv"))
  write.table(
    classified, out_tsv, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA"
  )
  
  plot_sj <- classified[classified$plot_included, , drop = FALSE]
  plot_sj <- plot_sj[order(plot_sj$unique_reads, decreasing = TRUE), , drop = FALSE]
  
  if (!nrow(plot_sj)) {
    warning(sample, ": no junctions pass min_junction_reads = ", min_junction_reads)
  }
  
  # Keep arcs from overlapping excessively by slightly limiting the largest number shown.
  max_reads <- if (nrow(plot_sj)) max(plot_sj$unique_reads) else 1
  max_span <- if (nrow(plot_sj)) max(plot_sj$span) else 1
  
  arcs <- list()
  labels <- list()
  if (nrow(plot_sj)) {
    for (i in seq_len(nrow(plot_sj))) {
      a <- arc_data(
        plot_sj$intron_start[i], plot_sj$intron_end[i],
        plot_sj$unique_reads[i], max_reads, max_span
      )
      a$junction <- plot_sj$junction[i]
      a$event <- plot_sj$event[i]
      a$reads <- plot_sj$unique_reads[i]
      arcs[[i]] <- a
      
      # Only label the more substantial junctions, avoiding a wall of tiny labels.
      if (plot_sj$unique_reads[i] >= min_junction_reads) {
        h <- 0.55 +
          1.15 * sqrt(plot_sj$span[i] / max(1, max_span)) +
          0.45 * sqrt(plot_sj$unique_reads[i] / max(1, max_reads))
        labels[[i]] <- data.frame(
          x = (plot_sj$intron_start[i] + plot_sj$intron_end[i]) / 2,
          y = 0.15 + h + 0.08,
          label = format(plot_sj$unique_reads[i], big.mark = ",", scientific = FALSE),
          event = plot_sj$event[i],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  arcs <- if (length(arcs)) do.call(rbind, arcs) else data.frame()
  labels <- if (length(labels)) do.call(rbind, labels) else data.frame()
  
  x_min <- max(1, min(exons$start) - 8)
  x_max <- max(exons$end) + 8
  
  # Coverage track.
  coverage <- get_coverage(
    bam_file = bam_file,
    seqname = seqname,
    start_pos = min(exons$start),
    end_pos = max(exons$end)
  )
  
  if (!is.null(coverage)) {
    ymax <- max(1, as.numeric(stats::quantile(coverage$coverage, coverage_quantile, na.rm = TRUE)))
    if (max(coverage$coverage, na.rm = TRUE) > ymax * 1.25) {
      # Give unusual coverage peaks a little headroom while avoiding one giant outlier dominating the plot.
      ymax <- max(ymax, as.numeric(stats::quantile(coverage$coverage, 0.999, na.rm = TRUE)))
    }
    
    coverage_plot <- ggplot2::ggplot(coverage, ggplot2::aes(x = position, y = coverage)) +
      ggplot2::geom_area(fill = "#B9C4CE", color = "#52616B", linewidth = 0.3) +
      ggplot2::geom_rect(
        data = exons,
        ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = ymax),
        inherit.aes = FALSE, fill = "white", alpha = 0.18,
        color = NA
      ) +
      ggplot2::coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, ymax), expand = FALSE) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::labs(x = NULL, y = "Coverage") +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(angle = 90, margin = ggplot2::margin(r = 6)),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.margin = ggplot2::margin(2, 6, 0, 6)
      )
  } else {
    coverage_plot <- ggplot2::ggplot() +
      ggplot2::annotate(
        "text", x = mean(c(x_min, x_max)), y = 0.5,
        label = "BAM not found: coverage track omitted",
        size = 3.5, color = "grey40"
      ) +
      ggplot2::xlim(x_min, x_max) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(2, 6, 0, 6))
  }
  
  # Main Sashimi panel.
  sashimi <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = data.frame(x = x_min, xend = x_max, y = 0.12, yend = 0.12),
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      color = "#8A949E", linewidth = 0.5
    ) +
    ggplot2::geom_path(
      data = arcs,
      ggplot2::aes(
        x = x, y = y, group = junction,
        color = event, linewidth = reads
      ),
      lineend = "round", show.legend = TRUE
    ) +
    ggplot2::geom_rect(
      data = exons,
      ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 0.32),
      fill = "#26394A", color = "#17232E", linewidth = 0.35
    ) +
    ggplot2::geom_text(
      data = exons,
      ggplot2::aes(
        x = (start + end) / 2,
        y = 0.16,
        label = paste0("E", exon_number)
      ),
      color = "white", fontface = "bold", size = 3.1
    )
  
  if (nrow(labels)) {
    sashimi <- sashimi + ggplot2::geom_text(
      data = labels,
      ggplot2::aes(x = x, y = y, label = label, color = event),
      fontface = "bold", size = 3.0, show.legend = FALSE
    )
  }
  
  arc_ymax <- if (nrow(arcs)) max(arcs$y, na.rm = TRUE) + 0.25 else 1
  sashimi <- sashimi +
    ggplot2::scale_color_manual(
      values = event_colors,
      breaks = names(event_colors),
      labels = event_labels,
      drop = FALSE,
      name = "Junction type"
    ) +
    ggplot2::scale_linewidth(
      range = c(0.45, 2.4), guide = "none",
      trans = "sqrt"
    ) +
    ggplot2::coord_cartesian(
      xlim = c(x_min, x_max), ylim = c(0, max(1.0, arc_ymax)), expand = FALSE
    ) +
    ggplot2::scale_x_continuous(breaks = pretty(c(x_min, x_max), n = 7), expand = c(0, 0)) +
    ggplot2::labs(x = paste0("Construct coordinate (", seqname, ")"), y = NULL) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = 9),
      axis.title.x = ggplot2::element_text(face = "bold", margin = ggplot2::margin(t = 6)),
      plot.margin = ggplot2::margin(0, 6, 4, 6)
    )
  
  # Title/header.
  header <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0, label = sample,
      hjust = 0, vjust = 0.5, fontface = "bold", size = 5.2
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 6, 1, 6))
  
  # Assemble without patchwork/gridExtra dependency.
  pdf_file <- file.path(out_dir, paste0(sample, "_sashimi.pdf"))
  png_file <- file.path(out_dir, paste0(sample, "_sashimi.png"))
  
  # Use ggplot2's built-in multipage support through grid graphics.
  grDevices::cairo_pdf(pdf_file, width = pdf_width, height = pdf_height)
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 3, ncol = 1,
    heights = grid::unit(c(0.55, 1.55, 4.3), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(header, vp = grid::viewport(layout.pos.row = 1))
  print(coverage_plot, vp = grid::viewport(layout.pos.row = 2))
  print(sashimi, vp = grid::viewport(layout.pos.row = 3))
  grid::upViewport()
  grDevices::dev.off()
  
  png_device <- if (capabilities("cairo")) grDevices::png else grDevices::png
  png(png_file, width = round(pdf_width * 600), height = round(pdf_height * 600), res = 600, type = if (capabilities("cairo")) "cairo" else "Xlib")
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 3, ncol = 1,
    heights = grid::unit(c(0.55, 1.55, 4.3), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(header, vp = grid::viewport(layout.pos.row = 1))
  print(coverage_plot, vp = grid::viewport(layout.pos.row = 2))
  print(sashimi, vp = grid::viewport(layout.pos.row = 3))
  grid::upViewport()
  dev.off()
  
  list(
    sample = sample,
    junctions = classified,
    pdf = pdf_file,
    png = png_file
  )
}

# -----------------------------
# Process samples
# -----------------------------
results <- vector("list", length(samples))
names(results) <- samples
log_rows <- vector("list", length(samples))

for (i in seq_along(samples)) {
  sample <- samples[[i]]
  gtf <- find_matching_file(gtf_files, sample)
  sj <- find_matching_file(sj_files, sample)
  bam <- find_matching_file(bam_files, sample)
  
  if (is.na(gtf) || is.na(sj)) {
    msg <- paste0(
      "Missing required input: ",
      if (is.na(gtf)) "GTF " else "",
      if (is.na(sj)) "SJ.out.tab" else ""
    )
    warning(sample, ": ", msg)
    log_rows[[i]] <- data.frame(sample = sample, status = "SKIPPED", message = msg, stringsAsFactors = FALSE)
    next
  }
  
  msg <- ""
  
  res <- tryCatch(
    plot_sample(sample, gtf, sj, bam, output_dir),
    error = function(e) {
      msg <<- conditionMessage(e)
      cat("FAILED:", sample, "\n")
      cat("ERROR:", msg, "\n")
      NULL
    }
  )
  
  if (is.null(res)) {
    log_rows[[i]] <- data.frame(
      sample = sample,
      status = "FAILED",
      message = msg,
      stringsAsFactors = FALSE
    )
  } else {
    results[[sample]] <- res
    log_rows[[i]] <- data.frame(
      sample = sample,
      status = "OK",
      message = "",
      stringsAsFactors = FALSE
    )
  }
}

log_df <- do.call(rbind, log_rows)
write.table(
  log_df, file.path(output_dir, "processing_log.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Combined summary and heatmap
# -----------------------------
good <- results[!vapply(results, is.null, logical(1))]
if (length(good)) {
  summary <- do.call(rbind, lapply(good, function(x) x$junctions))
  rownames(summary) <- NULL
  write.table(
    summary, file.path(output_dir, "all_mutants_junction_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # Build labels from event + coordinates, so cryptic junctions are retained explicitly.
  summary$junction_label <- paste0(
    summary$event, ": ", summary$intron_start, "-", summary$intron_end
  )
  
  # Top junctions per mutant, retaining all junctions above threshold where practical.
  keep <- summary$unique_reads >= min_junction_reads
  heat <- summary[keep, , drop = FALSE]
  if (nrow(heat)) {
    # Keep the globally strongest junction labels, with a reasonable ceiling for readability.
    totals <- tapply(heat$unique_reads, heat$junction_label, sum, na.rm = TRUE)
    top_labels <- names(sort(totals, decreasing = TRUE))[seq_len(min(30L, length(totals)))]
    heat <- heat[heat$junction_label %in% top_labels, , drop = FALSE]
    
    # Use a simple deterministic pseudo-normalization for display only: log10(reads + 1).
    heat$value <- log10(heat$unique_reads + 1)
    heat$sample <- factor(heat$sample, levels = samples[samples %in% unique(heat$sample)])
    
    p_heat <- ggplot2::ggplot(
      heat,
      ggplot2::aes(x = junction_label, y = sample, fill = value)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(name = "log10\nunique reads + 1") +
      ggplot2::labs(
        x = NULL, y = NULL,
        title = "Splice junction usage across mutants"
      ) +
      ggplot2::theme_classic(base_size = 9) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1),
        axis.text.y = ggplot2::element_text(size = 8),
        axis.title = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold"),
        legend.title = ggplot2::element_text(face = "bold")
      )
    
    heat_h <- max(5, 0.25 * length(unique(heat$sample)) + 2)
    ggplot2::ggsave(
      file.path(output_dir, "all_mutants_junction_heatmap.pdf"),
      p_heat, width = 11, height = heat_h, device = grDevices::cairo_pdf
    )
    ggplot2::ggsave(
      file.path(output_dir, "all_mutants_junction_heatmap.png"),
      p_heat, width = 11, height = heat_h, dpi = 600
    )
  }
}

message("\nFinished.")
message("Output directory: ", output_dir)
message("Processing log: ", file.path(output_dir, "processing_log.tsv"))