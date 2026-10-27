#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_figure2.R <processed_file_dir> <output_dir>", call. = FALSE)
}
processed_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
source_dirs <- source_dir
write_source_table <- function(data, filename) {
  for (directory in unique(source_dirs)) {
    write.table(data, file.path(directory, filename), sep = "\t",
                quote = FALSE, row.names = FALSE)
  }
}

states <- c("mESC", "neuron")
genotypes <- c("WT", "SF", "KO")
replicates <- c("TR1", "TR2")
conditions <- expand.grid(state = states, genotype = genotypes,
                          stringsAsFactors = FALSE)
conditions <- conditions[order(match(conditions$state, states),
                               match(conditions$genotype, genotypes)), ]
conditions$condition <- paste(conditions$state, conditions$genotype, sep = "_")

interaction_path <- function(state, genotype, replicate) {
  file.path(processed_dir, "final_interactions",
            sprintf("%s_%s_CTCF_%s_merged.CTCF.q30.iPET3.bothAnchors.bedpe",
                    state, genotype, replicate))
}
peak_path <- function(state, genotype, replicate) {
  file.path(processed_dir, "peaks_q1e-30",
            sprintf("%s_%s_CTCF_%s_merged_peaks.q1e-30.narrowPeak",
                    state, genotype, replicate))
}
read_interactions <- function(path) {
  x <- read.delim(path, header = FALSE, sep = "\t", comment.char = "",
                  quote = "", stringsAsFactors = FALSE,
                  colClasses = c("character", "integer", "integer",
                                 "character", "integer", "integer", "integer"))
  if (ncol(x) != 7L) stop("Expected 7 BEDPE columns: ", path)
  names(x) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "iPET")
  x$mid1 <- floor((as.double(x$start1) + x$end1) / 2)
  x$mid2 <- floor((as.double(x$start2) + x$end2) / 2)
  swap <- x$chr1 > x$chr2 | (x$chr1 == x$chr2 & x$mid1 > x$mid2)
  if (any(swap)) {
    for (pair in list(c("chr1", "chr2"), c("start1", "start2"),
                      c("end1", "end2"), c("mid1", "mid2"))) {
      tmp <- x[[pair[[1]]]][swap]
      x[[pair[[1]]]][swap] <- x[[pair[[2]]]][swap]
      x[[pair[[2]]]][swap] <- tmp
    }
  }
  x
}
read_peaks <- function(path) {
  x <- read.delim(path, header = FALSE, sep = "\t", comment.char = "",
                  quote = "", stringsAsFactors = FALSE,
                  colClasses = c("character", "integer", "integer", "character",
                                 "integer", "character", "numeric", "numeric",
                                 "numeric", "integer"))
  if (ncol(x) < 3L) stop("Expected at least 3 narrowPeak columns: ", path)
  GRanges(seqnames = x[[1]], ranges = IRanges(start = x[[2]] + 1L, end = x[[3]]))
}

# Fig. 2a: calculate replicate-specific mean iPET in 10-kb bins, then average
# the replicate means within each genotype and cell state, as described.
bin_rows <- list()
for (i in seq_len(nrow(conditions))) {
  state <- conditions$state[[i]]
  genotype <- conditions$genotype[[i]]
  condition <- conditions$condition[[i]]
  for (replicate in replicates) {
    path <- interaction_path(state, genotype, replicate)
    if (!file.exists(path)) stop("Missing interaction file: ", path)
    x <- read_interactions(path)
    span <- abs(x$mid2 - x$mid1)
    keep <- span >= 8000 & span <= 1000000
    bin_start <- floor(span[keep] / 10000) * 10000
    sums <- tapply(x$iPET[keep], bin_start, sum)
    counts <- tapply(x$iPET[keep], bin_start, length)
    means <- sums / counts
    bin_rows[[length(bin_rows) + 1L]] <- data.frame(
      state = state, genotype = genotype, condition = condition,
      replicate = replicate, bin_start_bp = as.integer(names(means)),
      bin_end_bp = as.integer(names(means)) + 10000L,
      interactions_in_bin = as.integer(counts), mean_iPET = as.numeric(means),
      stringsAsFactors = FALSE
    )
  }
}
bin_by_replicate <- do.call(rbind, bin_rows)
bin_by_replicate <- bin_by_replicate[order(bin_by_replicate$condition,
                                           bin_by_replicate$bin_start_bp,
                                           bin_by_replicate$replicate), ]
write_source_table(bin_by_replicate,
                   "fig2a_mean_iPET_by_replicate_and_10kb_bin.tsv")

group_key <- paste(bin_by_replicate$condition, bin_by_replicate$bin_start_bp, sep = "\r")
grouped <- split(seq_len(nrow(bin_by_replicate)), group_key)
mean_rows <- lapply(grouped, function(ix) {
  d <- bin_by_replicate[ix, ]
  data.frame(state = d$state[[1]], genotype = d$genotype[[1]],
             condition = d$condition[[1]], bin_start_bp = d$bin_start_bp[[1]],
             bin_end_bp = d$bin_end_bp[[1]],
             n_replicates_present = length(unique(d$replicate)),
             mean_iPET = mean(d$mean_iPET), stringsAsFactors = FALSE)
})
mean_by_condition <- do.call(rbind, mean_rows)
mean_by_condition <- mean_by_condition[order(mean_by_condition$condition,
                                               mean_by_condition$bin_start_bp), ]
write_source_table(mean_by_condition,
                   "fig2a_mean_iPET_by_condition_and_10kb_bin.tsv")

loess_rows <- list()
for (condition in unique(mean_by_condition$condition)) {
  d <- mean_by_condition[mean_by_condition$condition == condition, ]
  fit <- loess(mean_iPET ~ bin_start_bp, data = d, span = 0.35,
               degree = 2, family = "gaussian", surface = "direct")
  grid <- seq(min(d$bin_start_bp), max(d$bin_start_bp), length.out = 250L)
  smooth <- predict(fit, newdata = data.frame(bin_start_bp = grid))
  parts <- strsplit(condition, "_", fixed = TRUE)[[1]]
  loess_rows[[length(loess_rows) + 1L]] <- data.frame(
    state = parts[[1]], genotype = parts[[2]], condition = condition,
    span_bp = grid, mean_iPET_loess = as.numeric(smooth),
    stringsAsFactors = FALSE
  )
}
loess_data <- do.call(rbind, loess_rows)
write_source_table(loess_data, "fig2a_loess_curve_coordinates.tsv")

# Reciprocal peak overlap is summarized by the median of the two directional
# fractions, with an overlap defined by at least one shared genomic base.
peak_rows <- list()
for (i in seq_len(nrow(conditions))) {
  state <- conditions$state[[i]]
  genotype <- conditions$genotype[[i]]
  condition <- conditions$condition[[i]]
  p1 <- read_peaks(peak_path(state, genotype, "TR1"))
  p2 <- read_peaks(peak_path(state, genotype, "TR2"))
  overlap12 <- sum(countOverlaps(p1, p2, minoverlap = 1L) > 0L)
  overlap21 <- sum(countOverlaps(p2, p1, minoverlap = 1L) > 0L)
  peak_rows[[length(peak_rows) + 1L]] <- data.frame(
    state = state, genotype = genotype, condition = condition,
    direction = c("Rep1_to_Rep2", "Rep2_to_Rep1"),
    query_peak_count = c(length(p1), length(p2)),
    query_peaks_with_overlap = c(overlap12, overlap21),
    reciprocal_overlap_percent = 100 * c(overlap12 / length(p1), overlap21 / length(p2)),
    stringsAsFactors = FALSE
  )
}
peak_overlap <- do.call(rbind, peak_rows)
peak_medians <- aggregate(reciprocal_overlap_percent ~ state + genotype + condition,
                          data = peak_overlap, FUN = median)
names(peak_medians)[[4]] <- "median_reciprocal_overlap_percent"
peak_overlap <- merge(peak_overlap, peak_medians,
                      by = c("state", "genotype", "condition"), sort = FALSE)
write_source_table(peak_overlap, "fig2b_peak_reciprocal_overlap.tsv")

# A paired interaction matches when both corresponding anchor midpoints lie
# within 5 kb. findOverlaps generates candidate pairs independently for each
# anchor; intersecting those pairs enforces the two-anchor rule.
directional_interaction_matches <- function(query, subject, tolerance = 5000L) {
  q1 <- GRanges(query$chr1,
                IRanges(pmax(1, query$mid1 - tolerance), query$mid1 + tolerance))
  s1 <- GRanges(subject$chr1, IRanges(subject$mid1, subject$mid1))
  q2 <- GRanges(query$chr2,
                IRanges(pmax(1, query$mid2 - tolerance), query$mid2 + tolerance))
  s2 <- GRanges(subject$chr2, IRanges(subject$mid2, subject$mid2))
  h1 <- findOverlaps(q1, s1, ignore.strand = TRUE)
  h2 <- findOverlaps(q2, s2, ignore.strand = TRUE)
  key1 <- as.double(queryHits(h1)) + (as.double(subjectHits(h1)) - 1) * nrow(query)
  key2 <- as.double(queryHits(h2)) + (as.double(subjectHits(h2)) - 1) * nrow(query)
  matched_keys <- intersect(key1, key2)
  if (length(matched_keys) == 0L) return(0L)
  matched_queries <- unique(as.integer(((matched_keys - 1) %% nrow(query)) + 1))
  length(matched_queries)
}
interaction_rows <- list()
for (i in seq_len(nrow(conditions))) {
  state <- conditions$state[[i]]
  genotype <- conditions$genotype[[i]]
  condition <- conditions$condition[[i]]
  x1 <- read_interactions(interaction_path(state, genotype, "TR1"))
  x2 <- read_interactions(interaction_path(state, genotype, "TR2"))
  match12 <- directional_interaction_matches(x1, x2)
  match21 <- directional_interaction_matches(x2, x1)
  interaction_rows[[length(interaction_rows) + 1L]] <- data.frame(
    state = state, genotype = genotype, condition = condition,
    direction = c("Rep1_to_Rep2", "Rep2_to_Rep1"),
    query_interaction_count = c(nrow(x1), nrow(x2)),
    query_interactions_with_match = c(match12, match21),
    reciprocal_overlap_percent = 100 * c(match12 / nrow(x1), match21 / nrow(x2)),
    stringsAsFactors = FALSE
  )
  message("Computed interaction overlap: ", condition)
}
interaction_overlap <- do.call(rbind, interaction_rows)
interaction_medians <- aggregate(reciprocal_overlap_percent ~ state + genotype + condition,
                                 data = interaction_overlap, FUN = median)
names(interaction_medians)[[4]] <- "median_reciprocal_overlap_percent"
interaction_overlap <- merge(interaction_overlap, interaction_medians,
                             by = c("state", "genotype", "condition"), sort = FALSE)
write_source_table(interaction_overlap,
                   "fig2c_interaction_reciprocal_overlap.tsv")

script_arg <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_arg), fixed = FALSE)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
rscript <- file.path(R.home("bin"), "Rscript")
for (letter in c("a", "b", "c")) {
  panel_script <- file.path(script_dir, paste0("plot_figure2", letter, ".R"))
  panel_args <- vapply(c(panel_script, source_dir, out_dir), shQuote, character(1))
  status <- system2(rscript, args = panel_args, stdout = "", stderr = "")
  if (!is.null(status) && status != 0L) {
    stop("Panel plotting failed for Figure 2", letter, call. = FALSE)
  }
}
message("Wrote independent Figure 2 panels and source tables to: ", normalizePath(out_dir))
