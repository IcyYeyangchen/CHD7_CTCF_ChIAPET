#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_figure2b.R <source_data_dir> <output_dir>", call. = FALSE)
}
source_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_arg), fixed = FALSE)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
source(file.path(script_dir, "figure2_plot_common.R"), local = TRUE)

peak_overlap <- read.delim(file.path(source_dir, "fig2b_peak_reciprocal_overlap.tsv"),
                           stringsAsFactors = FALSE)
p_b <- make_figure2_overlap_plot(peak_overlap, "Reciprocal peak overlap (%)", "b")
save_figure2_panel(p_b, "b", 6.35, 3.30, out_dir)
