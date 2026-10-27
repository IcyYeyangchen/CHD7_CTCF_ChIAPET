#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_figure2c.R <source_data_dir> <output_dir>", call. = FALSE)
}
source_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_arg), fixed = FALSE)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
source(file.path(script_dir, "figure2_plot_common.R"), local = TRUE)

interaction_overlap <- read.delim(
  file.path(source_dir, "fig2c_interaction_reciprocal_overlap.tsv"),
  stringsAsFactors = FALSE
)
p_c <- make_figure2_overlap_plot(interaction_overlap,
                                 "Reciprocal interaction overlap (%)", "c")
save_figure2_panel(p_c, "c", 6.35, 3.30, out_dir)
