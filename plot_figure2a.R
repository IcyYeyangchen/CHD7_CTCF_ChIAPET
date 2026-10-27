#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript plot_figure2a.R <source_data_dir> <output_dir>", call. = FALSE)
}
source_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
script_arg <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_path <- gsub("~\\+~", " ", sub("^--file=", "", script_arg), fixed = FALSE)
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
source(file.path(script_dir, "figure2_plot_common.R"), local = TRUE)

loess_data <- read.delim(file.path(source_dir, "fig2a_loess_curve_coordinates.tsv"),
                         stringsAsFactors = FALSE)
loess_data$color_group <- factor(loess_data$condition,
                                 levels = names(figure2_palette))
p_a <- ggplot(loess_data, aes(x = span_bp / 1000, y = mean_iPET_loess,
                              colour = color_group, group = condition)) +
  geom_line(linewidth = 0.75, na.rm = TRUE) +
  facet_wrap(~state, ncol = 1, labeller = as_labeller(figure2_facet_labels)) +
  scale_colour_manual(values = figure2_palette, breaks = names(figure2_palette),
                      labels = c("WT (mESC)", "WT (neuron)", "SF (mESC)",
                                 "SF (neuron)", "KO (mESC)", "KO (neuron)")) +
  scale_x_continuous(limits = c(0, 1000), breaks = c(0, 250, 500, 750, 1000),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Interaction span (kb)", y = "Mean iPET support", colour = NULL,
       tag = "a") +
  figure2_theme + theme(legend.position = "right",
                        legend.text = element_text(size = 8),
                        plot.tag = element_text(face = "bold", size = 13),
                        plot.tag.position = c(0, 1.02),
                        plot.margin = margin(14, 5, 4, 5),
                        panel.spacing = unit(0.5, "lines"))

save_figure2_panel(p_a, "a", 3.15, 6.35, out_dir)
