# CTCF ChIA-PET processing and Figure 2 reproduction

This repository contains the custom scripts used to process the 12 CTCF in situ ChIA-PET libraries and to reproduce Figure 2 of:

`CTCF-associated chromatin interaction maps across Chd7 genotypes during mouse neuronal differentiation`

The package is organized as follows:

- `pipeline/`: Slurm/CPU processing and post-processing scripts.
- `source_data/`: tab-separated source tables used for the displayed Figure 2 panels.
- `plot_figure2.R`, `plot_figure2a.R`, `plot_figure2b.R`, `plot_figure2c.R`, `figure2_plot_common.R`: Figure 2 source-table generation and plotting.
- `sample_manifest.tsv`: the 12 CTCF samples included in the manuscript, with expected peak and interaction file names.

## Figure 2

The panel scripts require R 4.6.0 or later and the packages `ggplot2`, `ragg`, `GenomicRanges`, and `IRanges`.

To regenerate Figure 2 source tables and panels from processed BEDPE and narrowPeak files:

```bash
Rscript plot_figure2.R /path/to/processed_data /path/to/figure2_output
```

The processed-data directory must contain `final_interactions/` and `peaks_q1e-30/` subdirectories with the filenames specified in `sample_manifest.tsv`. The script expects mm10/GRCm38 coordinates and BED-style 0-based, half-open intervals.

The checked-in `source_data/` files are the exact tables used for the submitted Figure 2. They allow the panels to be redrawn without access to the large processed genomic files. A full run writes regenerated tables to `<figure2_output>/source_data/`.

## Peak and interaction post-processing

`pipeline/filter_q30_peaks.sh` filters each MACS2 narrowPeak file at column 9 >= 30, corresponding to q <= 1e-30. `pipeline/make_final_bedpe.sh` retains cis interaction clusters with iPET >= 3 and requires both anchors to overlap a q30 CTCF peak.

Both scripts accept explicit input/output directories and the manifest, so they do not depend on `/public/groups/...` paths. They process only the 12 CTCF samples in this dataset.

## Raw ChIA-PET processing

`pipeline/CPUpipeline_slurm_zyf.sh` is the lab's Slurm wrapper around ChIA-PET Utilities (CPU), Trim Galore, MACS2, samtools, bedtools, and UCSC utilities. It requires the CPU executable, mm10 reference FASTA/index, container image, and site-specific annotation paths. These external dependencies are not redistributed here. Before execution, copy `pipeline/config.example.sh` to a local configuration file and replace all site-specific paths.

The script is included to document the executed processing workflow; the deposited processed files and the software versions/parameters in the manuscript remain the authoritative reusable outputs.

## Data and security

No FASTQ files, access credentials, SSH keys, or private keys are included in this repository. In particular, do not upload the `aspsub_rsa` file from the working directory. If that key has been exposed outside the local machine, revoke and replace it before submission.
