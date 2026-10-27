# Copy this file to a private site-specific config before running the CPU wrapper.
# Do not commit credentials or local cluster paths.

R1suf="_R1_001.fastq.gz"
R2suf="_R2_001.fastq.gz"
genome="mm10"
extbp=500
selfbp=8000
macsq=1e-5
organism="mm"
peakext=250
PET=3
NTHREAD=16
memMb=64000
minMapScore=20
read_id_run_id="CHD7_CTCF"
sif="/path/to/cpu_3d.sif"
conda_activate="/path/to/conda/bin/activate"
trim_env="ATAC"
python_env_activate="/path/to/conda/envs/chiapet-py27/bin/activate"

# Site-specific paths. Replace every value before use.
fasta="/path/to/mm10.fa"
cpuprog="/path/to/cpu"
chiasigprog="/path/to/ChiaSig"
annot="/path/to/mm10.Genes.itxAnnotation.bed"
annotPerl="/path/to/itxAnnotation.pl"
extractsummary="/path/to/extract_summary_cpustat.sh"
Rdumpitx="/path/to/itx_from_annotated.R"
Rdumpitx_bedpe="/path/to/bedpe_itx_from_annotated.R"
addpath="/path/to/helper/scripts"
bin_dir="/path/to/bin"
BedGraphToBigWig="/path/to/bedGraphToBigWig"
BedSort="/path/to/bedSort"
