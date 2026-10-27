# Pipeline inputs

The main CPU wrapper takes two arguments:

```bash
bash CPUpipeline_slurm.sh config.sh fastq_list.tsv
```

`fastq_list.tsv` contains two tab-separated columns: the FASTQ prefix (without the R1/R2 suffix) and the output run identifier. The configuration must define the mm10 reference, CPU executable, helper scripts, and the values reported in the manuscript: Trim Galore quality 20, linker threshold 18, mapping score threshold 20, 8-kb self-ligation cutoff, MACS2 q-value 1e-5, 250-bp extension, and iPET threshold 3.
