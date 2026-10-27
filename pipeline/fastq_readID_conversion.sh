#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> <output_prefix> <run_id>" >&2
  exit 2
fi

r1=$1
r2=$2
prefix=$3
run_id=$4

convert_read() {
  local input=$1
  local mate=$2
  zcat "$input" | awk -v run_id="$run_id" -v mate="$mate" '
    NR % 4 == 1 { printf "@%s:%010d %d:N:0:0\\n", run_id, int((NR-1)/4)+1, mate; next }
    NR % 4 == 3 { print "+"; next }
    { print }
  '
}

convert_read "$r1" 1 | gzip -c > "${prefix}_R1.fastq.gz"
convert_read "$r2" 2 | gzip -c > "${prefix}_R2.fastq.gz"
