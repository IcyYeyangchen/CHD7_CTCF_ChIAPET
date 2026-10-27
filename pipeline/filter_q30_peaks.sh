#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <raw_peak_dir> <output_dir> <sample_manifest.tsv> <summary.tsv>" >&2
  exit 2
fi

raw_dir=$(cd "$1" && pwd)
out_dir=$2
manifest=$(cd "$(dirname "$3")" && pwd)/$(basename "$3")
summary=$4

mkdir -p "$out_dir" "$(dirname "$summary")"
printf 'Sample\tRaw_q1e-5\tQ30_peaks\tPct_retained\n' > "$summary"

tail -n +2 "$manifest" | while IFS=$'\t' read -r sample state genotype replicate peak_input peak_output interaction_input; do
  input="${raw_dir}/${peak_input}"
  output="${out_dir}/${peak_output}"
  [[ -s "$input" ]] || { echo "ERROR: missing input $input" >&2; exit 1; }

  awk 'BEGIN { OFS="\t" } $9 >= 30 { print }' "$input" > "$output"
  raw=$(wc -l < "$input")
  q30=$(wc -l < "$output")
  pct=$(awk -v n="$q30" -v d="$raw" 'BEGIN { if (d == 0) print "0.00"; else printf "%.2f", 100*n/d }')
  printf '%s\t%s\t%s\t%s\n' "$sample" "$raw" "$q30" "$pct" >> "$summary"
done

echo "Wrote q30 peak files to $out_dir"
echo "Wrote summary to $summary"
