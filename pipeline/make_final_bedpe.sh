#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <interaction_dir> <q30_peak_dir> <output_dir> <sample_manifest.tsv> <summary.tsv>" >&2
  exit 2
fi

interaction_dir=$(cd "$1" && pwd)
peak_dir=$(cd "$2" && pwd)
out_dir=$3
manifest=$(cd "$(dirname "$4")" && pwd)/$(basename "$4")
summary=$5

command -v bedtools >/dev/null || { echo "ERROR: bedtools is required" >&2; exit 1; }
mkdir -p "$out_dir" "$(dirname "$summary")"
printf 'Sample\tAll_cis_clusters\tiPET3_cis\tLeft_q30_supported\tRight_q30_supported\tBoth_q30_supported\tFinal_both_anchor_q30_iPET3\n' > "$summary"

tail -n +2 "$manifest" | while IFS=$'\t' read -r sample state genotype replicate peak_input peak_output interaction_input; do
  interaction="${interaction_dir}/${interaction_input}"
  peak="${peak_dir}/${peak_output}"
  sample_out="${out_dir}/${sample}"
  mkdir -p "$sample_out"
  [[ -s "$interaction" ]] || { echo "ERROR: missing interaction file $interaction" >&2; exit 1; }
  [[ -s "$peak" ]] || { echo "ERROR: missing peak file $peak" >&2; exit 1; }

  all_cis=$(wc -l < "$interaction")
  ipet3="${sample_out}/${sample}.cis.iPET3.bedpe"
  awk 'BEGIN { OFS="\t" } $7 >= 3 { print }' "$interaction" > "$ipet3"
  ipet3_n=$(wc -l < "$ipet3")

  left="${sample_out}/${sample}.left.bed"
  right="${sample_out}/${sample}.right.bed"
  awk 'BEGIN { OFS="\t" } { print $1,$2,$3,NR }' "$ipet3" > "$left"
  awk 'BEGIN { OFS="\t" } { print $4,$5,$6,NR }' "$ipet3" > "$right"

  left_hit="${sample_out}/${sample}.left.hit.bed"
  right_hit="${sample_out}/${sample}.right.hit.bed"
  bedtools intersect -u -a "$left" -b "$peak" > "$left_hit"
  bedtools intersect -u -a "$right" -b "$peak" > "$right_hit"
  left_n=$(wc -l < "$left_hit")
  right_n=$(wc -l < "$right_hit")

  left_ids="${sample_out}/${sample}.left.ids"
  right_ids="${sample_out}/${sample}.right.ids"
  both_ids="${sample_out}/${sample}.both.ids"
  cut -f4 "$left_hit" | sort -n -u > "$left_ids"
  cut -f4 "$right_hit" | sort -n -u > "$right_ids"
  awk 'NR == FNR { keep[$1] = 1; next } ($1 in keep) { print $1 }' \
    "$left_ids" "$right_ids" | sort -n -u > "$both_ids"
  both_n=$(wc -l < "$both_ids")

  final="${sample_out}/${sample}.q30.iPET3.bothAnchors.bedpe"
  awk 'NR == FNR { keep[$1] = 1; next } (FNR in keep) { print }' "$both_ids" "$ipet3" > "$final"
  final_n=$(wc -l < "$final")
  [[ "$final_n" -eq "$both_n" ]] || { echo "ERROR: count mismatch for $sample" >&2; exit 1; }

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$all_cis" "$ipet3_n" "$left_n" "$right_n" "$both_n" "$final_n" >> "$summary"
  rm -f "$ipet3" "$left" "$right" "$left_hit" "$right_hit" "$left_ids" "$right_ids" "$both_ids"
done

echo "Wrote final BEDPE files to $out_dir"
echo "Wrote summary to $summary"
