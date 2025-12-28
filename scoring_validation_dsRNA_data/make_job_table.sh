#!/usr/bin/env bash
set -euo pipefail

IN_TSV="kulkarni_amplicons.tsv"
OUT_TSV="siren_jobs.tsv"

# Output columns:
# idx  FBgn  rnai_length  amplicon  amp_type  gene_symbol
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{next}
{
  gene=$1; amp=$2; seq=$3; amptype=$4; fbgn=$5;
  gsub(/[ \r\n\t]/,"",seq);
  len=length(seq);
  if(len<1){next}
  print fbgn, len, amp, amptype, gene
}' "${IN_TSV}" > "${OUT_TSV}"

echo "Wrote ${OUT_TSV}"
echo "N jobs: $(wc -l < "${OUT_TSV}")"
