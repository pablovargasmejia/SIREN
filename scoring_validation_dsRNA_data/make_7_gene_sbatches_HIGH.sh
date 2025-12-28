#!/usr/bin/env bash
set -euo pipefail

# =======================
# CONFIG
# =======================
BASE="/fs/scratch/PAS1755/Pablo/SIREN/siren17/kulkarni"
OUTDIR="${BASE}/sbatches"
SLURMDIR="${BASE}/slurm"
RESULTS="${BASE}/results"

# IMPORTANT: use UNCOMPRESSED FASTA (SIREN can't read .gz with Biopython SeqIO.parse)
TRANSCRIPTS_FA="${BASE}/dmel-all-transcript-r6.66.fasta"

CONDA_ENV="siren_env"

mkdir -p "${OUTDIR}" "${SLURMDIR}" "${RESULTS}"

# quick sanity
if [[ ! -s "${TRANSCRIPTS_FA}" ]]; then
  echo "ERROR: Missing ${TRANSCRIPTS_FA}"
  echo "Fix: gunzip -c dmel-all-transcript-r6.66.fasta.gz > dmel-all-transcript-r6.66.fasta"
  exit 1
fi

write_job () {
  local FBGN="$1"
  local GENESYM="$2"
  local A1="$3" L1="$4"
  local A2="$5" L2="$6"
  local A3="$7" L3="$8"

  local JOBFILE="${OUTDIR}/sirenbench_Dmel_r6.66_${GENESYM}_${FBGN}.sbatch"

  cat > "${JOBFILE}" <<EOF
#!/bin/bash
#SBATCH --job-name=sirenbench_Dmel_r6.66_${GENESYM}_${FBGN}_b001
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --time=08:00:00
#SBATCH --account=PAS1755
#SBATCH --output=${SLURMDIR}/%x-%j.out
#SBATCH --error=${SLURMDIR}/%x-%j.err

set -euo pipefail

# ===== ENV =====
set +u
module load miniconda3
conda activate ${CONDA_ENV}
set -u

# ===== INPUTS =====
TRANSCRIPTS_FA="${TRANSCRIPTS_FA}"
RESULTS="${RESULTS}"

FBGN="${FBGN}"
GENESYM="${GENESYM}"

AMPS=("${A1}" "${A2}" "${A3}")
TYPES=("Amp1" "Amp2" "Amp3")
LENS=(${L1} ${L2} ${L3})

echo "=== JOB START ==="
echo "TRANSCRIPTS_FA: \${TRANSCRIPTS_FA}"
echo "RESULTS:        \${RESULTS}"
echo "FBGN:           \${FBGN}"
echo "GENESYM:        \${GENESYM}"
echo

# ===== Pick UNIQUE target key for SIREN (-g) =====
# If parent=FBgn... matches >1 transcript, choose first FBtr to make -g unique.
N_MATCH=\$(grep -c "parent=\${FBGN}" "\${TRANSCRIPTS_FA}" || true)

TARGET_KEY="\${FBGN}"
TARGET_KIND="FBgn"

if [[ "\${N_MATCH}" -gt 1 ]]; then
  FBTR=\$(grep -m1 "parent=\${FBGN}" "\${TRANSCRIPTS_FA}" | sed -n 's/.*ID=\\(FBtr[0-9]\\+\\).*/\\1/p')
  if [[ -n "\${FBTR}" ]]; then
    TARGET_KEY="\${FBTR}"
    TARGET_KIND="FBtr"
  fi
fi

echo "Target selection:"
echo "  parent hits : \${N_MATCH}"
echo "  -g key      : \${TARGET_KEY} (\${TARGET_KIND})"
echo

# ===== Run SIREN 3 times with rnai_length from experimental amplicon lengths =====
for i in 0 1 2; do
  AMP="\${AMPS[\$i]}"
  TYPE="\${TYPES[\$i]}"
  LEN="\${LENS[\$i]}"

  RUN_OUTDIR="\${RESULTS}/\${GENESYM}_\${FBGN}/\${TYPE}_\${AMP}_rnai_len\${LEN}"
  mkdir -p "\${RUN_OUTDIR}"

  echo "Running SIREN:"
  echo "  gene_symbol : \${GENESYM}"
  echo "  flybase_id  : \${FBGN}"
  echo "  length_tag  : \${TYPE} (\${AMP})"
  echo "  rnai_length : \${LEN}"
  echo "  outdir      : \${RUN_OUTDIR}"
  echo

  # ===== HIGH sensitivity (valid: high, medium) =====
  SIREN \\
    -T "\${TRANSCRIPTS_FA}" \\
    -g "\${TARGET_KEY}" \\
    -t "\${SLURM_CPUS_PER_TASK}" \\
    -r "\${LEN}" \\
    -S high \\
    -o "\${RUN_OUTDIR}" \\
    -g_o

  echo "Done: \${RUN_OUTDIR}"
  echo "---------------------------------------------"
done

echo "=== JOB END ==="
EOF

  chmod +x "${JOBFILE}"
  echo "Wrote ${JOBFILE}"
}

# ===== Your 7 genes =====
write_job "FBgn0029957" "CG12155" "DRSC17826" 516 "DRSC31035" 491 "DRSC31036" 202
write_job "FBgn0050421" "CG30421" "DRSC04532" 499 "DRSC30735" 498 "DRSC30736" 214
write_job "FBgn0052791" "CG32791" "DRSC17892" 509 "DRSC31039" 455 "DRSC31040" 353
write_job "FBgn0263929" "CG3563"  "DRSC13053" 351 "DRSC30936" 263 "DRSC30937" 335
write_job "FBgn0013733" "shot"    "DRSC05459" 474 "DRSC30757" 340 "DRSC30758" 479
write_job "FBgn0025800" "Smox"    "DRSC18716" 506 "DRSC31057" 342 "DRSC31058" 362
write_job "FBgn0024277" "trio"    "DRSC08527" 469 "DRSC30819" 431 "DRSC30820" 580

echo
echo "All sbatches created in: ${OUTDIR}"
echo "Slurm logs will be in : ${SLURMDIR}"
echo "Results will be in   : ${RESULTS}"