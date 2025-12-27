#!/usr/bin/env python3
"""
siren_bench_random100.py

Randomly sample transcripts from one or more FASTA files and generate Slurm
sbatch scripts where EACH sbatch processes multiple transcripts.

Requested behavior (Pablo):
- Sample 100 transcripts from EACH FASTA in db/
- Group into sbatches of 5 transcripts
- Inside each sbatch: run those 5 transcripts in medium AND the same 5 in high
  (i.e., 10 SIREN runs per sbatch)
- Use: --rnai_length 300, --threads 96
- Save outputs into per-transcript folders suffixed with *_medium and *_high
- Submit jobs “10 to 10”: we generate a helper submit script that submits 10
  sbatches at a time (chunked submission).

Usage (from /fs/scratch/PAS1755/Pablo/SIREN/siren17):
  # 1) Generate sbatches + manifest + submit helper
  python siren_bench_random100.py --dbdir db

  # 2) Submit chunk 1 (first 10 sbatches)
  bash sbatch_bench_5tx/submit_10.sh 1

  # 3) Submit chunk 2 (next 10 sbatches)
  bash sbatch_bench_5tx/submit_10.sh 2
"""

from __future__ import annotations

import argparse
import os
import random
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple


FASTA_EXT_RE = re.compile(r"\.(fa|fasta|fna)(\.gz)?$", re.IGNORECASE)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def fasta_ids(fasta_path: Path) -> List[str]:
    """
    Extract transcript IDs from a FASTA by taking the first token after '>'.
    This is typically what you want for --gene partial header matching.
    """
    ids: List[str] = []
    opener = open
    if fasta_path.suffix == ".gz":
        import gzip  # local import
        opener = gzip.open  # type: ignore

    with opener(fasta_path, "rt") as fh:  # type: ignore[arg-type]
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    continue
                tx_id = header.split()[0]
                ids.append(tx_id)
    # de-duplicate while preserving order
    seen = set()
    out = []
    for x in ids:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def chunked(seq: Sequence[str], n: int) -> List[List[str]]:
    return [list(seq[i:i+n]) for i in range(0, len(seq), n)]


def sanitize_tag(s: str) -> str:
    # for filenames / job names
    s = s.replace(" ", "_")
    s = s.replace("/", "_")
    s = s.replace(":", "_")
    s = s.replace("|", "_")
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    return s.strip("_") or "x"


@dataclass
class SlurmConfig:
    partition: str = "cpu"
    account: str = "PAS1755"
    cpus: int = 96
    time_limit: str = "08:00:00"
    slurm_out: str = "/fs/scratch/PAS1755/Pablo/SIREN/siren17/slurm-%x-%j.out"
    slurm_err: str = "/fs/scratch/PAS1755/Pablo/SIREN/siren17/slurm-%x-%j.err"
    conda_module: str = "miniconda3"
    conda_env: str = "siren_env"


def build_sbatch(
    sbatch_path: Path,
    job_name: str,
    fasta_path: Path,
    transcripts: Sequence[str],
    results_dir: Path,
    rnai_length: int,
    threads: int,
    slurm: SlurmConfig,
):
    """
    Build one sbatch that runs:
      - the 5 transcripts in medium
      - then the same 5 in high
    Each transcript writes to its own folder: <dbtag>__<tx>_medium / _high
    """
    # embed transcripts in bash array (quoted safely)
    tx_array = " ".join([f'"{t}"' for t in transcripts])

    # db tag for output folder prefix
    dbtag = sanitize_tag(FASTA_EXT_RE.sub("", fasta_path.name))

    content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={slurm.partition}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={slurm.cpus}
#SBATCH --time={slurm.time_limit}
#SBATCH --account={slurm.account}
#SBATCH --output={slurm.slurm_out}
#SBATCH --error={slurm.slurm_err}

set -euo pipefail

# ===== ENV =====
set +u
module load {slurm.conda_module}
conda activate {slurm.conda_env}
set -u

TARGETS="{fasta_path}"
THREADS={threads}
RNAI_LEN={rnai_length}
BASE_OUT="{results_dir}"

# 5 transcripts for this job:
TRANSCRIPTS=({tx_array})

echo "[ $(date) ] Starting: $SLURM_JOB_NAME"
echo "[ $(date) ] Targets: $TARGETS"
echo "[ $(date) ] Transcripts: ${{TRANSCRIPTS[*]}}"
echo "[ $(date) ] Threads: $THREADS  RNAi_len: $RNAI_LEN"
echo

run_one () {{
  local tx="$1"
  local sens="$2"

  # sanitize tx for folder name
  local txsafe
  txsafe="$(echo "$tx" | sed 's/[^A-Za-z0-9._-]/_/g')"

  local outdir="${{BASE_OUT}}/{dbtag}__${{txsafe}}_${{sens}}"
  mkdir -p "$outdir"

  echo "[ $(date) ] SIREN: tx=$tx sens=$sens outdir=$outdir"
  SIREN \\
    --targets "$TARGETS" \\
    --gene "$tx" \\
    --threads "$THREADS" \\
    --sensitivity "$sens" \\
    --rnai_length "$RNAI_LEN" \\
    --outdir "$outdir"
  echo
}}

# Run MEDIUM first (requested), then HIGH for the same 5
for tx in "${{TRANSCRIPTS[@]}}"; do
  run_one "$tx" "medium"
done

for tx in "${{TRANSCRIPTS[@]}}"; do
  run_one "$tx" "high"
done

echo "[ $(date) ] Done: $SLURM_JOB_NAME"
"""
    sbatch_path.write_text(content)


def write_submit_helper(scripts_dir: Path, chunk_size: int = 10) -> Path:
    """
    Write a submit_10.sh helper that submits sbatches in chunks (10 to 10).
    """
    helper = scripts_dir / "submit_10.sh"
    helper.write_text(f"""#!/bin/bash
set -euo pipefail

# Submit sbatch files in chunks of {chunk_size}.
# Usage:
#   bash {scripts_dir.name}/submit_10.sh 1   # submits 1..{chunk_size}
#   bash {scripts_dir.name}/submit_10.sh 2   # submits {chunk_size+1}..{2*chunk_size}
#   bash {scripts_dir.name}/submit_10.sh 3   # ...

CHUNK="${{1:-1}}"
if ! [[ "$CHUNK" =~ ^[0-9]+$ ]]; then
  echo "ERROR: chunk must be an integer (got: $CHUNK)" >&2
  exit 2
fi

mapfile -t FILES < <(ls -1 "{scripts_dir}"/*.sbatch 2>/dev/null | sort)
TOTAL="${{#FILES[@]}}"

if [[ "$TOTAL" -eq 0 ]]; then
  echo "No .sbatch files found in {scripts_dir}" >&2
  exit 1
fi

START=$(( (CHUNK - 1) * {chunk_size} ))
END=$(( START + {chunk_size} - 1 ))

echo "Found $TOTAL sbatch files."
echo "Submitting chunk $CHUNK: indices $START..$END"
echo

for i in $(seq "$START" "$END"); do
  if [[ "$i" -ge "$TOTAL" ]]; then
    break
  fi
  f="${{FILES[$i]}}"
  echo "sbatch $f"
  sbatch "$f"
done
""")
    helper.chmod(0o755)
    return helper


def parse_args():
    ap = argparse.ArgumentParser(
        description="Sample transcripts from FASTA(s) and generate multi-transcript SIREN sbatches (5 tx/job; medium then high)."
    )
    ap.add_argument("--dbdir", default="db", help="Directory containing FASTA files (default: db)")
    ap.add_argument(
        "--fasta",
        nargs="+",
        default=[
            "Homo_sapiens.GRCh38.cdna.all.fa",
            "TAIR10_cdna_20101214.fasta",
            "Phyca11_transcripts.fasta",
        ],
        help="FASTA filenames within --dbdir (default: the 3 benchmark FASTAs)",
    )
    ap.add_argument("--n", type=int, default=100, help="Random transcripts per FASTA (default: 100)")
    ap.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")

    ap.add_argument("--tx-per-job", type=int, default=5, help="Transcripts per sbatch (default: 5)")
    ap.add_argument("--chunk-size", type=int, default=10, help="Submission chunk size for submit_10.sh (default: 10)")

    ap.add_argument("--scripts-dir", default="sbatch_bench_5tx", help="Where to write sbatch scripts (default: sbatch_bench_5tx)")
    ap.add_argument("--results-dir", default="bench_results", help="Base results folder (default: bench_results)")
    ap.add_argument("--manifest", default="bench_manifest.tsv", help="TSV manifest of sampled transcripts (default: bench_manifest.tsv)")

    ap.add_argument("--threads", type=int, default=96, help="SIREN --threads (default: 96)")
    ap.add_argument("--rnai-length", type=int, default=300, help="SIREN --rnai_length (default: 300)")

    # Slurm header knobs (match your template)
    ap.add_argument("--partition", default="cpu", help="Slurm partition (default: cpu)")
    ap.add_argument("--account", default="PAS1755", help="Slurm account (default: PAS1755)")
    ap.add_argument("--cpus", type=int, default=96, help="--cpus-per-task (default: 96)")
    ap.add_argument("--time-limit", default="08:00:00", help="Slurm time (default: 08:00:00)")
    ap.add_argument("--slurm-out", default="/fs/scratch/PAS1755/Pablo/SIREN/siren17/slurm-%x-%j.out", help="Slurm --output path")
    ap.add_argument("--slurm-err", default="/fs/scratch/PAS1755/Pablo/SIREN/siren17/slurm-%x-%j.err", help="Slurm --error path")
    ap.add_argument("--conda-module", default="miniconda3", help="Module to load (default: miniconda3)")
    ap.add_argument("--conda-env", default="siren_env", help="Conda env name (default: siren_env)")

    ap.add_argument(
        "--submit",
        action="store_true",
        help="Submit ALL generated sbatches immediately (NOT chunked). For chunked, use submit_10.sh.",
    )
    return ap.parse_args()


def main():
    args = parse_args()

    project_dir = Path.cwd()
    dbdir = (project_dir / args.dbdir).resolve()
    scripts_dir = (project_dir / args.scripts_dir).resolve()
    results_dir = (project_dir / args.results_dir).resolve()
    manifest_path = (project_dir / args.manifest).resolve()

    scripts_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    slurm = SlurmConfig(
        partition=args.partition,
        account=args.account,
        cpus=args.cpus,
        time_limit=args.time_limit,
        slurm_out=args.slurm_out,
        slurm_err=args.slurm_err,
        conda_module=args.conda_module,
        conda_env=args.conda_env,
    )

    rng = random.Random(args.seed)

    rows: List[Tuple[str, str, int, str, str]] = []
    sbatch_paths: List[Path] = []

    global_batch = 0

    for fasta_name in args.fasta:
        fasta_path = (dbdir / fasta_name).resolve()
        if not fasta_path.exists():
            eprint(f"ERROR: FASTA not found: {fasta_path}")
            sys.exit(1)

        all_ids = fasta_ids(fasta_path)
        if len(all_ids) == 0:
            eprint(f"ERROR: No FASTA headers found in: {fasta_path}")
            sys.exit(1)

        n = min(args.n, len(all_ids))

        # make sampling stable per FASTA but still controlled by --seed:
        # derive a per-file seed from the global seed and fasta name
        file_seed = (args.seed, fasta_name)
        local_rng = random.Random(str(file_seed))
        sampled = local_rng.sample(all_ids, n)

        batches = chunked(sampled, args.tx_per_job)

        dbtag = sanitize_tag(FASTA_EXT_RE.sub("", fasta_path.name))

        for i, txs in enumerate(batches, start=1):
            global_batch += 1

            job_name = sanitize_tag(f"sirenbench_{dbtag}_b{i:03d}")
            sbatch_file = scripts_dir / f"{global_batch:03d}__{dbtag}__b{i:03d}.sbatch"

            build_sbatch(
                sbatch_path=sbatch_file,
                job_name=job_name,
                fasta_path=fasta_path,
                transcripts=txs,
                results_dir=results_dir,
                rnai_length=args.rnai_length,
                threads=args.threads,
                slurm=slurm,
            )

            sbatch_paths.append(sbatch_file)
            rows.append((dbtag, fasta_path.name, i, ",".join(txs), sbatch_file.name))

    # write manifest
    with open(manifest_path, "w") as out:
        out.write("dbtag\tfasta\tbatch_in_db\ttranscripts\tsbatch\n")
        for dbtag, fasta, bidx, txs, sb in rows:
            out.write(f"{dbtag}\t{fasta}\t{bidx}\t{txs}\t{sb}\n")

    helper = write_submit_helper(scripts_dir, chunk_size=args.chunk_size)

    print("Generated:")
    print(f"  sbatch scripts dir: {scripts_dir}")
    print(f"  results dir:        {results_dir}")
    print(f"  manifest:           {manifest_path}")
    print(f"  submit helper:      {helper}")
    print(f"  total sbatches:     {len(sbatch_paths)}  (each runs {args.tx_per_job} tx in medium + high)")

    if args.submit:
        submitted = 0
        for sb in sbatch_paths:
            res = subprocess.run(["sbatch", str(sb)], capture_output=True, text=True)
            if res.returncode == 0:
                submitted += 1
                print(res.stdout.strip())
            else:
                sys.stderr.write(f"[sbatch error] {sb.name}\n{res.stderr}\n")
        print(f"Submitted jobs: {submitted}/{len(sbatch_paths)}")


if __name__ == "__main__":
    main()
