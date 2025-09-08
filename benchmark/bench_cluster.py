#!/usr/bin/env python3
import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write(
        "ERROR: Biopython is required. On Cardinal run:\n"
        "  module load miniconda3\n"
        "  conda create -n siren_env -c conda-forge -c bioconda python=3.11 "
        "biopython pandas matplotlib numpy tqdm primer3-py rnahybrid\n"
        "  conda activate siren_env\n"
    )
    raise

def parse_args():
    p = argparse.ArgumentParser(
        description="Launch per-transcript SIREN runs as parallel Slurm jobs (OSC Cardinal)."
    )
    p.add_argument("--fasta", required=True, help="Input transcripts FASTA")
    p.add_argument("--outdir", required=True, help="Results root folder")
    p.add_argument("--threads", type=int, default=48, help="Threads for SIREN (and cpus-per-task)")
    p.add_argument("--cpus-per-task", type=int, default=48, help="Slurm --cpus-per-task (should match --threads)")
    p.add_argument("--account", default="PAS1755", help="Slurm account")
    p.add_argument("--time", dest="time_limit", default="12:00:00", help="Slurm time limit (HH:MM:SS)")
    p.add_argument("--partition", default=None, help="Optional Slurm partition")
    p.add_argument("--conda-module", default="miniconda3", help="OSC module providing conda")
    p.add_argument("--conda-env", default="siren_env", help="Conda environment name")
    p.add_argument("--siren_script", default="siren_masterV.py", help="Path to siren_masterV.py")
    p.add_argument(
        "--bins",
        default="500,1000,2000,4000,6000,8000",
        help="Comma-separated target lengths to sample (pick closest sequence to each bin)",
    )
    p.add_argument("--dry-run", action="store_true", help="Write sbatch files but do not submit")
    p.add_argument("--job-name-prefix", default="siren", help="Prefix for Slurm job names")
    # forward any unknown args directly to siren_masterV.py
    args, unknown = p.parse_known_args()
    args.forward_flags = unknown
    return args

def load_fasta_lengths(fasta_path: Path):
    lengths = {}
    with open(fasta_path, "r") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            lengths[rec.id] = len(rec.seq)
    return lengths

def pick_closest_by_bins(lengths_dict, bins):
    """
    For each bin length, pick 2 transcripts whose length is closest.
    Ensures unique picks overall.
    """
    items = list(lengths_dict.items())  # (id, length)
    chosen = []
    used_ids = set()
    for b in bins:
        sorted_by_proximity = sorted(items, key=lambda kv: abs(kv[1] - b))
        picked_count = 0
        for tid, L in sorted_by_proximity:
            if tid in used_ids:
                continue
            chosen.append((b, tid, L))
            used_ids.add(tid)
            picked_count += 1
            if picked_count >= 2:   # <<< hardcoded 2 picks per bin
                break
    return chosen

def ensure_dirs(path: Path):
    path.mkdir(parents=True, exist_ok=True)

def build_python_run_line(args, transcript_id, job_outdir: Path):
    """
    Build the command to run SIREN using host Python (conda env).
    RNAhybrid is expected to be available in the same conda env PATH.
    """
    cmd = [
        "python", args.siren_script,
        "--targets", str(args.fasta),
        "--gene", str(transcript_id),
        "--threads", str(args.threads),
        "--outdir", str(job_outdir),
    ]
    # forward any extra SIREN flags unchanged (e.g. -m windowed -k 9 -w 40 -H 2)
    cmd += args.forward_flags
    return " ".join(shlex.quote(x) for x in cmd)

def write_and_submit_sbatch(*, job_name, cpus, time_limit, account, partition,
                            out_log, err_log, conda_module, conda_env,
                            runtime_tsv, python_run_line, transcript_id, tx_len, sbatch_path):
    header = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --time={time_limit}
#SBATCH --account={account}
#SBATCH --output={out_log}
#SBATCH --error={err_log}
"""
    if partition:
        header += f"#SBATCH --partition={partition}\n"

    body = f"""
set -euo pipefail

# --- Load Conda and activate env on OSC Cardinal ---
module load {conda_module}
source activate {conda_env}

# Tie numerical libs to Slurm allocation (best-effort)
export OMP_NUM_THREADS={cpus}
export MKL_NUM_THREADS={cpus}
export OPENBLAS_NUM_THREADS={cpus}
export NUMEXPR_MAX_THREADS={cpus}

START_TS=$(date +%s)

echo "[ $(date) ] Starting {job_name}"
echo "Workdir: $(pwd)"
echo "Python: $(which python)"
echo "RNAhybrid: $(command -v RNAhybrid || true)"
echo "Command: {python_run_line}"

# Run SIREN
{python_run_line}

END_TS=$(date +%s)
ELAPSED=$((END_TS - START_TS))

# Append to shared runtime TSV with a file lock (one line per job)
RUNTIME_TSV="{runtime_tsv}"
mkdir -p "$(dirname "$RUNTIME_TSV")"
(
  flock -w 30 9
  printf "%s\\t%s\\t%s\\n" "{transcript_id}" "{tx_len}" "$ELAPSED" >> "$RUNTIME_TSV"
) 9>>"$RUNTIME_TSV"

echo "[ $(date) ] Done {job_name} in $ELAPSED s"
"""
    with open(sbatch_path, "w") as fh:
        fh.write(header + body)

    # submit
    return subprocess.run(["sbatch", sbatch_path], capture_output=True, text=True)

def main():
    args = parse_args()
    fasta = Path(args.fasta).resolve()
    outdir = Path(args.outdir).resolve()
    logs_dir = outdir / "logs"
    ensure_dirs(outdir)
    ensure_dirs(logs_dir)

    # Parse bins
    try:
        bin_targets = [int(x.strip()) for x in args.bins.split(",") if x.strip()]
    except Exception:
        sys.exit(f"Invalid --bins value: {args.bins}")

    # Load fasta lengths
    print(f">>> Scanning: {fasta.name}")
    lengths = load_fasta_lengths(fasta)

    # Pick closest per bin
    picks = pick_closest_by_bins(lengths, bin_targets)
    for b, tid, L in picks:
        print(f"Bin {b} nt: 1 selected from {fasta.name}\n  - {tid} ({L} nt)")

    runtime_tsv = outdir / "siren_runtime.tsv"
    if not runtime_tsv.exists():
        with open(runtime_tsv, "w") as fh:
            fh.write("transcript_id\tlength_nt\truntime_s\n")

    submitted = []
    for _, tid, L in picks:
        job_name = f"{args.job_name_prefix}_{tid}"
        job_outdir = outdir / tid
        ensure_dirs(job_outdir)

        out_log = logs_dir / f"{job_name}.out"
        err_log = logs_dir / f"{job_name}.err"
        sbatch_path = logs_dir / f"{job_name}.sbatch"

        python_run_line = build_python_run_line(args, tid, job_outdir)

        if args.dry_run:
            write_and_submit_sbatch(
                job_name=job_name,
                cpus=args.cpus_per_task,
                time_limit=args.time_limit,
                account=args.account,
                partition=args.partition,
                out_log=str(out_log),
                err_log=str(err_log),
                conda_module=args.conda_module,
                conda_env=args.conda_env,
                runtime_tsv=str(runtime_tsv),
                python_run_line=python_run_line,
                transcript_id=tid,
                tx_len=L,
                sbatch_path=str(sbatch_path),
            )
            print(f"[dry-run] Wrote: {sbatch_path}")
        else:
            res = write_and_submit_sbatch(
                job_name=job_name,
                cpus=args.cpus_per_task,
                time_limit=args.time_limit,
                account=args.account,
                partition=args.partition,
                out_log=str(out_log),
                err_log=str(err_log),
                conda_module=args.conda_module,
                conda_env=args.conda_env,
                runtime_tsv=str(runtime_tsv),
                python_run_line=python_run_line,
                transcript_id=tid,
                tx_len=L,
                sbatch_path=str(sbatch_path),
            )
            if res.returncode == 0:
                print(res.stdout.strip())
                submitted.append(job_name)
            else:
                sys.stderr.write(f"[sbatch error] {job_name}\n{res.stderr}\n")

    print(f"Submitted jobs: {len(submitted)}")
    for j in submitted:
        print(f"  - {j}")

if __name__ == "__main__":
    main()