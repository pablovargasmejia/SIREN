#!/usr/bin/env python3
import argparse
import subprocess
import sys


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run RNAhybrid on one targets shard vs siRNAs and write stdout+stderr to a file."
    )
    ap.add_argument("--targets_shard", required=True, help="FASTA of targets (one shard)")
    ap.add_argument("--sirnas", required=True, help="FASTA of siRNA queries")
    ap.add_argument("--out", required=True, help="Output file path (stdout+stderr will be written here)")

    # IMPORTANT: use REMAINDER so options like -e -25 -v 0 ... are captured
    ap.add_argument(
        "--opts",
        nargs=argparse.REMAINDER,
        default=[],
        help="Extra RNAhybrid options passed as-is. Must be the LAST argument, e.g. "
             "--opts -e -25 -v 0 -u 0 -f 2,7 -p 0.01 -d 0.5,0.1 -m 60000",
    )

    args = ap.parse_args()

    cmd = ["RNAhybrid", "-t", args.targets_shard, "-q", args.sirnas] + (args.opts or [])

    with open(args.out, "w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)

    # Make failures visible to Snakemake:
    if proc.returncode != 0:
        return proc.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())