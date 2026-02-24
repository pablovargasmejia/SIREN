#!/usr/bin/env python3
import argparse, os
from Bio import SeqIO

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--num_shards", type=int, required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    recs = list(SeqIO.parse(args.input, "fasta"))
    if not recs:
        raise SystemExit("Empty FASTA")

    lengths = [len(r.seq) for r in recs]
    order = sorted(range(len(recs)), key=lambda i: lengths[i], reverse=True)
    k = max(1, min(args.num_shards, len(recs)))

    buckets = [[] for _ in range(k)]
    loads = [0]*k
    for i in order:
        j = min(range(k), key=lambda x: loads[x])
        buckets[j].append(recs[i])
        loads[j] += lengths[i]

    for idx, chunk in enumerate(buckets, 1):
        if not chunk:
            continue
        out = os.path.join(args.outdir, f"sequences_split_{idx}.fa")
        SeqIO.write(chunk, out, "fasta")

if __name__ == "__main__":
    main()
