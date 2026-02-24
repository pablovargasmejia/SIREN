#!/usr/bin/env python3
import argparse, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--sirna_size", type=int, default=21)
    ap.add_argument("--sensitivity", choices=["high","medium"], default="high")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    other = os.path.join(args.outdir, "other_files")
    os.makedirs(other, exist_ok=True)

    records = list(SeqIO.parse(args.targets, "fasta"))
    matches = [r for r in records if args.gene in r.id]
    if len(matches) == 0:
        raise SystemExit(f"Gene '{args.gene}' not found in {args.targets}")
    if len(matches) > 1:
        raise SystemExit("More than one match found for gene; use a more specific string.\n" +
                         "\n".join([m.id for m in matches[:20]]))

    target = matches[0]
    off_targets = [r for r in records if r.id != target.id]

    target_fa = os.path.join(other, "target.fa")
    SeqIO.write(target, target_fa, "fasta")

    seqs_fa = os.path.join(other, "sequences.fa")
    SeqIO.write(off_targets, seqs_fa, "fasta")

    step = 1 if args.sensitivity == "high" else 2
    sirnas = []
    seq = target.seq
    for i in range(0, len(seq) - args.sirna_size + 1, step):
        sirna = seq[i:i+args.sirna_size]
        sirna_r = sirna.reverse_complement().transcribe()
        start = i+1
        end = start + args.sirna_size - 1
        name = f"sirna_{start}-{end}"
        sirnas.append(SeqRecord(sirna, id=name, description=""))
        sirnas.append(SeqRecord(sirna_r, id=f"{name}_r", description=""))

    sirnas_fa = os.path.join(args.outdir, "sirnas.fa")
    SeqIO.write(sirnas, sirnas_fa, "fasta")

if __name__ == "__main__":
    main()
