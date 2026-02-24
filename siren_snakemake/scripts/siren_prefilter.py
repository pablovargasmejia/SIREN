#!/usr/bin/env python3
"""
siren_prefilter.py
------------------
K-mer–based prefilter to shrink the RNAhybrid search database.

Two modes:
1) set (default)  : similarity of DISTINCT k-mer sets (dice/jaccard/containment_*)
2) windowed       : keep if there are >= N reverse-complement seed k-mer hits
                    within any W-bp window (local density of seed matches)

The windowed mode is tailored for RNAi: it retains sequences that share a dense
cluster of short exact reverse-complement seeds (e.g., k=9) with the gene,
which approximates local similarity without full alignment.
"""
import argparse
import os
import sys
from typing import Iterator, Tuple, List, Set
from Bio import SeqIO


DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def parse_args():
    p = argparse.ArgumentParser(
        prog="siren_prefilter",
        description="Alignment-free k-mer prefilter for SIREN."
    )
    p.add_argument("--targets", required=True, help="FASTA with all organism sequences.")
    p.add_argument("--gene", required=True, help="Gene ID or substring to select the target gene from --targets.")
    p.add_argument("--outdir", default=".", help="Output directory (default: current dir).")

    # Modes
    p.add_argument("--mode", choices=["set", "windowed"], default="set",
                   help="set: DISTINCT k-mer set similarity; windowed: local seed-density (default: set)")

    # Set-similarity params
    p.add_argument("--k", type=int, default=9, help="k-mer length for set mode (default: 9).")
    p.add_argument("--threshold", type=float, default=0.40, help="Similarity threshold (default: 0.40).")
    p.add_argument("--metric", choices=["dice", "jaccard", "containment_gene", "containment_target"],
                   default="dice", help="Similarity metric on DISTINCT k-mers (default: dice).")
    p.add_argument("--min_common_kmers", type=int, default=0,
                   help="Require at least this many shared DISTINCT k-mers (default: 0).")

    # Windowed seed-density params
    p.add_argument("--seed_k", type=int, default=9, help="Seed k-mer length for windowed mode (default: 9).")
    p.add_argument("--window_size", type=int, default=40, help="Window size in bp (default: 40).")
    p.add_argument("--min_window_hits", type=int, default=2, help="Minimum seed hits in any window (default: 2).")

    # Misc
    p.add_argument("--strand", choices=["rc", "fwd", "both"], default="rc",
                   help="Which gene strand to seed from (default: rc for complementarity).")
    p.add_argument("--case_sensitive", action="store_true", help="Match gene ID substring case-sensitively (default: off).")
    p.add_argument("--write_log", action="store_true", help="Write a TSV log with per-record metrics.")
    return p.parse_args()

def fasta_iter(path: str):
    for rec in SeqIO.parse(path, "fasta"):
        yield rec.id, rec.description, str(rec.seq)

def rc(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]

def clean_dna(seq: str) -> str:
    s = seq.upper()
    return "".join(ch if ch in "ACGTN" else "N" for ch in s)

def kmerset(seq: str, k: int) -> Set[str]:
    S = set()
    n = len(seq)
    if n < k:
        return S
    for i in range(n - k + 1):
        kmer = seq[i:i+k]
        if "N" in kmer:
            continue
        S.add(kmer)
    return S

def dice_similarity(A: Set[str], B: Set[str]) -> float:
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    inter = len(A & B)
    return (2.0 * inter) / (len(A) + len(B))

def jaccard_similarity(A: Set[str], B: Set[str]) -> float:
    if not A and not B:
        return 1.0
    if not A or not B:
        return 0.0
    inter = len(A & B)
    union = len(A | B)
    return inter / union if union else 0.0

def containment_similarity(A: Set[str], B: Set[str], denom: str) -> float:
    inter = len(A & B)
    if denom == "A":
        return inter / len(A) if A else 0.0
    else:  # denom == "B"
        return inter / len(B) if B else 0.0

def choose_gene_record(records, query: str, case_sensitive: bool=False):
    if not case_sensitive:
        q = query.lower()
        def low(s): return s.lower()
    else:
        q = query
        def low(s): return s

    for rid, desc, seq in records:
        if low(rid) == q:
            return rid, desc, seq
    for rid, desc, seq in records:
        if low(rid).startswith(q):
            return rid, desc, seq
    for rid, desc, seq in records:
        if q in low(rid):
            return rid, desc, seq
    for rid, desc, seq in records:
        if q in low(desc):
            return rid, desc, seq
    raise ValueError(f"Gene '{query}' not found in targets FASTA.")

def run_set_mode(args, records):
    if args.k < 6:
        print(f"[prefilter] WARNING: k={args.k} is very small (4^k={4**args.k}). "
              f"Expect similarity saturation; increase k to 8–11 for meaningful filtering.",
              file=sys.stderr)

    gene_rid, gene_desc, _ = choose_gene_record(records, args.gene, args.case_sensitive)
    gene_seq = None
    for rid, desc, seq in fasta_iter(args.targets):
        if rid == gene_rid:
            gene_seq = clean_dna(seq)
            break
    if gene_seq is None:
        print(f"ERROR: Chosen gene id '{gene_rid}' not found on second pass.", file=sys.stderr)
        sys.exit(3)

    k = args.k
    gene_set_f = kmerset(gene_seq, k) if args.strand in ("fwd","both") else set()
    gene_set_r = kmerset(rc(gene_seq), k) if args.strand in ("rc","both") else set()

    out_fa = os.path.join(args.outdir, "targets_prefiltered.fa")
    log_path = os.path.join(args.outdir, "prefilter_log.tsv") if args.write_log else None
    kept = 0
    total = 0

    def sim_metric(A, B):
        if args.metric == "dice":
            return dice_similarity(A, B)
        elif args.metric == "jaccard":
            return jaccard_similarity(A, B)
        elif args.metric == "containment_gene":
            return containment_similarity(A, B, "A")
        elif args.metric == "containment_target":
            return containment_similarity(A, B, "B")
        else:
            return dice_similarity(A, B)

    with open(out_fa, "w", encoding="utf-8") as out_fh:
        log_fh = open(log_path, "w", encoding="utf-8") if log_path else None
        if log_fh:
            log_fh.write("record_id\tsize_bp\tsim_forward\tsim_reverse\tmax_sim\tcommon_kmers\n")
        for rid, desc, seq in fasta_iter(args.targets):
            total += 1
            s = clean_dna(seq)
            S = kmerset(s, k)
            if not S:
                sim_f = 0.0
                sim_r = 0.0
                inter_f = inter_r = 0
            else:
                inter_f = len(S & gene_set_f) if gene_set_f else 0
                inter_r = len(S & gene_set_r) if gene_set_r else 0
                sim_f = sim_metric(S, gene_set_f) if gene_set_f else 0.0
                sim_r = sim_metric(S, gene_set_r) if gene_set_r else 0.0
            sim = max(sim_f, sim_r)
            common = max(inter_f, inter_r)
            if log_fh:
                log_fh.write(f"{rid}\t{len(s)}\t{sim_f:.4f}\t{sim_r:.4f}\t{sim:.4f}\t{common}\n")
            if (sim >= args.threshold and common >= args.min_common_kmers) or rid == gene_rid:
                kept += 1
                out_fh.write(f">{desc}\n")
                for i in range(0, len(seq), 60):
                    out_fh.write(seq[i:i+60] + "\n")

    print(f"[prefilter] {kept}/{total} records kept (mode=set, k={k}, threshold={args.threshold:.2f}, metric={args.metric}, min_common_kmers={args.min_common_kmers}).", file=sys.stderr)
    print(out_fa)

def run_windowed_mode(args, records):
    gene_rid, gene_desc, _ = choose_gene_record(records, args.gene, args.case_sensitive)
    gene_seq = None
    for rid, desc, seq in fasta_iter(args.targets):
        if rid == gene_rid:
            gene_seq = clean_dna(seq)
            break
    if gene_seq is None:
        print(f"ERROR: Chosen gene id '{gene_rid}' not found on second pass.", file=sys.stderr)
        sys.exit(3)

    seeds = set()
    if args.strand in ("rc", "both"):
        seeds |= kmerset(rc(gene_seq), args.seed_k)
    if args.strand in ("fwd", "both"):
        seeds |= kmerset(gene_seq, args.seed_k)

    out_fa = os.path.join(args.outdir, "targets_prefiltered.fa")
    log_path = os.path.join(args.outdir, "prefilter_log.tsv") if args.write_log else None
    kept = 0
    total = 0

    with open(out_fa, "w", encoding="utf-8") as out_fh:
        log_fh = open(log_path, "w", encoding="utf-8") if log_path else None
        if log_fh:
            log_fh.write("record_id\tsize_bp\tmax_window_hits\n")
        for rid, desc, seq in fasta_iter(args.targets):
            total += 1
            s = clean_dna(seq)
            # positions of exact seed matches in target
            k = args.seed_k
            pos = []
            if seeds and len(s) >= k:
                n = len(s)
                for i in range(n - k + 1):
                    kmer = s[i:i+k]
                    if "N" in kmer:
                        continue
                    if kmer in seeds:
                        pos.append(i)
            # compute local density
            left = 0
            best = 0
            for right in range(len(pos)):
                while pos[right] - pos[left] > args.window_size:
                    left += 1
                span = right - left + 1
                if span > best:
                    best = span
            if log_fh:
                log_fh.write(f"{rid}\t{len(s)}\t{best}\n")
            if best >= args.min_window_hits or rid == gene_rid:
                kept += 1
                out_fh.write(f">{desc}\n")
                for i in range(0, len(seq), 60):
                    out_fh.write(seq[i:i+60] + "\n")

    print(f"[prefilter] {kept}/{total} records kept (mode=windowed, seed_k={args.seed_k}, window_size={args.window_size}, min_window_hits={args.min_window_hits}).", file=sys.stderr)
    print(out_fa)

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    headers = []
    for rid, desc, seq in fasta_iter(args.targets):
        headers.append((rid, desc, ""))

    if not headers:
        print("ERROR: No records found in --targets.", file=sys.stderr)
        sys.exit(2)

    if args.mode == "set":
        run_set_mode(args, headers)
    else:
        run_windowed_mode(args, headers)

if __name__ == "__main__":
    main()
