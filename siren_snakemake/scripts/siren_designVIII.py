#!/usr/bin/env python
import pandas as pd
import argparse
from collections import defaultdict
import primer3
import csv
from multiprocessing import Pool
import math
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate RNAi sequences, score them, design primers with primer3, and output expected amplicon size."
    )
    parser.add_argument('--target_path', required=True, help="Path to the target.fa file.")
    parser.add_argument('--off_targets_summary_path', required=True, help="Path to the off_targets_summary.tsv file.")
    parser.add_argument('--rnai_seq_length', type=int, default=200, help="Base RNAi sequence length (default: 200).")
    parser.add_argument('--out', default="rna_sequences_with_scores_and_primers.tsv",
                        help="Output file name (TSV) for RNAi sequences, scores, primers, and amplicon sizes.")
    parser.add_argument('--threads', type=int, default=6, help="Number of threads (processes) to use for parallel scoring.")
    return parser.parse_args()

from Bio import SeqIO

def read_fasta(filepath):
    records = list(SeqIO.parse(filepath, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {filepath}")
    return str(records[0].seq)

def parse_sirna_name(sirna_name):
    sirna_name_clean = sirna_name.replace("sirna_", "").replace("_r", "")
    return [int(pos) for pos in sirna_name_clean.split("-")]

def generate_rnai_sequences(target_sequence, base_length):
    lengths = list(range(base_length - 50, base_length + 51, 50))
    rna_sequences = []
    for length in lengths:
        for start_pos in range(0, len(target_sequence) - length + 1, 4):
            rna_seq = target_sequence[start_pos:start_pos+length]
            rna_name = f"RNAi_seq{len(rna_sequences)+1}_{start_pos+1}-{start_pos+length}_len{length}"
            rna_sequences.append((rna_name, rna_seq, length))
    return rna_sequences

def score_chunk(args):
    sub_df, off_targets_summary = args
    penalized_sirnas = {i: set() for i in sub_df.index}
    for _, row in off_targets_summary.iterrows():
        sirna_names = row["siRNA names"].split(", ")
        off_target_counts = defaultdict(int)
        for sirna in sirna_names:
            sirna_start, sirna_end = parse_sirna_name(sirna)
            for i, row in sub_df.iterrows():
                rna_seq_start, rna_seq_end = [int(pos) for pos in row["Name"].split("_")[-2].split("-")]
                if sirna_start >= rna_seq_start and sirna_end <= rna_seq_end:
                    if sirna not in penalized_sirnas[i]:
                        sub_df.at[i, "Score"] -= 0.1
                        penalized_sirnas[i].add(sirna)
                    off_target_counts[i] += 1
        for i, count in off_target_counts.items():
            if count > 1:
                sub_df.at[i, "Score"] += (count - 1) * -30
    return sub_df

def score_rnai_sequences(rna_sequences, off_targets_summary, base_rnai_length, threads=6):
    rnai_df = pd.DataFrame(rna_sequences, columns=["Name", "Sequence", "Length"])
    rnai_df["Score"] = 0.0
    chunk_size = math.ceil(len(rnai_df) / threads)

    sub_dfs = [rnai_df.iloc[i:i+chunk_size].copy() for i in range(0, len(rnai_df), chunk_size)]
    task_args = [(chunk, off_targets_summary) for chunk in sub_dfs]

    with Pool(processes=threads) as pool:
        results = pool.map(score_chunk, task_args)

    return pd.concat(results)

def design_primers_for_sequence(seq_id, sequence):
    seq_args = {
        'SEQUENCE_ID': seq_id,
        'SEQUENCE_TEMPLATE': sequence,
    }
    global_args = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 17,
        'PRIMER_MAX_SIZE': 23,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_PRODUCT_SIZE_RANGE': [[len(sequence)-20, len(sequence)+20]],
        'PRIMER_NUM_RETURN': 1
    }
    
    primers = primer3.bindings.design_primers(seq_args, global_args)
    left_primer = primers.get('PRIMER_LEFT_0_SEQUENCE', "")
    right_primer = primers.get('PRIMER_RIGHT_0_SEQUENCE', "")
    product_size = primers.get('PRIMER_PAIR_0_PRODUCT_SIZE', "")
    
    if not left_primer or not right_primer:
        left_primer = "no suitable primers found"
        right_primer = "no suitable primers found"
        product_size = ""
        
    return left_primer, right_primer, product_size

def main():
    args = parse_args()
    
    target_seq = read_fasta(args.target_path)
    rna_sequences = generate_rnai_sequences(target_seq, args.rnai_seq_length)
    off_targets_summary = pd.read_csv(args.off_targets_summary_path, sep="\t")
    
    rnai_df = score_rnai_sequences(rna_sequences, off_targets_summary, args.rnai_seq_length, threads=args.threads)
    
    rnai_length_df = rnai_df[rnai_df["Length"] == args.rnai_seq_length].sort_values(by="Score", ascending=False)
    rnai_other_lengths_df = rnai_df[rnai_df["Length"] != args.rnai_seq_length].sort_values(by="Score", ascending=False)
    
    rnai_top_length_df = rnai_length_df.head(25)
    rnai_top_other_lengths_df = rnai_other_lengths_df.groupby("Length").head(10)
    
    final_rnai_df = pd.concat([rnai_top_length_df, rnai_top_other_lengths_df])
    
    left_primers, right_primers, amplicon_sizes = [], [], []
    for idx, row in final_rnai_df.iterrows():
        seq_id = row["Name"]
        seq = row["Sequence"]
        left, right, prod_size = design_primers_for_sequence(seq_id, seq)
        left_primers.append(left)
        right_primers.append(right)
        amplicon_sizes.append(prod_size)
    
    final_rnai_df["Primer_Left"] = left_primers
    final_rnai_df["Primer_Right"] = right_primers
    final_rnai_df["Expected_Amplicon_Size"] = amplicon_sizes
    
    final_rnai_df.to_csv(args.out, sep="\t", index=False)
    print(f"RNAi sequences, scores, primers, and expected amplicon saved to '{args.out}'")

if __name__ == "__main__":
    main()