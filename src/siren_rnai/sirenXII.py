#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import subprocess
from multiprocessing import Pool
import csv
from tqdm import tqdm
import shutil

script_has_run = False

def check_rnahybrid():
    try:
        subprocess.run(["RNAhybrid", "-h"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Error: RNAhybrid is not installed or not found in PATH.")
        print("Please install it using: ")
        print("    mamba install -c bioconda rnahybrid")
        print("#or")
        print("    conda install -c bioconda rnahybrid")
        return False
    return True

def generate_sirnas(gene, targets, sirna_size=21, sensitivity="low", outdir="siren_results"):
    other_files_dir = os.path.join(outdir, "other_files")
    if not os.path.exists(other_files_dir):
        os.makedirs(other_files_dir)

    target_seq = None
    off_target_records = []
    matching_records = []

    for record in SeqIO.parse(targets, "fasta"):
        if gene in record.id:
            matching_records.append(record)

    if len(matching_records) == 0:
        print(f"Error: Gene '{gene}' not found in the targets database. Check the name and try again with an existent gene. Quote the gene name if has special characters")
        return None, None, None
    elif len(matching_records) > 1:
        print(f"Error: More than one match found for gene '{gene}'. Please provide a more specific gene name.")
        for match in matching_records:
            print(f" - {match.id}")
        return None, None, None
    else:
        target_seq = matching_records[0]

    for record in SeqIO.parse(targets, "fasta"):
        if record.id != target_seq.id:
            off_target_records.append(record)

    target_file = os.path.join(other_files_dir, "target.fa")
    SeqIO.write(target_seq, target_file, "fasta")

    off_targets_file = os.path.join(other_files_dir, "sequences.fa")
    SeqIO.write(off_target_records, off_targets_file, "fasta")

    sirnas = []
    seq_len = len(target_seq.seq)

    step = 1 if sensitivity == "high" else 4 if sensitivity == "medium" else 8

    for i in range(0, seq_len - sirna_size + 1, step):
        sirna = target_seq.seq[i:i + sirna_size]
        sirna_r = sirna.reverse_complement().transcribe()
        start_pos = i + 1
        end_pos = start_pos + sirna_size - 1
        sirna_name = f"sirna_{start_pos}-{end_pos}"
        sirna_r_name = f"{sirna_name}_r"
        sirnas.append(SeqRecord(sirna, id=sirna_name, description=""))
        sirnas.append(SeqRecord(sirna_r, id=sirna_r_name, description=""))

    sirnas_file = os.path.join(outdir, "sirnas.fa")
    SeqIO.write(sirnas, sirnas_file, "fasta")

    print("")
    print("SIREN: Suite for Intelligent RNAi design and Evaluation of Nucleotide sequences")
    print("")
    print("Generating possible siRNAs from gene")
    return off_targets_file, sirnas_file, target_file

def run_rnahybrid(off_target_file, sirna_file, options, outdir, thread_id, results_file=None,sirna_size=21, min_align_length=None):
    output_file = os.path.join(outdir, f"rnahybrid_output_{thread_id}.txt")
    cmd = ["RNAhybrid", "-t", off_target_file, "-q", sirna_file] + options
    with open(output_file, "w") as f_out:
        subprocess.run(cmd, stdout=f_out, stderr=subprocess.STDOUT)
    if results_file:
        parse_rnahybrid_results(output_file, sirna_size=sirna_size, out_file=results_file, min_align_length=min_align_length)
        os.rename(output_file, results_file)
        return results_file
    return output_file

def split_fasta(input_file, num_splits):
    records = list(SeqIO.parse(input_file, "fasta"))
    split_size = len(records) // num_splits
    split_files = []
    for i in range(num_splits):
        split_file = os.path.join(os.path.dirname(input_file), f"sequences_split_{i + 1}.fa")
        SeqIO.write(records[i*split_size:(i+1)*split_size], split_file, "fasta")
        split_files.append(split_file)
    return split_files

def parse_rnahybrid_results(input_file, sirna_size, out_file, min_align_length=None):
    if not min_align_length:
        min_align_length = sirna_size - 4
    with open(input_file, "r") as infile, open(out_file, "w") as outfile:
        block = []
        inside_block = False
        for line in infile:
            if line.startswith("target:"):
                if block:
                    aligned_sirna_sequence = block[11].replace(" ", "").strip()
                    if len(aligned_sirna_sequence) >= min_align_length:
                        outfile.writelines(block)
                    block = []
                inside_block = True
            if inside_block:
                block.append(line)
        if block:
            aligned_sirna_sequence = block[11].replace(" ", "").strip()
            if len(aligned_sirna_sequence) >= min_align_length:
                outfile.writelines(block)

def generate_off_target_tsv(input_file, tsv_output):
    off_target_data = {}
    target = None
    sirna = None
    with open(input_file, "r") as infile:
        for line in infile:
            line = line.strip()
            if line.startswith("target:"):
                target = line.split(": ")[1].strip()
                if target not in off_target_data:
                    off_target_data[target] = {"count": 0, "sirnas": []}
            elif line.startswith("miRNA :"):
                sirna = line.split(": ")[1].strip()
                if target and sirna:
                    off_target_data[target]["count"] += 1
                    off_target_data[target]["sirnas"].append(sirna)
            elif line == "":
                target = None
                sirna = None
    sorted_data = sorted(off_target_data.items(), key=lambda x: x[1]["count"], reverse=True)
    with open(tsv_output, "w", newline='') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(["Off target", "siRNA number", "siRNA names"])
        for target, data in sorted_data:
            writer.writerow([target, data["count"], ", ".join(data["sirnas"])])
    print(f"Total potential off targets found: {len(off_target_data)}")

def process_rnahybrid_task(args):
    f, sirnas_file, options, outdir, off_target_splits = args
    thread_id = off_target_splits.index(f)
    return run_rnahybrid(f, sirnas_file, options, outdir, thread_id)

def clean_up_temp_files(directory):
    files_to_remove = [f for f in os.listdir(directory) if f.startswith("sequences_split_")]
    for file in files_to_remove:
        os.remove(os.path.join(directory, file))

if __name__ == "__main__":
    if not script_has_run:
        script_has_run = True
        parser = argparse.ArgumentParser(description="Generate siRNAs and evaluate potential off-targets")
        parser.add_argument("--targets", required=True)
        parser.add_argument("--gene", required=True)
        parser.add_argument("--sirna_size", type=int, default=21)
        parser.add_argument("--threads", type=int, default=6)
        parser.add_argument("--outdir", default="siren_results")
        parser.add_argument("--sensitivity", choices=["high", "medium", "low"], default="low")
        parser.add_argument("--min_align_length", type=int)
        parser.add_argument("--rnahybrid_options", nargs='+', default=["-e", "-25", "-v", "0", "-u", "0", "-f", "2,7", "-p", "0.01", "-d", "0.5,0.1", "-m", "60000"])
        args = parser.parse_args()

        if check_rnahybrid():
            off_targets_file, sirnas_file, target_file = generate_sirnas(args.gene, args.targets, args.sirna_size, args.sensitivity, args.outdir)
            if off_targets_file and sirnas_file:
                num_splits = args.threads * 8
                off_target_splits = split_fasta(off_targets_file, num_splits)
                print("Finding off targets... (this step may take a while, try to add more threads with --threads option or reduce --sensitivity)")
                combined_output = os.path.join(args.outdir, "other_files", "all_targets.txt")
                with Pool(args.threads) as pool:
                    results = list(tqdm(pool.imap_unordered(process_rnahybrid_task, [(f, sirnas_file, args.rnahybrid_options, args.outdir, off_target_splits) for f in off_target_splits]), total=len(off_target_splits)))
                with open(combined_output, "w") as outfile:
                    for result in results:
                        if result:
                            with open(result, "r") as infile:
                                outfile.write(infile.read())
                            os.remove(result)
                clean_up_temp_files(os.path.join(args.outdir, "other_files"))
                off_targets_results_file = os.path.join(args.outdir, "off_targets_results.txt")
                min_align_length = args.min_align_length if args.min_align_length else args.sirna_size - 4
                parse_rnahybrid_results(combined_output, args.sirna_size, off_targets_results_file, min_align_length)
                tsv_output = os.path.join(args.outdir, "off_targets_summary.tsv")
                generate_off_target_tsv(off_targets_results_file, tsv_output)
                print("Finding target silencing efficiency...")
                target_results_file = os.path.join(args.outdir, "target_results.txt")
                temp_target_file = os.path.join(args.outdir, "temp_target_results.txt")
                result = run_rnahybrid(target_file, sirnas_file, args.rnahybrid_options, args.outdir, 0, temp_target_file, sirna_size=args.sirna_size, min_align_length=min_align_length)
                if not os.path.exists(temp_target_file):
                    raise FileNotFoundError(f"The file {temp_target_file} was not created.")
                os.rename(temp_target_file, target_results_file)
                other_files_dir = os.path.join(args.outdir, "other_files")
                if not os.path.exists(other_files_dir):
                    os.makedirs(other_files_dir)
                files_to_move = [
                    os.path.join(args.outdir, "sirnas.fa"),
                    os.path.join(args.outdir, "target_results.txt"),
                    os.path.join(args.outdir, "off_targets_results.txt")
                ]
                for file_path in files_to_move:
                    if os.path.exists(file_path):
                        dest_file = os.path.join(other_files_dir, os.path.basename(file_path))
                        os.replace(file_path, dest_file)

