#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import subprocess
from multiprocessing import Pool
import csv
import tempfile
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil
import time

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

def generate_sirnas(gene, targets, sirna_size=21, sensitivity="medium", outdir="siren_results"):
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

    step = 1 if sensitivity == "high" else 2 if sensitivity == "medium" else 2

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

def run_rnahybrid(off_target_file, sirna_file, options, outdir, thread_id,
                  results_file=None, sirna_size=21, min_align_length=None):

    output_file = os.path.join(outdir, f"rnahybrid_output_{thread_id}.txt")

    # normalize options
    opts = options[:] if options else []   # make a copy or empty list
    if opts and opts[0] == "--":
        opts = opts[1:]  # strip a lone '--' if user typed it
    if not opts:  # if still empty, fall back to defaults
        opts = ["-e", "-25", "-v", "0", "-u", "0", "-f", "2,7", "-p", "0.01", "-d", "0.5,0.1", "-m", "60000"]

    # build command
    cmd = ["RNAhybrid", "-t", off_target_file, "-q", sirna_file] + opts

    # run RNAhybrid
    with open(output_file, "w") as f_out:
        subprocess.run(cmd, stdout=f_out, stderr=subprocess.STDOUT)

    # parse results if requested
    # parse results if requested
    if results_file:
        parse_rnahybrid_results(output_file,
                            sirna_size=sirna_size,
                            out_file=results_file,
                            min_align_length=min_align_length)
        return results_file
    else:
        return output_file




def run_rnahybrid_capture(off_target_file, sirna_file, options):
    """Run RNAhybrid and return its stdout as text (no temp files)."""
    # normalize options (mirror run_rnahybrid)
    opts = list(options) if options else []
    if opts and opts[0] == "--":
        opts = opts[1:]
    if not opts:
        # Same defaults as your run_rnahybrid()
        opts = ['-e', '-25', '-v', '0', '-u', '0', '-f', '2,7', '-p', '0.01', '-d', '0.5,0.1', '-m', '60000']

    cmd = ['RNAhybrid', '-t', off_target_file, '-q', sirna_file] + opts
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False
    )
    return proc.stdout or ""

def split_fasta(input_file, num_splits, tmp_dir=None):
    """Split FASTA into ~balanced shards by total sequence length.
    Writes shards into tmp_dir (defaults to /tmp), returns list of file paths.
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    if not records:
        return []
    # choose tmp_dir in RAM (/tmp) if not provided
    if tmp_dir is None:
        tmp_dir = tempfile.mkdtemp(prefix="siren_splits_")
    os.makedirs(tmp_dir, exist_ok=True)

    lengths = [len(r.seq) for r in records]
    order = sorted(range(len(records)), key=lambda i: lengths[i], reverse=True)
    k = max(1, min(num_splits, len(records)))

    buckets = [[] for _ in range(k)]
    bucket_loads = [0] * k

    for i in order:
        j = min(range(k), key=lambda idx: bucket_loads[idx])
        buckets[j].append(records[i])
        bucket_loads[j] += lengths[i]

    split_files = []
    for idx, chunk in enumerate(buckets, 1):
        if not chunk:
            continue
        split_file = os.path.join(tmp_dir, f"sequences_split_{idx}.fa")
        SeqIO.write(chunk, split_file, "fasta")
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
        parser.add_argument("--shard_factor", type=int, default=12,
                           help="Multiply factor for shards vs threads (default: 32).")
        parser.add_argument("--outdir", default="siren_results")
        parser.add_argument("--sensitivity", choices=["high", "medium", "low"], default="low")
        parser.add_argument("--min_align_length", type=int)
        parser.add_argument("--rnahybrid_options",
        nargs=argparse.REMAINDER,
        default=None,
        help=("Optional flags forwarded directly to RNAhybrid. Place this flag LAST. "
          "Example: --rnahybrid_options -s 3utr_human -e -30 -c -f 2,7 -p 0.001 -d 0.5,0.1 -m 60000")
        )
        args = parser.parse_args()

        if check_rnahybrid():
            t0 = time.time()
            off_targets_file, sirnas_file, target_file = generate_sirnas(args.gene, args.targets, args.sirna_size, args.sensitivity, args.outdir)
            print(f"[TIMING] generate_sirnas: {time.time()-t0:.1f}s")
            if off_targets_file and sirnas_file:
                # Compute number of shards: limit by number of records; aim for threads * shard_factor
                total_records = sum(1 for _ in SeqIO.parse(off_targets_file, "fasta"))
                target_factor = max(8, args.shard_factor)
                num_splits = max(min(total_records, args.threads * target_factor), args.threads)
                
                # Write splits into RAM-backed tmp to reduce disk I/O
                tmp_dir = tempfile.mkdtemp(prefix="siren_splits_")
                off_target_splits = split_fasta(off_targets_file, num_splits, tmp_dir=tmp_dir)
                
                print("Finding off targets... (this step may take a while, try to add more threads with --threads)")
                combined_output = os.path.join(args.outdir, "other_files", "all_targets.txt")
                os.makedirs(os.path.dirname(combined_output), exist_ok=True)

                if not (os.path.exists(combined_output) and os.path.getsize(combined_output) > 0):
                    t_rna = time.time()
                    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex, open(combined_output, "w") as out_all:
                        futs = {ex.submit(run_rnahybrid_capture, sp, sirnas_file, args.rnahybrid_options): sp for sp in off_target_splits}
                        for fut in tqdm(as_completed(futs), total=len(futs)):
                            out = fut.result()
                            if out:
                                out_all.write(out)
                    print(f"✔ RNAhybrid stage done in {time.time()-t_rna:.1f}s")
                else:
                    print("Skipping RNAhybrid: found existing all_targets.txt")

                # cleanup split files and tmp dir
                for _p in off_target_splits:
                    try:
                        os.remove(_p)
                    except Exception:
                        pass
                try:
                    os.rmdir(tmp_dir)
                except Exception:
                    pass

                clean_up_temp_files(os.path.join(args.outdir, "other_files"))
                off_targets_results_file = os.path.join(args.outdir, "off_targets_results.txt")
                min_align_length = args.min_align_length if args.min_align_length else args.sirna_size - 4
                t1 = time.time()
                parse_rnahybrid_results(combined_output, args.sirna_size, off_targets_results_file, min_align_length)
                print(f"[TIMING] parse_offtargets: {time.time()-t1:.1f}s")
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