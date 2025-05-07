#!/usr/bin/env python
import os
import subprocess
import argparse

def run_sirenXII(targets, gene, threads, sensitivity, outdir, sirna_size, min_align_length):
    sirenXII_path = os.path.join(os.path.dirname(__file__), 'sirenXII.py')
    sirenXII_cmd = [
        "python", sirenXII_path,
        "--targets", targets,
        "--gene", gene,
        "--threads", str(threads),
        "--sensitivity", sensitivity,
        "--outdir", outdir,
        "--sirna_size", str(sirna_size),
        "--min_align_length", str(min_align_length)
    ]
    subprocess.run(sirenXII_cmd, check=True)
    return os.path.join(outdir, "other_files", "target.fa"), os.path.join(outdir, "off_targets_summary.tsv")

def run_siren_plotIV(fasta, tsv, outdir):
    plot_output = os.path.join(outdir, "Off_targets_across_the_gene.png")
    siren_plotIV_path = os.path.join(os.path.dirname(__file__), 'siren_plotIV.py')
    siren_plot_cmd = [
        "python", siren_plotIV_path,
        "--fasta", fasta,
        "--input", tsv,
        "--out", plot_output
    ]
    subprocess.run(siren_plot_cmd, check=True)

def run_siren_designVII(fasta, tsv, rnai_length, outdir, threads):
    rnai_tsv_output = os.path.join(outdir, "rna_sequences_with_scores.tsv")
    graph_output = os.path.join(outdir, "rna_sequences_plot.png")
    siren_designVII_path = os.path.join(os.path.dirname(__file__), 'siren_designVIII.py')
    siren_design_cmd = [
        "python", siren_designVII_path,
        "--target_path", fasta,
        "--off_targets_summary_path", tsv,
        "--rnai_seq_length", str(rnai_length),
        "--threads", str(threads),
        "--out", rnai_tsv_output
    ]
    subprocess.run(siren_design_cmd, check=True)

def main():
    parser = argparse.ArgumentParser(
        prog="SIREN",
        description=r"""

 ____ ___ ____  _____ _   _ 
/ ___|_ _|  _ \| ____| \ | |
\___ \| || |_) |  _| |  \| |
 ___) | ||  _ <| |___| |\  |
|____/___|_| \_\_____|_| \_|: Suite for Intelligent RNAi design and Evaluation of Nucleotide sequences.

SIREN is a comprehensive toolset for designing RNA interference (RNAi) sequences to silence specific genes while minimizing off-target effects.

The workflow consists of three main steps:
  1. siRNA Generation: Generates all possible siRNAs from a target FASTA file and evaluates potential off-target interactions with RNAhybrid.
  2. Off-target Visualization: Creates a plot showing the distribution of siRNA hits and off-target events along the gene.
  3. RNAi Design & Primer Design: Scores RNAi sequences based on off-target penalties, designs primers using Primer3, and reports the expected amplicon size.

Usage:
  SIREN --targets <FASTA file> --gene <gene_name> [--threads <number>] [--sensitivity {high,medium,low}] [--rnai_length <length>] [--outdir <output_directory>]

Example for Arabidopsis:
  SIREN --targets TAIR10_cdna.fasta --gene AT1G50920 --threads 12 --rnai_length 300 --outdir results_AT1G50920
""",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--targets", required=True, help="FASTA file containing organism cDNA sequences.")
    parser.add_argument("--gene", required=True, help="Gene name or partial FASTA gene header to identify the target gene.")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for parallel processing (default: 8).")
    parser.add_argument("--sensitivity", choices=["high", "medium", "low"], default="low",
                        help="Sensitivity level for siRNA generation (default: low).")
    parser.add_argument("--rnai_length", type=int, default=200, help="Base RNAi sequence length (default: 200).")
    parser.add_argument("--outdir", default="siren_results", help="Directory to store output files (default: siren_results).")
    parser.add_argument("--sirna_size", type=int, default=21, help="Length of siRNAs (default: 21).")

    parser.add_argument("--min_align_length", type=int, help="Minimum alignment length for off-target detection (default: sirna_size - 4).")
    args = parser.parse_args()
    if args.min_align_length is None:
        args.min_align_length = args.sirna_size - 4


    # Create output directory if it doesn't exist.
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Run the three steps in sequence.
    target_fa, off_targets_summary_tsv = run_sirenXII(args.targets, args.gene, args.threads, args.sensitivity, args.outdir, args.sirna_size, args.min_align_length)
    run_siren_plotIV(target_fa, off_targets_summary_tsv, args.outdir)
    run_siren_designVII(target_fa, off_targets_summary_tsv, args.rnai_length, args.outdir, args.threads)

if __name__ == "__main__":
    main()

