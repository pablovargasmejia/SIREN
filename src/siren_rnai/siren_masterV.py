#!/usr/bin/env python
import os
import subprocess
import argparse
import sys

def run_sirenXII(targets, gene, threads, sensitivity, outdir, sirna_size, min_align_length, rnahybrid_options=None):    
    sirenXII_path = os.path.join(os.path.dirname(__file__), 'sirenXII.py')
    sirenXII_cmd = [sys.executable, sirenXII_path,
        "--targets", targets,
        "--gene", gene,
        "--threads", str(threads),
        "--sensitivity", sensitivity,
        "--outdir", outdir,
        "--sirna_size", str(sirna_size),
        "--min_align_length", str(min_align_length)
    ]
    if rnahybrid_options:
        if rnahybrid_options and rnahybrid_options[0] == "--":
            rnahybrid_options = rnahybrid_options[1:]
        if rnahybrid_options:
            sirenXII_cmd.extend(["--rnahybrid_options"] + rnahybrid_options)
    subprocess.run(sirenXII_cmd, check=True)
    return os.path.join(outdir, "other_files", "target.fa"), os.path.join(outdir, "off_targets_summary.tsv")

def run_siren_plotIV(fasta, tsv, outdir):
    plot_output = os.path.join(outdir, "Off_targets_across_the_gene.png")
    siren_plotIV_path = os.path.join(os.path.dirname(__file__), 'siren_plotIV.py')
    siren_plot_cmd = [sys.executable, siren_plotIV_path,
        "--fasta", fasta,
        "--input", tsv,
        "--out", plot_output
    ]
    subprocess.run(siren_plot_cmd, check=True)

def run_siren_designVII(fasta, tsv, rnai_length, outdir, threads):
    rnai_tsv_output = os.path.join(outdir, "rna_sequences_with_scores.tsv")
    graph_output = os.path.join(outdir, "rna_sequences_plot.png")
    siren_designVII_path = os.path.join(os.path.dirname(__file__), 'siren_designVIII.py')
    siren_design_cmd = [sys.executable, siren_designVII_path,
        "--target_path", fasta,
        "--off_targets_summary_path", tsv,
        "--rnai_seq_length", str(rnai_length),
        "--threads", str(threads),
        "--out", rnai_tsv_output
    ]
    subprocess.run(siren_design_cmd, check=True)


def run_prefilter(targets, gene, outdir, mode, strand, seed_k, window_size, min_window_hits, write_log):
    prefilter_path = os.path.join(os.path.dirname(__file__), 'siren_prefilter.py')
    filtered_path = os.path.join(outdir if outdir else ".", "targets_prefiltered.fa")
    if not os.path.exists(prefilter_path):
        return targets
    cmd = [sys.executable, prefilter_path,
           "--targets", targets,
           "--gene", gene,
           "--outdir", outdir if outdir else ".",
           "--strand", strand]
    if mode == "windowed":
        cmd += ["--mode", "windowed",
                "--seed_k", str(seed_k),
                "--window_size", str(window_size),
                "--min_window_hits", str(min_window_hits)]
    else:
        cmd += ["--mode", "set"]
    if write_log:
        cmd.append("--write_log")
    try:
        subprocess.run(cmd, check=True)
        return filtered_path if os.path.exists(filtered_path) else targets
    except Exception:
        return targets


def main():
    parser = argparse.ArgumentParser(
    prog="SIREN",
    description=r"""
    
     ____ ___ ____  _____ _   _ 
    / ___|_ _|  _ \| ____| \ | |
    \___ \| || |_) |  _| |  \| |
     ___) | ||  _ <| |___| |\  |
    |____/___|_| \_\_____|_| \_|: Suite for Intelligent RNAi design and Evaluation of Nucleotide sequences.
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Prefilter options
    parser.add_argument("-m", "--prefilter_mode", choices=["set","windowed"], default="windowed",
                    help="Prefilter mode (default: windowed)")
    parser.add_argument("-s", "--prefilter_strand", choices=["rc","fwd","both"], default="rc",
                    help="Strand to seed from (default: rc)")
    parser.add_argument("-k", "--prefilter_seed_k", type=int, default=9,
                    help="Seed k-mer length for windowed mode (default: 9)")
    parser.add_argument("-w", "--prefilter_window_size", type=int, default=40,
                    help="Window size in bp for windowed mode (default: 40)")
    parser.add_argument("-H", "--prefilter_min_window_hits", type=int, default=2,
                    help="Minimum seed hits in any window (default: 2)")
    parser.add_argument("-L", "--prefilter_write_log", dest="prefilter_write_log", action="store_true",
                    help="Write prefilter_log.tsv (default: on)")
    parser.add_argument("-N", "--no_prefilter_write_log", dest="prefilter_write_log", action="store_false",
                    help="Disable prefilter log")
    parser.set_defaults(prefilter_write_log=True)
    parser.add_argument("-X", "--no_prefilter", action="store_true",
                    help="Skip the prefilter step and run the full database")

    # Core workflow
    parser.add_argument("-T", "--targets", required=True, help="FASTA file containing organism cDNA sequences.")
    parser.add_argument("-g", "--gene", required=True, help="Gene name or partial FASTA gene header to identify the target gene.")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads for parallel processing (default: 8).")
    parser.add_argument("-S", "--sensitivity", choices=["high", "medium"], default="medium",
                        help="Sensitivity level for siRNA generation (default: medium).")
    parser.add_argument("-r", "--rnai_length", type=int, default=200, help="Base RNAi sequence length (default: 200).")
    parser.add_argument("-o", "--outdir", default="siren_results", help="Directory to store output files (default: siren_results).")
    parser.add_argument("-z", "--sirna_size", type=int, default=21, help="Length of siRNAs (default: 21).")
    parser.add_argument("-a", "--min_align_length", type=int, help="Minimum alignment length for off-target detection (default: sirna_size - 4).")
    parser.add_argument("-g_o", "--graphical_output", action="store_true",
    help="Run siren_plotIV to produce graphical off-target plots")
    parser.add_argument("-R", "--rnahybrid_options",
        nargs=argparse.REMAINDER,
        default=None,
        help=("Optional pass-through flags for RNAhybrid. Place this flag LAST. "
              "Example: -R -e -25 -v 0 -u 0 -f 2,7 -p 0.01 -d 0.5,0.1 -m 60000")
    )

    args = parser.parse_args()
    if args.min_align_length is None:
        args.min_align_length = args.sirna_size - 4

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    targets_for_run = args.targets
    if not args.no_prefilter:
        targets_for_run = run_prefilter(
            args.targets, args.gene, args.outdir,
            args.prefilter_mode, args.prefilter_strand,
            args.prefilter_seed_k, args.prefilter_window_size,
            args.prefilter_min_window_hits, args.prefilter_write_log
        )

    target_fa, off_targets_summary_tsv = run_sirenXII(
        targets_for_run, args.gene, args.threads, args.sensitivity,
        args.outdir, args.sirna_size, args.min_align_length,
        args.rnahybrid_options
    )
    if args.graphical_output:
        run_siren_plotIV(target_fa, off_targets_summary_tsv, args.outdir)
    run_siren_designVII(target_fa, off_targets_summary_tsv, args.rnai_length, args.outdir, args.threads)

if __name__ == "__main__":
    main()