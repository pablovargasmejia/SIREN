# SIREN: Suite for Intelligent RNAi Design and Evaluation of Nucleotide Sequences

SIREN is a comprehensive toolset for designing RNA interference (RNAi) sequences to silence specific genes while minimizing off-target effects. It integrates siRNA generation, off-target evaluation, off-target visualization, and RNAi sequence plus primer design into a streamlined workflow.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
  - [Conda / Mamba](#conda--mamba)
  - [Pip (PyPI)](#pip-pypi)
  - [Docker](#docker)
- [Requirements](#requirements)
- [Usage](#usage)
  - [Example](#example)
  - [Options](#options)
- [Pipeline Overview](#pipeline-overview)
  - [Prefilter alignment-free k-mer screen](#prefilter-alignment-free-k-mer-screen)
  - [siRNA Generation and Off-target Evaluation](#sirna-generation-and-off-target-evaluation)
  - [Off-target Visualization](#off-target-visualization)
  - [RNAi Selection and Primer Design](#rnai-selection-and-primer-design)
- [Snakemake workflow (optional; cluster-friendly RNAhybrid parallelization)](#snakemake-workflow-optional-cluster-friendly-rnahybrid-parallelization)
- [License](#license)
- [Citations](#citations)

## Features

- **siRNA Generation:** Automatically extracts the target gene from a multi-FASTA file and generates all possible siRNAs.
- **Off-target Evaluation:** Uses RNAhybrid to assess potential off-target interactions.
- **Off-target Visualization:** Creates a plot showing the distribution of siRNAs and off-target events along the gene.
- **RNAi Sequence and Primer Design:** Generates RNAi sequences of various lengths, scores them based on off-target penalties, and designs primers with Primer3 while reporting expected amplicon sizes.

## Installation

### Conda / Mamba

SIREN calls **RNAhybrid** as an external executable. The most reliable way to install is via Bioconda:

```bash
mamba install bioconda::siren-rnai
# or
conda install -c conda-forge -c bioconda siren-rnai
```


### Pip (PyPI)

```bash
mamba install -c conda-forge -c bioconda rnahybrid
pip install siren-rnai
```

#### Apple Silicon Installation

If you're on a Mac with Apple Silicon, this is a safe and compatible setup:

```bash
# 1) Create and activate a new environment with Python 3.12
mamba create -n osx64_env python=3.12.9 -y
mamba activate osx64_env

# 2) Install RNAhybrid from bioconda
mamba install -c conda-forge -c bioconda rnahybrid -y

# 3) Install SIREN
pip install siren-rnai
```

### Docker

A Dockerfile is provided to run SIREN reproducibly without installing dependencies on the host.

**Build:**

```bash
docker build -t siren-rnai:latest .
```

**Run:**

```bash
docker run --rm -it \
  -v "$PWD":/work \
  -w /work \
  siren-rnai:latest \
  SIREN --targets your_db.fa --gene YOUR_GENE --outdir siren_results --threads 12
```

Notes:
- For large databases, mount fast storage for best performance.

## Requirements

- **Python 3.x**
- **Mamba/Conda:** Recommended for installing RNAhybrid from Bioconda.
- **RNAhybrid:** Evaluates off-target interactions (external executable).
- **Primer3:** Required for primer design.
- **BioPython:** For sequence processing.
- **Matplotlib:** For generating visualizations.
- **Additional Python libraries:** Pandas, argparse, csv, tqdm, etc.

## Usage

Run SIREN using the command-line interface:

```bash
SIREN --targets <FASTA file> --gene <gene_name> [--threads <number>] [--sensitivity {high,medium}] [--rnai_length <length>] [--outdir <output_directory>] [--min_align_length <length>]
```

### Example

```bash
SIREN --targets TAIR10_cdna.fasta --gene AT1G50920 --threads 12 --rnai_length 300 --outdir results_AT1G50920
```

This runs the complete SIREN pipeline for the gene `AT1G50920` from the provided FASTA file, using 12 threads and a base RNAi length of 300 nucleotides, storing results in `results_AT1G50920`.

### Options

**Required:**
- `--targets <FASTA>`: FASTA containing organism cDNA sequences.
- `--gene <STRING>`: Gene name or a substring of the FASTA header to select the target.

**Common options:**
- `--threads <INT>`: Parallelism for heavy steps (default: 8).
- `--sensitivity {high,medium}`: Pipeline sensitivity (default: `high`).
- `--rnai_length <INT>`: RNAi region length used downstream (default: 200).
- `--sirna_size <INT>`: siRNA length (default: 21).
- `--min_align_length <INT>`: Minimum alignment length filter for off-target detection (optional).
- `--outdir <DIR>`: Output directory (default: `siren_results`).

**Prefilter controls (alignment-free k-mer screen; optional):**
- `-X, --no_prefilter`: Skip the prefilter step and run on the full database.
- `-m, --prefilter_mode {set,windowed}`: Prefilter mode (default: `windowed`).
- `-s, --prefilter_strand {rc,fwd,both}`: Strand used for seeding (default: `rc`).
- `-k, --prefilter_seed_k <INT>`: Seed k-mer length for windowed mode (default: 9).
- `-w, --prefilter_window_size <INT>`: Window size for density criterion (default: 40).
- `-H, --prefilter_min_window_hits <INT>`: Minimum seed hits in a window (default: 2).
- `-L, --prefilter_write_log` / `-N, --no_prefilter_write_log`: Toggle writing a TSV log (`prefilter_log.tsv`).

**Visualization:**
- `-g_o, --graphical_output`: Also run `siren_plotIV.py` to produce the off-target plot.

**Pass-through to RNAhybrid (placed last):**
- `-R, --rnahybrid_options ...`: Any extra flags forwarded directly to `RNAhybrid`.
  - Example:
    ```bash
    -R -e -25 -v 0 -u 0 -f 2,7 -p 0.01 -d 0.5,0.1 -m 60000
    ```


## Pipeline Overview

### Prefilter alignment-free k-mer screen

The `siren_prefilter.py` module:

- **Purpose:** Rapidly shrink the RNAhybrid search space by keeping only sequences likely to be similar to the target **without full alignments**.
- **Two modes (from the code):**
  - **`set` mode:** Compares distinct k-mer sets between each candidate sequence and the target using similarity measures.
  - **`windowed` mode (pipeline default):** Retains a sequence if there are ≥ H reverse-complement seed hits (exact k-mers) within any W-bp window. Pipeline defaults: **k=9**, **W=40**, **H=2**.
- **Strand control:** Seeds can be taken from **`rc`**, **`fwd`**, or **`both`** (pipeline default: **`rc`**).
- **Logging (optional):** When enabled, writes a per-record metrics table to **`prefilter_log.tsv`**.
- **Output:** Writes the filtered database to **`targets_prefiltered.fa`** (in the specified output directory).
- **Integration:** Enabled by default; can be skipped with **`--no_prefilter`**.

### siRNA Generation and Off-target Evaluation

The `sirenXII.py` module:
- **Target Extraction:** Searches the provided FASTA for the specified gene and extracts a unique target sequence.
- **siRNA Generation:** Generates siRNAs (typically 21 nucleotides) using a sensitivity-dependent step size.
- **Off-target Evaluation:** Evaluates off-target interactions via RNAhybrid on sequences not matching the target.
- **Parallel Processing:** Splits off-target data into chunks for parallel processing.
- **Output Files:** Produces files such as `target.fa` and `off_targets_summary.tsv` for downstream steps.

### Off-target Visualization

The `siren_plotIV.py` module:
- **Data Parsing:** Reads the target FASTA and the off-target summary TSV.
- **Aggregation:** Computes the distribution of siRNAs and off-target events along the gene.
- **Plot Generation:** Uses Matplotlib to create a plot with:
  - A red line for the count of siRNAs with off-target events.
  - A blue line for the count of off-target events per nucleotide position.
- **Output:** Saves the plot (e.g., `Off_targets_across_the_gene.png`).

### RNAi Selection and Primer Design

The `siren_designVIII.py` module:
- **RNAi Sequence Generation:** Creates RNAi sequences with lengths from (base length - 50) to (base length + 100) in steps of 50.
- **Scoring:** RNAi sequences are penalized for containing siRNAs with off-target potential. Each unique siRNA contributing to off-targets reduces the score slightly. If multiple off-targets are caused by the same siRNA within a given RNAi sequence, a stronger penalty is applied.
- **Primer Design:** Utilizes Primer3 to design primer pairs and calculates expected amplicon sizes.
- **Output:** Generates a TSV file (`rna_sequences_with_scores_and_primers.tsv`) with RNAi sequences, scores, primer details, and expected amplicon sizes.


## Snakemake workflow (optional; cluster-friendly RNAhybrid parallelization)

A Snakemake workflow is available in `siren_snakemake/` to parallelize the **RNAhybrid** stage by splitting the target database into shards and running shards concurrently (local or SLURM).

```bash
# 1) Create a small driver env for Snakemake (once)
conda create -n smk -c conda-forge -y snakemake "conda>=24.7.1"
conda activate smk

# 2) Run the workflow
cd siren_snakemake
# edit config.yaml (targets, gene, outdir, num_shards, rnahybrid_options)
snakemake --cores 32 --use-conda -p

# If a previous run was interrupted (stale lock)
snakemake --unlock
```



## License

SIREN is released under the GPLv3 license.

## Citations

If you use **SIREN** in your research, please cite:

- **SIREN** – Vargas Mejía, P., & Vega-Arreguín, J. C. (2025). SIREN: Suite for Intelligent RNAi Design and Evaluation of Nucleotide Sequences. bioRxiv. https://doi.org/10.1101/2025.05.26.656188.
- **RNAhybrid** – Rehmsmeier, M., Steffen, P., Höchsmann, M., & Giegerich, R. (2004). Fast and effective prediction of microRNA/target duplexes. RNA, 10(10), 1507–1517.

For any issues, feature requests, or further questions, please open an issue on GitHub. Happy RNAi designing!

