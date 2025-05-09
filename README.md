
# SIREN: Suite for Intelligent RNAi Design and Evaluation of Nucleotide Sequences

SIREN is a comprehensive toolset for designing RNA interference (RNAi) sequences to silence specific genes while minimizing off‑target effects. It integrates siRNA generation, off‑target evaluation, off‑target visualization, and RNAi sequence plus primer design into a streamlined workflow.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Requirements](#requirements)
- [Usage](#usage)
- [Pipeline Overview](#pipeline-overview)
  - [siRNA Generation and Off‑target Evaluation](#sirna-generation-and-offtarget-evaluation)
  - [Off‑target Visualization](#offtarget-visualization)
  - [RNAi Selection and Primer Design](#rnai-selection-and-primer-design)
- [License](#license)
- [Citations](#citations)

## Features

- **siRNA Generation:** Automatically extracts the target gene from a multi‑FASTA file and generates all possible siRNAs.
- **Off‑target Evaluation:** Uses RNAhybrid to assess potential off‑target interactions.
- **Off‑target Visualization:** Creates a plot showing the distribution of siRNAs and off‑target events along the gene.
- **RNAi Sequence and Primer Design:** Generates RNAi sequences of various lengths, scores them based on off‑target penalties, and designs primers with Primer3 while reporting expected amplicon sizes.

## Installation

First install RNAhybrid which is available on Bioconda. Install it using mamba:

```bash
    mamba install bioconda::rnahybrid
    # or
    conda install bioconda::rnahybrid
```
SIREN is available on PyPi. Install it using pip:

```bash
    pip install siren-rnai
```

This command installs SIREN along with all required dependencies.

### Apple Silicon Installation

If you're on a Mac with Apple Silicon, follow these steps to install SIREN in a clean and compatible environment:

```bash
# 1. Create and activate a new environment with Python 3.12
mamba create -n osx64_env python=3.12.9 -y
mamba activate osx64_env

# 2. Install RNAhybrid from bioconda
mamba install -c bioconda rnahybrid -y

# 3. Install SIREN
pip install siren-rnai
```

You can now use `SIREN` from the command line inside this environment.

## Requirements

- **Python 3.x**
- **pip**: For installing Python packages when needed.
- **Mamba/Conda:** For installation from Bioconda.
- **RNAhybrid:** Evaluates off‑target interactions.
- **Primer3:** Required for primer design.
- **BioPython:** For sequence processing.
- **Matplotlib:** For generating visualizations.
- **Additional Python libraries:** Pandas, argparse, csv, tqdm, etc.

## Usage

Run SIREN using the command-line interface as shown below. The options are as follows:

```bash
SIREN --targets <FASTA file> --gene <gene_name> [--threads <number>] [--sensitivity {high,medium,low}] [--rnai_length <length>] [--outdir <output_directory>] [--min_align_length <length>]
```

### Options:

- **`--targets <FASTA file>`**: Path to the FASTA file containing organism cDNA sequences.
- **`--gene <gene_name>`**: Gene name (or partial FASTA header) used to uniquely identify the target gene.
- **`--threads <number>`**: Number of threads for parallel processing (default: 6).
- **`--sensitivity {high,medium,low}`**: Sensitivity level for siRNA generation; controls the step size for candidate siRNA generation (default: low).
- **`--rnai_length <length>`**: Base RNAi sequence length to guide sequence generation (default: 200).
- **`--outdir <output_directory>`**: Directory where all output files will be stored (default: `siren_results`).
- **`--min_align_length <length>`**: (Optional) Minimum alignment length for off‑target detection.

### Example:

```bash
SIREN --targets TAIR10_cdna.fasta --gene AT1G50920 --threads 12 --rnai_length 300 --outdir results_AT1G50920
```

This command runs the complete SIREN pipeline for the gene `AT1G50920` from the provided Arabidopsis cds FASTA file, using 12 threads and a base RNAi length of 300 nucleotides, storing results in the `results_AT1G50920` directory.

## Pipeline Overview

### siRNA Generation and Off‑target Evaluation

The `sirenXII.py` module:
- **Target Extraction:** Searches the provided FASTA for the specified gene and extracts a unique target sequence.
- **siRNA Generation:** Generates siRNAs (typically 21 nucleotides) using a sensitivity-dependent step size.
- **Off‑target Evaluation:** Evaluates off‑target interactions via RNAhybrid on sequences not matching the target.
- **Parallel Processing:** Splits off‑target data into chunks for parallel processing.
- **Output Files:** Produces files such as `target.fa` and `off_targets_summary.tsv` for downstream steps.

### Off‑target Visualization

The `siren_plotIV.py` module:
- **Data Parsing:** Reads the target FASTA and the off‑target summary TSV.
- **Aggregation:** Computes the distribution of siRNAs and off‑target events along the gene.
- **Plot Generation:** Uses Matplotlib to create a plot with:
  - A red line for the count of siRNAs with off‑target events.
  - A blue line for the count of off‑target events per nucleotide position.
- **Output:** Saves the plot (e.g., `Off_targets_across_the_gene.png`).

### RNAi Selection and Primer Design

The `siren_designVIII.py` module:
- **RNAi Sequence Generation:** Creates RNAi sequences with lengths from (base length - 50) to (base length + 100) in steps of 50.
- **Scoring:** RNAi sequences are penalized for containing siRNAs with off-target potential. Each unique siRNA contributing to off-targets reduces the score slightly. However, if multiple off-targets are caused by the same siRNA within a given RNAi sequence, a strong penalty is applied. This discourages designs that repeatedly include problematic siRNAs, improving overall targeting specificity.
- **Primer Design:** Utilizes Primer3 to design primer pairs and calculates expected amplicon sizes.
- **Output:** Generates a TSV file (`rna_sequences_with_scores_and_primers.tsv`) with RNAi sequences, scores, primer details, and expected amplicon sizes.

## License

SIREN is released under the [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html) license.

## Citations

If you use **SIREN** in your research or projects, please cite the following tools and resources:

- **SIREN** – please cite this GitHub repository.
- **RNAhybrid** – Rehmsmeier, M., Steffen, P., Höchsmann, M., & Giegerich, R. (2004). Fast and effective prediction of microRNA/target duplexes. *RNA*, 10(10), 1507–1517.


For any issues, feature requests, or further questions, please open an issue on GitHub. Happy RNAi designing!
