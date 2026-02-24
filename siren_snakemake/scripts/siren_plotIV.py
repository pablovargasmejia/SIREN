import matplotlib.pyplot as plt
from collections import defaultdict
import argparse
import os

def parse_fasta(fasta_path):
    """Parse the target.fa file to get the length of the gene."""
    with open(fasta_path, 'r') as f:
        lines = f.readlines()
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence, len(sequence)

def parse_off_targets_summary(tsv_path):
    """Parse the off_targets_summary.tsv file and extract unique siRNA positions and off-target positions."""
    siRNA_positions = set()
    off_target_counts = defaultdict(int)
    
    with open(tsv_path, 'r') as f:
        # Skip the header
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            siRNA_names = parts[2].split(', ')
            for sirna in siRNA_names:
                # Extract the start and end positions
                sirna_name = sirna.split('_')[1].split('-')
                start = int(sirna_name[0])
                end = int(sirna_name[1])
                # Add each unique siRNA position to the set
                siRNA_positions.add((start, end))
                # Count the off-targets at each position
                for pos in range(start, end + 1):
                    off_target_counts[pos] += 1
    return siRNA_positions, off_target_counts

def calculate_sirna_counts(siRNA_positions, gene_length):
    """Calculate the siRNA count for each position in the gene."""
    siRNA_counts = defaultdict(int)
    
    for start, end in siRNA_positions:
        for pos in range(start, end + 1):
            siRNA_counts[pos] += 1
            
    return siRNA_counts

def plot_sirna_distribution(siRNA_counts, off_target_counts, gene_length, output_name):
    """Plot the distribution of siRNAs and off-targets across the gene."""
    x = list(range(1, gene_length + 1))
    y_sirnas = [siRNA_counts.get(pos, 0) for pos in x]
    y_off_targets = [off_target_counts.get(pos, 0) for pos in x]

    # Plot siRNAs with off-targets (red line)
    plt.figure(figsize=(10, 6))
    plt.plot(x, y_sirnas, label='siRNAs with off-targets', color='red')

    # Plot off-target counts (blue line)
    plt.plot(x, y_off_targets, label='Off-targets per position', color='blue')

    plt.xlabel('Nucleotide position in the target gene')
    plt.ylabel('Count')
    plt.title('siRNAs and Off-targets across the gene')
    plt.grid(True)
    plt.legend()
    plt.savefig(output_name)

def main():
    parser = argparse.ArgumentParser(description="Parse target.fa and off_targets_summary.tsv and plot siRNA off-target distribution.")
    parser.add_argument('--fasta', required=True, help="Path to the target.fa file.")
    parser.add_argument('--input', required=True, help="Path to the off_targets_summary.tsv file.")
    parser.add_argument('--out', required=True, help="Output plot file name (e.g., Off_targets_across_the_gene.png).")

    args = parser.parse_args()

    # Parse files
    sequence, gene_length = parse_fasta(args.fasta)
    siRNA_positions, off_target_counts = parse_off_targets_summary(args.input)

    # Calculate siRNA counts per position
    siRNA_counts = calculate_sirna_counts(siRNA_positions, gene_length)

    # Save the plot to the specified output file
    plot_sirna_distribution(siRNA_counts, off_target_counts, gene_length, args.out)

if __name__ == "__main__":
    main()

