#!/usr/bin/env python3

### Simple python script written with assistance from Claude to parse a gff annotation file and output statistics relating to gene density
### Simplifying the process is that any overlapping alignments are simply removed, but this number is tracked. It's around 5% of genes in the handful of test cases I've checked

import sys
import argparse
from collections import defaultdict
import os
import statistics

def parse_gff(gff_file):
    """Parse GFF file and extract gene coordinates by chromosome."""
    genes_by_chr = defaultdict(list)
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2].lower() != 'gene':
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            genes_by_chr[chrom].append((start, end, strand))
    return genes_by_chr

def calculate_intergenic_distances(genes_by_chr):
    """Calculate intergenic distances for each chromosome."""
    distances = []
    overlap_count = 0
    for chrom, genes in genes_by_chr.items():
        # Sort genes by start position
        genes.sort()
        # Calculate distances between adjacent genes
        for i in range(len(genes) - 1):
            current_gene_end = genes[i][1]
            next_gene_start = genes[i+1][0]
            # Check for overlap
            if current_gene_end >= next_gene_start:
                overlap_count += 1
                continue
            # Calculate distance
            distance = next_gene_start - current_gene_end - 1
            distances.append(distance)
    return distances, overlap_count

def get_summary_stats(distances, overlap_count, gff_file):
    """Calculate summary statistics and return as a dictionary."""
    total_pairs = len(distances) + overlap_count
    
    stats = {
        'gff_file': os.path.basename(gff_file),
        'total_gene_pairs': total_pairs,
        'overlapping_genes': overlap_count,
        'valid_distances': len(distances),
        'mean_distance': statistics.mean(distances) if distances else 0,
        'median_distance': statistics.median(distances) if distances else 0,
        'min_distance': min(distances) if distances else 0,
        'max_distance': max(distances) if distances else 0,
        'std_dev': statistics.stdev(distances) if len(distances) > 1 else 0
    }
    return stats

def write_summary_header(output_file):
    """Write header line to the summary file."""
    headers = [
        'gff_file',
        'total_gene_pairs',
        'overlapping_genes',
        'valid_distances',
        'mean_distance',
        'median_distance',
        'min_distance',
        'max_distance',
        'std_dev'
    ]
    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')

def append_summary_stats(stats, output_file):
    """Append statistics to the summary file."""
    values = [
        str(stats['gff_file']),
        str(stats['total_gene_pairs']),
        str(stats['overlapping_genes']),
        str(stats['valid_distances']),
        f"{stats['mean_distance']:.2f}",
        f"{stats['median_distance']:.2f}",
        str(stats['min_distance']),
        str(stats['max_distance']),
        f"{stats['std_dev']:.2f}"
    ]
    with open(output_file, 'a') as f:
        f.write('\t'.join(values) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Calculate intergenic distances from GFF file')
    parser.add_argument('gff_file', help='Input GFF file')
    parser.add_argument('--summary', '-s', help='Summary output file (tab-delimited)', required=True)
    args = parser.parse_args()

    # Create summary file with headers if it doesn't exist
    if not os.path.exists(args.summary):
        write_summary_header(args.summary)

    # Parse GFF and calculate distances
    genes_by_chr = parse_gff(args.gff_file)
    distances, overlap_count = calculate_intergenic_distances(genes_by_chr)

    # Calculate and append summary statistics
    stats = get_summary_stats(distances, overlap_count, args.gff_file)
    append_summary_stats(stats, args.summary)

if __name__ == '__main__':
    main()
