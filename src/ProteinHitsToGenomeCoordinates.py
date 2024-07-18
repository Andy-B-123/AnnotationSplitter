### ProteinCoordinatesToGenome
### Take regions of a protein and translate them back to the genome coordinates
### Things to catch:
###     - Codons on splice sites aren't neatly definable on genome coordinates
###     - Protein hits always positive direction, genome coordiantes can be negative strand
###     - Sometimes I want to round things so that a single amino acid doesn't make the region 'jump' to the next exon 

from tqdm import tqdm
from cogent3 import load_annotations
from cogent3 import load_seq
from cogent3.core.annotation_db import GffAnnotationDb
import pandas as pd
import os
import matplotlib.pyplot as plt

### Nicely colour output for bed file
def bits_to_rgb(bits, min_bits, max_bits):
    # Cap the bits value at the max_bits value
    capped_bits = min(bits, max_bits)
    norm = plt.Normalize(min_bits, max_bits)
    cmap = plt.get_cmap('coolwarm')  # Use a color map from blue to red
    rgba = cmap(norm(capped_bits))
    return f'{int(rgba[0]*255)},{int(rgba[1]*255)},{int(rgba[2]*255)}'

### Get the nucleotide position
# Function to extract nucleotide positions based on qstart and qend
def extract_nucleotide_positions(row):
    nucleotide_to_aa_mapping = row['nucleotide_to_aa_mapping']
     
    # Check strand:
    if row['strand'] == "+":
        # Get nucleotide positions for qstart and qend amino acids
        start_positions = nucleotide_to_aa_mapping.get(row['qstart'] , [(None, None)])
        end_positions = nucleotide_to_aa_mapping.get(row['qend'], [(None, None)])
        # Extract the relevant nucleotide positions and exon numbers
        start_nucleotide, start_exon = start_positions[0]
        end_nucleotide, end_exon = end_positions[-1]
    else:
        # Get nucleotide positions for qstart and qend amino acids
        start_positions = nucleotide_to_aa_mapping.get(row['qstart'] , [(None, None)])
        end_positions = nucleotide_to_aa_mapping.get(row['qend'], [(None, None)])        
        # Extract the relevant nucleotide positions and exon numbers
        start_nucleotide, start_exon = start_positions[0]
        end_nucleotide, end_exon = end_positions[-1]
        # Swap em because they are negative, this feels dumb but works?
        start_nucleotide, end_nucleotide = end_nucleotide, start_nucleotide
    return pd.Series([start_nucleotide, end_nucleotide, start_exon, end_exon])

def process_mmseqs_to_genome(mmseqs_output_file,annotation_df,output_directory):
    ### Read in the MMSeqs results to a DataFrame
    mmseqs_dataframe = pd.read_csv(mmseqs_output_file, sep = '\t')

    ### Add 'parent' information for each hit:
    mmseqs_dataframe
    slim_df_metadata = annotation_df[['name', 'strand','seqid','parent_id','nucleotide_to_aa_mapping']]
    mmseqs_dataframe_meta = pd.merge(mmseqs_dataframe, slim_df_metadata, left_on='query', right_on='name', how='left')

    # Apply the function to each row and create new columns
    mmseqs_dataframe_meta[['genome_start', 'genome_end', 'exon_start', 'exon_end']] = mmseqs_dataframe_meta.apply(extract_nucleotide_positions, axis=1)

    # Convert the relevant columns to integers
    mmseqs_dataframe_meta['genome_start'] = mmseqs_dataframe_meta['genome_start'].astype(int)
    mmseqs_dataframe_meta['genome_end'] = mmseqs_dataframe_meta['genome_end'].astype(int)
    mmseqs_dataframe_meta['bits'] = mmseqs_dataframe_meta['bits'].astype(int)
    mmseqs_dataframe_meta['thickStart'] = mmseqs_dataframe_meta['genome_start'].astype(int)  # Default thickStart
    mmseqs_dataframe_meta['thickEnd'] = mmseqs_dataframe_meta['genome_end'].astype(int)  # Default thickEnd

    mmseqs_dataframe_meta['itemRgb'] = mmseqs_dataframe_meta['tcov'].apply(bits_to_rgb, args=(0, 1.0))
    mmseqs_dataframe_meta_sort = mmseqs_dataframe_meta.sort_values(by=['seqid', 'genome_start'])
    mmseqs_dataframe_meta_sort['unique_id'] = mmseqs_dataframe_meta_sort['parent_id'] + "_" + mmseqs_dataframe_meta_sort['target']
    # Specify the columns to write
    columns_to_write = ['seqid', 'genome_start', 'genome_end', 'unique_id', 'tcov' ,'strand','thickStart',"thickEnd","itemRgb"]

    # Write the specified columns to a BED file
    output_file = os.path.join(output_directory,"filtered_proteins.mmseqs.relativeToGenome.bed")

    mmseqs_dataframe_meta_sort.to_csv(output_file, sep='\t', columns=columns_to_write, header=False, index=False)

    print(f"Data written to {output_file}")