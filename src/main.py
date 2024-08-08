import argparse
import os
from icecream import ic
import pandas as pd
import warnings
from AnnotationParsing import *
from MMSeqs_runner import *
from ProteinHitsToGenomeCoordinates import *
from Clustering_trials import *

def create_and_check_output_folder(output_folder):
    # Create the output folder and any necessary parent directories
    os.makedirs(output_folder, exist_ok=True)
    if os.listdir(output_folder):
        print(f"Warning: The output folder '{output_folder}' is not empty.")

def run_pipeline(fasta_path, gff_path, output_folder, database_path, mmseqs_path, threads, mmseqs_params,create_bed_file):

    create_and_check_output_folder(output_folder)


    gff_db = load_gff(gff_path)
    # Suppress specific FutureWarning
    warnings.filterwarnings("ignore", category=FutureWarning, module='gffpandas')
    gff_stats = gff_db.stats_dic()
    annot_num_genes = gff_stats.get('Counted_feature_types').get('gene')
    annot_num_trans = gff_stats.get('Counted_feature_types').get('mRNA')
    print(f'Number of genes in annotation: {annot_num_genes}')
    print(f'Number of transcripts in annotation: {annot_num_trans}')
    print(f'Average number transcripts per gene: {round(annot_num_trans/annot_num_genes,2)}')

    genes_df = gff_db.filter_feature_of_type(['gene'])
    genes_df_to_columns = genes_df.attributes_to_columns()
    column_name_list = genes_df_to_columns.columns.tolist()
    biotype_index = None
    biotype_column = None
    protein_coding_df = None
    for index, column_header in enumerate(column_name_list):
        if 'biotype' in column_header:
            biotype_index = index
            biotype_column = column_header
    if biotype_index == None:
        print("Hmmm no biotype present in attributes of 'genes', going to assume all are protein coding!")
        protein_coding_df = genes_df_to_columns.reset_index(drop=True)
    else:
        protein_coding_df = genes_df_to_columns[genes_df_to_columns[biotype_column] == 'protein_coding'].reset_index(drop=True)
    print(f'Number of protein coding genes: {len(protein_coding_df)}')

    ### Processing to get longest isoform for each gene
    print("Working on getting the longest isoform for each gene from the annotation file.")
    transcript_df = gff_db.filter_feature_of_type(['mRNA']).attributes_to_columns()
    transcript_df['length'] = transcript_df['end'] - transcript_df['start']
    largest_isoforms = transcript_df.loc[transcript_df.groupby('Parent')['length'].idxmax()]
    largest_isoforms = largest_isoforms[['Parent', 'ID', 'seq_id', 'start', 'end', 'length']].reset_index(drop = True)

    ### Extract the DNA and protein sequences from the annotation with the fasta file
    print("Get the protein sequences for each isoform.")
    CDS_df = gff_db.filter_feature_of_type(['CDS']).attributes_to_columns()
    filtered_CDS_df = CDS_df[CDS_df['Parent'].isin(largest_isoforms['ID'])].reset_index(drop = True)
    filtered_CDS_df_data = filtered_CDS_df.iloc[:, 8:].drop_duplicates()

    grouped = filtered_CDS_df.sort_values(by='start').groupby(['Parent','strand', 'seq_id'], group_keys=False).apply(
        lambda x: pd.Series({
            'spans': [(row['start']-1, (row['end'])) for idx, row in x.iterrows()]
        })
    ).reset_index()
    merged_df = pd.merge(grouped, filtered_CDS_df_data, on='Parent', how='left')
    CDS_df_with_sequences = extract_dna_sequences(merged_df, fasta_path)
    CDS_df_with_proteins = translate_to_protein(CDS_df_with_sequences)
    CDS_df_with_proteins['nucleotide_to_aa_mapping'] = CDS_df_with_proteins.apply(lambda row: map_nucleotides_to_amino_acids(row['spans'], row['strand']), axis=1)
    CDS_df_with_proteins['ID'] = CDS_df_with_proteins['ID'].str.replace('cds-', '')
    CDS_df_with_proteins['Parent'] = CDS_df_with_proteins['Parent'].str.replace('rna-', '')

    ### Write a nice summary dataframe for later
    output_csv_file = f'{output_folder}/protein_sequences.csv'
    excluded_columns = ['nucleotide_to_aa_mapping']
    CDS_df_with_proteins.drop(columns = excluded_columns).to_csv(output_csv_file, sep="\t", index=False)
    print(f"DataFrame written to {output_csv_file}")

    ### Write proteins to file
    output_fasta_file = f'{output_folder}/protein_sequences.fasta'
    write_protein_sequences_to_fasta(CDS_df_with_proteins[CDS_df_with_proteins['is_partial'] != True], output_fasta_file)
    print(f"Protein sequences written to {output_fasta_file}")

    ### Check that MMSeqs is installed and that the SwissProt database exists
    mmseqs_path = check_mmseqs_existence(mmseqs_path)
    database_path = check_database(database_path)

    ### Run MMSeqs with the SwissProt database on the generated protein fasta file
    run_mmseqs(output_fasta_file, database_path, output_folder, mmseqs_path, threads, mmseqs_params)

    print("Finished running MMSeqs.")

    os.environ['MPLCONFIGDIR'] = output_folder + "/plot_runtime"
    ### Process MMSeqs output for clusters
    mmseqs_output_file = output_folder + '/filtered_proteins.mmseqs.out'
    df_of_bad_genes = convert_mmseqs_output(mmseqs_output_file,output_folder)

    ### Map hits to genome, make a nice bed file for viewing
    if create_bed_file:
        process_mmseqs_to_genome(mmseqs_output_file, df_of_bad_genes, CDS_df_with_proteins, output_folder)

    ### Extract incorrect genes to fasta to make my life easier
    list_of_bad_genes = df_of_bad_genes['query'].unique().tolist()
    output_incorrect_fasta_file = output_folder + '/incorrect_protein_sequences.fasta'
    write_protein_sequences_to_fasta(CDS_df_with_proteins[CDS_df_with_proteins['ID'].isin(list_of_bad_genes)], output_incorrect_fasta_file)
    print(f"Protein sequences written to {output_incorrect_fasta_file}")

def main():
    parser = argparse.ArgumentParser(description='Process a genome and annotation file and try and identify mis-annotated genes.')
    parser.add_argument('--fasta_path', type=str, required=True, help='Path to the FASTA file')
    parser.add_argument('--gff_path', type=str, required=True, help='Path to the GFF file')
    parser.add_argument('--output_folder', type=str, required=True, help='Path to where the output data should be stored.')
    parser.add_argument('--database_path', type=str, required=True, help='Path to where the database has been downloaded (or needs to be downloaded).')
    parser.add_argument('--mmseqs_path', type=str, default=None, help='Path to the mmseqs executable. If not provided, the system PATH will be used.')
    parser.add_argument('--mmseqs_params', type=str, default="-s 7.5", help='Paramaters used instead of the default ones ("-s 7.5", high sensitivty). Does not change output format parameters.')
    parser.add_argument('--threads', type=int, default=16, help='Number of threads to run with mmseqs, defaults to 16.')
    parser.add_argument('--create_bed_file', action='store_true', help='If set, create a BED file.')

    args = parser.parse_args()
    run_pipeline(args.fasta_path, args.gff_path, args.output_folder, args.database_path, args.mmseqs_path, args.threads,args.mmseqs_params,args.create_bed_file)

if __name__ == "__main__":
    main()
