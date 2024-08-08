import argparse
import os
from icecream import ic
import pandas as pd
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

    gff_db = load_or_create_gff_db(gff_path)
    ic(gff_db)
    print(f'Number of genes in annotation: {gff_db.num_matches(biotype="gene")}')
    print(f'Number of transcripts in annotation: {gff_db.num_matches(biotype="mRNA")}')
    print(f'Average number transcripts per gene: {round(gff_db.num_matches(biotype="mRNA")/gff_db.num_matches(biotype="gene"),2)}')

    genes_list = list(gff_db.get_records_matching(biotype="gene"))
    protein_coding_gene_list = [gene["attributes"].split(";")[0].split("=")[1] for gene in genes_list if "biotype=protein_coding" in gene["attributes"]]
    print(f'Number of protein coding genes: {len(protein_coding_gene_list)}')

    ### Processing to get longest isoform for each gene
    print("Working on getting the longest isoform for each gene from the annotation file.")
    transcript_list = list(gff_db.get_records_matching(biotype="mRNA"))
    transcript_df = pd.DataFrame(transcript_list)
    transcript_df['gene_id'] = transcript_df['attributes'].apply(extract_gene_id)
    transcript_df['length'] = transcript_df['stop'] - transcript_df['start']
    largest_isoforms = transcript_df.loc[transcript_df.groupby('gene_id')['length'].idxmax()]
    largest_isoforms = largest_isoforms[['gene_id', 'name', 'seqid', 'start', 'stop', 'length']].reset_index()

    ### Extract the DNA and protein sequences from the annotation with the fasta file
    print("Get the protein sequences for each isoform.")
    CDS_list = list(gff_db.get_records_matching(biotype="CDS"))
    CDS_df = pd.DataFrame(CDS_list)
    filtered_CDS_df = CDS_df[CDS_df['parent_id'].isin(largest_isoforms['name'])].reset_index()
    CDS_df_with_sequences = extract_dna_sequences(filtered_CDS_df, fasta_path)
    CDS_df_with_proteins = translate_to_protein(CDS_df_with_sequences)
    CDS_df_with_proteins['nucleotide_to_aa_mapping'] = CDS_df_with_proteins.apply(lambda row: map_nucleotides_to_amino_acids(row['spans'], row['strand']), axis=1)

    ### Write a nice summary dataframe for later
    output_csv_file = f'{output_folder}/protein_sequences.csv'
    CDS_df_with_proteins.to_csv(output_csv_file, columns=['seqid', 'source', 'start', 'biotype','start','stop', 'strand', 'attributes', 'name', 'parent_id','protein_sequence', 'is_partial','nucleotide_to_aa_mapping'], index=False)
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

    ### Process MMSeqs output for clusters
    mmseqs_output_file = f'{output_folder}/filtered_proteins.mmseqs.out'
    df_of_bad_genes = convert_mmseqs_output(mmseqs_output_file,output_folder)

    ### Map hits to genome, make a nice bed file for viewing
    if create_bed_file:
        process_mmseqs_to_genome(mmseqs_output_file, df_of_bad_genes, CDS_df_with_proteins, output_folder)

def main():
    parser = argparse.ArgumentParser(description='Process a genome and annotation file and try and identify mis-annotated genes.')
    parser.add_argument('--fasta_path', type=str, required=True, help='Path to the FASTA file')
    parser.add_argument('--gff_path', type=str, required=True, help='Path to the GFF file')
    parser.add_argument('--output_folder', type=str, required=True, help='Path to where the output data should be stored.')
    parser.add_argument('--database_path', type=str, required=True, help='Path to where the database has been downloaded (or needs to be downloaded).')
    parser.add_argument('--mmseqs_path', type=str, default=None, help='Path to the mmseqs executable. If not provided, the system PATH will be used.')
    parser.add_argument('--mmseqs_params', type=str, default="--min-aln-len 100 -e 1.000E-010", help='Paramaters used instead of the default ones (min-aln-len 100 and e-filter -10). Does not change output format parameters.')
    parser.add_argument('--threads', type=int, default=16, help='Number of threads to run with mmseqs, defaults to 16.')
    parser.add_argument('--create_bed_file', action='store_true', help='If set, create a BED file.')

    args = parser.parse_args()
    run_pipeline(args.fasta_path, args.gff_path, args.output_folder, args.database_path, args.mmseqs_path, args.threads,args.mmseqs_params,args.create_bed_file)

if __name__ == "__main__":
    main()
