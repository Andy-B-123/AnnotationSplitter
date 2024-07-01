import argparse
from icecream import ic
import pandas as pd
import AnnotationParsing
import MMSeqs_runner

fasta_path = "U:/Chap4/Sf9_denovo_transcriptome/references/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna"
gff_path = "U:/Chap4/Sf9_denovo_transcriptome/references/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff"

#def main(fasta_path, gff_path):
gff_db = load_or_create_gff_db(gff_path)
ic(gff_db)
print(f'Number of genes in annotation: {gff_db.num_matches(biotype="gene")}')
print(f'Number of transcripts in annotation: {gff_db.num_matches(biotype="mRNA")}')
print(f'Average number transcripts per gene: {round(gff_db.num_matches(biotype="mRNA")/gff_db.num_matches(biotype="gene"),2)}')

genes_list = list(gff_db.get_records_matching(biotype="gene"))
protein_coding_gene_list = [gene["attributes"].split(";")[0].split("=")[1] for gene in genes_list if "gene_biotype=protein_coding" in gene["attributes"]]
print(f'Number of protein coding genes: {len(protein_coding_gene_list)}')

### Processing to get longest isoform for each gene
transcript_list = list(gff_db.get_records_matching(biotype="mRNA"))
transcript_df = pd.DataFrame(transcript_list)
transcript_df['gene_id'] = transcript_df['attributes'].apply(extract_gene_id)
transcript_df['length'] = transcript_df['stop'] - transcript_df['start']
largest_isoforms = transcript_df.loc[transcript_df.groupby('gene_id')['length'].idxmax()]
largest_isoforms = largest_isoforms[['gene_id', 'name', 'seqid', 'start', 'stop', 'length']].reset_index()

### Extract the DNA and protein sequences from the annotation with the fasta file
CDS_list = list(gff_db.get_records_matching(biotype="CDS"))
CDS_df = pd.DataFrame(CDS_list)
filtered_CDS_df = CDS_df[CDS_df['parent_id'].isin(largest_isoforms['name'])].reset_index()
CDS_df_with_sequences = extract_dna_sequences(filtered_CDS_df, fasta_path)
CDS_df_with_proteins = translate_to_protein(CDS_df_with_sequences)

### Write a nice summary dataframe for later
output_csv_file = 'protein_sequences.csv'
CDS_df_with_proteins.to_csv(output_csv_file, columns=['seqid', 'source', 'start', 'biotype','start','stop', 'strand', 'attributes', 'dna_sequence', 'name', 'parent_id','protein_sequence', 'is_partial'], index=False)
print(f"DataFrame written to {output_csv_file}")

### Write proteins to file
output_fasta_file = 'protein_sequences.fasta'
write_protein_sequences_to_fasta(CDS_df_with_proteins[CDS_df_with_proteins['is_partial'] != True], output_fasta_file)
print(f"Protein sequences written to {output_fasta_file}")

### Check that MMSeqs is installed and that the SwissProt database exists
mmseqs_path = check_mmseqs_existence(mmseqs_path)
database_path = check_database(database_path)

### Run MMSeqs with the SwissProt database on the generated protein fasta file
run_mmseqs(output_fasta_file, database_path, output_directory, mmseqs_path, threads)

print("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a genome and annotation file and try and identify mis-annotated genes.')
    parser.add_argument('fasta_path', type=str, help='Path to the FASTA file')
    parser.add_argument('gff_path', type=str, help='Path to the GFF file')
    parser.add_argument('--output_folder', type=str, help='Path to where the output data should be stored.')
    parser.add_argument('--database_path', type=str, help='Path to where the database has been downloaded (or needs to be downloaded).')
    parser.add_argument('--mmseqs_path', type=str, default=None, help='Path to the mmseqs executable. If not provided, the system PATH will be used.')
    parser.add_argument('--threads', type=int, default=16, help='Number of threads to run with mmseqs, defaults to 16.')

    args = parser.parse_args()
    main(args.fasta_path, args.gff_path, args.output_folder, args.database_path, args.mmseqs_path, args.threads)