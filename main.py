import argparse
from icecream import ic
import pandas as pd
from tqdm import tqdm
from AnnotationParsing import (
    load_or_create_gff_db, extract_gene_id, extract_dna_sequences,
    translate_to_protein, write_protein_sequences_to_fasta
)

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

transcript_list = list(gff_db.get_records_matching(biotype="mRNA"))
transcript_df = pd.DataFrame(transcript_list)
transcript_df['gene_id'] = transcript_df['attributes'].apply(extract_gene_id)
transcript_df['length'] = transcript_df['stop'] - transcript_df['start']
largest_isoforms = transcript_df.loc[transcript_df.groupby('gene_id')['length'].idxmax()]
largest_isoforms = largest_isoforms[['gene_id', 'name', 'seqid', 'start', 'stop', 'length']].reset_index()

CDS_list = list(gff_db.get_records_matching(biotype="CDS"))
CDS_df = pd.DataFrame(CDS_list)
filtered_CDS_df = CDS_df[CDS_df['parent_id'].isin(largest_isoforms['name'])].reset_index()

CDS_df_with_sequences = extract_dna_sequences(filtered_CDS_df, fasta_path)
CDS_df_with_proteins = translate_to_protein(CDS_df_with_sequences)

output_csv_file = 'protein_sequences.csv'
CDS_df_with_proteins.to_csv(output_csv_file, columns=['seqid', 'source', 'start', 'biotype','start','stop', 'strand', 'attributes', 'dna_sequence', 'name', 'parent_id','protein_sequence', 'is_partial'], index=False)
print(f"DataFrame written to {output_csv_file}")

output_fasta_file = 'protein_sequences.fasta'
write_protein_sequences_to_fasta(CDS_df_with_proteins[CDS_df_with_proteins['is_partial'] != True], output_fasta_file)
print(f"Protein sequences written to {output_fasta_file}")

print("Done")

""" if __name__ == "__main__":
    # Set these variables for direct assignment
    fasta_path = "U:/Chap4/Sf9_denovo_transcriptome/references/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.fna"
    gff_path = "U:/Chap4/Sf9_denovo_transcriptome/references/GCF_023101765.2_AGI-APGP_CSIRO_Sfru_2.0_genomic.gff"

    if fasta_path and gff_path:
        main(fasta_path, gff_path)
    else:
        parser = argparse.ArgumentParser(description='Process genome and annotation files.')
        parser.add_argument('fasta_path', type=str, help='Path to the FASTA file')
        parser.add_argument('gff_path', type=str, help='Path to the GFF file')
        args = parser.parse_args()
        main(args.fasta_path, args.gff_path)
 """