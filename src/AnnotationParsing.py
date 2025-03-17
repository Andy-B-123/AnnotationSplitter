import argparse
import os
import subprocess
import gffpandas.gffpandas as gffpd
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from pyfaidx import Fasta
from icecream import ic

pd.set_option('display.max_columns', None)

def summarize_annotations(gff_file, label):
    gff_db = gffpd.read_gff3(gff_file)

    # Suppress specific FutureWarning
    warnings.filterwarnings("ignore", category=FutureWarning, module='gffpandas')
    genes_df_to_columns = gff_db.filter_feature_of_type(['gene']).attributes_to_columns()

    trans_df_to_columns = gff_db.filter_feature_of_type(['transcript']).attributes_to_columns()
    mrna_df_to_columns = gff_db.filter_feature_of_type(['mRNA']).attributes_to_columns()

    # Check the number of records in each
    trans_records = len(trans_df_to_columns)
    mrna_records = len(mrna_df_to_columns)

    if trans_records > 0 and mrna_records > 0:
        if trans_records > mrna_records:
            ic(f"Transcript has more records: {trans_records} > {mrna_records}")
        else:
            ic(f"mRNA has more records: {mrna_records} > {trans_records}")
            trans_df_to_columns = mrna_df_to_columns
    elif trans_records > 0:
        ic(f"Only transcript records present, count: {trans_records}")
    elif mrna_records > 0:
        ic(f"Only mRNA records present, count: {mrna_records}")
        trans_df_to_columns = mrna_df_to_columns
    else:
        ic("No transcript or mRNA records found.")

    if label == "Reference":
        column_name_list = genes_df_to_columns.columns.tolist()
        biotype_index, biotype_column, protein_coding_df = None, None, None
        for index, column_header in enumerate(column_name_list):
            if 'biotype' in column_header:
                biotype_index = index
                biotype_column = column_header
        if biotype_index == None:
            print("Hmmm no biotype present in attributes of 'genes', going to assume all are protein coding!")
            protein_coding_df = genes_df_to_columns.reset_index(drop=True)
        else:
            protein_coding_df = genes_df_to_columns[genes_df_to_columns[biotype_column] == 'protein_coding'].reset_index(drop=True)
        protein_coding_df_trans = trans_df_to_columns[trans_df_to_columns['Parent'].isin(protein_coding_df['ID'])]
        print(f'Number of protein coding genes in {label} annotation: {len(protein_coding_df)}')
        print(f'Number of transcripts in annotation: {len(protein_coding_df_trans)}')
        print(f'Average number transcripts per gene: {round(((len(protein_coding_df_trans))/(len(protein_coding_df))),2)}')
    elif label == "Helixer":
        print("Helixer only outputs single, longest isoform per gene")
        protein_coding_df = genes_df_to_columns.reset_index(drop=True)
        protein_coding_df_trans = trans_df_to_columns[trans_df_to_columns['Parent'].isin(protein_coding_df['ID'])]
        print(f'Number of protein coding genes in {label} annotation: {len(protein_coding_df)}')
        print(f'Number of transcripts in annotation: {len(protein_coding_df_trans)}')
        print(f'Average number transcripts per gene: {round(((len(protein_coding_df_trans))/(len(protein_coding_df))),2)}')
    else:
        print("Something very wrong, neither Helixer nor Reference provided...")
    print("")
    return gff_db, protein_coding_df, protein_coding_df_trans

def reverse_complement(sequence):
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n',
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 
        'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', '-': '-', 
        'r': 'y', 'y': 'r', 's': 's', 'w': 'w', 'k': 'm', 'm': 'k', 
        'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b'
    }
    
    return ''.join(complement.get(base, base) for base in reversed(sequence))

def extract_dna_sequences(CDS_df, fasta_file):
    fasta = Fasta(fasta_file)
    def get_sequence(row):
        seqid = row['seq_id']
        spans = row['spans']
        strand = row['strand']
        sequence = ''
        for span in spans:
            start, end = span
            sequence += fasta[seqid][start:end].seq
        if strand == '-':
            sequence = reverse_complement(sequence)
        return sequence
    CDS_df['dna_sequence'] = CDS_df.apply(get_sequence, axis=1)
    fasta.close()
    return CDS_df

def translate_to_protein(CDS_df):
    CDS_df['attributes'] = CDS_df['attributes'].fillna("unknown")
    def get_protein_sequence(row):
        dna_sequence = row['dna_sequence']
        protein_metadata = row['attributes']
        if isinstance(protein_metadata, float):
            print(f"protein_metadata is a float for row: {row.name}, value: {protein_metadata}")
            print(protein_metadata)
        is_partial = len(dna_sequence) % 3 != 0
        if is_partial or "partial=true" in protein_metadata:
            trimmed_sequence = dna_sequence[:-(len(dna_sequence) % 3)]
        else:
            trimmed_sequence = dna_sequence
        my_seq = Seq(trimmed_sequence)
        protein_sequence = my_seq.translate()
        return str(protein_sequence), is_partial
    CDS_df[['protein_sequence', 'is_partial']] = CDS_df.apply(
        lambda row: pd.Series(get_protein_sequence(row)), axis=1)
    return CDS_df

def clean_cds_id(id_string):
    # Split the string at the last occurrence of ".CDS"
    cleaned_id = id_string.rsplit(".CDS", 1)[0] + ".CDS"
    return cleaned_id

def write_protein_sequences_to_fasta(CDS_df, output_file):
    records = []
    for index, row in CDS_df.iterrows():
        header = row['Parent']
        transcript_id = row['ID'].replace("cds-", "").replace("transcript:", "")
        description = f"transcriptID={transcript_id}"
        sequence = row['protein_sequence']
        record = SeqRecord(Seq(sequence), id=header, description=description)
        records.append(record)
    with open(output_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

def extract_protein_sequences(annot_db, annot_db_trans, fasta_file, output_fasta_path):
    ### Processing to get longest isoform for each gene, redudant for Helixer
    print("Working on getting the longest isoform for each gene from the annotation file.")
    transcript_df = annot_db_trans.copy()
    transcript_df['length'] = transcript_df['end'] - transcript_df['start']
    largest_isoforms = transcript_df.loc[transcript_df.groupby('Parent')['length'].idxmax()]
    largest_isoforms = largest_isoforms[['Parent', 'ID', 'seq_id', 'start', 'end', 'length']].reset_index(drop = True)
    
    ### Extract the DNA and protein sequences from the annotation with the fasta file
    print("Get the protein sequences for each isoform.")
    CDS_df = annot_db.filter_feature_of_type(['CDS']).attributes_to_columns()
    filtered_CDS_df = CDS_df[CDS_df['Parent'].isin(largest_isoforms['ID'])].reset_index(drop = True)
    
    if "helixer_proteins.fasta" in output_fasta_path:
        filtered_CDS_df["ID"] = filtered_CDS_df["ID"].apply(clean_cds_id)
        filtered_CDS_df["attributes"] = filtered_CDS_df["attributes"].apply(clean_cds_id)
    
    filtered_CDS_df_data = filtered_CDS_df.iloc[:, 8:].drop_duplicates()
    filtered_CDS_df['seq_id'] = filtered_CDS_df['seq_id'].astype(str)
    
    # Modified groupby operation
    def get_spans(group):
        return pd.Series({
            'spans': list(zip(group['start'] - 1, group['end']))
        })

    grouped = (filtered_CDS_df.sort_values(by='start')
              .groupby(['Parent', 'strand', 'seq_id'])
              .apply(get_spans)
              .reset_index())
    
    merged_df = pd.merge(grouped, filtered_CDS_df_data, on='Parent', how='left')
    CDS_df_with_sequences = extract_dna_sequences(merged_df, fasta_file)
    CDS_df_with_proteins = translate_to_protein(CDS_df_with_sequences)
    write_protein_sequences_to_fasta(CDS_df_with_proteins, output_fasta_path)


