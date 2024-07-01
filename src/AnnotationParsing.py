import os
import pandas as pd
from cogent3 import load_annotations
from cogent3.core.annotation_db import GffAnnotationDb
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def load_or_create_gff_db(gff_path):
    db_path = f"{gff_path}.db"
    if os.path.exists(db_path):
        print(f"Loading existing database from {db_path}")
        gff_db = GffAnnotationDb(source=db_path)
    else:
        print(f"Creating new database from {gff_path}")
        gff_db = load_annotations(path=gff_path)
        gff_db.write(db_path)
    return gff_db

def extract_gene_id(attributes):
    for attribute in attributes.split(';'):
        if attribute.startswith('Parent='):
            return attribute.split('=')[1]
    return None

def reverse_complement(sequence):
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'
    }
    return ''.join(complement[base] for base in reversed(sequence))

def extract_dna_sequences(CDS_df, fasta_file):
    fasta = Fasta(fasta_file)
    def get_sequence(row):
        seqid = row['seqid']
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
    def get_protein_sequence(row):
        dna_sequence = row['dna_sequence']
        protein_metadata = row['attributes']
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

def write_protein_sequences_to_fasta(CDS_df, output_file):
    records = []
    for index, row in CDS_df.iterrows():
        header = row['name'].replace('cds-', '')
        gene_id = extract_gene_id(row['attributes'])
        transcript_id = row['parent_id'].replace('rna-', '')
        description = f"geneID={gene_id} transcriptID={transcript_id}"
        sequence = row['protein_sequence']
        record = SeqRecord(Seq(sequence), id=header, description=description)
        records.append(record)
    with open(output_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

