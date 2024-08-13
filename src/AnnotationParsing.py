import os
import pandas as pd
import gffpandas.gffpandas as gffpd
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def load_gff(gff_path):
    annotation = gffpd.read_gff3(gff_path)
    return annotation

def extract_gene_id(attributes):
    for attribute in attributes.split(';'):
        if attribute.startswith('Parent='):
            return attribute.split('=')[1]
    return None

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

def map_nucleotides_to_amino_acids(CDS_spans, strand='+'):
    # Initialize the mapping dictionary
    nucleotide_to_aa_exon = {}
    nucleotide_index = 0

    # Adjust the order of CDS spans based on the strand
    if strand == '-':
        CDS_spans = [(end, start) for start, end in CDS_spans[::-1]]
    for exon_number, (start, end) in enumerate(CDS_spans, start=1):
        if strand == '+':
            range_nucleotides = range(start, end + 1)
        else:  # strand == '-'
            range_nucleotides = range(start, end - 1, -1)
        for nucleotide in range_nucleotides:
            # Calculate the corresponding amino acid position
            amino_acid_position = nucleotide_index // 3 + 1
            nucleotide_to_aa_exon[nucleotide] = (amino_acid_position, exon_number)
            nucleotide_index += 1

    # Create reverse mapping for amino acid positions to nucleotide positions and exon numbers
    reverse_mapping = {}
    for nucleotide, (amino_acid_pos, exon_num) in nucleotide_to_aa_exon.items():
        if amino_acid_pos not in reverse_mapping:
            reverse_mapping[amino_acid_pos] = []
        reverse_mapping[amino_acid_pos].append((nucleotide, exon_num))
    return reverse_mapping

def write_protein_sequences_to_fasta(CDS_df, output_file):
    records = []
    gene_column_exists = 'gene' in CDS_df.columns

    for index, row in CDS_df.iterrows():
        header = row['ID'].replace('cds-', '')
        if gene_column_exists and row['gene']:
            gene_id = row['gene']
        else:
            gene_id = "NotAvailable"
        transcript_id = row['Parent'].replace('rna-', '')
        description = f"geneID={gene_id} transcriptID={transcript_id}"
        sequence = row['protein_sequence']
        record = SeqRecord(Seq(sequence), id=header, description=description)
        records.append(record)
    with open(output_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

def write_protein_sequences_to_fasta_from_mmseqs(mmseqs_file, output_file):
    mmseqs_df = pd.read_csv(mmseqs_file, sep ='\t')
    records = []
    for index, row in mmseqs_df.iterrows():
        header = row['target'].replace('cds-', '')
        gene_id = row['query']
        description = f"geneID={gene_id} organism={row['taxname']}"
        sequence = row['tseq']
        record = SeqRecord(Seq(sequence), id=header, description=description)
        records.append(record)
    with open(output_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')
