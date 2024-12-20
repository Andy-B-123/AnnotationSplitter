import argparse
import os
import subprocess
import warnings
from icecream import ic
from Bio import SeqIO
from AnnotationParsing import *
from program_runners import *
from MMSeqs_analysis import *


def parse_arguments():
    parser = argparse.ArgumentParser(description="GeneSplitter comparing reference annotations with Helixer generated annotations.")
    parser.add_argument("reference_gff", type=str, help="Path to the Reference GFF file")
    parser.add_argument("helixer_gff", type=str, help="Path to the Helixer GFF file")
    parser.add_argument("reference_fasta", type=str, help="Path to the Reference FASTA file")
    parser.add_argument("output_directory", type=str, help="Directory to store the output files")
    parser.add_argument("--mmseqs_db", type=str, required=True, help="Path to the MMseqs2  SwissProt database")
    parser.add_argument("--mmseqs_path", type=str, default=None, help="Path to the MMseqs2 executable, defaults to being present in the system path if not provided.")
    parser.add_argument("--mmseqs_threads", type=str, default=16, help="Threads to use for mmseqs.")
    return parser.parse_args()

def check_mmseqs_path(provided_path=None):
    """
    Checks if the 'mmseqs' program exists in the provided path or the current system path.
    Returns the path to the 'mmseqs' executable to be used.
    """
    if provided_path:
        if os.path.isfile(provided_path) and os.access(provided_path, os.X_OK):
            print(f"'mmseqs' found at provided path: {provided_path}")
            return provided_path
        else:
            raise FileNotFoundError(f"Provided 'mmseqs' path is not valid or not executable: {provided_path}")
    
    # Check if mmseqs is in the system path
    mmseqs_path = shutil.which("mmseqs")
    if mmseqs_path:
        print(f"'mmseqs' found in system path: {mmseqs_path}")
        return mmseqs_path
    else:
        raise FileNotFoundError("'mmseqs' not found in system path. Please provide a valid path using --mmseqs_path.")

def main():
    print("Starting process\n")
    args = parse_arguments()

    # Create the output directory if it doesn't exist
    os.makedirs(args.output_directory, exist_ok=True)
    # Run gffcompare
    gffcompare_output = os.path.join(args.output_directory, "gffcompare_output")
    if not os.path.exists(gffcompare_output + ".tracking"): 
        run_gffcompare(args.reference_gff, args.helixer_gff, args.output_directory)
    else:
        print("gffcompare output already exists, skipping...")

    # Process Reference Annotation
    ref_protein_output = os.path.join(args.output_directory, "reference_proteins.fasta")
    if not os.path.exists(ref_protein_output): 
        # Extract protein sequences for Reference
        print("For the Reference annotation:")
        annot_db_ref, annot_db_ref_gene_df, annot_db_ref_trans_df = summarize_annotations(args.reference_gff, "Reference")
        extract_protein_sequences(annot_db_ref,annot_db_ref_trans_df, args.reference_fasta, ref_protein_output)
    else:
        print("Reference protein output already exists, skipping...")

    # Process Helixer Annotation
    helixer_protein_output = os.path.join(args.output_directory, "helixer_proteins.fasta")
    if not os.path.exists(helixer_protein_output): 
        # Extract protein sequences for Helixer
        print("For the Helixer annotation:")
        annot_db_helixer, annot_db_helixer_gene_df, annot_db_helixer_trans_df = summarize_annotations(args.helixer_gff, "Helixer")
        extract_protein_sequences(annot_db_helixer,annot_db_helixer_trans_df, args.reference_fasta, helixer_protein_output)
    else:
        print("Helixer protein output already exists, skipping...")

    # Run mmseqs for each dataset:
    mmseqs_path = check_mmseqs_path(args.mmseqs_path)
    helixer_protein_output_mmseqs = os.path.join(args.output_directory, "helixer_proteins.mmseqs.out")
    reference_protein_output_mmseqs = os.path.join(args.output_directory, "reference_proteins.mmseqs.out")
    if not os.path.exists(helixer_protein_output_mmseqs) or not os.path.exists(reference_protein_output_mmseqs): 
        run_mmseqs(ref_protein_output, args.mmseqs_db, args.output_directory, mmseqs_path, threads=args.mmseqs_threads)
        run_mmseqs(helixer_protein_output, args.mmseqs_db, args.output_directory, mmseqs_path, threads=args.mmseqs_threads)
    else:
        print("Both mmseqs outputs exist, skipping...")

    # Process mmseqs data:
    filtered_hits_RefSeq, filtered_hits_Helixer = process_mmseqs_data(gffcompare_output, reference_protein_output_mmseqs, helixer_protein_output_mmseqs, args.output_directory)
    if filtered_hits_RefSeq is None:
        print("No detected chimeras, exiting!")
        exit()

    # Output summary tables and protein sequences:
    # Reference:
    reference_protein_output_mmseqs_summary= os.path.join(args.output_directory, "PotentialChimeras.RefSeqs.Summary.csv")
    idx = filtered_hits_RefSeq.groupby(['query'])['appended_cluster_count'].idxmax()
    filtered_hits_summary_table_refseq = filtered_hits_RefSeq.loc[idx].reset_index()[['query','Reference_ID','gene','appended_cluster_count']]
    filtered_hits_summary_table_refseq['link_raw'] = "https://www.ncbi.nlm.nih.gov/gene/" + filtered_hits_summary_table_refseq["gene"].str.replace('gene-LOC', '')
    filtered_hits_summary_table_refseq['link_clickable'] = '=HYPERLINK("' + filtered_hits_summary_table_refseq['link_raw'] + '")'
    filtered_hits_summary_table_refseq.to_csv(reference_protein_output_mmseqs_summary, index=False)

    ref_protein_output_potentialChimeras = os.path.join(args.output_directory, "reference_proteins.potential_chimeras.fasta")
    with open(ref_protein_output_potentialChimeras, 'w') as output_handle:
        for record in SeqIO.parse(ref_protein_output, 'fasta'):
            if record.id in filtered_hits_summary_table_refseq['query'].tolist():
                SeqIO.write(record, output_handle, 'fasta')
    # Helixer
    helixer_protein_output_mmseqs_summary = os.path.join(args.output_directory, "PotentialChimeras.Helixer.Summary.csv")
    filtered_hits_Helixer_summary = filtered_hits_Helixer.drop_duplicates(subset=['query', 'Reference_ID', 'Loc_ID']).sort_values(by='Reference_ID')[['query', 'Reference_ID', 'Loc_ID']].reset_index(drop= True)
    filtered_hits_Helixer_summary.to_csv(helixer_protein_output_mmseqs_summary, index=False)
        
    helixer_protein_output_potentialChimeras = os.path.join(args.output_directory, "helixer_proteins.potential_chimeras.fasta")
    with open(helixer_protein_output_potentialChimeras, 'w') as output_handle:
        for record in SeqIO.parse(helixer_protein_output, 'fasta'):
            if record.id in filtered_hits_summary_table_refseq['query'].tolist():
                SeqIO.write(record, output_handle, 'fasta')
    

if __name__ == "__main__":
    main()


""" 

output_directory = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/Dre"
reference_gff = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/raw_data_downloads/Danio_rerio.GRCz10.80.longest.gff"
reference_fasta = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/raw_data_downloads/Danio_rerio.GRCz10.dna.toplevel.fa"
helixer_gff = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/outputs/Danio_rerio.Helixer.gff"

output_directory = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/Sma"
reference_gff = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/raw_data_downloads/Strigamia_maritima.Smar1.26.longest.gff3"
reference_fasta = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/raw_data_downloads/Strigamia_maritima.Smar1.26.dna.toplevel.fa"
helixer_gff = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/outputs/Strigamia_maritima.Helixer.gff"

mmseqs_path = "/home/bac050/mmseqs/bin/mmseqs"
mmseqs_db = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/JustTryRegularMMSeqs/database_SwissProt/mmseqs_db_SProt_protostomia"
mmseqs_threads = 16

print("Starting process\n")

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)
# Run gffcompare
gffcompare_output = os.path.join(output_directory, "gffcompare_output")
if not os.path.exists(gffcompare_output + ".tracking"): 
    run_gffcompare(reference_gff, helixer_gff, output_directory)
else:
    print("gffcompare output already exists, skipping...")

# Process Reference Annotation
ref_protein_output = os.path.join(output_directory, "reference_proteins.fasta")
if not os.path.exists(ref_protein_output): 
    # Extract protein sequences for Reference
    print("For the Reference annotation:")
    annot_db_ref, annot_db_ref_gene_df, annot_db_ref_trans_df = summarize_annotations(reference_gff, "Reference")
    extract_protein_sequences(annot_db_ref,annot_db_ref_trans_df, reference_fasta, ref_protein_output)
else:
    print("Reference protein output already exists, skipping...")

# Process Helixer Annotation
helixer_protein_output = os.path.join(output_directory, "helixer_proteins.fasta")
if not os.path.exists(helixer_protein_output): 
    # Extract protein sequences for Helixer
    print("For the Helixer annotation:")
    annot_db_helixer, annot_db_helixer_gene_df, annot_db_helixer_trans_df = summarize_annotations(helixer_gff, "Helixer")
    extract_protein_sequences(annot_db_helixer,annot_db_helixer_trans_df, reference_fasta, helixer_protein_output)
else:
    print("Helixer protein output already exists, skipping...")

# Run mmseqs for each dataset:
check_mmseqs_path(mmseqs_path)
helixer_protein_output_mmseqs = os.path.join(output_directory, "helixer_proteins.mmseqs.out")
reference_protein_output_mmseqs = os.path.join(output_directory, "reference_proteins.mmseqs.out")
if not os.path.exists(helixer_protein_output_mmseqs) or not os.path.exists(reference_protein_output_mmseqs): 
    run_mmseqs(ref_protein_output, mmseqs_db, output_directory, mmseqs_path, threads=mmseqs_threads)
    run_mmseqs(helixer_protein_output, mmseqs_db, output_directory, mmseqs_path, threads=mmseqs_threads)
else:
    print("Both mmseqs outputs exist, skipping...")

# Process mmseqs data:
filtered_hits_RefSeq, filtered_hits_Helixer = process_mmseqs_data(gffcompare_output, reference_protein_output_mmseqs, helixer_protein_output_mmseqs, output_directory)
if filtered_hits_RefSeq is None:
    print("No detected chimeras, exiting!")
    exit()

# Output summary tables and protein sequences:
# Reference:
reference_protein_output_mmseqs_summary= os.path.join(output_directory, "PotentialChimeras.RefSeqs.Summary.csv")
idx = filtered_hits_RefSeq.groupby(['query'])['appended_cluster_count'].idxmax()
filtered_hits_summary_table_refseq = filtered_hits_RefSeq.loc[idx].reset_index()[['query','Reference_ID','gene','appended_cluster_count']]
filtered_hits_summary_table_refseq['link_raw'] = "https://www.ncbi.nlm.nih.gov/gene/" + filtered_hits_summary_table_refseq["gene"].str.replace('gene-LOC', '')
filtered_hits_summary_table_refseq['link_clickable'] = '=HYPERLINK("' + filtered_hits_summary_table_refseq['link_raw'] + '")'
filtered_hits_summary_table_refseq.to_csv(reference_protein_output_mmseqs_summary, index=False)

ref_protein_output_potentialChimeras = os.path.join(output_directory, "reference_proteins.potential_chimeras.fasta")
with open(ref_protein_output_potentialChimeras, 'w') as output_handle:
    for record in SeqIO.parse(ref_protein_output, 'fasta'):
        if record.id in filtered_hits_summary_table_refseq['query'].tolist():
            SeqIO.write(record, output_handle, 'fasta')
# Helixer
helixer_protein_output_mmseqs_summary = os.path.join(output_directory, "PotentialChimeras.Helixer.Summary.csv")
filtered_hits_Helixer_summary = filtered_hits_Helixer.drop_duplicates(subset=['query', 'Reference_ID', 'Loc_ID']).sort_values(by='Reference_ID')[['query', 'Reference_ID', 'Loc_ID']].reset_index(drop= True)
filtered_hits_Helixer_summary.to_csv(helixer_protein_output_mmseqs_summary, index=False)
    
helixer_protein_output_potentialChimeras = os.path.join(output_directory, "helixer_proteins.potential_chimeras.fasta")
with open(helixer_protein_output_potentialChimeras, 'w') as output_handle:
    for record in SeqIO.parse(helixer_protein_output, 'fasta'):
        if record.id in filtered_hits_summary_table_refseq['query'].tolist():
            SeqIO.write(record, output_handle, 'fasta')
 """
