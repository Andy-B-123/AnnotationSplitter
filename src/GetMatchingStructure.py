import pandas as pd

# Load the data (equivalent to fread)
file_path = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/mmseqs.UniProtTrEMBL"
confirmed_chimerics_to_TrEMBL = pd.read_csv(file_path, sep="\t")

# Filter the data (equivalent to filter in R)
filtered_data = confirmed_chimerics_to_TrEMBL[
    (confirmed_chimerics_to_TrEMBL['fident'] > 0.95) & 
    (confirmed_chimerics_to_TrEMBL['qcov'] > 0.95) & 
    (confirmed_chimerics_to_TrEMBL['tcov'] > 0.95)
]
filtered_data.to_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/mmseqs.UniProtTrEMBL.filtered.csv")
# Display the filtered data
print(filtered_data)

import os
import requests
from Bio import PDB

# Function to download AlphaFold structure
def download_alphafold_structure(accession, output_dir="structures"):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Save the structure to a file
        file_path = os.path.join(output_dir, f"{accession}.pdb")
        with open(file_path, 'wb') as f:
            f.write(response.content)
        
        return file_path
    else:
        print(f"Failed to download structure for {accession}. Status code: {response.status_code}")
        return None

# Function to parse the PDB file and extract pLDDT for each amino acid
def extract_plddt_per_residue(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    
    structure = parser.get_structure('structure', pdb_file)
    plddt_per_residue = []

    for model in structure:
        for chain in model:l
            for residue in chain:
                # Extract pLDDT score (stored in B-factor field for each atom)
                residue_plddt = []
                for atom in residue:
                    residue_plddt.append(atom.bfactor)
                
                # Take the average pLDDT for the residue (across atoms)
                avg_plddt_residue = sum(residue_plddt) / len(residue_plddt) if residue_plddt else None
                plddt_per_residue.append(avg_plddt_residue)
    
    return plddt_per_residue

# Function to save pLDDT per residue to a file
def save_plddt_to_file(accession, plddt_per_residue, output_dir="plddt_scores"):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Save pLDDT scores to a file
    output_file = os.path.join(output_dir, f"{accession}_plddt.txt")
    with open(output_file, 'w') as f:
        for i, plddt in enumerate(plddt_per_residue, start=1):
            f.write(f"Residue {i}: {plddt}\n")

# Example usage
accessions = ["A0A7M7MSF0"]  # Add your list of 200 accessions here
failed_accessions = []

for accession in filtered_data['target'].unique().tolist():
    time.sleep(2)
    print(f"Processing {accession}...")
    pdb_file = download_alphafold_structure(accession, output_dir="/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/ReferenceLookups")
    
    if pdb_file:
        plddt_residues = extract_plddt_per_residue(pdb_file)
        save_plddt_to_file(accession, plddt_residues, output_dir="/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/ReferenceLookups")
        print(f"Saved pLDDT scores for {accession}")
    else:
        print(f"Failed to process {accession}")
        failed_accessions.append(accession)  