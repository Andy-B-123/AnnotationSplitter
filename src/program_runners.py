import os
import shutil
import subprocess
import requests
import tarfile
from pathlib import Path

def run_gffcompare(reference_gff, helixer_gff, output_directory):
    print("Running gffcompare on input annotations")
    output_prefix = os.path.join(output_directory, "gffcompare_output")
    cmd = f"gffcompare -r {reference_gff} {helixer_gff} -o {output_prefix}"
    subprocess.run(cmd, shell=True, check=True)
    print(f"gffcompare run completed. Output files prefixed with {output_prefix}")

def check_mmseqs_path_fine(mmseqs_path):
    """
    Checks if the 'mmseqs' program exists in the current system path or a provided path.
    If not, use the file in the 'program_binaries' directory of the repo.
    Returns the path to the 'mmseqs' executable to be used.
    """
    mmseqs_path = shutil.which("mmseqs")
    
    if mmseqs_path:
        print(f"'mmseqs' found in system path: {mmseqs_path}")
        return True

def run_mmseqs(protein_fasta_file, database_path, output_directory, mmseqs_path, threads=16, mmseqs_params="-s 7.5", include_seq=False):
    if check_mmseqs_path_fine(mmseqs_path):
        mmseqs_params_output = r"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,pident,qcov,tcov,qlen,tlen"
        if include_seq == True:
            mmseqs_params_output = r"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,pident,qcov,tcov,qlen,tlen,tseq,taxname,taxlineage"
        mmseqs_params_input = mmseqs_params
        mmseqs_params_input_list = mmseqs_params_input.split()
        base_output = protein_fasta_file.replace(".fasta","")
        output_file = f"{base_output}.mmseqs.out"
        if include_seq == True:
            output_file = f"{base_output}.mmseqs.database.hits"
        tmp_dir = output_directory + "/mmseqs.tmp"

        command = [
            str(mmseqs_path), "easy-search", str(protein_fasta_file), str(database_path), str(output_file), str(tmp_dir),
            "--format-mode", "4", "--format-output", mmseqs_params_output, "--threads", str(threads)] + mmseqs_params_input_list 
        print(command)
        print(f"Executing command: {' '.join(command)}")
        try:
            subprocess.run(command, check=True)
            print(f"Command executed successfully: {' '.join(command)}")
        except subprocess.CalledProcessError as e:
            print(f"Command failed with error: {e}")

