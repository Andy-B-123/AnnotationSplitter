import os
import shutil
import subprocess
import requests
import tarfile
from pathlib import Path

def check_database(prefix):
    """
    Checks if the database files with the specified prefix exist and if the .version file contains 'Swiss-Prot'. 
    If not, downloads them.

    Args:
        prefix (str): The prefix of the database files.
    """
    required_files = [
        f"{prefix}.version",
        f"{prefix}.source"
    ]

    # Check if the required files exist
    missing_files = [file for file in required_files if not os.path.exists(file)]

    if not missing_files:
        # Check if the .version file contains 'Swiss-Prot'
        with open(f"{prefix}.version", 'r') as version_file:
            content = version_file.read()
            if 'Swiss-Prot' in content:
                print("All required database files are present and the .version file has been verified as containing 'Swiss-Prot'.")
                return prefix
            else:
                print("The .version file does not contain 'Swiss-Prot'. Downloading the database files again...")
    else:
        print("Some required database files are missing. Downloading...")
    # Call download_database if necessary
    return download_database(prefix)


def download_database(prefix):
    """
    Downloads the database using the provided mmseqs command.

    Args:
        prefix (str): The prefix of the database files.
    """
    output_db_path = prefix
    tmp_db_download = "./tmp_db_download"

    # Ensure the temporary download directory exists
    if not os.path.exists(tmp_db_download):
        os.makedirs(tmp_db_download)

    # Construct the mmseqs command
    command = f"mmseqs databases UniProtKB/Swiss-Prot {output_db_path} {tmp_db_download}"

    try:
        # Execute the command
        subprocess.run(command, shell=True, check=True)
        print("Database downloaded and extracted successfully.")
        return(prefix)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while downloading the database: {e}")

def check_mmseqs_existence(mmseqs_path):
    """
    Checks if the 'mmseqs' program exists in the current system path or a provided path.
    If not, downloads it from a specified URL.

    Args:
        mmseqs_path (str): The path to check for the 'mmseqs' program. If None, checks the system PATH.
    """
    # URL to download mmseqs if it doesn't exist
    download_url = "https://github.com/soedinglab/MMseqs2/releases/download/latest/mmseqs-linux.tar.gz"
    download_path = "./mmseqs-linux.tar.gz"
    extract_path = "./mmseqs"

    # Check if mmseqs exists in the provided path
    if mmseqs_path:
        if os.path.isfile(mmseqs_path) and os.access(mmseqs_path, os.X_OK):
            print(f"'mmseqs' found at {mmseqs_path}.")
            return mmseqs_path
    else:
        # Check if mmseqs exists in the system PATH
        mmseqs_path = shutil.which("mmseqs")
        if mmseqs_path:
            print(f"'mmseqs' found in system PATH at {mmseqs_path}.")
            return mmseqs_path

    # If mmseqs is not found, download it
    print("'mmseqs' not found. Downloading from the provided URL...")
    response = requests.get(download_url, stream=True)
    if response.status_code == 200:
        with open(download_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
        print("Download complete. Extracting files...")

        # Extract the tar.gz file
        with tarfile.open(download_path, "r:gz") as tar:
            tar.extractall(path=extract_path)
        print("Extraction complete.")

        mmseqs_executable = os.path.join(extract_path, "mmseqs")
        if os.path.isfile(mmseqs_executable) and os.access(mmseqs_executable, os.X_OK):
            print(f"'mmseqs' is now available at {mmseqs_executable}.")
            return mmseqs_executable
        else:
            print("Extraction failed or 'mmseqs' is not executable.")
            return None
    else:
        print("Failed to download 'mmseqs'. Please check the URL and try again.")
        return None

def run_mmseqs(protein_fasta_file, database_path, output_directory, mmseqs_path, threads=16, mmseqs_params="-s 7.5", include_seq=False):
    mmseqs_params_output = r"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,pident,qcov,tcov,qlen,tlen"
    if include_seq == True:
        mmseqs_params_output = r"query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,pident,qcov,tcov,qlen,tlen,tseq,taxname,taxlineage"
    mmseqs_params_input = mmseqs_params
    mmseqs_params_input_list = mmseqs_params_input.split()
    output_file = output_directory + "/filtered_proteins.mmseqs.out"
    if include_seq == True:
        output_file = output_directory + "/incorrect_protein_sequences.Database_hits.out"
    tmp_dir = output_directory + "/filtered_proteins.mmseqs.tmp"

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

