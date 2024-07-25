import webbrowser

def open_genbank_pages(file_path):
    base_url = "https://www.ncbi.nlm.nih.gov/gene/?term="
    
    with open(file_path, 'r') as file:
        sequences = file.readlines()
    
    for sequence in sequences:
        sequence = sequence.strip()  # Remove any surrounding whitespace or newlines
        if sequence:  # Ensure the sequence is not empty
            url = base_url + sequence
            webbrowser.open(url)

# Specify the path to your file containing protein sequences
file_path = r"U:\AnnotationCheckerWithStructure\Development\TCA_New\results.final.list.order"
open_genbank_pages(file_path)