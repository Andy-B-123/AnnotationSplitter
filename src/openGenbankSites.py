import os

def create_hyperlinks_file(file_path, output_file):
    base_url = "https://www.ncbi.nlm.nih.gov/gene/?term="
    
    with open(file_path, 'r') as file:
        sequences = file.readlines()
    
    with open(output_file, 'w') as f:
        f.write('<html><body>\n')
        for sequence in sequences:
            sequence = sequence.strip()  # Remove any surrounding whitespace or newlines
            if sequence:  # Ensure the sequence is not empty
                url = base_url + sequence
                f.write(f'<a href="{url}">{sequence}</a><br>\n')
        f.write('</body></html>\n')
    print(f"Hyperlinks file created: {output_file}")

# Specify the path to your file containing protein sequences and the output folder
file_path = r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\results\\ESF2_results.list"
output_folder = r"U:\AnnotationCheckerWithStructure\Development\\check_other_organisms\\results\\GenBankLinks.html"

create_hyperlinks_file(file_path, output_folder)