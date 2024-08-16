import os
import pandas as pd

def create_hyperlinks_file(file_path, output_file):
    results_df = pd.read_csv(file_path, sep = '\t')

    base_url = "https://www.ncbi.nlm.nih.gov/gene/?term="
    
    with open(output_file, 'w') as f:
        f.write('<html><body>\n')
        for sequence in results_df['query'].unique():
            sequence = sequence.strip()  # Remove any surrounding whitespace or newlines
            if sequence:  # Ensure the sequence is not empty
                url = base_url + sequence
                f.write(f'<a href="{url}">{sequence}</a><br>\n')
        f.write('</body></html>\n')
    print(f"Hyperlinks file created: {output_file}")

# Specify the path to your file containing protein sequences and the output folder
#file_path = r"U:\\AnnotationCheckerWithStructure\\Development\\Redux4.UncharacterisedOnly\\filtered_proteins.mmseqs.out.finalResults.tsv"
#output_folder = r"U:\\AnnotationCheckerWithStructure\\Development\\Redux4.UncharacterisedOnly\\filtered_proteins.mmseqs.out.finalResults.links.html"
#create_hyperlinks_file(file_path, output_folder)
