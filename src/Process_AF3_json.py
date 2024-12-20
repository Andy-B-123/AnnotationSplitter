### Visualise AlphaFold3 json output

import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Example function to load data from a JSON file
def load_alphafold3_output(json_file):
    with open(json_file, 'r') as f:
        return json.load(f)

def print_json_structure(data, indent=0):
    """
    Recursively prints the structure of a JSON object without printing all the data.
    """
    indent_str = '  ' * indent
    if isinstance(data, dict):
        for key, value in data.items():
            print(f"{indent_str}{key}: {type(value).__name__}")
            if isinstance(value, (dict, list)):
                print_json_structure(value, indent + 1)
    elif isinstance(data, list):
        print(f"{indent_str}List of {len(data)} items:")
        if data:
            # Print structure of first item to give an idea
            print_json_structure(data[0], indent + 1)

# Plotting functions
def plot_plddts(atom_plddts):
    plt.figure(figsize=(10, 5))
    plt.plot(atom_plddts, marker='o', linestyle='-', color='b')
    plt.title('Per-residue Confidence (pLDDT)')
    plt.xlabel('Residue Index')
    plt.ylabel('pLDDT Score')
    plt.grid(True)
    plt.show()

def plot_contact_probs(contact_probs):
    plt.figure(figsize=(10, 8))
    sns.heatmap(contact_probs, cmap='coolwarm', cbar_kws={'label': 'Contact Probability'})
    plt.title('Contact Probability Heatmap')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.show()

def plot_pae(pae):
    plt.figure(figsize=(10, 8))
    sns.heatmap(pae, cmap='coolwarm', cbar_kws={'label': 'PAE (Å)'})
    plt.title('Predicted Alignment Error (PAE) Heatmap')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.show()

input_json_path = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/AlphaFold3Results/fold_2024_09_18_14_28_xp_067170518_full_data_0.json"
data_json = load_alphafold3_output(input_json_path)
print_json_structure(data_json)

atom_plddts = data_json.get('atom_plddts', [])
contact_probs = data_json.get('contact_probs', [])
pae = data_json.get('pae', [])

contact_probs = np.array(contact_probs) if contact_probs else []
pae = np.array(pae) if pae else []

plot_plddts(atom_plddts)
plot_contact_probs(contact_probs)
plot_pae(pae)




import zipfile
import json
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Function to extract the 'full_data_0.json' file from a zip archive
def extract_json_from_zip(zip_path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Look for the file ending with 'full_data_0.json'
        for file in zip_ref.namelist():
            if file.endswith('full_data_0.json'):
                with zip_ref.open(file) as json_file:
                    return json.load(json_file), file  # Return JSON data and filename

# Function to save PAE data and associated filenames to a JSON file
def save_pae_data(pae_data_list, filenames, output_file="pae_data.json"):
    data_to_save = {
        "files": filenames,
        "paes": pae_data_list
    }
    # Save as JSON
    with open(output_file, 'w') as f:
        json.dump(data_to_save, f, indent=2)
    print(f"PAE data saved to {output_file}")

def plot_multiple_paes(pae_data_list, filenames, ncols=4):
    n = len(pae_data_list)
    nrows = (n + ncols - 1) // ncols  # Calculate number of rows required

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 3, nrows * 3))
    axes = axes.flatten()  # Flatten axes for easier iteration

    # Plot each PAE heatmap
    for i, (pae, filename) in enumerate(zip(pae_data_list, filenames)):
        heatmap = sns.heatmap(pae, cmap='coolwarm', cbar=False, ax=axes[i])
        axes[i].set_title(filename, fontsize=12)
        axes[i].axis('off')  # Turn off axis for minimal plot

    # Turn off any extra axes (in case the grid has more spaces than heatmaps)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    # Add a single color bar at the bottom of the grid
    cbar_ax = fig.add_axes([0.3, 0.05, 0.4, 0.02])  # Position [left, bottom, width, height]
    colorbar = plt.colorbar(heatmap.get_children()[0], cax=cbar_ax, orientation='horizontal')
    
    # Set color bar label for confidence levels
    colorbar.set_label('Confidence (Blue: High, Red: Low)', fontsize=10)
    cbar_ax.xaxis.set_ticks_position('bottom')  # Ensure ticks are at the bottom

    # Adjust layout to make space for the color bar at the bottom
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])
    plt.show()

# Main function to process all zip files in a directory, save and plot the PAEs
def process_zip_files(directory, output_file="pae_data.json"):
    pae_data_list = []
    filenames = []

    # Iterate over all zip files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.zip'):
            zip_path = os.path.join(directory, filename)
            json_data, json_filename = extract_json_from_zip(zip_path)
            
            # Extract the PAE data from the JSON
            pae = np.array(json_data.get('pae', []))
            
            # Append the PAE data and filename to the lists
            if pae.size > 0:
                pae_data_list.append(pae.tolist())  # Convert numpy array to list for JSON saving
                filenames.append(json_filename)
    
    # Save the PAE data and associated filenames to a file
    if pae_data_list:
        save_pae_data(pae_data_list, filenames, output_file)

    # Plot all PAEs if any data was collected
    if pae_data_list:
        plot_multiple_paes(pae_data_list, filenames)

directory = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/AlphaFold3Results/"  # Replace with your directory path
output_file = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/AlphaFold3Results/pae_data.json"  # Replace with your desired output filename
process_zip_files(directory, output_file)


import zipfile
import json
import os
import numpy as np

# Function to extract the 'full_data_0.json' file from a zip archive
def extract_json_from_zip(zip_path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Look for the file ending with 'full_data_0.json'
        for file in zip_ref.namelist():
            if file.endswith('full_data_0.json'):
                with zip_ref.open(file) as json_file:
                    return json.load(json_file)  # Return JSON data and filename

# Load data from a directory once
def load_pae_data(directory):
    pae_data_list = []
    filenames = []

    # Iterate over all zip files in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.zip'):
            zip_path = os.path.join(directory, filename)
            json_data = extract_json_from_zip(zip_path)
            
            # Extract the PAE data from the JSON
            pae = np.array(json_data.get('pae', []))
            
            # Append the PAE data and filename to the lists
            if pae.size > 0:
                pae_data_list.append(pae)
                filenames.append(filename)
    
    return pae_data_list, filenames

# Example usage for loading
directory = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/check_structure/AlphaFold3Results/Subset"
pae_data_list, filenames = load_pae_data(directory)

import re

# Pair filenames with their associated data
paired = list(zip(filenames, pae_data_list))

# Extract the leading number from each filename and sort based on it
sorted_pairs = sorted(paired, key=lambda x: int(re.search(r'^\d+', x[0]).group()))

# Unzip the sorted pairs back into separate lists
sorted_filenames, sorted_pae_data_list = zip(*sorted_pairs)

# Convert to lists (optional, if further modification is needed)
sorted_filenames = list(sorted_filenames)
sorted_pae_data_list = list(sorted_pae_data_list)

# Output the results
print("Sorted Filenames:", sorted_filenames)
print("Sorted Data List:", sorted_pae_data_list)
# Apply transformations to each filename
tidied_filenames = [
    filename.replace('.zip', '')               # Remove ".zip"
            .replace('_1', '.1')              # Change "_1" to ".1"
            .replace('_fold_rna_', '. ')       # Replace "_fold_rna-" with a single space
            .replace('xm', 'XM')              # Make "xm" uppercase (without underscores)
    for filename in sorted_filenames
]

# Output the tidied filenames
tidied_filenames


import matplotlib.pyplot as plt
import seaborn as sns

def plot_multiple_paes(pae_data_list, filenames, ncols=4, save_as_png=False, output_file="pae_heatmaps.png"):
    n = len(pae_data_list)
    nrows = (n + ncols - 1) // ncols  # Calculate number of rows required

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 3, nrows * 3))
    axes = axes.flatten()  # Flatten axes for easier iteration

    # Plot each PAE heatmap using 'cividis' colormap (colorblind-friendly)
    for i, (pae, filename) in enumerate(zip(pae_data_list, filenames)):
        heatmap = sns.heatmap(pae, cmap='YlOrRd_r', cbar=False, ax=axes[i])
        axes[i].set_title(filename, fontsize=12)
        axes[i].axis('off')  # Turn off axis for minimal plot

    # Turn off any extra axes (in case the grid has more spaces than heatmaps)
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')

    # Add a single color bar at the bottom of the grid
    cbar_ax = fig.add_axes([0.3, 0.05, 0.4, 0.02])  # Position [left, bottom, width, height]
    colorbar = plt.colorbar(heatmap.get_children()[0], cax=cbar_ax, orientation='horizontal')

    # Set custom ticks and labels for color bar
    colorbar.set_label('Expected Position Error (0: Low, 30: High)', fontsize=10)  # Custom label

    # Adjust layout to make space for the color bar at the bottom
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])

    # Option to save the plot to a PNG file
    if save_as_png:
        plt.savefig(output_file, format="png", dpi=300)
        print(f"Plot saved as {output_file}")
    else:
        plt.show()

# Example usage for plotting
plot_multiple_paes(sorted_pae_data_list, tidied_filenames, ncols=3, save_as_png=True, output_file="output_heatmaps.png")

