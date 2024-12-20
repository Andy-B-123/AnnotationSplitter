import os
import pandas as pd

# Path to the main directory containing all subdirectories
base_dir = '/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/'
dataframes_results = []
# Loop through each folder in the base directory
for folder in os.listdir(base_dir):
    print(folder)
    folder_path = os.path.join(base_dir, folder)
    
    # Check if it is a directory
    if os.path.isdir(folder_path):
        # Define the file path
        file_path = os.path.join(folder_path, 'PotentialChimeras.RefSeqs.Summary.Update.csv')
        
        # Check if the file exists
        if os.path.isfile(file_path):
            # Read the file into a DataFrame
            df = pd.read_csv(file_path)  # Adjust 'sep' and 'header' as needed
            # Optionally, add a column for the folder name to track source
            df['source_folder'] = folder
            # Append the DataFrame to the list
            dataframes_results.append(df)
        else:
            print(f"File not found in {folder_path}")

# Concatenate all DataFrames in the list
all_data_results = pd.concat(dataframes_results, ignore_index=True)

# Display the combined DataFrame
print(all_data_results)

confirmed_chimerics = all_data_results.query('CorrectChimericLabelAssigned == "Yes"').copy().reset_index(drop=True)
confirmed_chimerics['chimeric_status'] = 1
unclear_chimerics = all_data_results.query('CorrectChimericLabelAssigned == "Unclear"').copy().reset_index(drop=True)
unclear_chimerics['chimeric_status'] = -1
false_chimerics = all_data_results.query('CorrectChimericLabelAssigned == "No"').copy().reset_index(drop=True)
false_chimerics['chimeric_status'] = 0
# Now, we will concatenate these DataFrames to create a label mapping DataFrame
label_df = pd.concat([confirmed_chimerics[['query', 'chimeric_status']], 
                      unclear_chimerics[['query', 'chimeric_status']], 
                      false_chimerics[['query', 'chimeric_status']]], 
                     ignore_index=True)
label_df


# Initialize an empty list to store DataFrames
dataframes = []
# Loop through each folder in the base directory
for folder in os.listdir(base_dir):
    print(folder)
    folder_path = os.path.join(base_dir, folder)
    
    # Check if it is a directory
    if os.path.isdir(folder_path):
        # Define the file path
        file_path = os.path.join(folder_path, 'reference_proteins.mmseqs.out')
        
        # Check if the file exists
        if os.path.isfile(file_path):
            # Read the file into a DataFrame
            df = pd.read_csv(file_path, sep='\t')  # Adjust 'sep' and 'header' as needed
            df = df.merge(label_df, on='query', how='left')
            # Optionally, add a column for the folder name to track source
            df['source_folder'] = folder
            # Append the DataFrame to the list
            dataframes.append(df)
        else:
            print(f"File not found in {folder_path}")

# Concatenate all DataFrames in the list
all_data = pd.concat(dataframes, ignore_index=True)

# Create relative start and end columns
all_data['relative_qstart'] = all_data['qstart'] / all_data['qlen']
all_data['relative_qend'] = all_data['qend'] / all_data['qlen']
all_data['relative_tstart'] = all_data['tstart'] / all_data['tlen']
all_data['relative_tend'] = all_data['tend'] / all_data['tlen']
all_data['alnlen_prop_target'] = all_data['alnlen'] / all_data['tlen']
all_data['alnlen_prop_query'] = all_data['alnlen'] / all_data['qlen']


# Display the combined DataFrame
print(all_data)
import pandas as pd
import pyranges as pr
chimera_raw_data = all_data[all_data['chimeric_status'] == 1].copy().reset_index(drop=True)

chimera_pyranges_df = chimera_raw_data.rename(columns={
    'query': 'Chromosome',
    'qstart': 'Start',
    'qend': 'End'
})




# Make sure 'Start' and 'End' are integers
chimera_pyranges_df['Start'] = chimera_pyranges_df['Start'].astype(int)
chimera_pyranges_df['End'] = chimera_pyranges_df['End'].astype(int)

# Now we can convert it to a PyRanges object
chimera_pr = pr.PyRanges(chimera_pyranges_df)
merged = chimera_pr.merge()
# Merge overlapping intervals

# Assign cluster numbers by grouping on merged intervals
merged_df = merged.as_df()
merged_df['cluster'] = range(1, len(merged_df) + 1)

# Join the cluster information back to the original dataframe
# Here we use an interval join to assign clusters to the original data
df_with_clusters = chimera_pr.join(pr.PyRanges(merged_df)).as_df()

# Print the result with clusters
print(df_with_clusters)


### Getting some features:
# Assuming your dataframe is named `df`
# Group by 'query' and 'source_folder', and then aggregate the mean, median, and std for 'qcov' and 'tcov'
# Assuming your main DataFrame is named `all_data`
all_data_filter = all_data[all_data['chimeric_status'] == 1].groupby(['query', 'source_folder']).agg(
    alnlen_prop_target_mean = ('alnlen_prop_target','mean'),
    alnlen_prop_target_median = ('alnlen_prop_target','mean'),
    alnlen_prop_target_std = ('alnlen_prop_target','mean'),
    alnlen_prop_query_mean = ('alnlen_prop_query','mean'),
    alnlen_prop_query_median = ('alnlen_prop_query','mean'),
    alnlen_prop_query_std = ('alnlen_prop_query','mean'),
    alnlen_mean = ('alnlen', 'mean'),
    alnlen_median =('alnlen','median')
).reset_index()

import matplotlib.pyplot as plt
import seaborn as sns

# Set the style for the plots
sns.set(style="whitegrid")

# List of columns to plot
columns_to_plot = ['alnlen_prop_target_mean', 'alnlen_prop_target_median', 'alnlen_prop_target_std', 'alnlen_prop_query_mean', 'alnlen_prop_query_median', 'alnlen_prop_query_std','alnlen_mean','alnlen_median']

# Define the size of the plot
plt.figure(figsize=(15, 10))

# Loop through each column and create a box plot for each, grouped by 'chimeric_status'
for i, column in enumerate(columns_to_plot, 1):
    plt.subplot(3, 3, i)  # create a 3x3 grid of subplots
    sns.boxplot(x='source_folder', y=column, data=all_data_filter)
    plt.title(f'{column} by Chimeric Status')
    plt.xlabel('Chimeric Status')
    plt.ylabel(column)

# Adjust layout for better spacing between plots
plt.tight_layout()
plt.show()

all_data_labelled = all_data.merge(label_df, on='query', how='left')


aggregated_df = all_data_labelled[all_data_labelled['source_folder'] == 'Daucus_carota'].groupby(['query', 'source_folder']).agg(
    # Coverage statistics
    qcov_mean=('qcov', 'mean'),
    qcov_median=('qcov', 'median'),
    qcov_stddev=('qcov', 'std'),
    tcov_mean=('tcov', 'mean'),
    tcov_median=('tcov', 'median'),
    tcov_stddev=('tcov', 'std'),
    # Count of hits
    counts=('qcov', 'size'),
    # Query start and end statistics
    qstart_mean=('qstart', 'mean'),
    qstart_median=('qstart', 'median'),
    qstart_stddev=('qstart', 'std'),
    qend_mean=('qend', 'mean'),
    qend_median=('qend', 'median'),
    qend_stddev=('qend', 'std'),
    # Target start and end statistics
    tstart_mean=('tstart', 'mean'),
    tstart_median=('tstart', 'median'),
    tstart_stddev=('tstart', 'std'),
    tend_mean=('tend', 'mean'),
    tend_median=('tend', 'median'),
    tend_stddev=('tend', 'std'),
    # Cluster-based features - measuring gaps between hits
    avg_gap_qstart=('qstart', lambda x: x.sort_values().diff().mean()),   # average gap between qstart positions
    avg_gap_qend=('qend', lambda x: x.sort_values().diff().mean()),       # average gap between qend positions
    avg_gap_tstart=('tstart', lambda x: x.sort_values().diff().mean()),   # average gap between tstart positions
    avg_gap_tend=('tend', lambda x: x.sort_values().diff().mean())        # average gap between tend positions
).reset_index()
aggregated_df


aggregated_df = aggregated_df.merge(label_df, on='query', how='left')
### Positive only:


# Display the resulting DataFrame
print(aggregated_df)

aggregated_df['chimeric_status'] = aggregated_df['chimeric_status'].fillna(0).copy()  # -2 for unknown status

import matplotlib.pyplot as plt
import seaborn as sns

# Set the style for the plots
sns.set(style="whitegrid")

# List of columns to plot
columns_to_plot = ['qcov_mean', 'qcov_median', 'qcov_stddev', 'tcov_mean', 'tcov_median', 'tcov_stddev', 'counts']

# Define the size of the plot
plt.figure(figsize=(15, 10))

# Loop through each column and create a box plot for each, grouped by 'chimeric_status'
for i, column in enumerate(columns_to_plot, 1):
    plt.subplot(3, 3, i)  # create a 3x3 grid of subplots
    sns.boxplot(x='chimeric_status', y=column, data=aggregated_df)
    plt.title(f'{column} by Chimeric Status')
    plt.xlabel('Chimeric Status')
    plt.ylabel(column)

# Adjust layout for better spacing between plots
plt.tight_layout()
plt.show()


import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix

# Assume your feature dataframe is named `df` and it includes a 'label' column
# with 1 for chimeric, 0 for non-chimeric, and -1 for unclear.

# Separate the labeled data
df_labeled = aggregated_df

# Define features and labels
X = df_labeled.drop(columns=['chimeric_status',"source_folder","query"])
y = df_labeled['chimeric_status']

# Split data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# Initialize the RandomForest model
model = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
model.fit(X_train, y_train)

# Make predictions on the test set
y_pred = model.predict(X_test)

# Evaluate the model
print("Classification Report:")
print(classification_report(y_test, y_pred))

print("Confusion Matrix:")
cm = confusion_matrix(y_test, y_pred)
# Display the confusion matrix
print(cm)
TN, FP, FN, TP = cm.ravel()

print(f'True Positives (TP): {TP}')
print(f'True Negatives (TN): {TN}')
print(f'False Positives (FP): {FP}')
print(f'False Negatives (FN): {FN}')




from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report
import pandas as pd
source_folders = aggregated_df['source_folder'].unique()

# Initialize a list to store the results
results = []

# Iterate over each source folder for leave-one-out cross-validation
for folder in source_folders:
    # Split the data: use all folders except the current one for training
    train_data = aggregated_df[aggregated_df['source_folder'] != folder]
    test_data = aggregated_df[aggregated_df['source_folder'] == folder]
    
    # Separate features and target labels for training and testing
    X_train = train_data.drop(columns=['query', 'source_folder', 'chimeric_status'])
    y_train = train_data['chimeric_status']
    X_test = test_data.drop(columns=['query', 'source_folder', 'chimeric_status'])
    y_test = test_data['chimeric_status']
    
    # Initialize the RandomForest model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    
    # Train the model on the training data
    model.fit(X_train, y_train)
    
    # Predict on the test data
    y_pred = model.predict(X_test)
    
    # Calculate accuracy and store the result
    accuracy = accuracy_score(y_test, y_pred)
    results.append({'source_folder': folder, 'accuracy': accuracy})
    
    # Print detailed classification report for each folder
    print(f"Results for source folder '{folder}':")
    print(classification_report(y_test, y_pred))
    print("-" * 50)

# Convert results to a DataFrame for easier analysis
results_df = pd.DataFrame(results)

# Display the final results
print("\nSummary of Accuracy for each Source Folder:")
print(results_df)