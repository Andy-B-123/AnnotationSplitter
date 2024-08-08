import pandas as pd
from icecream import ic
from tqdm import tqdm
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import seaborn as sns
import matplotlib.colors as mcolors
import os

def split_by_row_count(data):
    # Count the number of rows per gene
    row_counts = data['query'].value_counts()
    # Filter genes based on the row counts
    two_rows = data[data['query'].isin(row_counts[row_counts == 2].index)]
    three_rows = data[data['query'].isin(row_counts[row_counts == 3].index)]
    four_or_more_rows = data[data['query'].isin(row_counts[row_counts >= 4].index)]
    return two_rows, three_rows, four_or_more_rows

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def plot_faceted_segments(df_raw, output_filename = 'temp.png'):
    if df_raw.empty:
        print("The DataFrame is empty. No plot will be generated.")
        return
    df = df_raw.sort_values(by='query', key=lambda x: x.map(natural_sort_key))

    # Normalize the thickness values for color mapping
    norm = mcolors.Normalize(vmin=df['thickness'].min(), vmax=df['thickness'].max())
    colormap = plt.get_cmap('viridis')

    # Extract unique queries
    unique_queries = df['query'].unique()
    num_queries = len(unique_queries)
    
    # Set up the FacetGrid
    col_wrap = min(10, num_queries)
    g = sns.FacetGrid(df, col='query', col_wrap=col_wrap, height=2, sharex=True, sharey=True)

    # Precompute color values to avoid recalculating in the loop
    df['color'] = df['thickness'].apply(lambda x: colormap(norm(x)))

    # Function to draw line segments
    def draw_segments(data, **kwargs):
        for _, row in data.iterrows():
            plt.plot([row['qstart_relative'], row['qend_relative']], 
                     [row['tstart_relative'], row['tend_relative']], 
                     color=row['color'], linewidth=1)
        plt.xlim(0, 1)
        plt.ylim(0, 1)

    # Apply the plotting function to each facet
    g.map_dataframe(draw_segments)

    # Remove excess formatting
    g.set_titles("{col_name}", size=8)
    plt.subplots_adjust(top=0.9)

    plt.savefig(output_filename, bbox_inches='tight')

def plot_query_counts_histogram(df):
    """
    Plots a histogram of the counts from the 'query_counts' series.
    
    Parameters:
    - query_counts (pd.Series): A pandas Series containing the counts of each query.
    """
    query_counts = df['query'].value_counts()
    plt.figure(figsize=(10, 6))
    plt.hist(query_counts, bins=range(1, query_counts.max() + 1), edgecolor='black')
    plt.xlabel('Number of Rows')
    plt.ylabel('Frequency')
    plt.title('Histogram of Number of Rows per Query')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()

def collapse_segments_initial(data, wiggle_room=0.05, collapse_direction='left'):
    if collapse_direction == 'left':
        data = data.sort_values(by=['qstart_relative', 'qend_relative'])
    else:
        data = data.sort_values(by=['qend_relative', 'qstart_relative'], ascending=[False, False])
    
    grouped = data.groupby('query')
    output_data = []

    for gene, gene_segments in tqdm(grouped):
        gene_segments = gene_segments.reset_index(drop=True)
        current = gene_segments.iloc[0]
        count = 1
        tstart_relative = current['tstart_relative']
        tend_relative = current['tend_relative']
        
        collapsed_segments = []
        for i in range(1, len(gene_segments)):
            row = gene_segments.iloc[i]
            if (current['qstart_relative'] - wiggle_room <= row['qstart_relative'] and 
                current['qend_relative'] + wiggle_room >= row['qend_relative']):
                count += 1
                tstart_relative = min(tstart_relative, row['tstart_relative'])
                tend_relative = max(tend_relative, row['tend_relative'])
            else:
                collapsed_segments.append({
                    'qstart_relative': current['qstart_relative'],
                    'qend_relative': current['qend_relative'],
                    'tstart_relative': tstart_relative,
                    'tend_relative': tend_relative,
                    'thickness': count
                })
                current = row
                count = 1
                tstart_relative = current['tstart_relative']
                tend_relative = current['tend_relative']
        
        collapsed_segments.append({
            'qstart_relative': current['qstart_relative'],
            'qend_relative': current['qend_relative'],
            'tstart_relative': tstart_relative,
            'tend_relative': tend_relative,
            'thickness': count
        })
        
        collapsed_df = pd.DataFrame(collapsed_segments)
        collapsed_df['query'] = gene
        output_data.append(collapsed_df)
    
    collapsed_data = pd.concat(output_data, ignore_index=True)
    return collapsed_data

def collapse_segments_subsequent(data, wiggle_room=0.05, collapse_direction='left'):
    if collapse_direction == 'left':
        data = data.sort_values(by=['qstart_relative', 'qend_relative'])
    else:
        data = data.sort_values(by=['qend_relative', 'qstart_relative'], ascending=[False, False])

    grouped = data.groupby('query')
    output_data = []

    for gene, gene_segments in tqdm(grouped):
        gene_segments = gene_segments.reset_index(drop=True)
        current = gene_segments.iloc[0]
        count = current['thickness']
        tstart_relative = current['tstart_relative']
        tend_relative = current['tend_relative']
        
        collapsed_segments = []
        for i in range(1, len(gene_segments)):
            row = gene_segments.iloc[i]
            if (current['qstart_relative'] - wiggle_room <= row['qstart_relative'] and 
                current['qend_relative'] + wiggle_room >= row['qend_relative']):
                count += row['thickness']
                tstart_relative = min(tstart_relative, row['tstart_relative'])
                tend_relative = max(tend_relative, row['tend_relative'])
            else:
                collapsed_segments.append({
                    'qstart_relative': current['qstart_relative'],
                    'qend_relative': current['qend_relative'],
                    'tstart_relative': tstart_relative,
                    'tend_relative': tend_relative,
                    'thickness': count
                })
                current = row
                count = current['thickness']
                tstart_relative = current['tstart_relative']
                tend_relative = current['tend_relative']
        
        collapsed_segments.append({
            'qstart_relative': current['qstart_relative'],
            'qend_relative': current['qend_relative'],
            'tstart_relative': tstart_relative,
            'tend_relative': tend_relative,
            'thickness': count
        })
        
        collapsed_df = pd.DataFrame(collapsed_segments)
        collapsed_df['query'] = gene
        output_data.append(collapsed_df)
    
    collapsed_data = pd.concat(output_data, ignore_index=True)
    return collapsed_data

### Remove any where the overlap is too high:
def calculate_overlap_percentage(df):
    overlap_percentages = []
    for gene, gene_segments in df.groupby('query'):
        gene_segments = gene_segments.sort_values(by='qstart_relative').reset_index(drop=True)
        for i in range(len(gene_segments) - 1):
            current_end = gene_segments.loc[i, 'qend_relative']
            next_start = gene_segments.loc[i + 1, 'qstart_relative']
            segment_length = current_end - gene_segments.loc[i, 'qstart_relative']
            
            if next_start <= current_end:
                overlap_length = current_end - next_start
                if segment_length > 0:  # To avoid division by zero
                    overlap_percentage = (overlap_length / segment_length) * 100
                else:
                    overlap_percentage = 0
            else:
                overlap_percentage = 0  # No overlap
            
            overlap_percentages.append(overlap_percentage)
        
        # Ensure the last segment has an overlap of 0 (as it cannot have a next segment)
        overlap_percentages.append(0)
    
    # Add the overlap percentages to the dataframe
    df['overlap_percentage'] = overlap_percentages
    return df


def convert_mmseqs_output(input_mmseqs_filepath, output_dir):
    ### Read dataframe and convert start position to relative positons
    input_mmseqs_full_results = pd.read_csv(input_mmseqs_filepath, sep='\t')
    input_mmseqs_full_results['qstart_relative'] = (input_mmseqs_full_results['qstart'] / input_mmseqs_full_results['qlen']).round(2)
    input_mmseqs_full_results['qend_relative'] = (input_mmseqs_full_results['qend'] / input_mmseqs_full_results['qlen']).round(2)
    input_mmseqs_full_results['tstart_relative'] = (input_mmseqs_full_results['tstart'] / input_mmseqs_full_results['tlen']).round(2)
    input_mmseqs_full_results['tend_relative'] = (input_mmseqs_full_results['tend'] / input_mmseqs_full_results['tlen']).round(2)
    # Basic stats
    total_number_of_genes_with_hits = input_mmseqs_full_results['query'].nunique()
    print(f"Total number of genes with hits in MMSeqs output: {total_number_of_genes_with_hits}" )

    collapsed_data = collapse_segments_initial(input_mmseqs_full_results, collapse_direction='left')
    #plot_query_counts_histogram(collapsed_data)
    #plot_faceted_segments(collapsed_data)

    collapsed_data_2 = collapse_segments_subsequent(collapsed_data, collapse_direction='right')
    #plot_query_counts_histogram(collapsed_data_2)
    #plot_faceted_segments(collapsed_data_2)

    ### initial filtering:
    counts_of_hits = collapsed_data_2['query'].value_counts()
    normal_genes_list = counts_of_hits[counts_of_hits == 1].index.tolist()
    print(f"Number of genes without indiciation of chimeric pattern: {len(normal_genes_list)}")
    potential_misannotated_list = counts_of_hits[counts_of_hits != 1].index.tolist()
    print(f"Number of potential chimerics to assess after simple collapsing: {len(potential_misannotated_list)}")

    potential_misannotated = collapsed_data_2[collapsed_data_2['query'].isin(potential_misannotated_list)]

    ### Fine-tuning:
    """ rand_list = random.sample(potential_misannotated['query'].unique().tolist(), 150)
    potential_misannotated_sample = potential_misannotated[potential_misannotated['query'].isin(rand_list)]
    plot_faceted_segments(potential_misannotated_sample) """

    ### Remove any where the hit is linear along both the query and target(indicates that it's not a mis-annotation)
    print("Removing false positives based on slope...")
    potential_misannotated = potential_misannotated.copy()
    potential_misannotated.loc[:,'slope'] = (potential_misannotated['tend_relative'] - potential_misannotated['tstart_relative']) / (potential_misannotated['qend_relative'] - potential_misannotated['qstart_relative'])
    filtered_sample = potential_misannotated.groupby('query').filter(lambda x: (x['slope'] >= 1.5).all())
    print(f"Remaining after first filter: {filtered_sample['query'].nunique()}")

    # Calculate overlap percentages and add to the dataframe
    print("Removing false positives based on overlap of hits...")
    filtered_sample_with_overlap = calculate_overlap_percentage(filtered_sample)
    threshold = 10
    filtered_genes_overlap = filtered_sample_with_overlap.groupby('query').filter(lambda x: (x['overlap_percentage'] <= threshold).all())
    print(f"Remaining after second filter: {filtered_genes_overlap['query'].nunique()}")    

    # threshold for target coverage for all hits for a gene
    print("Removing low coverage target hits...")
    filtered_genes_overlap.loc[:,'tlength_relative'] = filtered_genes_overlap['tend_relative'] - filtered_genes_overlap['tstart_relative']
    filtered_genes_overlap.loc[:,'qlength_relative'] = filtered_genes_overlap['qend_relative'] - filtered_genes_overlap['qstart_relative']
    ### Removing any where the hits are very far apart
    query_distance_threshold = 0.5
    filtered_genes_overlap_sorted = filtered_genes_overlap.sort_values(by=['query', 'qstart_relative'])
    filtered_genes_overlap_sorted['distance'] = (filtered_genes_overlap_sorted.groupby('query')['qstart_relative'].shift(-1) - filtered_genes_overlap_sorted['qend_relative'])
    low_confidence_genes_query_distance = filtered_genes_overlap_sorted[filtered_genes_overlap_sorted['distance'] >= query_distance_threshold]['query'].unique().tolist()
    ### Removing genes where any have low coverage hits
    target_threshold = 0.70
    low_confidence_genes = filtered_genes_overlap[filtered_genes_overlap['tlength_relative'] <= target_threshold]['query'].unique().tolist()
    ### Combining both filters
    query_and_target_bad_genes = low_confidence_genes_query_distance + low_confidence_genes
    filtered_genes_overlap_high_confidence = filtered_genes_overlap[~filtered_genes_overlap['query'].isin(query_and_target_bad_genes)]

    print(f"Remaining after third filter: {filtered_genes_overlap_high_confidence['query'].nunique()}")

    print("Making some pretty plots...")
    basename = os.path.basename(input_mmseqs_filepath)
    two_rows_df, three_rows_df, four_or_more_rows_df = split_by_row_count(filtered_genes_overlap_high_confidence)
    plot_faceted_segments(two_rows_df, f'{output_dir}/{basename}.Chimerics-2mer.png')
    plot_faceted_segments(three_rows_df, f'{output_dir}/{basename}.Chimerics-3mer.png')
    plot_faceted_segments(four_or_more_rows_df, f'{output_dir}/{basename}.Chimerics-4OrMore-mer.png')
    filtered_genes_overlap_high_confidence.to_csv(f'{output_dir}/{basename}.finalResults.tsv', sep ='\t')

    return filtered_genes_overlap_high_confidence



#input_mmseqs_file = r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\ESF.mmseqs.out"
#chimerics_df = convert_mmseqs_output(input_mmseqs_file, r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\results")
#chimerics_df.to_csv(r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\results\\ESF2_results.tsv", sep ='\t')
#iput_mmseqs_file = r"U:\\AnnotationCheckerWithStructure\\Development\\BSF_test\\output_testParams.2.minAlign.out"
#iput_mmseqs_file = r"U:\\AnnotationCheckerWithStructure\\Development\\BSF_test\\output_testParams.2.minAlign.sensitivty.out"
#iput_mmseqs_file = r"U:\\AnnotationCheckerWithStructure\\Development\\BSF_test\\output_testParams.2.minAlign.sensitivty.altAli.out"




""" total_alignment_length = []
for gene, group in filtered_genes_overlap.groupby('query'):
    total_length = group['qlength_relative'].sum()
    total_alignment_length.append(total_length)

# Assuming total_alignment_length contains the calculated lengths
plt.figure(figsize=(10, 6))
plt.hist(total_alignment_length, bins=30, edgecolor='black')
plt.xlabel('Total Alignment Length')
plt.ylabel('Frequency')
plt.title('Histogram of Total Alignment Lengths')
plt.grid(axis='y', alpha=0.75)
plt.show()
# Assuming BSF_data_check is a pandas DataFrame containing the data
collapsed_data = []
for gene in BSF_data_check['query'].unique():
    chunked_data = BSF_data_check[BSF_data_check['query'] == gene]
    collapsed_chunk = collapse_segments_initial(chunked_data, geneName=gene)
    collapsed_data.append(collapsed_chunk)

collapsed_data = pd.concat(collapsed_data, ignore_index=True)

plot_faceted_segments(collapsed_data)

def full_collapse(data, wiggle_room=0.05):
    Continuously collapse overlapping segments until no more segments can be collapsed.
    while True:
        collapsed_df = collapse_segments_initial(data, wiggle_room)
        # Check if collapsing has stabilized
        if len(collapsed_df) == len(data):
            break
        data = collapsed_df

    return collapsed_df

collapsed_data_again = []
for gene in BSF_data_check['query'].unique():
    chunked_data = BSF_data_check[BSF_data_check['query'] == gene]
    collapsed_chunk = full_collapse(chunked_data, geneName=gene)
    collapsed_data.append(collapsed_chunk)

collapsed_data_again = pd.concat(collapsed_data_again, ignore_index=True)

print(collapsed_data)

plot_faceted_segments(collapsed_data)


viridis = cm.get_cmap('viridis', 12)

BSF_data_check = BSF_data_all[BSF_data_all['query'] == "XP_037905664.1"]
# Plot the line segments
plt.figure(figsize=(10, 6))
for index, row in BSF_data_check.iterrows():
    plt.plot([row['qstart_relative'], row['qend_relative']],
             [row['tstart_relative'], row['tend_relative']],
             marker='o')

plt.xlabel('Query (relative)')
plt.ylabel('Target (relative)')
plt.title('Relative Positions of Query and Target Segments')
plt.grid(True)
plt.show()



collapsed_data = collapse_segments(BSF_data_check)
# Set up the Viridis color map
norm = mcolors.Normalize(vmin=collapsed_data['thickness'].min(), vmax=collapsed_data['thickness'].max())
cmap = plt.cm.viridis
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Plot the collapsed segments with Viridis colormap
plt.figure(figsize=(10, 6))
for index, row in collapsed_data.iterrows():
    plt.plot([row['qstart_relative'], row['qend_relative']],
             [row['tstart_min'], row['tend_max']],
             marker='o', color=cmap(norm(row['thickness'])))

plt.colorbar(sm, label='Thickness')
plt.xlabel('Query (relative)')
plt.ylabel('Target (relative)')
plt.title('Collapsed Segments with Relative Positions (Colored by Thickness)')
plt.grid(True)
plt.show()


cluster_threshold=0.17
Z = linkage(BSF_data_check[["qmidpoint_relative"]], method='ward')
cluster_output = fcluster(Z, cluster_threshold, criterion='distance')
max(cluster_output)
plt.figure(figsize=(10, 7))
dendrogram(Z)
plt.show()

def hierarchial_clustering(df, column, cluster_threshold):
    results = []

    for gene in tqdm(df['query'].unique()):
        gene_data = df[df['query'] == gene].copy()
        Z = linkage(gene_data[[column]], method='centroid')
        cluster_output = fcluster(Z, cluster_threshold, criterion='distance')
        gene_data.loc[:, 'cluster'] = cluster_output
        gene_data.loc[:, 'k'] = max(cluster_output)

        results.append(gene_data)    
    return pd.concat(results)


Z = linkage(data_scaled, method='centroid')
dendrogram(Z)
cluster_sizes = [1,2,3,4,5,6,7,8]
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.preprocessing import StandardScaler

points_1 = np.random.normal(loc=0.1, scale=0.01, size=10)
points_2 = np.random.normal(loc=0.2, scale=0.01, size=10)
points_3 = np.random.normal(loc=0.8, scale=0.01, size=100)

points_concat = np.concatenate([points_1,points_2,points_3])
points_concat_clip = np.clip(points_concat, 0, 1)
data_scaled = points_concat_clip.reshape(-1, 1)
data_scaled = np.tile(data_scaled, (1,1))

plt.figure(figsize=(10, 5))
plt.plot(data_scaled, 'o', markersize=5, color='blue', label='Data Points')
plt.xlabel('Index')
plt.ylabel('Value')
plt.title('Plot of Scaled Data Points')
plt.legend()
plt.grid(True)
plt.show()


# Example data: each row represents a sample with start and end points
data_scaled = BSF_data_check[["qmidpoint_relative"]]

# Define clustering methods and distances
clustering_methods = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']
distance_thresholds = np.arange(0.05, 0.3, 0.025)

# Prepare a DataFrame to store the results
results = []
data = data_scaled
# Perform clustering and evaluate silhouette scores
for method in clustering_methods:
    for threshold in distance_thresholds:
        # Perform hierarchical clustering
        Z = linkage(data, method=method)
        clusters = fcluster(Z, t=threshold, criterion='distance')
        
        # Calculate silhouette score
        num_clusters = len(set(clusters))
        if 2 <= num_clusters < len(data):  # Silhouette score requires at least 2 clusters and at most n_samples - 1
            silhouette_avg = silhouette_score(data, clusters)
        else:
            silhouette_avg = -1  # Assign a default low score if the number of clusters is out of valid range
        
        # Store the results
        results.append({
            'Method': method,
            'Threshold': threshold,
            'Silhouette Score': silhouette_avg,
            'Number of Clusters': num_clusters
        })

# Convert results to DataFrame
results_df = pd.DataFrame(results)

# Visualize the silhouette scores
plt.figure(figsize=(14, 10))
for method in clustering_methods:
    method_results = results_df[results_df['Method'] == method]
    plt.plot(method_results['Threshold'], method_results['Silhouette Score'], label=method)

plt.xlabel('Distance Threshold')
plt.ylabel('Silhouette Score')
plt.title('Silhouette Scores for Different Clustering Methods and Thresholds')
plt.legend()
plt.show()

# Visualize the number of clusters
plt.figure(figsize=(14, 10))
for method in clustering_methods:
    method_results = results_df[results_df['Method'] == method]
    plt.plot(method_results['Threshold'], method_results['Number of Clusters'], label=method)

plt.xlabel('Distance Threshold')
plt.ylabel('Number of Clusters')
plt.title('Number of Clusters for Different Clustering Methods and Thresholds')
plt.ylim(0, 8)  # Set the y-axis limit to a maximum of 5
plt.legend()
plt.show()



### ALright, check the BSF ones:
known_results = r"U:\\AnnotationCheckerWithStructure\\Development\\BSF_test\\results.final.tsv"
known_results_df = pd.read_csv(known_results, sep ='\t')
known_results_list = known_results_df['query'].unique().tolist()
BSF_data_check = BSF_data_all[BSF_data_all['query'].isin(known_results_list)]
BSF_data_check_long = BSF_data_check
BSF_data_check = pd.concat([BSF_data_check,BSF_data_check_long])
distance_thresholds = np.arange(0.02, 0.25, 0.01)

# Function to perform clustering and plot results
def cluster_and_plot(data, gene, method='centroid', thresholds=distance_thresholds):
    try:
        print(gene)
        print(data)
        gene_data = []
        for threshold in thresholds:
            Z = linkage(data[['qmidpoint_relative']], method=method)
            clusters = fcluster(Z, t=threshold, criterion='distance')
            num_clusters = len(set(clusters))
            gene_data.append([threshold, num_clusters])
        print(gene_data)
        x_values = [item[0] for item in gene_data]
        y_values = [item[1] for item in gene_data]
        # Create the plot
        plt.figure(figsize=(10, 6))
        plt.plot(x_values, y_values, marker='o', linestyle='-', color='b')
        # Adding labels and title
        plt.axhline(y=2, color='r', linestyle='--')  # Add horizontal line at y=2
        plt.xlabel('Threshold')
        plt.ylabel('Number of Clusters')
        plt.title(gene)
        plt.grid(True)
        plt.show(block=True)        
    except KeyboardInterrupt:
        print("Interrupted by user. Exiting the function.")
        return 'c'

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

BSF_data_check_sort = BSF_data_check.sort_values(by='query', key=lambda x: x.map(natural_sort_key))

# Iterate over each query in the dataframe
for gene in BSF_data_check_sort['query'].unique():
    data = BSF_data_check_sort[BSF_data_check_sort['query'] == gene]
    cluster_and_plot(data, gene)
    key = input("Press 'c' to cancel or any other key to continue: ").strip().lower()
    if key == 'c':
        print("Exiting loop.")
        break


### OMG - intervals! I should collapse the overlapping segments to handle the annoying problem  of mis-clustering! Try this before hand
### Main thoughts - I really only want to sum them if they are complete overlaps, not partial! And then I do clustering.
### Try it.


BSF_data_check

known_results = r"U:\\AnnotationCheckerWithStructure\\Development\\BSF_test\\results.final.tsv"
known_results_df = pd.read_csv(known_results, sep ='\t')
known_results_list = known_results_df['query'].unique().tolist()
len(known_results_list)
BSF_data_check = BSF_data_all[BSF_data_all['query'].isin(known_results_list)]

import random
rand_list = random.sample(BSF_data_all['query'].unique().tolist(), 150)
BSF_data_rand = BSF_data_all[BSF_data_all['query'].isin(rand_list)] 








BSF_data_check_sort = pd.DataFrame(data)

# Ensure unique indices
BSF_data_check_sort = BSF_data_check_sort.reset_index(drop=True)

def collapse_segments_with_count(df):
    collapsed_data = []

    # Group by 'query'
    for gene in df['query'].unique():
        gene_data = df[df['query'] == gene]

        # Create an IntervalTree for each group
        tree = IntervalTree(Interval(row.qstart_relative, row.qend_relative, index)
                            for index, row in group.iterrows())
        
        # Sort the intervals by start point, then by end point
        sorted_intervals = sorted(tree, key=lambda x: (x.begin, -x.end))
        
        collapsed_intervals = []
        max_end = -float('inf')
        encapsulated_count = 0
        
        for interval in sorted_intervals:
            if interval.end > max_end:
                collapsed_intervals.append((interval.begin, interval.end, encapsulated_count, interval.data))
                max_end = interval.end
                encapsulated_count = 0  # Reset count for the next largest segment
            else:
                encapsulated_count += 1  # Increment count for each encapsulated segment
        
        # Append the collapsed intervals along with the original data
        for start, end, count, index in collapsed_intervals:
            row = group.iloc[index].copy()  # Access the row using iloc to ensure correct index handling
            row['encapsulated_count'] = count
            collapsed_data.append(row)
    
    # Create a new DataFrame from the collapsed data
    collapsed_df = pd.DataFrame(collapsed_data).reset_index(drop=True)
    return collapsed_df

# Process the DataFrame
collapsed_df = collapse_segments_with_count(BSF_data_check_sort)
print(collapsed_df) """