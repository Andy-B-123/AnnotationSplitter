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
from plotnine import *
from matplotlib import cm
from pathlib import Path


def split_by_row_count(data):
    # Count the number of rows per gene
    row_counts = data['query'].value_counts()
    # Filter genes based on the row counts
    two_rows = data[data['query'].isin(row_counts[row_counts == 2].index)]
    three_rows = data[data['query'].isin(row_counts[row_counts == 3].index)]
    four_or_more_rows = data[data['query'].isin(row_counts[row_counts >= 4].index)]

    # Function to assign cluster numbers based on qstart
    def assign_cluster(df):
        df = df.sort_values(by=['query', 'qstart_relative']).copy()
        df['cluster'] = df.groupby('query').cumcount() + 1
        return df
    
    # Assign clusters to each subset
    two_rows = assign_cluster(two_rows)
    three_rows = assign_cluster(three_rows)
    four_or_more_rows = assign_cluster(four_or_more_rows)

    two_rows['Type'] = "Summary"
    three_rows['Type'] = "Summary"
    four_or_more_rows['Type'] = "Summary"
    
    return two_rows, three_rows, four_or_more_rows

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def plot_faceted_segments(df_raw, output_dir, output_filename = 'temp.png'):
    os.environ['MPLCONFIGDIR'] = output_dir + "/plot_runtime"
    print(f"Setting plotting temp folder to: {os.environ['MPLCONFIGDIR']}")
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

def create_plot(df, output_dir, output_filepath):
    # Ensure the output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    os.environ['MPLCONFIGDIR'] = output_dir + "/plot_runtime"
    print(f"Setting plotting temp folder to: {os.environ['MPLCONFIGDIR']}")
    if df.empty:
        print("The DataFrame is empty. No plot will be generated.")
        return
    # Generate the plot
    plot = (
        ggplot() +
        geom_segment(df[df['Type'] == "Raw"], 
                     aes(x="qstart_relative", xend="qend_relative", 
                         y="tstart_relative", yend="tend_relative", 
                         size="thickness"), 
                     alpha=0.25) +
        geom_segment(df[df['Type'] == "Summary"], 
                     aes(x="qstart_relative", xend="qend_relative", 
                         y="tstart_relative", yend="tend_relative", 
                         size="thickness", color='cluster')) +
        scale_color_brewer(type='qual', palette='Set1') +
        theme_bw() +
        facet_wrap("query") +
        xlab("Query position") +
        ylab("Target hit position")
    )
    
    # Save the plot to the specified file path
    full_output_path = os.path.join(output_dir, output_filepath)
    plot.save(filename=full_output_path, height=12, width=12, verbose = False)
    
    print(f"Plot saved to {full_output_path}")

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

#input_mmseqs_filepath = r"U:\\AnnotationCheckerWithStructure\\Development\\Redux2.ESF\\filtered_proteins.mmseqs.out"
#output_dir=r"U:\\AnnotationCheckerWithStructure\\Development\\Redux2.ESF\\new_plots"
#chimerics_df = convert_mmseqs_output(input_mmseqs_file, r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\results")
#chimerics_df.to_csv(r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\results\\ESF2_results.tsv", sep ='\t')
def convert_mmseqs_output(input_mmseqs_filepath, output_dir):
    ### Read dataframe and convert start position to relative positons
    basename = os.path.basename(input_mmseqs_filepath)
    input_mmseqs_full_results_raw = pd.read_csv(input_mmseqs_filepath, sep='\t')
    input_mmseqs_full_results_raw['qstart_relative'] = (input_mmseqs_full_results_raw['qstart'] / input_mmseqs_full_results_raw['qlen']).round(2)
    input_mmseqs_full_results_raw['qend_relative'] = (input_mmseqs_full_results_raw['qend'] / input_mmseqs_full_results_raw['qlen']).round(2)
    input_mmseqs_full_results_raw['tstart_relative'] = (input_mmseqs_full_results_raw['tstart'] / input_mmseqs_full_results_raw['tlen']).round(2)
    input_mmseqs_full_results_raw['tend_relative'] = (input_mmseqs_full_results_raw['tend'] / input_mmseqs_full_results_raw['tlen']).round(2)
    columns_to_keep = ['query','qstart_relative','qend_relative','tstart_relative','tend_relative','tlen','qlen']
    input_mmseqs_full_results = input_mmseqs_full_results_raw[columns_to_keep]

    # Basic stats
    total_number_of_genes_with_hits = input_mmseqs_full_results['query'].nunique()
    print(f"Total number of genes with hits in MMSeqs output: {total_number_of_genes_with_hits}" )
    print("Processing hits to remove redundant hits...")
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
    potential_misannotated.to_csv(f'{output_dir}/{basename}.filter_1.slope.tsv', sep ='\t')
    filtered_sample = potential_misannotated.groupby('query').filter(lambda x: (x['slope'] >= 1.5).all())
    print(f"Remaining after first filter: {filtered_sample['query'].nunique()}")

    # Calculate overlap percentages and add to the dataframe
    print("Removing false positives based on overlap of hits...")
    filtered_sample_with_overlap = calculate_overlap_percentage(filtered_sample)
    threshold = 10
    filtered_genes_overlap = filtered_sample_with_overlap.groupby('query').filter(lambda x: (x['overlap_percentage'] <= threshold).all())
    filtered_sample_with_overlap.to_csv(f'{output_dir}/{basename}.filter_2.overlap.tsv', sep ='\t')
    print(f"Remaining after second filter: {filtered_genes_overlap['query'].nunique()}")    

    # threshold for target coverage for all hits for a gene
    print("Removing low coverage target hits...")
    filtered_genes_overlap.loc[:,'tlength_relative'] = filtered_genes_overlap['tend_relative'] - filtered_genes_overlap['tstart_relative']
    filtered_genes_overlap.loc[:,'qlength_relative'] = filtered_genes_overlap['qend_relative'] - filtered_genes_overlap['qstart_relative']
    ### Removing any where the hits are very far apart
    query_distance_threshold = 0.5
    filtered_genes_overlap_sorted = filtered_genes_overlap.sort_values(by=['query', 'qstart_relative'])
    filtered_genes_overlap_sorted['distance'] = (filtered_genes_overlap_sorted.groupby('query')['qstart_relative'].shift(-1) - filtered_genes_overlap_sorted['qend_relative'])
    filtered_genes_overlap_sorted.to_csv(f'{output_dir}/{basename}.filter_3.distance.tsv', sep ='\t')
    low_confidence_genes_query_distance = filtered_genes_overlap_sorted[filtered_genes_overlap_sorted['distance'] >= query_distance_threshold]['query'].unique().tolist()
    ### Removing genes where any have low coverage hits
    target_threshold = 0.70
    low_confidence_genes = filtered_genes_overlap[filtered_genes_overlap['tlength_relative'] <= target_threshold]['query'].unique().tolist()
    ### Combining both filters
    query_and_target_bad_genes = low_confidence_genes_query_distance + low_confidence_genes
    filtered_genes_overlap_high_confidence = filtered_genes_overlap[~filtered_genes_overlap['query'].isin(query_and_target_bad_genes)]
    print(f"Remaining after third filter: {filtered_genes_overlap_high_confidence['query'].nunique()}")

    print("Making some pretty plots...")
    two_rows_df, three_rows_df, four_or_more_rows_df = split_by_row_count(filtered_genes_overlap_high_confidence)
    all_results_incorrect = input_mmseqs_full_results_raw[input_mmseqs_full_results_raw['query'].isin(filtered_genes_overlap_high_confidence['query'].unique())]
    all_results_incorrect = all_results_incorrect.copy()
    all_results_incorrect.loc[:,'Type'] = "Raw"
    all_results_incorrect.loc[:,'thickness'] = 1
    all_results_incorrect.loc[:,'cluster'] = None
    required_plot_columns = ["query",'qstart_relative','qend_relative','tstart_relative','tend_relative','thickness', 'cluster', "Type"]

    two_rows_df_raw = all_results_incorrect[all_results_incorrect['query'].isin(two_rows_df['query'])]
    two_rows_df_all = pd.concat([two_rows_df[required_plot_columns],two_rows_df_raw[required_plot_columns]])

    three_rows_df_raw = all_results_incorrect[all_results_incorrect['query'].isin(three_rows_df['query'])]
    three_rows_df_all = pd.concat([three_rows_df[required_plot_columns],three_rows_df_raw[required_plot_columns]])

    four_or_more_rows_df_raw = all_results_incorrect[all_results_incorrect['query'].isin(four_or_more_rows_df['query'])]
    four_or_more_rows_df_all = pd.concat([four_or_more_rows_df[required_plot_columns],four_or_more_rows_df_raw[required_plot_columns]])

    create_plot(two_rows_df_all,output_dir,f'{output_dir}/{basename}.Chimerics-2mer.png')
    create_plot(three_rows_df_all,output_dir,f'{output_dir}/{basename}.Chimerics-3mer.png')
    create_plot(four_or_more_rows_df_all,output_dir,f'{output_dir}/{basename}.Chimerics-4mer.png')

    filtered_genes_overlap_high_confidence.to_csv(f'{output_dir}/{basename}.finalResults.tsv', sep ='\t')

    return filtered_genes_overlap_high_confidence

#input_mmseqs_file = r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\ESF.mmseqs.out"
#chimerics_df = convert_mmseqs_output(input_mmseqs_file, r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\results")
#chimerics_df.to_csv(r"U:\\AnnotationCheckerWithStructure\\Development\\check_other_organisms\\results\\ESF2_results.tsv", sep ='\t')
