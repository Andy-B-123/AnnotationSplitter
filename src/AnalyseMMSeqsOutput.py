import pandas as pd
from icecream import ic
from tqdm import tqdm
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re
### Analysis workflow based on the Tribolium R Markdown document

def flag_large_target_spans(df):
    # Filter for spans between 15% and 85% of the query region
    df_filtered = df[(df['qstart_relative'] < 0.15) & (df['qend_relative'] > 0.85)]
    
    # Flag queries with two or more hits in the specified range
    query_flags = df_filtered.groupby('query').size().reset_index(name='count')
    query_flags['LargeTargetSpan'] = query_flags['count'] >= 1

    # Merge with original dataframe to get the 'query' column
    result_df = pd.merge(df[['query']].drop_duplicates(), query_flags[['query', 'LargeTargetSpan']], on='query', how='left')
    result_df['LargeTargetSpan'] = result_df['LargeTargetSpan'].fillna(False)
    result_df = result_df.infer_objects(copy=False)
    
    return result_df

def filter_low_hits(df):
    # Group by 'query' and count the number of hits
    grouped_df = df.groupby('query').size().reset_index(name='count_hits')
    
    # Filter groups with less than 10 hits
    filtered_df = grouped_df[grouped_df['count_hits'] < 10][['query']].drop_duplicates().reset_index(drop=True).copy()
    filtered_df['LessThan10Hits'] = True
    ic(filtered_df)
    result_df = pd.merge(df[['query']].drop_duplicates(), filtered_df, on='query', how='left')
    result_df['LessThan10Hits'] = result_df['LessThan10Hits'].fillna(False)
    return result_df

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

def check_distinct_cluster(df, cluster_distance_threshold=0.20):
    # Step 1: Group by 'query' and 'cluster' to calculate the mean of 'qmidpoint_relative'
    mean_midpoint_cluster_df = (
        df
        .groupby(['query', 'cluster'])['qmidpoint_relative']
        .mean()
        .reset_index(name='mean_midpoint_cluster')
    )
    # Step 2: Arrange the grouped data by 'query' and 'mean_midpoint_cluster'
    mean_midpoint_cluster_df = mean_midpoint_cluster_df.sort_values(by=['query', 'mean_midpoint_cluster'])
    # Step 3: Calculate distance and determine DistinctCluster
    def calculate_distance_and_proper_cluster(sub_df):
        sub_df['distance'] = sub_df['mean_midpoint_cluster'].diff().shift(-1)
        num_clusters = sub_df['cluster'].max() + 1
        proper_cluster = not (sub_df['distance'] < cluster_distance_threshold).any()
        return pd.Series({
            'num_clusters': num_clusters,
            'DistinctCluster': proper_cluster
        })
    clustering_results = (
        mean_midpoint_cluster_df
        .groupby('query')
        .apply(calculate_distance_and_proper_cluster)
        .reset_index()
    )
    # Merge the 'DistinctCluster' column back to the original dataframe
    merged_df = df.merge(clustering_results[['query', 'DistinctCluster']], on='query', how='left')
    return merged_df

def check_target_coverage(final_results_clusters, target_threshold=0.80):
    # Filter for DistinctCluster
    filtered_df = final_results_clusters[final_results_clusters['DistinctCluster']]
    
    # Step 1: Group by 'query' and 'cluster' to get the max 'tcov'
    max_target_coverage_df = (
        filtered_df
        .groupby(['query', 'cluster'])['tcov']
        .max()
        .reset_index(name='max_target_coverage_cluster')
    )

    # Step 2: Pivot wider
    pivot_df = max_target_coverage_df.pivot(index='query', columns='cluster', values='max_target_coverage_cluster').reset_index()
    
    # Step 3: Check for low coverage clusters
    pivot_df['CheckForLowCoverageCluster'] = pivot_df.apply(
        lambda row: any(row[1:] < target_threshold), axis=1
    )
    
    # Step 4: Summarize results
    final_results_clusters_annotated = pivot_df[['query', 'CheckForLowCoverageCluster']]

    # Step 5: Summarize GoodOrBad and filter for queries with no low coverage clusters
    good_or_bad_df = (
        final_results_clusters_annotated
        .groupby('query')['CheckForLowCoverageCluster']
        .sum()
        .reset_index(name='GoodOrBad')
    )
    
    high_likelihood_chimeras = good_or_bad_df[good_or_bad_df['GoodOrBad'] == 0]['query'].tolist()
    
    return high_likelihood_chimeras

def check_overlap_between_clusters(final_results_clusters, overlap_threshold=0.25):
    # Step 1: Group by 'query' and 'cluster' to calculate necessary metrics
    summary_df = (
        final_results_clusters
        .groupby(['query', 'cluster'])
        .agg(
            lowest_query_point=('qstart_relative', 'min'),
            highest_query_point=('qend_relative', 'max')
        )
        .reset_index()
    )
    summary_df['total_length'] = summary_df['highest_query_point'] - summary_df['lowest_query_point']

    # Step 2: Arrange by 'query' and 'lowest_query_point'
    summary_df = summary_df.sort_values(by=['query', 'lowest_query_point'])

    # Step 3: Calculate overlap amount and proportion
    def calculate_overlap(sub_df):
        sub_df['overlap_amount'] = sub_df['highest_query_point'] - sub_df['lowest_query_point'].shift(-1)
        sub_df['overlap_prop'] = sub_df['overlap_amount'] / sub_df['total_length'].shift(-1)
        return sub_df

    overlap_df = summary_df.groupby('query').apply(calculate_overlap).reset_index(drop=True)

    # Step 4: Filter based on overlap proportion
    filtered_overlap_df = overlap_df.groupby('query').filter(lambda x: all(x['overlap_prop'].fillna(0) <= overlap_threshold))
    
    # Step 5: Get queries with high likelihood chimerics
    higher_likelihood_chimerics = filtered_overlap_df['query'].unique().tolist()
    
    return higher_likelihood_chimerics

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]

def generate_faceted_plot(df, output_filepath):
    print("Making a nice plot of final candidate chimeras.")
    # Ensure the dataframe is copied to avoid modifications on the original
    df_copy = df.copy()
    
    # Sort the dataframe by 'query' column using natural order
    df_copy['query'] = df_copy['query'].astype(str)
    df_copy = df_copy.sort_values(by='query', key=lambda x: x.map(natural_sort_key))
    
    # Create a colormap
    max_clusters = df_copy['k'].max()
    cmap = plt.get_cmap('viridis')
    colors = cmap(np.linspace(0, 1, max_clusters))
    color_dict = {i: colors[i-1] for i in range(1, max_clusters+1)}

    # Plot with FacetGrid
    g = sns.FacetGrid(df_copy, col="query", col_wrap=10, height=2, aspect=1)

    def plot_facets(data, **kwargs):
        for index, row in data.iterrows():
            cluster_color = color_dict[row['cluster']]  # Get the color for the cluster
            plt.plot([row['qstart_relative'], row['qend_relative']], [row['tstart_relative'], row['tend_relative']], 
                     marker='o', color=cluster_color)
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xticks([])
        plt.yticks([])

    g.map_dataframe(plot_facets)

    # Set facet titles to gene names only
    for ax, title in zip(g.axes.flatten(), g.col_names):
        ax.set_title(title)

    # Remove the x and y labels
    for ax in g.axes.flatten():
        ax.set_xlabel('')
        ax.set_ylabel('')

    # Save the plot to the provided filepath
    plt.savefig(output_filepath)
    plt.close()

def convert_mmseqs_output(input_mmseqs_filepath, output_dir):
    ### Read dataframe and convert start position to relative positons
    input_mmseqs_full_results = pd.read_csv(input_mmseqs_filepath, sep='\t')
    input_mmseqs_full_results['qstart_relative'] = input_mmseqs_full_results['qstart'] / input_mmseqs_full_results['qlen']
    input_mmseqs_full_results['qend_relative'] = input_mmseqs_full_results['qend'] / input_mmseqs_full_results['qlen']
    input_mmseqs_full_results['qmidpoint_relative'] = (input_mmseqs_full_results['qstart_relative'] + input_mmseqs_full_results['qend_relative']) / 2
    input_mmseqs_full_results['tstart_relative'] = input_mmseqs_full_results['tstart'] / input_mmseqs_full_results['tlen']
    input_mmseqs_full_results['tend_relative'] = input_mmseqs_full_results['tend'] / input_mmseqs_full_results['tlen']
    # Basic stats
    total_number_of_genes_with_hits = input_mmseqs_full_results['query'].nunique()
    print(f"Total number of genes with hits in MMSeqs output: {total_number_of_genes_with_hits}" )
    # Some pre-filtering critera
    large_span_check = flag_large_target_spans(input_mmseqs_full_results)
    low_hits_check = filter_low_hits(input_mmseqs_full_results)

    merged_checks = pd.merge(large_span_check,low_hits_check)
    true_counts = merged_checks[['LargeTargetSpan', 'LessThan10Hits']].sum()
    genes_to_check = merged_checks[(merged_checks['LargeTargetSpan'] == False) & (merged_checks['LessThan10Hits'] == False)].reset_index(drop=True)
    print(f"Number of genes with good match in target database: {true_counts.iloc[0]}. Not analysing these." )
    print(f"Number of genes with less than 10 database hits: {true_counts.iloc[1]}. Not analysing these." )
    print(f"Remaining genes to analyse for chimeric pattern: {genes_to_check.shape[0]}. \nStarting clustering...")
    # Start clustering
    input_mmseqs_full_results_postChecks = pd.merge(input_mmseqs_full_results,genes_to_check[['query']])
    input_mmseqs_full_results_postChecks_clustered = hierarchial_clustering(input_mmseqs_full_results_postChecks, 'qmidpoint_relative', 0.17)
    input_mmseqs_full_results_postChecks_clustered = input_mmseqs_full_results_postChecks_clustered[input_mmseqs_full_results_postChecks_clustered['k'] != 1]
    input_mmseqs_full_results_postChecks_clustered.drop_duplicates(subset=['query'])['k'].value_counts().reset_index(name='count')

    # Filter clusters check 1
    input_mmseqs_full_results_postChecks_clustered_distance = check_distinct_cluster(input_mmseqs_full_results_postChecks_clustered, cluster_distance_threshold=0.15)
    input_mmseqs_full_results_postChecks_clustered_distance.drop_duplicates(subset=['query'])[['k',"DistinctCluster"]].value_counts().reset_index(name='count')

    # Filter clusters check 2
    high_coverage_list = check_target_coverage(input_mmseqs_full_results_postChecks_clustered_distance, target_threshold=0.80)
    high_coverage_df = pd.DataFrame(high_coverage_list, columns = ['query'])
    high_coverage_df['HighTargetCoverage'] = True
    len(high_coverage_list)
    input_mmseqs_full_results_postChecks_clustered_distance_cov = pd.merge(
        input_mmseqs_full_results_postChecks_clustered_distance,
        high_coverage_df,
        how='left'
    )
    input_mmseqs_full_results_postChecks_clustered_distance_cov['HighTargetCoverage'].fillna('false', inplace=True)
    input_mmseqs_full_results_postChecks_clustered_distance_cov.drop_duplicates(subset=['query'])[['k',"HighTargetCoverage"]].value_counts().reset_index(name='count')

    # Filter clusters check 3
    overlap_cluster_check = check_overlap_between_clusters(input_mmseqs_full_results_postChecks_clustered_distance_cov,  overlap_threshold=0.30)
    overlap_cluster_check_df = pd.DataFrame(overlap_cluster_check, columns = ['query'])
    overlap_cluster_check_df['NoClusterOverlap'] = True
    input_mmseqs_full_results_postChecks_clustered_distance_cov_overlap = pd.merge(
        input_mmseqs_full_results_postChecks_clustered_distance_cov,
        overlap_cluster_check_df,
        how='left'
    )

    input_mmseqs_full_results_postChecks_clustered_distance_cov_overlap['NoClusterOverlap'].fillna(False, inplace=True)
    input_mmseqs_full_results_postChecks_clustered_distance_cov_overlap.drop_duplicates(subset=['query'])[['k',"NoClusterOverlap"]].value_counts().reset_index(name='count')

    # Filter the DataFrame where all conditions are True
    filtered_df = input_mmseqs_full_results_postChecks_clustered_distance_cov_overlap[
        (input_mmseqs_full_results_postChecks_clustered_distance_cov_overlap['DistinctCluster'] == True) &
        (input_mmseqs_full_results_postChecks_clustered_distance_cov_overlap['HighTargetCoverage'] == True) &
        (input_mmseqs_full_results_postChecks_clustered_distance_cov_overlap['NoClusterOverlap'] == True)
    ]

    # Select only the 'query' column and drop duplicates
    result_queries = filtered_df[['query','k']].drop_duplicates()
    # Display the result
    print(result_queries)
    # Save output to tsv file
    output_results_fiile = f'{output_dir}/results.final.tsv'
    filtered_df.to_csv(output_results_fiile, sep='\t', index=False)
    # Save plot to png file
    if result_queries['query'].nunique() < 150:
        output_plot_file = f'{output_dir}/results.plot.png'
        generate_faceted_plot(filtered_df, output_plot_file)