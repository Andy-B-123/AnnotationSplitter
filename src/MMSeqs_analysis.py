import pandas as pd
import re
import pyranges as pr
import os
import plotnine as plotnine
from plotnine import ggplot, aes, geom_segment, facet_wrap, scale_color_gradient, theme_bw, theme, element_blank,element_text, scale_x_continuous, scale_y_continuous, labs,scale_color_brewer

#gffcompare_base = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/Pocillopora_verrucosa/gffcompare_output"
#reference_mmseqs_output = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/Pocillopora_verrucosa/reference_proteins.mmseqs.out"
#helixer_mmseqs_output = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/Pocillopora_verrucosa/helixer_proteins.mmseqs.out"
#output_directory = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/Pocillopora_verrucosa/"

def process_mmseqs_data(gffcompare_base, reference_mmseqs_output, helixer_mmseqs_output, output_directory):
    ### Process gffcompare file:
    print("Starting analysis of MMSeqs results...")
    SpodRefLongestCompareHelixer = pd.read_csv(gffcompare_base + ".tracking", delimiter="\t", header=None)
    SpodRefLongestCompareHelixer.columns = ["Query_transfrags", "Query_loci", "Reference_ID", "class", "Helixer"]
    filtered_df = SpodRefLongestCompareHelixer[SpodRefLongestCompareHelixer["Reference_ID"] != "-"]
    grouped_df = filtered_df.groupby(["Reference_ID", "class"]).size().reset_index(name='number_cluster')
    # Fixed version
    totals = (grouped_df.groupby("Reference_ID", group_keys=False)
          .apply(lambda x: x.assign(total_number=x['number_cluster'].sum()))
          .reset_index(drop=True))
    contingency_table = pd.crosstab(totals['total_number'], totals['class'])
    print("Distribution of gffcompare hits across Reference and Helixer annotations:")
    print(contingency_table)

    ### Make linker tables:
    # First part: Creating linker_table for Helixer
    linker_table = (
        SpodRefLongestCompareHelixer
        .groupby("Reference_ID")
        .filter(lambda x: len(x) >= 2)  # Filtering groups with counts >= 2
        .copy()
    )
    linker_table[['gene', 'query']] = linker_table['Helixer'].str.split('|', expand=True).iloc[:, :2]
    linker_table = linker_table.drop(columns=[ "Query_transfrags",  "class"])

    # Second part: Creating linker_table_RefSeq
    linker_table_RefSeq = (
        SpodRefLongestCompareHelixer
        .groupby("Reference_ID")
        .filter(lambda x: len(x) >= 2)  # Filtering groups with counts >= 2
        .copy()
    )

    linker_table_RefSeq[['gene', 'query']] = linker_table_RefSeq['Reference_ID'].str.split('|', expand=True).iloc[:, :2]
    linker_table_RefSeq = linker_table_RefSeq.drop(columns=["Helixer", "Query_transfrags", "class"])
    linker_table_RefSeq = linker_table_RefSeq.drop_duplicates()

    ### Read in mmseqs output and filter:
    # Load the data
    Helixer_output = pd.read_csv(helixer_mmseqs_output, delimiter="\t")
    Helixer_output['source'] = "Helixer"

    RefSeq_output = pd.read_csv(reference_mmseqs_output, delimiter="\t")
    RefSeq_output['source'] = "RefSeq"

    # Merge the output data with the linker tables
    Helixer_output_matches_only = pd.merge(Helixer_output, linker_table, how='inner')
    RefSeq_output_matches_only = pd.merge(RefSeq_output, linker_table_RefSeq, how='inner')

    # Combine the results
    all_results = pd.concat([Helixer_output_matches_only, RefSeq_output_matches_only],ignore_index=True)
    check_valid_results = all_results['source'].value_counts().shape[0]
    print("Number of hits for the Reference and Helixer annotations:")
    print(all_results['source'].value_counts())
    if check_valid_results != 2:
        print("Something wrong with the mmseqs output! Exiting...")
        return None, None

    ### First filter - difference in distributions
    # Grouping by Reference_ID and source
    summary_stats = all_results.groupby(['source','Reference_ID' ]).agg(
        count=('Reference_ID', 'size'),  # Counting the number of occurrences
        qcov_median=('qcov', 'median'),
        qcov_q1=('qcov', lambda x: x.quantile(0.25)),
        qcov_q3=('qcov', lambda x: x.quantile(0.75)),
        qcov_iqr=('qcov', lambda x: x.quantile(0.75) - x.quantile(0.25)),
        tcov_median=('tcov', 'median'),
        tcov_q1=('tcov', lambda x: x.quantile(0.25)),
        tcov_q3=('tcov', lambda x: x.quantile(0.75)),
        tcov_iqr=('tcov', lambda x: x.quantile(0.75) - x.quantile(0.25))
    ).reset_index()

    # Separate Helixer and RefSeq data
    helixer_data = summary_stats[summary_stats['source'] == 'Helixer'].copy()
    refseq_data = summary_stats[summary_stats['source'] == 'RefSeq'].copy()

    # Join the data on Reference_ID
    overlap_data = pd.merge(
        helixer_data,
        refseq_data,
        on='Reference_ID',
        suffixes=('_helixer', '_refseq')
    )

    # Calculate the IQR overlap
    overlap_data['iqr_overlap'] = (
        (overlap_data['qcov_q3_helixer'] >= overlap_data['qcov_q1_refseq']) &
        (overlap_data['qcov_q3_refseq'] >= overlap_data['qcov_q1_helixer'])
    )

    # Filter distinct queries based on the conditions
    distinct_queries = overlap_data[
        (~overlap_data['iqr_overlap']) &
        (overlap_data['qcov_median_refseq'] < 0.6) &
        (overlap_data['qcov_median_helixer'] > 0.7)
    ]
    num_passing_filter_1 = len(distinct_queries["Reference_ID"])
    if num_passing_filter_1 == 0:
        print("No valid hits after first filter, exiting...")
        return None, None
    else:
        print("Filter 1 - Number of annotations with distinct hit profiles:")
        print(num_passing_filter_1)
        output_filter1_path = os.path.join(output_directory, "mmseqs_filter1.csv")
        distinct_queries.to_csv(output_filter1_path)


    ### Second filter - Matching sub-hits for RefSeq and Helixer annotations:
    # Load the data
    gffcompare_file_loci = gffcompare_base + ".loci"
    #gffcompare_file_loci = r"U:\\AnnotationCheckerWithStructure\\Helixer_time\\TestSpodCompare\\gffcompare_output.loci"
    compare_loci = pd.read_csv(gffcompare_file_loci, delimiter="\t", header=None)

    # Rename columns
    compare_loci.columns = ["Query_loci", "reference_location", "reference_annotation", "query_annotation"]

    # Extract seq_id, strand, start, and end from reference_location
    compare_loci[['seq_id', 'strand', 'start', 'end']] = compare_loci['reference_location'].str.extract(r'^([^\[]+)\[(\+|\-)\](\d+)-(\d+)$')

    # Ensure that start and end are converted to integers
    compare_loci['start'] = compare_loci['start'].astype(int)
    compare_loci['end'] = compare_loci['end'].astype(int)

    # Filter based on distinct_queries' Reference_ID
    pattern = "|".join(map(re.escape, distinct_queries['Reference_ID']))
    compare_loci_formatted_distinct = compare_loci[compare_loci['reference_annotation'].str.contains(pattern, regex=True)]

    # Filter all_results based on Reference_ID in distinct_queries
    all_results_distinct = all_results[all_results['Reference_ID'].isin(distinct_queries['Reference_ID'])].copy()

    # Separate Reference_ID into Loc_ID and Reference_ID
    all_results_distinct[['Loc_ID', 'Reference_ID']] = all_results_distinct['Reference_ID'].str.split('|', expand=True)

    # Group by Reference_ID and target, then summarize counts
    grouped_counts = (
        all_results_distinct.groupby(['Reference_ID', 'target'])
        .agg(counts=('source', pd.Series.nunique))
        .reset_index()
    )

    # Group by Reference_ID and calculate median hits per target
    median_hits = (
        grouped_counts.groupby('Reference_ID')
        .agg(median_hits_per_target=('counts', 'median'))
        .reset_index()
    )

    # Filter where median_hits_per_target is 2
    passed_hit_filter = median_hits[median_hits['median_hits_per_target'] == 2]
    if len(passed_hit_filter) == 0:
        print("No valid hits after second filter, exiting...")
        return None, None
    else:
        print("Filter 2 - Number of annotations with matching hits in RefSeq and Helixer annotations:")
        print(len(passed_hit_filter))
        output_filter2_path = os.path.join(output_directory, "mmseqs_filter2.csv")
        passed_hit_filter.to_csv(output_filter2_path)

    ### Cluster of RefSeq hits:
    # Assuming all_results_distinct is already created
    shrink_hit_factor = 0.10  # Shrink hits by this amount to make overlaps clearer
    # Filter for RefSeq and calculate shrunk coordinates
    RefIDLens = (
        all_results_distinct[all_results_distinct['source'] == "RefSeq"]
        .loc[:, ['Reference_ID', 'qlen']]
        .drop_duplicates()
    )
    Reference_data_hits = (
        all_results_distinct[all_results_distinct['source'] == "RefSeq"]
        .copy()
    )
    Reference_data_hits['qstart_shrunk'] = Reference_data_hits['qstart'] + (Reference_data_hits['alnlen'] * shrink_hit_factor).astype('int32')
    Reference_data_hits['qend_shrunk'] = Reference_data_hits['qend'] - (Reference_data_hits['alnlen'] * shrink_hit_factor).astype('int32')
    Reference_data_hits['alnlen_shrunk'] = Reference_data_hits['qend_shrunk'] - Reference_data_hits['qstart_shrunk']
    Reference_data_hits['alnlen_prop'] = Reference_data_hits['alnlen_shrunk'] / Reference_data_hits['qlen']
    # Filter based on alnlen_shrunk and alnlen_prop
    Reference_data_hits = Reference_data_hits[
        (Reference_data_hits['alnlen_shrunk'] > 1) & 
        (Reference_data_hits['alnlen_prop'] > 0.05)
    ]
    # Create a PyRanges object ensuring types are consistent
    gr = pr.PyRanges(
        chromosomes=Reference_data_hits['Reference_ID'].astype(str),
        starts=Reference_data_hits['qstart_shrunk'].astype('int32'),
        ends=Reference_data_hits['qend_shrunk'].astype('int32')
    )
    # Reduce the ranges (equivalent to GRanges::reduce)
    ### URGH, this isn't supported on Windows, only works on Linux/Mac... Work around running this and saving it:
    output_cluster_path = os.path.join(output_directory, "cluster_results.csv")
    reduced_gr = gr.merge()
    reduced_gr.to_csv(output_cluster_path)  #Workaround for this not working on Windows............
    reduced_gr = pd.read_csv(output_cluster_path)
    clusters = reduced_gr.rename(columns={'Chromosome': 'Reference_ID', 'Start': 'qstart_shrunk', 'End': 'qend_shrunk'})
    clusters['width'] = clusters['qend_shrunk'] - clusters['qstart_shrunk']
    cluster_results = (
        clusters.groupby('Reference_ID')
        .filter(lambda x: len(x) >= 2)
        .merge(RefIDLens, how='left', on='Reference_ID')
        .assign(proportion_covered=lambda x: x['width'] / x['qlen'])
        .groupby('Reference_ID', as_index=False)
        .agg(total_covered=('proportion_covered', 'sum'))
        .query('total_covered > 0.50')
        .sort_values(by='total_covered')
    )

    if len(cluster_results) == 0:
        print("No valid hits after third filter, exiting...")
        return None, None
    else:
        print("Filter 3 - Number of annotations with distinct clusters for Reference annotation hits:")
        print(len(cluster_results))
        output_filter3_path = os.path.join(output_directory, "mmseqs_filter3.csv")
        cluster_results.to_csv(output_filter3_path)

    clusters_valid = pd.merge(clusters,cluster_results, how = 'inner')
    clusters_valid['cluster_count'] = clusters_valid.groupby('Reference_ID').cumcount() + 1

    # Filter to keep rows where the value is "RefSeq" and are present in the valid clusters DF
    filtered_df = all_results_distinct[all_results_distinct['source'] == "RefSeq"]
    filtered_df = Reference_data_hits[Reference_data_hits['Reference_ID'].isin(clusters_valid['Reference_ID'])]
    # Assign each hit back to it's identified cluster for the RefSeq hits:
    results_clustered = []
    for target in clusters_valid['Reference_ID'].unique():
        cluster_count_list = []
        # Filter filtered_df to match the Reference_ID
        matching_rows = filtered_df[filtered_df['Reference_ID'] == target].copy()
        # Iterate over each matching row
        for _, filtered_row in matching_rows.iterrows():
            for _, cluster_row in clusters_valid[clusters_valid['Reference_ID'] == target].iterrows():
                # Check if the start and end in filtered_df are within the qstart_shrunk and qend_shrunk range in clusters_valid
                if cluster_row['qstart_shrunk'] <= (filtered_row['qstart_shrunk']) <= cluster_row['qend_shrunk'] and \
                cluster_row['qstart_shrunk'] <= (filtered_row['qend_shrunk']) <= cluster_row['qend_shrunk']:
                    # If they are, append the value from the 'cluster_count' column
                    cluster_count_list.append(cluster_row['cluster_count'])
        matching_rows.loc[:,'appended_cluster_count'] = cluster_count_list
        results_clustered.append(matching_rows)

    final_df_RefSeq = pd.concat(results_clustered, ignore_index=True)
    final_df_Helixer = all_results_distinct[all_results_distinct['source'] == "Helixer"]
    final_df_Helixer = final_df_Helixer[final_df_Helixer['Reference_ID'].isin(final_df_RefSeq['Reference_ID'])]
    # Final sanity check: Number of gene annotations in RefSeq and Helixer should match:
    if len(final_df_RefSeq['Reference_ID'].unique()) != len(final_df_Helixer['Reference_ID'].unique()):
        print("Something funky happening after passing third filter, aborting...")
        return None, None

    print("Plotting...")

    final_df_RefSeq['qstart_relative'] = final_df_RefSeq['qstart'] / final_df_RefSeq['qlen']
    final_df_RefSeq['qend_relative'] = final_df_RefSeq['qend'] / final_df_RefSeq['qlen']
    final_df_RefSeq['tstart_relative'] = final_df_RefSeq['tstart'] / final_df_RefSeq['tlen']
    final_df_RefSeq['tend_relative'] = final_df_RefSeq['tend'] / final_df_RefSeq['tlen']
    # Create the plot
    plot = (ggplot(final_df_RefSeq)
            + aes(x='qstart_relative', xend='qend_relative', y='tstart_relative', yend='tend_relative', color='factor(appended_cluster_count)')
            + geom_segment(alpha = 0.8)
            + facet_wrap('~Reference_ID')
            + theme_bw()
            + scale_x_continuous(breaks=[0, 0.5, 1])  # Set x-axis breaks
            + scale_y_continuous(breaks=[0, 0.5, 1])  # Set y-axis breaks
            + theme(strip_text=element_text(size=4))
            + scale_color_brewer(type='qual', palette='Set2')
            + labs( x = "Query length",
            y = "Target Length",
            color = 'Cluster'))

    output_plot_path = os.path.join(output_directory, "cluster_results.plot.svg")

    plot.save(output_plot_path)
    final_df_RefSeq.to_csv(os.path.join(output_directory, "PotentialChimeras.RefSeqs.csv"))
    final_df_Helixer.to_csv(os.path.join(output_directory, "PotentialChimeras.Helixer.csv"))
    return final_df_RefSeq, final_df_Helixer

#""" gffcompare_file = r"/scratch3/bac050/AnnotationCheckerWithStructure/Helixer_time/TestSpodCompare/gffcompare_output.tracking"
#    reference_mmseqs_output = r"/scratch3/bac050/AnnotationCheckerWithStructure/Helixer_time/TestSpodCompare/reference_proteins.mmseqs.out"
#    helixer_mmseqs_output = r"/scratch3/bac050/AnnotationCheckerWithStructure/Helixer_time/TestSpodCompare/helixer_proteins.mmseqs.out"
#    
#    gffcompare_base = r"U:\\AnnotationCheckerWithStructure\\Helixer_time\\TestSpodCompare\\gffcompare_output"
#    reference_mmseqs_output = r"U:\\AnnotationCheckerWithStructure\\Helixer_time\\TestSpodCompare\\reference_proteins.mmseqs.out"
#    helixer_mmseqs_output = r"U:\\AnnotationCheckerWithStructure\\Helixer_time\\TestSpodCompare\\helixer_proteins.mmseqs.out" 
#"""

#process_mmseqs_data(gffcompare_base,reference_mmseqs_output,helixer_mmseqs_output,r"U:\\AnnotationCheckerWithStructure\\Helixer_time\\TestSpodCompare\\")
