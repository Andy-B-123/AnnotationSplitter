from unipressed import UniprotkbClient
import pandas as pd
import time
from tqdm import tqdm

import_mmseqs_output_path = "/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/summary_stats/AllPotentialChimeras.mmseqs.slim.out"
mmseqs_df_raw = pd.read_csv(import_mmseqs_output_path, sep="\t")
confirmed_chimeras_list = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/summary_stats/cluster_results_summary_plot_data.csv")
mmseqs_df_raw

matching_rows = pd.merge(mmseqs_df_raw, confirmed_chimeras_list, left_on='query', right_on='cluster_rep', how='inner')

def fetch_recommended_names(accessions):
    results = {}
    # Fetch information for the provided accessions
    time.sleep(1)  # Pause for 1 second before each request to avoid hitting rate limits
    response = UniprotkbClient.fetch_many(accessions)  # Fetch for one accession at a time
    for entry in response:
        primary_accession = entry.get('primaryAccession', 'N/A')
        recommended_name = entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'N/A')
        # Store the primary accession and recommended name in a dictionary
        results[primary_accession] = recommended_name
    return results

# Load your mmseqs_df
mmseqs_df = matching_rows
mmseqs_df = mmseqs_df.groupby('query').head(10)
mmseqs_df.reset_index(drop=True,inplace=True)
# Extract the unique target values from the dataframe
queries = mmseqs_df['query'].unique().tolist()

name_results = pd.DataFrame(columns=['Accession', 'Recommended_Name'])
for query in tqdm(queries):
    targets_to_check = mmseqs_df[mmseqs_df['query'] == query]['target'].to_list()
    target_names = fetch_recommended_names(targets_to_check)
    target_names_df = pd.DataFrame(list(target_names.items()), columns=['Accession', 'Recommended_Name'])
    name_results = pd.concat([name_results,target_names_df])
name_results.reset_index(drop=True,inplace=True)

mmseqs_df_names = pd.concat([mmseqs_df,name_results],axis = 1)

# Save or display the updated dataframe with recommended names
mmseqs_df_names.to_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/summary_stats/cluster_results_summary_plot_data.names.csv", sep="\t", index=False)
print(mmseqs_df.head(n=20))
