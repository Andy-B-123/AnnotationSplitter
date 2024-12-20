

# Get all my data together!

import pandas as pd

# Load BLAST results
blast_df_correct = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllSamples.Processed.Confirmed.query.list.MMSeqsRaw.out", sep ='\t')
blast_df_correct['Chimera'] = 1
blast_df_incorrect = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllSamples.Processed.Incorrect.query.list.MMSeqsRaw.out", sep ='\t')
blast_df_incorrect['Chimera'] = 0
blast_df_noise = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllGenes.list.filtered.random10000.MatchInHelixer.Reference.MMSeqsRaw.out", sep ='\t')
blast_df_noise['Chimera'] = 0
blast_df = pd.concat([blast_df_correct,blast_df_incorrect,blast_df_noise])

### Feature engineering time!
# Get linkers between Ref and Helixer annotations:
RefHelixerMatch_correct = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllSamples.Processed.Confirmed.query.MatchInHelixer.csv", header=None)
RefHelixerMatch_correct.columns = ["query","Ref_query","Ref_gene"]
RefHelixerMatch_incorrect = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllSamples.Processed.Incorrect.query.MatchInHelixer.csv", header=None)
RefHelixerMatch_incorrect.columns = ["query","Ref_query","Ref_gene"]
RefHelixerMatch_noise = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllGenes.list.filtered.random10000.MatchInHelixer.tsv", header=None, sep = '\t')
RefHelixerMatch_noise.columns = ["Ref_query","query"]

# Counts of Helixer annotations per reference annotation:
all_linkers = pd.concat([RefHelixerMatch_correct,RefHelixerMatch_incorrect,RefHelixerMatch_noise])
all_linkers = all_linkers.drop_duplicates()
all_linkers_summary = all_linkers.groupby('Ref_query').agg({
    'query' : 'count'
}).reset_index()
all_linkers_summary.columns = ['query','count_helixer_annotations']
# Proportion of coverage of Helixer annotation?
blast_df_correct_Helixer = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllSamples.Processed.Confirmed.query.MatchInHelixer.query.list.MMSeqsRaw.out", sep ='\t')
blast_df_correct_Helixer['Chimera'] = 1
blast_df_incorrect_Helixer = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllSamples.Processed.Incorrect.query.MatchInHelixer.query.list.MMSeqsRaw.out", sep ='\t')
blast_df_incorrect_Helixer['Chimera'] = 0
blast_df_noise_Helixer = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/pulling_data_together_for_classifier/AllGenes.list.filtered.random10000.MatchInHelixer.Helixer.MMSeqsRaw.out", sep ='\t')
blast_df_noise_Helixer['Chimera'] = 0

blast_df_Helixer = pd.concat([blast_df_correct_Helixer,blast_df_incorrect_Helixer,blast_df_noise_Helixer])
blast_df_Helixer_summary = blast_df_Helixer.groupby('query').agg({
    "qcov" : ['mean', 'median','std'],
    'tcov' : ['mean', 'median','std'],
    'qlen' : ['mean', 'median','std'],  # Include std for qlen
    "tlen" : ['mean', 'median','std']
}).reset_index()

blast_df_Helixer_summary.columns = ['_'.join(col).strip() for col in blast_df_Helixer_summary.columns.values]
blast_df_Helixer_summary.rename(columns={'query_': 'query'}, inplace=True)
blast_df_Helixer_summary['ratio_qlen_tlen'] = blast_df_Helixer_summary['qlen_mean'] / blast_df_Helixer_summary['tlen_mean']
blast_df_Helixer_summary.columns = ['Helixer_' + col if col != 'query' else col for col in blast_df_Helixer_summary.columns]
blast_df_Helixer_summary_meta = pd.merge(blast_df_Helixer_summary, all_linkers, how = 'left', on = 'query')
blast_df_Helixer_summary_meta_RefSummary = blast_df_Helixer_summary_meta.groupby('Ref_query').agg({
    'Helixer_ratio_qlen_tlen' : ['mean','median','std'],
    "Helixer_qcov_mean" : ['mean',"median","std"],
    "Helixer_tcov_mean" : ['mean','median','std'],
}).reset_index()
blast_df_Helixer_summary_meta_RefSummary = blast_df_Helixer_summary_meta_RefSummary.fillna(0) # Fill any missing values and process back to DF:
blast_df_Helixer_summary_meta_RefSummary.columns = ['_'.join(col).strip() for col in blast_df_Helixer_summary_meta_RefSummary.columns.values]
blast_df_Helixer_summary_meta_RefSummary.rename(columns={'Ref_query_': 'Ref_query'}, inplace=True)

# Process the reference MMSeqs data:
blast_summary = blast_df.groupby('query').agg({
    "qcov" : ['mean', 'median','std'],
    'tcov' : ['mean', 'median','std'],
    'qlen' : ['mean', 'median','std'],  # Include std for qlen
    "tlen" : ['mean', 'median','std'],
    "target" : ['count']
     
}).reset_index()
blast_summary.columns = ['_'.join(col).strip() for col in blast_summary.columns.values]
blast_summary.rename(columns={'query_': 'query'}, inplace=True)
blast_summary = blast_summary.fillna(0)

### Combine all of the feature tables!
model_data = pd.merge(blast_summary,all_linkers_summary)
model_data = pd.merge(model_data,blast_df_Helixer_summary_meta_RefSummary, how = 'left', left_on='query',right_on='Ref_query')
model_data = pd.merge(model_data,blast_df[['query',"Chimera"]].drop_duplicates(), how = 'left')

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

# Assume blast_summary has a column 'is_chimeric' (1 for chimeric, 0 for non-chimeric)
X = model_data.drop(columns=['query', 'Chimera',"Ref_query"])
y = model_data['Chimera'].astype(int)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train a Random Forest classifier
clf = RandomForestClassifier(n_estimators=1000, random_state=42, )
clf.fit(X_train, y_train)

# Predict on the test set
y_pred = clf.predict(X_test)

from sklearn.metrics import confusion_matrix
# Evaluate the classifier
print(classification_report(y_test, y_pred))
# Compute the confusion matrix
cm = confusion_matrix(y_test, y_pred)

# Display the confusion matrix
print(cm)
# Boolean mask for false positives
# Step 1: Identify the False Positives
false_positive_mask = (y_pred == 1) & (y_test == 0)
# Step 2: Get the indices of the false positives
false_positive_indices = X_test[false_positive_mask].index  # Get the index of the false positives
# Step 3: Retrieve the original rows from the original dataframe
false_positives_from_original = model_data.loc[false_positive_indices]



from sklearn.tree import export_text

# Print the rules for the first tree
tree_rules = export_text(clf.estimators_[0], feature_names=list(X.columns))
print(tree_rules)

y_pred = clf.predict(X_test)  # Get the hard predictions (0 or 1)

# Compute the confusion matrix
cm = confusion_matrix(y_test, y_pred)

# Display the confusion matrix
print(cm)
TN, FP, FN, TP = cm.ravel()

print(f'True Positives (TP): {TP}')
print(f'True Negatives (TN): {TN}')
print(f'False Positives (FP): {FP}')
print(f'False Negatives (FN): {FN}')

from sklearn.metrics import precision_recall_curve
y_proba = clf.predict_proba(X_test)[:, 1]  # Probabilities for the positive class (chimeric genes)

# Calculate precision, recall, and thresholds
precision, recall, thresholds = precision_recall_curve(y_test, y_proba)

# Plot precision and recall against the threshold
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
plt.plot(thresholds, precision[:-1], 'b-', label='Precision')
plt.plot(thresholds, recall[:-1], 'g-', label='Recall')
plt.xlabel('Threshold')
plt.ylabel('Precision / Recall')
plt.title('Precision and Recall vs Threshold')
plt.legend()
plt.grid(True)
plt.show()





from sklearn.metrics import roc_curve

# Calculate FPR, TPR, and threshold values
fpr, tpr, thresholds = roc_curve(y_test, y_proba)

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, marker='.', label='Random Forest')
plt.plot([0, 1], [0, 1], linestyle='--', label='No Skill')  # Line for a random classifier
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate (Recall)')
plt.title('ROC Curve')
plt.legend()
plt.grid(True)
plt.show()

from sklearn.metrics import roc_auc_score

# Calculate AUC
auc = roc_auc_score(y_test, y_proba)
print(f'AUC: {auc:.4f}')

from sklearn.metrics import precision_recall_curve

# Calculate precision, recall, and threshold values
precision, recall, thresholds = precision_recall_curve(y_test, y_proba)

import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
plt.plot(recall, precision, marker='.', label='Random Forest')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend()
plt.grid(True)
plt.show()

from sklearn.metrics import average_precision_score

# Calculate average precision score
avg_precision = average_precision_score(y_test, y_proba)
print(f'Average Precision Score: {avg_precision:.4f}')

from treeinterpreter import treeinterpreter as ti

# Get the contributions for a specific prediction (e.g., the first test sample)
prediction, bias, contributions = ti.predict(clf, X_test[:1])

# Print out feature contributions for the prediction
for feature, contribution in zip(X.columns, contributions[0]):
    print(f"{feature}: {contribution}")


### Running this on one of the Nature samples:
resolved_chimeras = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/Supplementary_Data_1_ResolvedChimericGenes.csv")
resolved_chimeras['Chimera_type'] = "resolved"
resolved_chimeras['Chimera'] = 1
unresolved_chimeras = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/Supplementary_Data_1_Unresolved_chimeric_genes.csv")
unresolved_chimeras['Chimera_type'] = "unresolved"
unresolved_chimeras['Chimera'] = 1

chimera_list = pd.concat([resolved_chimeras,unresolved_chimeras])
chimera_list.groupby('species').agg({
    'chimeric_geneID' : "count"
})

### Looking at Sma to start, 299 chimeras.
"Tca"
linker = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/Tca/gffcompare_output.tracking", sep ='\t', header = None)
blast_hits_helixer = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/Tca/helixer_proteins.mmseqs.out", sep ='\t')
blast_hits_ref = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/ReCheckNaturePaperData/Tca/reference_proteins.mmseqs.out", sep = '\t')

linker.columns = ["Query_transfrags", "Query_loci", "Reference_ID", "class", "Helixer"]
linker_filtered = linker[linker["Reference_ID"] != "-"]
linker_filtered[["Ref_gene","Ref_query"]] = linker_filtered['Reference_ID'].str.split("|", expand=True).loc[:,:2]
linker_filtered[['gene', 'query']] = linker_filtered['Helixer'].str.split('|', expand=True).iloc[:, :2]
linker_filtered_grouped = linker_filtered.groupby("Ref_query").size().reset_index(name='count_helixer_annotations')

blast_df_Helixer_summary = blast_hits_helixer.groupby('query').agg({
    "qcov" : ['mean', 'median','std'],
    'tcov' : ['mean', 'median','std'],
    'qlen' : ['mean', 'median','std'],  # Include std for qlen
    "tlen" : ['mean', 'median','std']
}).reset_index()

blast_df_Helixer_summary.columns = ['_'.join(col).strip() for col in blast_df_Helixer_summary.columns.values]
blast_df_Helixer_summary.rename(columns={'query_': 'query'}, inplace=True)
blast_df_Helixer_summary['ratio_qlen_tlen'] = blast_df_Helixer_summary['qlen_mean'] / blast_df_Helixer_summary['tlen_mean']
blast_df_Helixer_summary.columns = ['Helixer_' + col if col != 'query' else col for col in blast_df_Helixer_summary.columns]
blast_df_Helixer_summary_meta = pd.merge(blast_df_Helixer_summary, linker_filtered, how = 'left', on = 'query')
blast_df_Helixer_summary_meta_RefSummary = blast_df_Helixer_summary_meta.groupby('Ref_query').agg({
    'Helixer_ratio_qlen_tlen' : ['mean','median','std'],
    "Helixer_qcov_mean" : ['mean',"median","std"],
    "Helixer_tcov_mean" : ['mean','median','std'],
}).reset_index()
blast_df_Helixer_summary_meta_RefSummary = blast_df_Helixer_summary_meta_RefSummary.fillna(0) # Fill any missing values and process back to DF:
blast_df_Helixer_summary_meta_RefSummary.columns = ['_'.join(col).strip() for col in blast_df_Helixer_summary_meta_RefSummary.columns.values]
blast_df_Helixer_summary_meta_RefSummary.rename(columns={'Ref_query_': 'Ref_query'}, inplace=True)

# Process the reference MMSeqs data:
blast_summary = blast_hits_ref.groupby('query').agg({
    "qcov" : ['mean', 'median','std'],
    'tcov' : ['mean', 'median','std'],
    'qlen' : ['mean', 'median','std'],  # Include std for qlen
    "tlen" : ['mean', 'median','std'],
    "target" : ['count']
}).reset_index()
blast_summary.columns = ['_'.join(col).strip() for col in blast_summary.columns.values]
blast_summary.rename(columns={'query_': 'query'}, inplace=True)
blast_summary = blast_summary.fillna(0)


### Combine all of the feature tables!
# Step 1: Create a function to check for substring matches efficiently
def find_substring_matches(main_df, lookup_df, main_col, lookup_col, value_col):
    # Create a dictionary for quick lookups based on the chimeric_geneID
    lookup_dict = pd.Series(lookup_df[value_col].values, index=lookup_df[lookup_col]).to_dict()

    # Step 2: Vectorized substring matching using str.contains
    match_series = main_df[main_col].apply(lambda x: next((lookup_dict[key] for key in lookup_dict if key in x), None))

    return match_series
model_data = pd.merge(blast_summary,linker_filtered_grouped, how = 'left', left_on = 'query', right_on='Ref_query')
model_data = model_data.drop(columns = "Ref_query")
model_data = pd.merge(model_data,blast_df_Helixer_summary_meta_RefSummary, how = 'left', left_on='query',right_on='Ref_query')
model_data['Chimera_gene'] = find_substring_matches(model_data, chimera_list, 'query', 'chimeric_geneID', 'chimeric_geneID')
model_data['Chimera'] = model_data['Chimera_gene'].apply(lambda x: 1 if pd.notna(x) else 0)
# Assume blast_summary has a column 'is_chimeric' (1 for chimeric, 0 for non-chimeric)
X_new_test = model_data.drop(columns=['query', 'Chimera',"Chimera_gene","Ref_query"])
y_new_test= model_data['Chimera'].astype(int)

# Step 2: Make predictions using the trained model
y_new_pred = clf.predict(X_new_test)

# Step 3: Evaluate the model performance
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

# Accuracy
accuracy = accuracy_score(y_new_test, y_new_pred)
print(f'Accuracy on new data: {accuracy:.4f}')

# Classification report for detailed precision, recall, and F1-score
print("Classification Report on New Data:")
print(classification_report(y_new_test, y_new_pred))

# Confusion Matrix
cm = confusion_matrix(y_new_test, y_new_pred)
print("Confusion Matrix on New Data:")
print(cm)
# Compute the confusion matrix
TN, FP, FN, TP = cm.ravel()

print(f'True Positives (TP): {TP}')
print(f'True Negatives (TN): {TN}')
print(f'False Positives (FP): {FP}')
print(f'False Negatives (FN): {FN}')

from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score

# Get predicted probabilities for the positive class
y_new_proba = clf.predict_proba(X_new_test)[:, 1]

# ROC Curve
fpr, tpr, thresholds = roc_curve(y_new_test, y_new_proba)
plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, marker='.', label='Random Forest')
plt.plot([0, 1], [0, 1], linestyle='--', label='No Skill')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve on New Data')
plt.legend()
plt.show()

# AUC score
auc = roc_auc_score(y_new_test, y_new_proba)
print(f'ROC AUC Score on new data: {auc:.4f}')

# Precision-Recall Curve
precision, recall, thresholds_pr = precision_recall_curve(y_new_test, y_new_proba)
plt.figure(figsize=(8, 6))
plt.plot(recall, precision, marker='.', label='Random Forest')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve on New Data')
plt.legend()
plt.show()




### Looking at individual samples with the model?
"Achroia_grisella"
linker = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/Achroia_grisella/gffcompare_output.tracking", sep ='\t', header = None)
blast_hits_helixer = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/Achroia_grisella/helixer_proteins.mmseqs.out", sep ='\t')
blast_hits_ref = pd.read_csv("/home/bac050/mounts/scratch3/AnnotationCheckerWithStructure/Helixer_time/AGI_genomes_trials/Achroia_grisella/reference_proteins.mmseqs.out", sep = '\t')

linker.columns = ["Query_transfrags", "Query_loci", "Reference_ID", "class", "Helixer"]
linker_filtered = linker[linker["Reference_ID"] != "-"]
linker_filtered[["Ref_gene","Ref_query"]] = linker_filtered['Reference_ID'].str.split("|", expand=True).loc[:,:2]
linker_filtered[['gene', 'query']] = linker_filtered['Helixer'].str.split('|', expand=True).iloc[:, :2]
linker_filtered_grouped = linker_filtered.groupby("Ref_query").size().reset_index(name='count_helixer_annotations')

blast_df_Helixer_summary = blast_hits_helixer.groupby('query').agg({
    "qcov" : ['mean', 'median','std'],
    'tcov' : ['mean', 'median','std'],
    'qlen' : ['mean', 'median','std'],  # Include std for qlen
    "tlen" : ['mean', 'median','std']
}).reset_index()

blast_df_Helixer_summary.columns = ['_'.join(col).strip() for col in blast_df_Helixer_summary.columns.values]
blast_df_Helixer_summary.rename(columns={'query_': 'query'}, inplace=True)
blast_df_Helixer_summary['ratio_qlen_tlen'] = blast_df_Helixer_summary['qlen_mean'] / blast_df_Helixer_summary['tlen_mean']
blast_df_Helixer_summary.columns = ['Helixer_' + col if col != 'query' else col for col in blast_df_Helixer_summary.columns]
blast_df_Helixer_summary_meta = pd.merge(blast_df_Helixer_summary, linker_filtered, how = 'left', on = 'query')
blast_df_Helixer_summary_meta_RefSummary = blast_df_Helixer_summary_meta.groupby('Ref_query').agg({
    'Helixer_ratio_qlen_tlen' : ['mean','median','std'],
    "Helixer_qcov_mean" : ['mean',"median","std"],
    "Helixer_tcov_mean" : ['mean','median','std'],
}).reset_index()
blast_df_Helixer_summary_meta_RefSummary = blast_df_Helixer_summary_meta_RefSummary.fillna(0) # Fill any missing values and process back to DF:
blast_df_Helixer_summary_meta_RefSummary.columns = ['_'.join(col).strip() for col in blast_df_Helixer_summary_meta_RefSummary.columns.values]
blast_df_Helixer_summary_meta_RefSummary.rename(columns={'Ref_query_': 'Ref_query'}, inplace=True)

# Process the reference MMSeqs data:
blast_summary = blast_hits_ref.groupby('query').agg({
    "qcov" : ['mean', 'median','std'],
    'tcov' : ['mean', 'median','std'],
    'qlen' : ['mean', 'median','std'],  # Include std for qlen
    "tlen" : ['mean', 'median','std'],
    "target" : ['count']
}).reset_index()
blast_summary.columns = ['_'.join(col).strip() for col in blast_summary.columns.values]
blast_summary.rename(columns={'query_': 'query'}, inplace=True)
blast_summary = blast_summary.fillna(0)


### Combine all of the feature tables!
# Step 1: Create a function to check for substring matches efficiently
def find_substring_matches(main_df, lookup_df, main_col, lookup_col, value_col):
    # Create a dictionary for quick lookups based on the chimeric_geneID
    lookup_dict = pd.Series(lookup_df[value_col].values, index=lookup_df[lookup_col]).to_dict()

    # Step 2: Vectorized substring matching using str.contains
    match_series = main_df[main_col].apply(lambda x: next((lookup_dict[key] for key in lookup_dict if key in x), None))

    return match_series
model_data = pd.merge(blast_summary,linker_filtered_grouped, how = 'left', left_on = 'query', right_on='Ref_query')
model_data = model_data.drop(columns = "Ref_query")
model_data = pd.merge(model_data,blast_df_Helixer_summary_meta_RefSummary, how = 'left', left_on='query',right_on='Ref_query')
model_data = model_data.drop(columns = "Ref_query")
base_df = pd.DataFrame(RefHelixerMatch_correct['Ref_query'].drop_duplicates().reset_index(drop=True))
base_df['Chimera'] = 1
model_data = pd.merge(model_data,base_df, how = 'left', left_on = 'query', right_on="Ref_query")

# Assume blast_summary has a column 'is_chimeric' (1 for chimeric, 0 for non-chimeric)
X_new_test = model_data.drop(columns=['query', 'Chimera',"Ref_query"])
y_new_test= model_data['Chimera'].fillna(0)

# Step 2: Make predictions using the trained model
y_new_pred = clf.predict(X_new_test)

# Step 3: Evaluate the model performance
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score

# Accuracy
accuracy = accuracy_score(y_new_test, y_new_pred)
print(f'Accuracy on new data: {accuracy:.4f}')

# Classification report for detailed precision, recall, and F1-score
print("Classification Report on New Data:")
print(classification_report(y_new_test, y_new_pred))

# Confusion Matrix
cm = confusion_matrix(y_new_test, y_new_pred)
print("Confusion Matrix on New Data:")
print(cm)
# Compute the confusion matrix
TN, FP, FN, TP = cm.ravel()

print(f'True Positives (TP): {TP}')
print(f'True Negatives (TN): {TN}')
print(f'False Positives (FP): {FP}')
print(f'False Negatives (FN): {FN}')

### List false-positives and check them!!!!!!!
# Identify false positives
false_positives = model_data[(y_new_test == 0) & (y_new_pred == 1)]

# Display the list of false positives, including the original columns for review
print("List of False Positives:")
print(false_positives)
