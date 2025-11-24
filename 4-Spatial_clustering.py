# Function: Spatial clustering of circRNA on enriched pathway.
# Programmer: Zhangli Yuan, Yi Shi, 20251124.

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load All_bins.xlsx
all_bin_data = pd.read_excel('All_bins.xlsx')

# Load diff_bin_co.csv
diff_data = pd.read_csv('diff_bin_co .csv')
diff_data = diff_data[['circRNA_ID', 'bin', 'X', 'Y', 'Z']]

# Load go_bin_co.csv
go_data = pd.read_csv('go_bin_co .csv')
go_data = go_data[['circRNA_ID', 'go_path', 'gene_correlations', 'bin', 'X', 'Y', 'Z']]


unique_paths = go_data['go_path'].unique()
optimal_k_values = {}
path_silhouette_scores = {}

for path in unique_paths:
    path_data = go_data[go_data['go_path'] == path][['X', 'Y', 'Z']].drop_duplicates()
    best_k = 0
    best_score = -1
    max_k = min(len(path_data), 10)  # Maximum k value does not exceed the number of data points and 10

    if max_k > 1:  # At least two clusters are needed
        for k in range(2, max_k + 1):
            try:
                kmeans = KMeans(n_clusters=k, random_state=42)
                labels = kmeans.fit_predict(path_data)
                score = silhouette_score(path_data, labels)

                if score > best_score:
                    best_k = k
                    best_score = score

            except Exception as e:
                print(f"Path {path}, k = {k} encountered an error: {e}")

        optimal_k_values[path] = best_k
        path_silhouette_scores[path] = best_score
        print(f"Path {path}: Optimal k value is {best_k}, silhouette coefficient is {best_score}")
    else:
        print(f"Path {path}: Too few data points, unable to cluster")

print("Optimal k values:", optimal_k_values)
print("Silhouette coefficients of paths:", path_silhouette_scores)


background_silhouettes_bin = {}
background_silhouettes_diff = {}

# Cluster and calculate the silhouette coefficients of the background dataset for all bins
background_bins_data = all_bin_data[['X', 'Y', 'Z']]
for k in range(2, min(len(background_bins_data), 11)):
    kmeans_bins = KMeans(n_clusters=k, random_state=42)
    labels_bins = kmeans_bins.fit_predict(background_bins_data)
    score_bins = silhouette_score(background_bins_data, labels_bins)
    background_silhouettes_bin[f'bin_circRNA_k={k}'] = score_bins

# Cluster and calculate the silhouette coefficients of the background dataset for all differentially expressed circRNAs
background_diff_data = diff_data[['X', 'Y', 'Z']]
for k in range(2, min(len(background_diff_data), 11)):
    kmeans_diff = KMeans(n_clusters=k, random_state=42)
    labels_diff = kmeans_diff.fit_predict(background_diff_data)
    score_diff = silhouette_score(background_diff_data, labels_diff)
    background_silhouettes_diff[f'diff_circRNA_k={k}'] = score_diff

print("Background dataset silhouette coefficients (bin): ", background_silhouettes_bin)
print("Background dataset silhouette coefficients (diff): ", background_silhouettes_diff)


# Perform final clustering for each path using the optimal k values, and calculate silhouette scores
final_silhouettes = {}
for path, k in optimal_k_values.items():
    path_data = go_data[go_data['go_path'] == path][['X', 'Y', 'Z']].drop_duplicates()
    if k > 1 and len(path_data) >= k:  # Ensure there are enough data points
        kmeans_path = KMeans(n_clusters=k, random_state=42)
        labels_path = kmeans_path.fit_predict(path_data)
        score_path = silhouette_score(path_data, labels_path)
        final_silhouettes[path] = score_path
        print(f"Path {path}: Optimal k value is {k}, silhouette score is {score_path}")
    else:
        print(f"Path {path}: Not enough data points, cannot calculate silhouette score")


silhouette_cutoff = 0.5
margin = 0.02

# Calculate the silhouette coefficient threshold for the background dataset
background_bins_threshold = silhouette_cutoff - margin * 2
background_diff_threshold = silhouette_cutoff - margin


# Select the background silhouette coefficients that meet the conditions from the previous calculations
valid_bins = all(value < background_bins_threshold for value in background_silhouettes_bin.values())
valid_diff = all(value < background_diff_threshold for value in background_silhouettes_bin.values())

# Initialize the list of paths that meet the conditions
qualified_paths = {}

# Check if there is valid background bins/diff (Change cutoff and margin!!)
print(valid_bins)

print(valid_diff)

final_path = {}

if valid_bins and valid_diff:
    for path, score in final_silhouettes.items():
        if score > silhouette_cutoff:
            # Get the top 5 circRNAs with the highest Spearman correlation coefficients in the path
            top_circRNAs = go_data[go_data['go_path'] == path].nlargest(5, 'gene_correlations')
            qualified_paths[path] = top_circRNAs[['circRNA_ID', 'gene_correlations']]
            final_path[path] = score
            print(f"Path {path} meets the conditions, its silhouette coefficient is {score}")
            print(top_circRNAs[['circRNA_ID', 'gene_correlations']])
else:
    print("The silhouette coefficients of the background dataset do not meet the threshold conditions, no paths are selected")


# Flatten the dictionaries into a list of tuples
data_tuples = ([(k, 'selected pathways', v) for k, v in final_path.items()] +
               [(k, 'bin_back', v) for k, v in background_silhouettes_bin.items()] +
               [(k, 'circRNA_back', v) for k, v in background_silhouettes_diff.items()])

df = pd.DataFrame(data_tuples, columns=['Pathway', 'Category', 'Silhouette Score'])

sns.violinplot(x='Category', y='Silhouette Score', data=df)

plt.tight_layout()

plt.show()